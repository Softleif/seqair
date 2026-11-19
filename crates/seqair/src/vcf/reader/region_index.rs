//! A VCF or BCF file's TBI or CSI index, its references keyed by the file's
//! contig table (`r[vcf_reader.reuse_io+2]`).

use super::error::OpenError;
use super::file::ContigTable;
use crate::bam::bgzf::VirtualOffset;
use crate::bam::csi_index::CsiIndex;
use crate::bam::index::{AlignmentIndex, BamIndex, Chunk, decompress_bgzf_file, tabix_names};
use crate::vcf::header::VcfHeader;
use seqair_types::Pos0;
use tracing::warn;

/// A parsed index and the contig-table tid of each of its references.
#[derive(Debug)]
pub(crate) struct VcfIndex {
    index: AlignmentIndex,
    /// The index's reference number of each contig-table tid that it lists.
    by_contig: Vec<Option<u32>>,
}

/// The sequence names a CSI's tabix aux block lists (`format`, `col_seq`,
/// `col_beg`, `col_end`, `meta`, `skip`, `l_nm`, names); `None` without one,
/// as a BCF's CSI has none.
fn tabix_aux_names(aux: &[u8]) -> Option<Vec<Vec<u8>>> {
    let l_nm = aux.get(24..28)?.first_chunk::<4>().map(|b| i32::from_le_bytes(*b))?;
    let names = aux.get(28..28usize.checked_add(usize::try_from(l_nm).ok()?)?)?;
    Some(tabix_names(names))
}

impl VcfIndex {
    // r[impl vcf_reader.file+2]
    // r[impl vcf_reader.reuse_io+2]
    /// Parse `bytes` (TBI or CSI, BGZF-compressed or not, told apart by
    /// content) and key its references by `contigs`, adding the contigs only
    /// the index names.
    pub(crate) fn parse(
        bytes: &[u8],
        header: &VcfHeader,
        contigs: &mut ContigTable,
    ) -> Result<Self, OpenError> {
        let decompressed;
        let data = if bytes.starts_with(&[0x1f, 0x8b]) {
            decompressed = decompress_bgzf_file(bytes)?;
            decompressed.as_slice()
        } else {
            bytes
        };
        let (index, names) = if data.starts_with(b"TBI\x01") {
            let (index, names) = BamIndex::from_tabix_named(data)?;
            (AlignmentIndex::Bai(index), Some(names))
        } else if data.starts_with(b"CSI\x01") {
            let index = CsiIndex::from_bytes(data)?;
            let names = tabix_aux_names(index.aux());
            (AlignmentIndex::Csi(index), names)
        } else {
            return Err(OpenError::UnknownIndexFormat);
        };
        let mut by_contig = Vec::new();
        let mut set = |tid: u32, reference: usize| {
            let (Ok(tid), Ok(reference)) = (usize::try_from(tid), u32::try_from(reference)) else {
                return;
            };
            if by_contig.len() <= tid {
                by_contig.resize(tid.saturating_add(1), None);
            }
            if let Some(slot) = by_contig.get_mut(tid) {
                *slot = Some(reference);
            }
        };
        match names {
            // A tabix index numbers the sequences it lists, by name.
            Some(names) => {
                for (reference, name) in names.iter().enumerate() {
                    match std::str::from_utf8(name) {
                        Ok(name) => set(contigs.add_indexed(header, name), reference),
                        Err(_) => warn!(
                            name = %name.escape_ascii(),
                            "the index names a sequence that is not UTF-8; it is not queried"
                        ),
                    }
                }
            }
            // A BCF's CSI numbers the header's contig dictionary.
            None => {
                for rid in 0..index.reference_stats().len() {
                    if let Some(id) = i32::try_from(rid).ok().and_then(|rid| contigs.by_rid(rid)) {
                        set(id.tid(), rid);
                    }
                }
            }
        }
        Ok(Self { index, by_contig })
    }

    /// The index's reference number for contig-table tid `tid`; `None` if
    /// the index does not list it (it has no records).
    pub(crate) fn reference(&self, tid: u32) -> Option<u32> {
        self.by_contig.get(usize::try_from(tid).ok()?).copied().flatten()
    }

    /// The contig-table tids the index lists.
    pub(crate) fn contigs(&self) -> impl Iterator<Item = u32> + '_ {
        (0u32..).zip(&self.by_contig).filter_map(|(tid, r)| r.map(|_| tid))
    }

    /// The chunks holding the records of reference `reference` that may
    /// overlap the 0-based half-open `[start, end)`.
    pub(crate) fn chunks(&self, reference: u32, start: i64, end: i64) -> Vec<Chunk> {
        // Positions stop at the binning's reach and at `Pos0::MAX`.
        let max = i64::try_from(self.index.max_position()).unwrap_or(i64::MAX);
        let last = end.min(max).saturating_sub(1).min(Pos0::MAX.as_i64());
        let (Ok(start), Ok(last)) = (Pos0::try_from(start), Pos0::try_from(last)) else {
            return Vec::new();
        };
        if start > last {
            return Vec::new();
        }
        self.index.query(reference, start, last)
    }

    /// Where `reference`'s positions fall in the file: leaf-window starts and
    /// the virtual offsets records near them start at, by position.
    pub(crate) fn leaf_offsets(&self, reference: u32) -> Vec<(u64, VirtualOffset)> {
        self.index.leaf_offsets(reference)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aux_names() {
        let mut aux = vec![0u8; 24];
        aux.extend_from_slice(&8i32.to_le_bytes());
        aux.extend_from_slice(b"c1\0chr2\0");
        assert_eq!(tabix_aux_names(&aux), Some(vec![b"c1".to_vec(), b"chr2".to_vec()]));
        assert_eq!(tabix_aux_names(&[]), None);
        // r[verify tabix.names]
        // An empty name numbers a reference too.
        let mut aux = vec![0u8; 24];
        aux.extend_from_slice(&9i32.to_le_bytes());
        aux.extend_from_slice(b"c1\0\0chr2\0");
        assert_eq!(tabix_aux_names(&aux), Some(vec![b"c1".to_vec(), vec![], b"chr2".to_vec()]));
        assert_eq!(tabix_aux_names(aux.get(..30).unwrap_or_default()), None);
    }
}
