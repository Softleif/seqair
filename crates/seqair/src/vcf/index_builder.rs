//! Single-pass index builder for TBI/CSI. Mirrors htslib's `hts_idx_push`
//! algorithm: accumulates bin/chunk/linear index data during writing.

use crate::bam::bgzf::VirtualOffset;
use seqair_types::SmolStr;
use std::collections::BTreeMap;

/// A chunk: contiguous range of virtual offsets in the compressed file.
#[derive(Debug, Clone, Copy)]
pub struct IndexChunk {
    pub begin: VirtualOffset,
    pub end: VirtualOffset,
}

#[derive(Debug, thiserror::Error)]
pub enum IndexError {
    #[error("input not sorted: tid={tid}, pos={pos}")]
    UnsortedInput { tid: i32, pos: u64 },

    #[error("I/O error writing index")]
    Io(#[from] std::io::Error),

    #[error("BGZF error writing index")]
    Bgzf(#[from] crate::bam::bgzf::BgzfError),

    #[error("finish() must be called before writing the index")]
    NotFinished,

    #[error("index field `{field}` count {value} exceeds i32 limit")]
    CountOverflow { field: &'static str, value: usize },
}

/// Per-reference accumulated index data.
#[derive(Debug)]
struct RefIndexBuilder {
    bins: BTreeMap<u32, Vec<IndexChunk>>,
    linear_index: Vec<VirtualOffset>,
    // Pseudo-bin metadata
    off_beg: VirtualOffset,
    off_end: VirtualOffset,
    n_mapped: u64,
    n_unmapped: u64,
}

impl RefIndexBuilder {
    fn new() -> Self {
        Self {
            bins: BTreeMap::new(),
            linear_index: Vec::new(),
            off_beg: VirtualOffset(u64::MAX),
            off_end: VirtualOffset(0),
            n_mapped: 0,
            n_unmapped: 0,
        }
    }
}

const UNSET: u64 = u64::MAX;

// r[impl index_builder.single_pass]
// r[impl index_builder.sort_validation]
// r[impl index_builder.offset_convention]
/// Incrementally builds a TBI or CSI index during writing.
///
/// Follows htslib's single-pass `hts_idx_push` algorithm. The offset
/// passed to `push()` is the virtual offset AFTER writing the record;
/// the previous call's offset is used as the start of the current record.
pub struct IndexBuilder {
    min_shift: u32,
    depth: u32,

    refs: Vec<RefIndexBuilder>,

    // State machine (mirrors htslib's hts_idx_t.z)
    last_off: VirtualOffset,
    save_off: VirtualOffset,
    last_bin: u32,
    save_bin: u32,
    last_tid: i32,
    save_tid: i32,
    last_coor: u64,
    finished: bool,
}

impl IndexBuilder {
    /// Create a new index builder.
    ///
    /// - `n_refs`: number of reference sequences
    /// - `min_shift`: minimum bin shift (14 for BAI/TBI)
    /// - `depth`: binning depth (5 for BAI/TBI)
    /// - `header_end_offset`: virtual offset after writing the file header
    pub fn new(
        n_refs: usize,
        min_shift: u32,
        depth: u32,
        header_end_offset: VirtualOffset,
    ) -> Self {
        let refs = (0..n_refs).map(|_| RefIndexBuilder::new()).collect();
        Self {
            min_shift,
            depth,
            refs,
            last_off: header_end_offset,
            save_off: header_end_offset,
            last_bin: u32::MAX,
            save_bin: u32::MAX,
            last_tid: -1,
            save_tid: -1,
            last_coor: 0,
            finished: false,
        }
    }

    /// Create a TBI-compatible builder (`min_shift=14`, depth=5).
    pub fn tbi(n_refs: usize, header_end_offset: VirtualOffset) -> Self {
        Self::new(n_refs, 14, 5, header_end_offset)
    }

    // r[impl index_builder.bai_constructor]
    /// Create a BAI-compatible builder (`min_shift=14`, depth=5).
    pub fn bai(n_refs: usize, header_end_offset: VirtualOffset) -> Self {
        Self::new(n_refs, 14, 5, header_end_offset)
    }

    // r[impl index_builder.single_pass]
    // r[impl index_builder.sort_validation]
    // r[impl index_builder.binning]
    // r[impl index_builder.chunk_flush]
    // r[impl index_builder.linear_index]
    /// Register a record that was just written.
    ///
    /// - `tid`: reference sequence index (0-based)
    /// - `beg`: 0-based start position
    /// - `end`: 0-based exclusive end position
    /// - `offset`: virtual offset AFTER writing the record
    pub fn push(
        &mut self,
        tid: i32,
        beg: u64,
        end: u64,
        offset: VirtualOffset,
    ) -> Result<(), IndexError> {
        // Validate sort order
        if tid < self.last_tid {
            return Err(IndexError::UnsortedInput { tid, pos: beg });
        }
        if tid == self.last_tid && beg < self.last_coor {
            return Err(IndexError::UnsortedInput { tid, pos: beg });
        }

        let bin = reg2bin(beg, end, self.min_shift, self.depth);

        if tid != self.last_tid {
            // New reference sequence
            if self.save_bin != u32::MAX {
                self.flush_chunk();
            }
            if self.last_tid >= 0 {
                self.flush_pseudo_bin();
            }
            // Initialize for new tid
            if let Some(r) = self.refs.get_mut(tid as usize) {
                r.off_beg = self.last_off;
            }
            self.save_off = self.last_off;
            self.save_bin = bin;
            self.save_tid = tid;
        } else if bin != self.last_bin {
            // Same tid, different bin
            self.flush_chunk();
            self.save_off = self.last_off;
            self.save_bin = bin;
        }

        // Update linear index
        if let Some(r) = self.refs.get_mut(tid as usize) {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "window index is beg/end >> min_shift (≥14); genomic coords fit well within usize on all supported targets"
            )]
            let beg_window = (beg >> self.min_shift) as usize;
            #[expect(
                clippy::cast_possible_truncation,
                reason = "window index is beg/end >> min_shift (≥14); genomic coords fit well within usize on all supported targets"
            )]
            let end_window = (end.saturating_sub(1) >> self.min_shift) as usize;
            if r.linear_index.len() <= end_window {
                r.linear_index.resize(end_window.saturating_add(1), VirtualOffset(UNSET));
            }
            for slot in r.linear_index.get_mut(beg_window..=end_window).unwrap_or(&mut []) {
                if slot.0 == UNSET {
                    *slot = self.last_off;
                }
            }
            r.n_mapped = r.n_mapped.saturating_add(1);
            if r.off_beg.0 == u64::MAX {
                r.off_beg = self.last_off;
            }
            r.off_end = offset;
        }

        self.last_bin = bin;
        self.last_tid = tid;
        self.last_coor = beg;
        self.last_off = offset;

        Ok(())
    }

    /// Register a placed-unmapped record (flag 0x4 set, `ref_id` >= 0).
    /// Same as `push()` but increments `n_unmapped` instead of `n_mapped` in the
    /// pseudo-bin, matching htslib's `hts_idx_push` behavior for unmapped reads.
    pub fn push_unmapped(
        &mut self,
        tid: i32,
        beg: u64,
        end: u64,
        offset: VirtualOffset,
    ) -> Result<(), IndexError> {
        self.push(tid, beg, end, offset)?;
        // push() incremented n_mapped; correct it: undo mapped, add unmapped
        if let Some(r) = self.refs.get_mut(tid as usize) {
            r.n_mapped = r.n_mapped.saturating_sub(1);
            r.n_unmapped = r.n_unmapped.saturating_add(1);
        }
        Ok(())
    }

    fn flush_chunk(&mut self) {
        if self.save_bin == u32::MAX || self.save_tid < 0 {
            return;
        }
        let chunk = IndexChunk { begin: self.save_off, end: self.last_off };
        if let Some(r) = self.refs.get_mut(self.save_tid as usize) {
            r.bins.entry(self.save_bin).or_default().push(chunk);
        }
    }

    // r[impl index_builder.pseudo_bin]
    fn flush_pseudo_bin(&mut self) {
        if self.last_tid < 0 {
            return;
        }
        let pseudo_bin_id = pseudo_bin(self.depth);
        if let Some(r) = self.refs.get_mut(self.last_tid as usize) {
            let chunks = vec![
                // Chunk 0: virtual offset range
                IndexChunk { begin: r.off_beg, end: r.off_end },
                // Chunk 1: mapped/unmapped counts packed as virtual offsets
                IndexChunk { begin: VirtualOffset(r.n_mapped), end: VirtualOffset(r.n_unmapped) },
            ];
            r.bins.insert(pseudo_bin_id, chunks);
        }
    }

    // r[impl index_builder.linear_backfill]
    /// Finalize the index after all records have been written.
    pub fn finish(&mut self, final_offset: VirtualOffset) -> Result<(), IndexError> {
        self.last_off = final_offset;
        if self.save_bin != u32::MAX {
            self.flush_chunk();
        }
        if self.last_tid >= 0 {
            self.flush_pseudo_bin();
        }

        // Backfill linear index gaps (right-to-left propagation)
        for r in &mut self.refs {
            backfill_linear_index(&mut r.linear_index);
        }
        self.finished = true;
        Ok(())
    }

    // r[impl index_builder.tbi_format]
    /// Write TBI format to a writer. The output is BGZF-compressed.
    /// Only references with actual records are included (matching bcftools behavior).
    pub fn write_tbi<W: std::io::Write>(
        &self,
        writer: W,
        contig_names: &[SmolStr],
    ) -> Result<(), IndexError> {
        use crate::bam::bgzf_writer::BgzfWriter;

        // r[impl index_builder.tbi_empty_refs]
        // Collect only references that have data (non-empty bins)
        let active_refs: Vec<(usize, &RefIndexBuilder)> =
            self.refs.iter().enumerate().filter(|(_, r)| !r.bins.is_empty()).collect();

        let mut bgzf = BgzfWriter::new(writer);
        let mut buf = Vec::new();

        // Magic
        buf.extend_from_slice(b"TBI\x01");
        // n_ref — only count references with data
        write_i32(&mut buf, count_i32(active_refs.len(), "n_ref")?);
        // format = 2 (VCF)
        write_i32(&mut buf, 2);
        // col_seq = 1, col_beg = 2, col_end = 0
        write_i32(&mut buf, 1);
        write_i32(&mut buf, 2);
        write_i32(&mut buf, 0);
        // meta = '#' = 35
        write_i32(&mut buf, 35);
        // skip = 0
        write_i32(&mut buf, 0);

        // Concatenated null-terminated sequence names — only for active refs
        let mut names_buf = Vec::new();
        for &(idx, _) in &active_refs {
            if let Some(name) = contig_names.get(idx) {
                names_buf.extend_from_slice(name.as_bytes());
                names_buf.push(0);
            }
        }
        write_i32(&mut buf, count_i32(names_buf.len(), "l_nm")?);
        buf.extend_from_slice(&names_buf);

        // Per-reference bin/chunk/linear data — only for active refs
        for &(_, r) in &active_refs {
            write_ref_index(&mut buf, r)?;
        }

        bgzf.write_all(&buf)?;
        bgzf.finish()?;
        Ok(())
    }

    // r[impl index_builder.bai_format]
    // r[impl index_builder.bai_all_refs]
    // r[impl index_builder.bai_write]
    /// Write BAI format to a writer. The output is uncompressed (not BGZF-wrapped).
    ///
    /// `n_refs` is the total number of reference sequences in the BAM header — BAI
    /// includes an entry for every reference, even those with no records (`n_bin=0`, `n_intv=0`).
    pub fn write_bai<W: std::io::Write>(
        &self,
        mut writer: W,
        n_refs: usize,
    ) -> Result<(), IndexError> {
        if !self.finished {
            return Err(IndexError::NotFinished);
        }
        let mut buf = Vec::new();

        // Magic
        buf.extend_from_slice(b"BAI\x01");
        // n_ref — always the full header count (BAI uses positional lookup by tid)
        write_i32(&mut buf, count_i32(n_refs, "n_ref")?);

        // Per-reference data: include all refs, even empty ones
        for tid in 0..n_refs {
            match self.refs.get(tid) {
                Some(r) if !r.bins.is_empty() => write_ref_index(&mut buf, r)?,
                _ => {
                    // Empty ref: n_bin=0, n_intv=0
                    write_i32(&mut buf, 0); // n_bin
                    write_i32(&mut buf, 0); // n_intv
                }
            }
        }

        writer.write_all(&buf).map_err(IndexError::Io)?;
        Ok(())
    }

    // r[impl csi.write]
    // r[impl csi.write_format]
    // r[impl csi.write_loffset]
    /// Write CSI format to a writer. The output is uncompressed.
    ///
    /// `n_refs` is the total number of reference sequences in the BAM header —
    /// CSI includes an entry for every reference, even those with no records.
    pub fn write_csi<W: std::io::Write>(
        &self,
        mut writer: W,
        n_refs: usize,
        aux: &[u8],
    ) -> Result<(), IndexError> {
        if !self.finished {
            return Err(IndexError::NotFinished);
        }
        let mut buf = Vec::new();

        // r[impl csi.write_format]
        buf.extend_from_slice(b"CSI\x01");
        write_i32(&mut buf, count_i32(self.min_shift as usize, "min_shift")?);
        write_i32(&mut buf, count_i32(self.depth as usize, "depth")?);
        write_i32(&mut buf, count_i32(aux.len(), "l_aux")?);
        buf.extend_from_slice(aux);
        write_i32(&mut buf, count_i32(n_refs, "n_ref")?);

        for tid in 0..n_refs {
            match self.refs.get(tid) {
                Some(r) if !r.bins.is_empty() => write_csi_ref_index(&mut buf, r)?,
                _ => {
                    write_i32(&mut buf, 0); // n_bin
                }
            }
        }

        writer.write_all(&buf).map_err(IndexError::Io)?;
        Ok(())
    }

    // r[impl csi.write_tabix_aux]
    /// Write CSI format with tabix metadata in the aux block (for VCF/BED/GFF).
    /// BGZF-compressed output.
    pub fn write_csi_tabix<W: std::io::Write>(
        &self,
        writer: W,
        contig_names: &[SmolStr],
    ) -> Result<(), IndexError> {
        use crate::bam::bgzf_writer::BgzfWriter;

        // Build tabix aux block
        let mut aux = Vec::new();
        // format = 2 (VCF)
        write_i32(&mut aux, 2);
        // col_seq = 1, col_beg = 2, col_end = 0
        write_i32(&mut aux, 1);
        write_i32(&mut aux, 2);
        write_i32(&mut aux, 0);
        // meta = '#' = 35
        write_i32(&mut aux, 35);
        // skip = 0
        write_i32(&mut aux, 0);

        // Collect only active references
        let active_refs: Vec<(usize, &RefIndexBuilder)> =
            self.refs.iter().enumerate().filter(|(_, r)| !r.bins.is_empty()).collect();

        let mut names_buf = Vec::new();
        for &(idx, _) in &active_refs {
            if let Some(name) = contig_names.get(idx) {
                names_buf.extend_from_slice(name.as_bytes());
                names_buf.push(0);
            }
        }
        write_i32(&mut aux, count_i32(names_buf.len(), "l_nm")?);
        aux.extend_from_slice(&names_buf);

        // Write into BGZF
        let mut bgzf = BgzfWriter::new(writer);
        let mut buf = Vec::new();

        buf.extend_from_slice(b"CSI\x01");
        write_i32(&mut buf, count_i32(self.min_shift as usize, "min_shift")?);
        write_i32(&mut buf, count_i32(self.depth as usize, "depth")?);
        write_i32(&mut buf, count_i32(aux.len(), "l_aux")?);
        buf.extend_from_slice(&aux);
        write_i32(&mut buf, count_i32(active_refs.len(), "n_ref")?);

        for &(_, r) in &active_refs {
            write_csi_ref_index(&mut buf, r)?;
        }

        bgzf.write_all(&buf)?;
        bgzf.finish()?;
        Ok(())
    }

    /// Number of references in the index.
    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }
}

/// Safe i32 cast for index field counts. These are structurally bounded
/// (max ~37451 bins for depth=5, linear index ≤ `genome_size` / 16KiB), but
/// we validate at the serialization boundary to catch corruption.
fn count_i32(n: usize, field: &'static str) -> Result<i32, IndexError> {
    i32::try_from(n).map_err(|_| IndexError::CountOverflow { field, value: n })
}

// r[impl csi.write_loffset]
/// Write one reference's bin data in CSI format (with per-bin loffset, no linear index).
fn write_csi_ref_index(buf: &mut Vec<u8>, r: &RefIndexBuilder) -> Result<(), IndexError> {
    write_i32(buf, count_i32(r.bins.len(), "n_bin")?);

    for (&bin_id, chunks) in &r.bins {
        write_u32(buf, bin_id);
        // loffset: virtual offset of the first record in this bin
        let loffset = chunks.first().map_or(0, |c| c.begin.0);
        write_u64(buf, loffset);
        write_i32(buf, count_i32(chunks.len(), "n_chunk")?);
        for chunk in chunks {
            write_u64(buf, chunk.begin.0);
            write_u64(buf, chunk.end.0);
        }
    }
    Ok(())
}

/// Write one reference's bin data in BAI format (no per-bin loffset, with linear index).
fn write_ref_index(buf: &mut Vec<u8>, r: &RefIndexBuilder) -> Result<(), IndexError> {
    write_i32(buf, count_i32(r.bins.len(), "n_bin")?);

    for (&bin_id, chunks) in &r.bins {
        write_u32(buf, bin_id);
        write_i32(buf, count_i32(chunks.len(), "n_chunk")?);
        for chunk in chunks {
            write_u64(buf, chunk.begin.0);
            write_u64(buf, chunk.end.0);
        }
    }

    write_i32(buf, count_i32(r.linear_index.len(), "n_intv")?);
    for &offset in &r.linear_index {
        let val = if offset.0 == UNSET { 0 } else { offset.0 };
        write_u64(buf, val);
    }
    Ok(())
}

// r[impl index_builder.binning]
/// Compute the bin for a genomic interval.
/// Matches htslib's `hts_reg2bin(beg, end, min_shift, n_lvls)`.
#[expect(
    clippy::cast_possible_truncation,
    reason = "bin values are bounded by the TBI/CSI spec (e.g. depth=5 → max bin 37449), fits in u32"
)]
pub fn reg2bin(beg: u64, end: u64, min_shift: u32, depth: u32) -> u32 {
    let end = end.saturating_sub(1);
    let mut s = min_shift;
    // Leaf-level offset: ((1 << depth*3) - 1) / 7
    let mut t = ((1u64 << depth.saturating_mul(3)).saturating_sub(1)).checked_div(7).unwrap_or(0);

    let mut l = depth;
    while l > 0 {
        if beg >> s == end >> s {
            return t.saturating_add(beg >> s) as u32;
        }
        // htslib: decrement l first, then update s and t with new l
        l = l.saturating_sub(1);
        s = s.saturating_add(3);
        t = t.saturating_sub(1u64 << l.saturating_mul(3));
    }
    0
}

/// Pseudo-bin ID for the given depth (one past the max valid bin).
/// For depth=5: 37450.
#[expect(
    clippy::cast_possible_truncation,
    reason = "pseudo-bin is one past the max valid bin, bounded by depth (e.g. depth=5 → 37450), fits in u32"
)]
fn pseudo_bin(depth: u32) -> u32 {
    // max_bin_id = ((1 << (depth+1)*3) - 1) / 7
    let bits = depth.saturating_add(1).saturating_mul(3);
    let max_id = ((1u64 << bits).saturating_sub(1)).checked_div(7).unwrap_or(0);
    max_id.saturating_add(1) as u32
}

// r[impl index_builder.linear_backfill]
fn backfill_linear_index(linear: &mut [VirtualOffset]) {
    if linear.is_empty() {
        return;
    }
    // Right-to-left: each unset slot gets the next valid value
    let mut last_valid = VirtualOffset(0);
    for slot in linear.iter_mut().rev() {
        if slot.0 == UNSET {
            *slot = last_valid;
        } else {
            last_valid = *slot;
        }
    }
}

fn write_i32(buf: &mut Vec<u8>, val: i32) {
    buf.extend_from_slice(&val.to_le_bytes());
}

fn write_u32(buf: &mut Vec<u8>, val: u32) {
    buf.extend_from_slice(&val.to_le_bytes());
}

fn write_u64(buf: &mut Vec<u8>, val: u64) {
    buf.extend_from_slice(&val.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify index_builder.binning]
    #[test]
    fn reg2bin_leaf_level() {
        // Position 0..16384 should be in bin 4681 (first leaf)
        assert_eq!(reg2bin(0, 16384, 14, 5), 4681);
        // Position 16384..32768 should be in bin 4682
        assert_eq!(reg2bin(16384, 32768, 14, 5), 4682);
    }

    // r[verify index_builder.binning]
    #[test]
    fn reg2bin_spanning_bins() {
        // Interval spanning two leaf bins → goes to parent
        assert_eq!(reg2bin(0, 32768, 14, 5), 585);
    }

    // r[verify index_builder.binning]
    #[test]
    fn reg2bin_root() {
        // Very large interval → root bin 0
        assert_eq!(reg2bin(0, 512_000_000, 14, 5), 0);
    }

    #[test]
    fn pseudo_bin_id() {
        assert_eq!(pseudo_bin(5), 37450);
    }

    // r[verify index_builder.linear_backfill]
    #[test]
    fn backfill_linear() {
        let mut linear = vec![
            VirtualOffset(100),
            VirtualOffset(UNSET),
            VirtualOffset(UNSET),
            VirtualOffset(300),
            VirtualOffset(UNSET),
        ];
        backfill_linear_index(&mut linear);
        assert_eq!(linear[0].0, 100);
        assert_eq!(linear[1].0, 300); // backfilled from slot 3
        assert_eq!(linear[2].0, 300);
        assert_eq!(linear[3].0, 300);
        assert_eq!(linear[4].0, 0); // no valid slot to the right → 0
    }

    // r[verify index_builder.sort_validation]
    #[test]
    fn rejects_unsorted_tid() {
        let mut builder = IndexBuilder::tbi(2, VirtualOffset(0));
        builder.push(1, 100, 200, VirtualOffset(100)).unwrap();
        let result = builder.push(0, 50, 100, VirtualOffset(200));
        assert!(result.is_err());
    }

    // r[verify index_builder.sort_validation]
    #[test]
    fn rejects_unsorted_pos() {
        let mut builder = IndexBuilder::tbi(1, VirtualOffset(0));
        builder.push(0, 200, 300, VirtualOffset(100)).unwrap();
        let result = builder.push(0, 100, 200, VirtualOffset(200));
        assert!(result.is_err());
    }

    // r[verify index_builder.single_pass]
    // r[verify index_builder.chunk_flush]
    #[test]
    fn basic_index_building() {
        let mut builder = IndexBuilder::tbi(1, VirtualOffset(0));
        // Push a few records in the same bin
        builder.push(0, 100, 200, VirtualOffset(1000)).unwrap();
        builder.push(0, 150, 250, VirtualOffset(2000)).unwrap();
        // Push record in a different bin
        builder.push(0, 20000, 20100, VirtualOffset(3000)).unwrap();
        builder.finish(VirtualOffset(4000)).unwrap();

        // Should have bins + pseudo-bin for ref 0
        let r = &builder.refs[0];
        assert!(!r.bins.is_empty());
        // Pseudo-bin should exist
        assert!(r.bins.contains_key(&37450));
    }

    // r[verify index_builder.tbi_format]
    #[test]
    fn write_tbi_produces_valid_bgzf() {
        let mut builder = IndexBuilder::tbi(1, VirtualOffset(0));
        builder.push(0, 100, 200, VirtualOffset(1000)).unwrap();
        builder.finish(VirtualOffset(2000)).unwrap();

        let names = vec![SmolStr::from("chr1")];
        let mut output = Vec::new();
        builder.write_tbi(&mut output, &names).unwrap();

        // Output should be BGZF-compressed — starts with gzip magic
        assert!(output.len() > 28);
        assert_eq!(output[0], 0x1f);
        assert_eq!(output[1], 0x8b);

        // Decompress and check TBI magic
        let mut reader = crate::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(output));
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();
        assert_eq!(&decompressed[..4], b"TBI\x01");
    }
}
