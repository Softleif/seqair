//! BAM file writer with optional BAI index co-production.
//!
//! [`BamWriter`] serializes [`OwnedBamRecord`] values
//! into BGZF-compressed BAM format and optionally co-produces a BAI index during writing.

use super::bgzf::BgzfError;
use super::bgzf_writer::BgzfWriter;
use super::header::{BamHeader, BamHeaderError, TargetInfoAccess};
use super::owned_record::{OwnedBamRecord, OwnedRecordError};
use crate::vcf::index_builder::{IndexBuilder, IndexError};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use thiserror::Error;

/// Maximum record size the writer will accept (matching reader's 2 MiB cap).
const MAX_RECORD_SIZE: usize = 2 * 1024 * 1024;

// Limits matching the reader (header.rs) — the writer should never produce
// data that the reader would reject.
const MAX_HEADER_TEXT: usize = 256 * 1024 * 1024; // 256 MiB (same as reader)
const MAX_REFERENCES: usize = 1_000_000; // same as reader
const MAX_REF_NAME: usize = 256 * 1024; // 256 KiB (same as reader)

// r[impl bam_writer.error_type]
/// Errors from BAM writing operations.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum BamWriteError {
    /// I/O error from the underlying stream.
    #[error("I/O error")]
    Io {
        #[from]
        source: io::Error,
    },

    /// BGZF compression error.
    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    /// Record serialization error.
    #[error("record error")]
    Record {
        #[from]
        source: OwnedRecordError,
    },

    /// Index building error (e.g. unsorted input).
    #[error("index error")]
    Index {
        #[from]
        source: IndexError,
    },

    /// Header serialization error.
    #[error("header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    /// Writer was poisoned after a previous write failure.
    #[error("writer poisoned after previous error")]
    Poisoned,

    /// Record exceeds the 2 MiB size limit.
    #[error("record size {size} exceeds {MAX_RECORD_SIZE} byte limit")]
    RecordTooLarge { size: usize },

    /// Mapped record (flags & 0x4 == 0) has `ref_id` == -1.
    #[error("mapped record has ref_id == -1 (structurally invalid)")]
    MappedWithoutReference,
}

// r[impl bam_writer.create_from_path]
// r[impl bam_writer.header_eager]
/// BAM file writer with optional BAI index co-production.
///
/// The header is written eagerly during construction. Records are written one
/// at a time via [`write`](Self::write). Call [`finish`](Self::finish) to flush
/// and retrieve the inner writer and optional index builder.
pub struct BamWriter<W: Write> {
    bgzf: BgzfWriter<W>,
    index: Option<IndexBuilder>,
    buf: Vec<u8>,
    poisoned: bool,
}

impl BamWriter<BufWriter<File>> {
    /// Create a BAM writer writing to a file path.
    ///
    /// Uses the default compression level (6). If `build_index` is true,
    /// a BAI index is co-produced during writing.
    pub fn from_path(
        path: &Path,
        header: &BamHeader,
        build_index: bool,
    ) -> Result<Self, BamWriteError> {
        Self::from_path_with_level(path, header, build_index, 6)
    }

    // r[impl bam_writer.compression_level]
    /// Create a BAM writer with a specific compression level.
    pub fn from_path_with_level(
        path: &Path,
        header: &BamHeader,
        build_index: bool,
        level: i32,
    ) -> Result<Self, BamWriteError> {
        let file = File::create(path)?;
        let inner = BufWriter::new(file);
        let bgzf = BgzfWriter::with_compression_level(inner, level);
        Self::new_inner(bgzf, header, build_index)
    }
}

// r[impl bam_writer.create_from_stdout]
impl BamWriter<io::Stdout> {
    /// Create a BAM writer writing to stdout. Index co-production is not supported.
    pub fn from_stdout(header: &BamHeader) -> Result<Self, BamWriteError> {
        let bgzf = BgzfWriter::new(io::stdout());
        Self::new_inner(bgzf, header, false)
    }
}

impl<W: Write> BamWriter<W> {
    /// Construct from any `BgzfWriter`. Writes the BAM header eagerly.
    pub fn new_inner(
        mut bgzf: BgzfWriter<W>,
        header: &BamHeader,
        build_index: bool,
    ) -> Result<Self, BamWriteError> {
        write_bam_header(&mut bgzf, header)?;

        // r[impl bam_writer.index_coproduction]
        let index = if build_index {
            let voff = bgzf.virtual_offset();
            Some(IndexBuilder::bai(header.target_count(), voff))
        } else {
            None
        };

        Ok(Self { bgzf, index, buf: Vec::with_capacity(512), poisoned: false })
    }

    // r[impl bam_writer.write_record]
    /// Write a single record to the BAM file.
    pub fn write(&mut self, record: &OwnedBamRecord) -> Result<(), BamWriteError> {
        // r[impl bam_writer.error_poisoning]
        if self.poisoned {
            return Err(BamWriteError::Poisoned);
        }

        match self.write_inner(record) {
            Ok(()) => Ok(()),
            Err(e) => {
                self.poisoned = true;
                Err(e)
            }
        }
    }

    fn write_inner(&mut self, record: &OwnedBamRecord) -> Result<(), BamWriteError> {
        // r[impl bam_writer.reuse_buffers]
        self.buf.clear();
        record.to_bam_bytes(&mut self.buf)?;

        // r[impl bam_writer.record_size_limit]
        if self.buf.len() > MAX_RECORD_SIZE {
            return Err(BamWriteError::RecordTooLarge { size: self.buf.len() });
        }

        // r[impl bam_writer.index_record_dispatch]
        // Validate indexability BEFORE writing to BGZF — a failed validation after writing
        // would leave the record in the stream but poison the writer.
        if self.index.is_some() {
            let is_unmapped = record.flags.is_unmapped();
            if record.ref_id == -1 && !is_unmapped {
                return Err(BamWriteError::MappedWithoutReference);
            }
        }

        // Safe: buf.len() <= MAX_RECORD_SIZE (2 MiB) < i32::MAX, checked above.
        debug_assert!(self.buf.len() <= MAX_RECORD_SIZE);
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "buf.len() ≤ MAX_RECORD_SIZE (2 MiB) < i32::MAX; debug_assert enforces invariant"
        )]
        let block_size = self.buf.len() as i32;
        let total = 4usize.saturating_add(self.buf.len());

        // r[impl bam_writer.flush_before_record]
        self.bgzf.flush_if_needed(total)?;

        // r[impl bam_writer.insertion_order]
        // Records are written immediately to the BGZF stream in call order.
        self.bgzf.write_all(&block_size.to_le_bytes())?;
        self.bgzf.write_all(&self.buf)?;

        // Push to index after writing (offset convention: offset AFTER the record)
        if let Some(ref mut index) = self.index {
            let is_unmapped = record.flags.is_unmapped();

            if record.ref_id != -1 {
                // r[impl bam_writer.index_sort_order]
                // Sort validation is delegated to IndexBuilder::push()
                let voff = self.bgzf.virtual_offset();
                let beg = if record.pos < 0 { 0u64 } else { record.pos as u64 };
                let end = if is_unmapped {
                    // Case 2: placed unmapped — beg = end = pos
                    beg
                } else {
                    // Case 1: mapped — use end_pos from CIGAR
                    let ep = record.end_pos();
                    if ep < 0 { beg } else { ep as u64 }
                };
                let end = end.max(beg.saturating_add(1));
                if is_unmapped {
                    // Placed unmapped: increment n_unmapped (not n_mapped) in pseudo-bin
                    index.push_unmapped(record.ref_id, beg, end, voff)?;
                } else {
                    index.push(record.ref_id, beg, end, voff)?;
                }
            }
            // Case 3: fully unmapped (ref_id == -1) — not pushed to index
        }

        Ok(())
    }

    // r[impl bam_writer.write_store_record]
    /// Write a record directly from a [`RecordStore`] without constructing an
    /// intermediate [`OwnedBamRecord`].
    ///
    /// Serializes the record's BAM binary representation into the writer's
    /// reusable buffer, re-encoding the sequence from `Base` to 4-bit packed
    /// format. CIGAR, quality, and aux data are already in BAM wire format
    /// in the store's slabs and are copied directly.
    pub fn write_store_record(
        &mut self,
        store: &super::record_store::RecordStore,
        idx: u32,
    ) -> Result<(), BamWriteError> {
        if self.poisoned {
            return Err(BamWriteError::Poisoned);
        }
        match self.write_store_record_inner(store, idx) {
            Ok(()) => Ok(()),
            Err(e) => {
                self.poisoned = true;
                Err(e)
            }
        }
    }

    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "all casts validated by bounds checks or BAM field limits"
    )]
    fn write_store_record_inner(
        &mut self,
        store: &super::record_store::RecordStore,
        idx: u32,
    ) -> Result<(), BamWriteError> {
        let rec = store.record(idx);
        let qname = store.qname(idx);
        let cigar = store.cigar(idx);
        let seq = store.seq(idx);
        let qual = store.qual(idx);
        let aux = store.aux(idx);

        // Validate field limits (matching OwnedBamRecord::to_bam_bytes).
        // n_cigar_ops is u16 (validated at push time), pos is Pos0 (always ≤ i32::MAX).
        if qname.len() > 254 {
            return Err(OwnedRecordError::QnameTooLong { len: qname.len() }.into());
        }
        if rec.seq_len > i32::MAX as u32 {
            return Err(OwnedRecordError::SeqLengthOverflow { len: rec.seq_len as usize }.into());
        }

        let l_read_name = qname.len() + 1; // +1 for NUL
        let n_cigar_op = rec.n_cigar_ops;

        // Compute BAI bin from pos and end_pos
        let beg = rec.pos.as_u64();
        let end = rec.end_pos.as_u64().max(beg.saturating_add(1));
        let bin = crate::vcf::index_builder::reg2bin(beg, end, 14, 5) as u16;

        // Build the 32-byte fixed header + variable data into self.buf
        self.buf.clear();

        // ref_id (4)
        self.buf.extend_from_slice(&rec.tid.to_le_bytes());
        // pos (4)
        self.buf.extend_from_slice(&rec.pos.as_i32().to_le_bytes());
        // bin_mq_nl (4): bin << 16 | mapq << 8 | l_read_name
        let bin_mq_nl = (u32::from(bin) << 16) | (u32::from(rec.mapq) << 8) | (l_read_name as u32);
        self.buf.extend_from_slice(&bin_mq_nl.to_le_bytes());
        // flag_nc (4): flags << 16 | n_cigar_op
        let flag_nc = (u32::from(rec.flags.raw()) << 16) | u32::from(n_cigar_op);
        self.buf.extend_from_slice(&flag_nc.to_le_bytes());
        // l_seq (4)
        self.buf.extend_from_slice(&(rec.seq_len as i32).to_le_bytes());
        // next_ref_id (4)
        self.buf.extend_from_slice(&rec.next_ref_id.to_le_bytes());
        // next_pos (4)
        self.buf.extend_from_slice(&rec.next_pos.to_le_bytes());
        // template_len (4)
        self.buf.extend_from_slice(&rec.template_len.to_le_bytes());

        // qname + NUL
        self.buf.extend_from_slice(qname);
        self.buf.push(0);

        // CIGAR (already packed u32 LE ops)
        self.buf.extend_from_slice(cigar);

        // Sequence: re-encode Base → 4-bit packed directly into self.buf.
        // Base is #[repr(u8)] with compile-time size_of::<Base>() == 1 assert;
        // transmuting &[Base] to &[u8] is sound.
        let seq_u8: &[u8] =
            unsafe { std::slice::from_raw_parts(seq.as_ptr().cast::<u8>(), seq.len()) };
        super::seq::encode_seq_into(seq_u8, &mut self.buf);

        // Quality scores
        if qual.is_empty() {
            self.buf.resize(self.buf.len() + rec.seq_len as usize, 0xFF);
        } else {
            self.buf.extend_from_slice(qual);
        }

        // Aux tags (raw BAM bytes)
        self.buf.extend_from_slice(aux);

        // Validate size
        if self.buf.len() > MAX_RECORD_SIZE {
            return Err(BamWriteError::RecordTooLarge { size: self.buf.len() });
        }

        // Validate indexability before writing
        if self.index.is_some() {
            let is_unmapped = rec.flags.is_unmapped();
            if rec.tid == -1 && !is_unmapped {
                return Err(BamWriteError::MappedWithoutReference);
            }
        }

        let block_size = self.buf.len() as i32;
        let total = 4usize.saturating_add(self.buf.len());

        self.bgzf.flush_if_needed(total)?;
        self.bgzf.write_all(&block_size.to_le_bytes())?;
        self.bgzf.write_all(&self.buf)?;

        // Index dispatch (same logic as write_inner)
        if let Some(ref mut index) = self.index {
            let is_unmapped = rec.flags.is_unmapped();
            if rec.tid != -1 {
                let voff = self.bgzf.virtual_offset();
                let idx_beg = beg;
                let idx_end = if is_unmapped {
                    // Placed unmapped: beg = end = pos
                    idx_beg
                } else {
                    // Mapped: use end_pos from CIGAR
                    rec.end_pos.as_u64()
                };
                let idx_end = idx_end.max(idx_beg.saturating_add(1));
                if is_unmapped {
                    index.push_unmapped(rec.tid, idx_beg, idx_end, voff)?;
                } else {
                    index.push(rec.tid, idx_beg, idx_end, voff)?;
                }
            }
        }

        Ok(())
    }

    // r[impl bam_writer.finish]
    // r[impl bam_writer.index_finish]
    /// Flush all data, finalize the index, write the BGZF EOF block, and return
    /// the inner writer and optional finished `IndexBuilder`.
    pub fn finish(mut self) -> Result<(W, Option<IndexBuilder>), BamWriteError> {
        if let Some(ref mut index) = self.index {
            let voff = self.bgzf.virtual_offset();
            index.finish(voff)?;
        }
        let inner = self.bgzf.finish()?;
        Ok((inner, self.index))
    }
}

// r[impl bam_writer.magic]
// r[impl bam_writer.header_text]
// r[impl bam_writer.header_references]
/// Write the BAM binary header (magic + text + references) to a BGZF stream.
fn write_bam_header<W: Write>(
    bgzf: &mut BgzfWriter<W>,
    header: &BamHeader,
) -> Result<(), BamWriteError> {
    let mut buf = Vec::new();

    // r[impl bam_writer.magic]
    buf.extend_from_slice(b"BAM\x01");

    // r[impl bam_writer.header_text]
    let text = header.header_text().as_bytes();
    if text.len() > MAX_HEADER_TEXT {
        return Err(BamHeaderError::FieldTooLarge {
            field: "l_text",
            value: text.len(),
            limit: MAX_HEADER_TEXT,
        }
        .into());
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "text.len() ≤ MAX_HEADER_TEXT (256 MiB) < i32::MAX; validated above"
    )]
    buf.extend_from_slice(&(text.len() as i32).to_le_bytes());
    buf.extend_from_slice(text);

    // r[impl bam_writer.header_references]
    if header.target_count() > MAX_REFERENCES {
        return Err(BamHeaderError::FieldTooLarge {
            field: "n_ref",
            value: header.target_count(),
            limit: MAX_REFERENCES,
        }
        .into());
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "target_count() ≤ MAX_REFERENCES (1M) < i32::MAX; validated above"
    )]
    buf.extend_from_slice(&(header.target_count() as i32).to_le_bytes());
    for target in header.targets() {
        let name = target.target_name().as_bytes();
        let name_with_nul = name.len().saturating_add(1);
        if name_with_nul > MAX_REF_NAME {
            return Err(BamHeaderError::FieldTooLarge {
                field: "l_name",
                value: name_with_nul,
                limit: MAX_REF_NAME,
            }
            .into());
        }
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "name_with_nul ≤ MAX_REF_NAME (256 KiB) < i32::MAX; validated above"
        )]
        buf.extend_from_slice(&(name_with_nul as i32).to_le_bytes());
        buf.extend_from_slice(name);
        buf.push(0); // NUL terminator
        // BAM stores l_ref as i32; reject contigs > i32::MAX (≈2.1 Gbp).
        // This is the BAM format limit — no reasonable contig exceeds this.
        let l_ref =
            i32::try_from(target.target_length()).map_err(|_| BamHeaderError::FieldTooLarge {
                field: "l_ref",
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "target_length() is u64; cast to usize for error reporting only (display value)"
                )]
                value: target.target_length() as usize,
                limit: i32::MAX as usize,
            })?;
        buf.extend_from_slice(&l_ref.to_le_bytes());
    }

    bgzf.write_all(&buf)?;
    Ok(())
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::super::aux_data::AuxData;
    use super::super::cigar::{CigarOp, CigarOpType};
    use super::*;
    use seqair_types::BamFlags;
    use seqair_types::Base;
    use seqair_types::BaseQuality;

    fn test_header() -> BamHeader {
        BamHeader::from_sam_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n",
        )
        .unwrap()
    }

    fn test_record(ref_id: i32, pos: i64, qname: &[u8]) -> OwnedBamRecord {
        OwnedBamRecord::builder(ref_id, pos, qname.to_vec())
            .mapq(30)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual([30, 31, 32, 33, 34].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap()
    }

    // r[verify bam_writer.test_empty_file]
    #[test]
    fn empty_bam_has_valid_structure() {
        let header = test_header();
        let mut output = Vec::new();

        {
            let bgzf = BgzfWriter::new(&mut output);
            let writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();
            writer.finish().unwrap();
        }

        // Should start with BGZF magic (gzip)
        assert!(output.len() > 28);
        assert_eq!(output[0], 0x1f);
        assert_eq!(output[1], 0x8b);

        // Decompress and check BAM magic
        let mut reader = super::super::bgzf::BgzfReader::from_reader(io::Cursor::new(&output));
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();
        assert_eq!(&decompressed[..4], b"BAM\x01");
    }

    // r[verify bam_writer.test_roundtrip]
    #[test]
    fn write_and_read_back_records() {
        let header = test_header();
        let mut output = Vec::new();

        {
            let bgzf = BgzfWriter::new(&mut output);
            let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();

            let rec1 = test_record(0, 100, b"read1");
            let rec2 = test_record(0, 200, b"read2");
            writer.write(&rec1).unwrap();
            writer.write(&rec2).unwrap();
            writer.finish().unwrap();
        }

        // Read back using BgzfReader + manual BAM record parsing
        let mut reader = super::super::bgzf::BgzfReader::from_reader(io::Cursor::new(&output));

        // Skip BAM header
        let parsed_header = BamHeader::parse(&mut reader).unwrap();
        assert_eq!(parsed_header.target_count(), 2);

        // Read first record
        let block_size = reader.read_i32().unwrap();
        assert!(block_size > 0);
        let mut rec_data = vec![0u8; block_size as usize];
        reader.read_exact_into(&mut rec_data).unwrap();
        let rec = super::super::record::BamRecord::decode(&rec_data).unwrap();
        assert_eq!(&*rec.qname, b"read1");
        assert_eq!(*rec.pos, 100);
        assert_eq!(rec.mapq, 30);

        // Read second record
        let block_size = reader.read_i32().unwrap();
        let mut rec_data = vec![0u8; block_size as usize];
        reader.read_exact_into(&mut rec_data).unwrap();
        let rec = super::super::record::BamRecord::decode(&rec_data).unwrap();
        assert_eq!(&*rec.qname, b"read2");
        assert_eq!(*rec.pos, 200);
    }

    // r[verify bam_writer.error_poisoning]
    #[test]
    fn poisoned_writer_rejects_subsequent_writes() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();

        // Force poison by writing a record with an oversized qname
        let bad_rec = OwnedBamRecord {
            ref_id: 0,
            pos: 0,
            mapq: 0,
            flags: BamFlags::empty(),
            next_ref_id: -1,
            next_pos: -1,
            template_len: 0,
            qname: vec![b'A'; 255], // too long
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            aux: AuxData::new(),
        };
        assert!(writer.write(&bad_rec).is_err());

        // Subsequent writes should fail with Poisoned
        let good_rec = test_record(0, 100, b"ok");
        let err = writer.write(&good_rec).unwrap_err();
        assert!(matches!(err, BamWriteError::Poisoned));
    }

    // r[verify bam_writer.index_coproduction]
    #[test]
    fn index_coproduction_produces_builder() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, true).unwrap();

        let rec = test_record(0, 100, b"r1");
        writer.write(&rec).unwrap();

        let (_inner, index) = writer.finish().unwrap();
        assert!(index.is_some());

        // Write the BAI to a buffer
        let index = index.unwrap();
        let mut bai_buf = Vec::new();
        index.write_bai(&mut bai_buf, header.target_count()).unwrap();

        // Verify BAI magic
        assert_eq!(&bai_buf[..4], b"BAI\x01");
        // Verify n_ref
        let n_ref = i32::from_le_bytes([bai_buf[4], bai_buf[5], bai_buf[6], bai_buf[7]]);
        assert_eq!(n_ref, 2);
    }

    // r[verify bam_writer.test_index_unsorted]
    #[test]
    fn unsorted_records_error_with_index() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, true).unwrap();

        let rec1 = test_record(0, 200, b"r1");
        let rec2 = test_record(0, 100, b"r2"); // out of order
        writer.write(&rec1).unwrap();
        let err = writer.write(&rec2).unwrap_err();
        assert!(matches!(err, BamWriteError::Index { .. }));
    }

    // r[verify bam_writer.index_record_dispatch]
    #[test]
    fn mapped_without_reference_errors() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, true).unwrap();

        // Mapped (flags & 0x4 == 0) but ref_id == -1
        let rec = OwnedBamRecord::builder(-1, 0, b"bad".to_vec())
            .flags(BamFlags::empty()) // mapped
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            .qual([30, 31, 32].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap();

        let err = writer.write(&rec).unwrap_err();
        assert!(matches!(err, BamWriteError::MappedWithoutReference));
    }

    #[test]
    fn fully_unmapped_not_indexed() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, true).unwrap();

        // Write a mapped record first (index needs at least one)
        let mapped = test_record(0, 100, b"mapped");
        writer.write(&mapped).unwrap();

        // Fully unmapped: ref_id=-1, flags=0x4
        let unmapped = OwnedBamRecord::builder(-1, -1, b"unmapped".to_vec())
            .flags(BamFlags::from(0x4))
            .seq(vec![Base::A])
            .build()
            .unwrap();
        // Should succeed (not pushed to index)
        writer.write(&unmapped).unwrap();

        let (_inner, index) = writer.finish().unwrap();
        assert!(index.is_some());
    }

    #[test]
    fn no_index_when_disabled() {
        let header = test_header();
        let mut output = Vec::new();

        let bgzf = BgzfWriter::new(&mut output);
        let writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();
        let (_inner, index) = writer.finish().unwrap();
        assert!(index.is_none());
    }

    // r[verify bam_writer.test_store_roundtrip]
    #[test]
    fn write_store_record_roundtrips() {
        use super::super::record_store::RecordStore;

        let header = test_header();

        // Build a RecordStore with two records
        let mut store = RecordStore::new();
        let rec1 = OwnedBamRecord::builder(0, 100, b"read1".to_vec())
            .mapq(30)
            .flags(BamFlags::from(0x63))
            .next_ref_id(1)
            .next_pos(500)
            .template_len(250)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual(vec![30, 31, 32, 33, 34])
            .aux({
                let mut a = AuxData::new();
                a.set_string(*b"RG", b"group1");
                a
            })
            .build()
            .unwrap();
        let rec2 = OwnedBamRecord::builder(0, 200, b"read2".to_vec())
            .mapq(60)
            .cigar(vec![
                CigarOp::new(CigarOpType::SoftClip, 2),
                CigarOp::new(CigarOpType::Match, 3),
            ])
            .seq(vec![Base::G, Base::T, Base::A, Base::C, Base::G])
            .qual(vec![35, 36, 37, 38, 39])
            .build()
            .unwrap();

        // Write via OwnedBamRecord to get them into the store
        let mut raw_buf = Vec::new();
        rec1.to_bam_bytes(&mut raw_buf).unwrap();
        store.push_raw(&raw_buf).unwrap();
        raw_buf.clear();
        rec2.to_bam_bytes(&mut raw_buf).unwrap();
        store.push_raw(&raw_buf).unwrap();

        // Write via write_store_record
        let mut output = Vec::new();
        {
            let bgzf = BgzfWriter::new(&mut output);
            let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();
            writer.write_store_record(&store, 0).unwrap();
            writer.write_store_record(&store, 1).unwrap();
            writer.finish().unwrap();
        }

        // Read back into a fresh store and compare
        let mut reader = super::super::bgzf::BgzfReader::from_reader(io::Cursor::new(&output));
        let _ = BamHeader::parse(&mut reader).unwrap();

        let mut store2 = RecordStore::new();
        for _ in 0..2 {
            let block_size = reader.read_i32().unwrap();
            let mut rec_data = vec![0u8; block_size as usize];
            reader.read_exact_into(&mut rec_data).unwrap();
            store2.push_raw(&rec_data).unwrap();
        }

        for i in 0..2u32 {
            let a = store.record(i);
            let b = store2.record(i);
            assert_eq!(*a.pos, *b.pos, "pos mismatch for record {i}");
            assert_eq!(a.flags, b.flags, "flags mismatch for record {i}");
            assert_eq!(a.mapq, b.mapq, "mapq mismatch for record {i}");
            assert_eq!(a.tid, b.tid, "tid mismatch for record {i}");
            assert_eq!(a.next_ref_id, b.next_ref_id, "next_ref_id mismatch for record {i}");
            assert_eq!(a.next_pos, b.next_pos, "next_pos mismatch for record {i}");
            assert_eq!(a.template_len, b.template_len, "template_len mismatch for record {i}");
            assert_eq!(a.seq_len, b.seq_len, "seq_len mismatch for record {i}");
            assert_eq!(a.n_cigar_ops, b.n_cigar_ops, "n_cigar_ops mismatch for record {i}");
            assert_eq!(store.qname(i), store2.qname(i), "qname mismatch for record {i}");
            assert_eq!(store.cigar(i), store2.cigar(i), "cigar mismatch for record {i}");
            assert_eq!(store.seq(i), store2.seq(i), "seq mismatch for record {i}");
            assert_eq!(store.qual(i), store2.qual(i), "qual mismatch for record {i}");
            assert_eq!(store.aux(i), store2.aux(i), "aux mismatch for record {i}");
        }
    }

    // r[verify bam_writer.write_store_record]
    #[test]
    fn write_store_record_with_index() {
        use super::super::record_store::RecordStore;

        let header = test_header();

        // Build a store with 3 sorted records (mapped, placed-unmapped, mapped)
        let mut store = RecordStore::new();
        let recs = [
            OwnedBamRecord::builder(0, 100, b"m1".to_vec())
                .mapq(30)
                .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
                .seq(vec![Base::A, Base::C, Base::G])
                .qual(vec![30, 31, 32])
                .build()
                .unwrap(),
            OwnedBamRecord::builder(0, 200, b"pu1".to_vec())
                .flags(BamFlags::from(0x4)) // unmapped but placed
                .seq(vec![Base::T, Base::A])
                .qual(vec![20, 21])
                .build()
                .unwrap(),
            OwnedBamRecord::builder(0, 300, b"m2".to_vec())
                .mapq(60)
                .cigar(vec![CigarOp::new(CigarOpType::Match, 4)])
                .seq(vec![Base::G, Base::T, Base::A, Base::C])
                .qual(vec![33, 34, 35, 36])
                .build()
                .unwrap(),
        ];
        let mut raw_buf = Vec::new();
        for rec in &recs {
            raw_buf.clear();
            rec.to_bam_bytes(&mut raw_buf).unwrap();
            store.push_raw(&raw_buf).unwrap();
        }

        let mut output = Vec::new();
        let bgzf = BgzfWriter::new(&mut output);
        let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, true).unwrap();
        for i in 0..store.len() as u32 {
            writer.write_store_record(&store, i).unwrap();
        }
        let (_inner, index) = writer.finish().unwrap();
        assert!(index.is_some());

        // Write the BAI and verify magic + n_ref
        let index = index.unwrap();
        let mut bai_buf = Vec::new();
        index.write_bai(&mut bai_buf, header.target_count()).unwrap();
        assert_eq!(&bai_buf[..4], b"BAI\x01");
        let n_ref = i32::from_le_bytes([bai_buf[4], bai_buf[5], bai_buf[6], bai_buf[7]]);
        assert_eq!(n_ref, 2);
    }

    // r[verify bam_writer.insertion_order]
    // r[verify bam_writer.record_spanning_blocks]
    #[test]
    fn records_written_in_insertion_order() {
        let header = test_header();
        let mut output = Vec::new();

        {
            let bgzf = BgzfWriter::new(&mut output);
            let mut writer = BamWriter::<&mut Vec<u8>>::new_inner(bgzf, &header, false).unwrap();

            for i in 0..5 {
                let name = format!("read{i}");
                let rec = test_record(0, i64::from(i) * 100, name.as_bytes());
                writer.write(&rec).unwrap();
            }
            writer.finish().unwrap();
        }

        // Read back and verify order
        let mut reader = super::super::bgzf::BgzfReader::from_reader(io::Cursor::new(&output));
        let _ = BamHeader::parse(&mut reader).unwrap();

        for i in 0..5 {
            let block_size = reader.read_i32().unwrap();
            let mut rec_data = vec![0u8; block_size as usize];
            reader.read_exact_into(&mut rec_data).unwrap();
            let rec = super::super::record::BamRecord::decode(&rec_data).unwrap();
            let expected = format!("read{i}");
            assert_eq!(&*rec.qname, expected.as_bytes());
        }
    }
}
