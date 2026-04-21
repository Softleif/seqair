//! BCF encoding primitives. [`BcfRecordEncoder`] is the internal BCF encoder
//! used by the unified [`Writer`](super::unified::Writer) via enum dispatch.

use super::bcf_encoding::*;
use super::error::VcfError;
use crate::io::{BgzfWriter, IndexBuilder, VirtualOffset};
use std::io::Write;

// ── BcfValue trait ──────────────────────────────────────────────────────

// r[impl bcf_encoder.bcf_value]
/// Trait for types that can be directly encoded into BCF typed values.
pub trait BcfValue: Copy {
    /// The BCF type code used for scalar encoding.
    const TYPE_CODE: u8;
    /// BCF type code for a specific scalar value. Defaults to `TYPE_CODE`.
    /// Overridden by i32 to select the smallest fitting type per `r[bcf_writer.smallest_int_type]`.
    fn scalar_type_code(self) -> u8 {
        Self::TYPE_CODE
    }
    /// Write this value's bytes into the buffer using the given type code.
    fn encode_bcf_as(self, buf: &mut Vec<u8>, type_code: u8);
    /// Write this value's bytes using the default type code.
    #[allow(dead_code, reason = "used by unified.rs via trait dispatch")]
    fn encode_bcf(self, buf: &mut Vec<u8>) {
        self.encode_bcf_as(buf, Self::TYPE_CODE)
    }
    /// Write the missing sentinel for this value type.
    #[allow(dead_code, reason = "used by unified.rs via trait dispatch")]
    fn encode_missing(buf: &mut Vec<u8>);
    // r[impl bcf_writer.end_of_vector]
    /// Write the end-of-vector sentinel for this value type.
    #[allow(dead_code, reason = "used by unified.rs via trait dispatch")]
    fn encode_end_of_vector(buf: &mut Vec<u8>);
}

// r[impl bcf_encoder.bcf_value_float]
impl BcfValue for f32 {
    const TYPE_CODE: u8 = BCF_BT_FLOAT;

    fn encode_bcf_as(self, buf: &mut Vec<u8>, _type_code: u8) {
        buf.extend_from_slice(&self.to_le_bytes());
    }

    fn encode_missing(buf: &mut Vec<u8>) {
        buf.extend_from_slice(&FLOAT_MISSING.to_le_bytes());
    }

    fn encode_end_of_vector(buf: &mut Vec<u8>) {
        buf.extend_from_slice(&FLOAT_END_OF_VECTOR.to_le_bytes());
    }
}

// r[impl bcf_encoder.bcf_value_int]
impl BcfValue for i32 {
    const TYPE_CODE: u8 = BCF_BT_INT32;

    /// Select smallest BCF int type for this value per `r[bcf_writer.smallest_int_type]`.
    fn scalar_type_code(self) -> u8 {
        if (INT8_MIN..=INT8_MAX).contains(&self) {
            BCF_BT_INT8
        } else if (INT16_MIN..=INT16_MAX).contains(&self) {
            BCF_BT_INT16
        } else {
            BCF_BT_INT32
        }
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "type_code is selected by scalar_type_code/smallest_int_type which verified value fits"
    )]
    fn encode_bcf_as(self, buf: &mut Vec<u8>, type_code: u8) {
        match type_code {
            BCF_BT_INT8 => buf.push(self as u8),
            BCF_BT_INT16 => buf.extend_from_slice(&(self as i16).to_le_bytes()),
            _ => buf.extend_from_slice(&self.to_le_bytes()),
        }
    }

    fn encode_bcf(self, buf: &mut Vec<u8>) {
        // Use smallest type for scalars
        self.encode_bcf_as(buf, self.scalar_type_code());
    }

    fn encode_missing(buf: &mut Vec<u8>) {
        buf.extend_from_slice(&INT32_MISSING.to_le_bytes());
    }

    fn encode_end_of_vector(buf: &mut Vec<u8>) {
        buf.extend_from_slice(&INT32_END_OF_VECTOR.to_le_bytes());
    }
}

/// Pre-resolved contig (chromosome) handle. Carries the tid.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ContigHandle(pub(crate) u32);

impl ContigHandle {
    #[allow(dead_code, reason = "used by unified.rs")]
    pub(crate) fn tid(self) -> u32 {
        self.0
    }
}

// ── Field deduplication tracker ─────────────────────────────────────────

// r[impl record_encoder.info_dedup]
// r[impl record_encoder.format_dedup]
/// Tracks written fields by `dict_idx` to detect duplicates within a record.
/// Reused across records: cleared per record, allocation retained.
#[derive(Default)]
pub(crate) struct FieldTracker {
    /// `(dict_idx, byte_offset_in_buffer)` in write order.
    entries: Vec<(u32, usize)>,
}

#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "idx is always from find() which guarantees it is in bounds; \
              delta is always ≤ buf.len() since it came from a drain range"
)]
impl FieldTracker {
    pub(crate) fn clear(&mut self) {
        self.entries.clear();
    }

    /// Returns `Some(index)` if `dict_idx` was already recorded.
    fn find(&self, dict_idx: u32) -> Option<usize> {
        self.entries.iter().position(|(id, _)| *id == dict_idx)
    }

    /// Byte range of the entry at `idx`. Uses `buf_len` as the end for the last entry.
    fn byte_range(&self, idx: usize, buf_len: usize) -> (usize, usize) {
        let start = self.entries[idx].1;
        let end = self.entries.get(idx + 1).map_or(buf_len, |e| e.1);
        (start, end)
    }

    /// Remove entry at `idx` and subtract `delta` from all subsequent offsets.
    fn remove_and_adjust(&mut self, idx: usize, delta: usize) {
        self.entries.remove(idx);
        for entry in &mut self.entries[idx..] {
            entry.1 -= delta;
        }
    }

    pub(crate) fn push(&mut self, dict_idx: u32, offset: usize) {
        self.entries.push((dict_idx, offset));
    }

    /// If `dict_idx` was previously written, remove its bytes from `buf` and
    /// return `true`. The caller should skip incrementing the field count.
    pub(crate) fn remove_duplicate(
        &mut self,
        buf: &mut Vec<u8>,
        dict_idx: u32,
        name: &str,
    ) -> bool {
        let Some(idx) = self.find(dict_idx) else {
            return false;
        };
        tracing::warn!(field = name, "field encoded twice; overwriting previous value");
        let (start, end) = self.byte_range(idx, buf.len());
        buf.drain(start..end);
        self.remove_and_adjust(idx, end - start);
        true
    }

    /// Like [`remove_duplicate`](Self::remove_duplicate) but also handles the
    /// VCF text `;` separator: when the first field is removed, the next field's
    /// leading `;` must be stripped.
    pub(crate) fn remove_duplicate_vcf(
        &mut self,
        buf: &mut Vec<u8>,
        info_count: &mut u16,
        dict_idx: u32,
        name: &str,
    ) -> bool {
        let Some(idx) = self.find(dict_idx) else {
            return false;
        };
        tracing::warn!(field = name, "field encoded twice; overwriting previous value");
        let (start, end) = self.byte_range(idx, buf.len());
        buf.drain(start..end);
        let delta = end - start;
        // First field has no `;` prefix, but the next field does — strip it.
        // The `;` is part of the next entry's tracked range, so removing it
        // does not increase the offset delta for subsequent entries.
        if idx == 0 && self.entries.len() > idx + 1 {
            let sep_pos = start; // after drain, `;` is now at `start`
            buf.drain(sep_pos..=sep_pos);
        }
        self.remove_and_adjust(idx, delta);
        *info_count = info_count.saturating_sub(1);
        true
    }
}

// ── BcfRecordEncoder ────────────────────────────────────────────────────

/// Zero-allocation BCF record encoder. Borrows the writer's buffers.
pub struct BcfRecordEncoder<'a> {
    pub(crate) shared_buf: &'a mut Vec<u8>,
    pub(crate) indiv_buf: &'a mut Vec<u8>,
    #[allow(dead_code, reason = "read via unified.rs do_emit_bcf through BgzfWrite dyn trait")]
    pub(crate) bgzf: &'a mut dyn BgzfWrite,
    #[allow(dead_code, reason = "read via unified.rs do_emit_bcf")]
    pub(crate) index: Option<&'a mut IndexBuilder>,
    // Record state
    pub(crate) n_allele: u16,
    pub(crate) n_alt: u16,
    pub(crate) n_info: u16,
    pub(crate) n_fmt: u8,
    pub(crate) n_sample: u32,
    pub(crate) tid: i32,
    pub(crate) pos_0based: i32,
    pub(crate) rlen: i32,
    // Dedup trackers (borrowed from WriterInner for cross-record reuse)
    pub(crate) info_tracker: &'a mut FieldTracker,
    pub(crate) fmt_tracker: &'a mut FieldTracker,
}

// r[impl bcf_writer.bgzf_blocks]
/// Trait to abstract over `BgzfWriter`<W> for different W types.
pub(crate) trait BgzfWrite {
    #[allow(dead_code, reason = "called through dyn BgzfWrite in unified.rs")]
    fn virtual_offset(&self) -> VirtualOffset;
    #[allow(dead_code, reason = "called through dyn BgzfWrite in unified.rs")]
    fn flush_if_needed(&mut self, upcoming: usize) -> Result<(), crate::io::BgzfError>;
    #[allow(dead_code, reason = "called through dyn BgzfWrite in unified.rs")]
    fn write_all(&mut self, data: &[u8]) -> Result<(), crate::io::BgzfError>;
}

impl<W: Write> BgzfWrite for BgzfWriter<W> {
    fn virtual_offset(&self) -> VirtualOffset {
        self.virtual_offset()
    }
    fn flush_if_needed(&mut self, upcoming: usize) -> Result<(), crate::io::BgzfError> {
        self.flush_if_needed(upcoming)
    }
    fn write_all(&mut self, data: &[u8]) -> Result<(), crate::io::BgzfError> {
        self.write_all(data)
    }
}

impl<'a> BcfRecordEncoder<'a> {
    // r[impl bcf_encoder.emit]
    /// Patch the header, write the record to BGZF, push to index.
    #[allow(dead_code, reason = "called by unified.rs do_emit_bcf")]
    pub(crate) fn emit_inner(&mut self) -> Result<(), VcfError> {
        // Patch n_info|n_allele and n_fmt|n_sample in the 24-byte fixed header
        // These are at offsets 16 and 20 in shared_buf
        let n_info_allele = (u32::from(self.n_allele) << 16) | u32::from(self.n_info);
        let n_fmt_sample = (u32::from(self.n_fmt) << 24) | self.n_sample;

        // Patch bytes 16..20 and 20..24 of shared_buf
        if let Some(dest) = self.shared_buf.get_mut(16..20) {
            dest.copy_from_slice(&n_info_allele.to_le_bytes());
        }
        if let Some(dest) = self.shared_buf.get_mut(20..24) {
            dest.copy_from_slice(&n_fmt_sample.to_le_bytes());
        }

        // r[impl bcf_writer.record_layout]
        let l_shared = u32::try_from(self.shared_buf.len()).map_err(|_| {
            VcfError::RecordTooLarge { section: "shared", size: self.shared_buf.len() }
        })?;
        let l_indiv = u32::try_from(self.indiv_buf.len()).map_err(|_| {
            VcfError::RecordTooLarge { section: "individual", size: self.indiv_buf.len() }
        })?;
        let total =
            8usize.saturating_add(self.shared_buf.len()).saturating_add(self.indiv_buf.len());

        self.bgzf.flush_if_needed(total)?;
        self.bgzf.write_all(&l_shared.to_le_bytes())?;
        self.bgzf.write_all(&l_indiv.to_le_bytes())?;
        self.bgzf.write_all(self.shared_buf)?;
        self.bgzf.write_all(self.indiv_buf)?;

        // Index co-production
        if let Some(ref mut index) = self.index {
            let beg = self.pos_0based as u64;
            let end = beg.saturating_add(self.rlen as u64);
            index.push(self.tid, beg, end, self.bgzf.virtual_offset())?;
        }

        Ok(())
    }
}

// BcfRecordEncoder is used internally by the unified Writer/RecordEncoder
// in unified.rs. No trait impls needed — the unified encoder delegates directly.

// ── Alleles integration ─────────────────────────────────────────────────

use super::alleles::Alleles;
use seqair_types::Pos1;

impl Alleles {
    // r[impl bcf_encoder.begin_record]
    // r[impl bcf_encoder.checked_casts]
    // r[impl bcf_writer.shared_variable]
    /// Begin a BCF record: write the 24-byte fixed header, ID, and alleles.
    /// Sets `n_allele/n_alt` on the encoder for downstream field validation.
    pub(crate) fn begin_record(
        &self,
        enc: &mut BcfRecordEncoder<'_>,
        contig: ContigHandle,
        pos: Pos1,
        qual: Option<f32>,
    ) -> Result<(), VcfError> {
        enc.shared_buf.clear();
        enc.indiv_buf.clear();
        enc.n_info = 0;
        enc.n_fmt = 0;
        // n_sample is set from the header via the struct literal in Writer::begin_record each call.

        enc.tid = i32::try_from(contig.0).map_err(|_| VcfError::ValueOverflow {
            field: "contig_tid",
            value: u64::from(contig.0),
            target_type: "i32",
        })?;
        // r[impl bcf_writer.coordinate_system]
        enc.pos_0based = pos.to_zero_based().as_i32();
        enc.rlen = i32::try_from(self.rlen()).map_err(|_| VcfError::ValueOverflow {
            field: "rlen",
            value: self.rlen() as u64,
            target_type: "i32",
        })?;
        enc.n_allele = u16::try_from(self.n_allele()).map_err(|_| VcfError::ValueOverflow {
            field: "n_allele",
            value: self.n_allele() as u64,
            target_type: "u16",
        })?;
        enc.n_alt = u16::try_from(self.n_allele().saturating_sub(1)).map_err(|_| {
            VcfError::ValueOverflow {
                field: "n_alt",
                value: self.n_allele().saturating_sub(1) as u64,
                target_type: "u16",
            }
        })?;

        // r[impl bcf_writer.qual_missing]
        let qual_bits = match qual {
            Some(q) => q.to_bits(),
            None => FLOAT_MISSING,
        };

        // r[impl bcf_writer.fixed_fields]
        // 24-byte fixed header (n_info/n_allele and n_fmt/n_sample patched at emit)
        enc.shared_buf.extend_from_slice(&enc.tid.to_le_bytes());
        enc.shared_buf.extend_from_slice(&enc.pos_0based.to_le_bytes());
        enc.shared_buf.extend_from_slice(&enc.rlen.to_le_bytes());
        enc.shared_buf.extend_from_slice(&qual_bits.to_le_bytes());
        // Placeholder for n_info|n_allele (patched at emit)
        enc.shared_buf.extend_from_slice(&0u32.to_le_bytes());
        // Placeholder for n_fmt|n_sample (patched at emit)
        enc.shared_buf.extend_from_slice(&0u32.to_le_bytes());

        // ID = "."
        encode_type_byte(enc.shared_buf, 1, BCF_BT_CHAR);
        enc.shared_buf.push(b'.');

        // REF allele (zero-alloc): compute length, write type header, then data.
        let ref_len = self.ref_byte_len();
        encode_type_byte(enc.shared_buf, ref_len, BCF_BT_CHAR);
        self.write_ref_into(enc.shared_buf);

        // ALT alleles (zero-alloc, one typed string each)
        match self {
            Alleles::Reference { .. } => {} // no ALT
            Alleles::Snv { alt_bases, .. } => {
                for b in alt_bases {
                    encode_type_byte(enc.shared_buf, 1, BCF_BT_CHAR);
                    enc.shared_buf.push(b.as_char() as u8);
                }
            }
            Alleles::Insertion { anchor, inserted } => {
                let alt_len = 1usize.saturating_add(inserted.len());
                encode_type_byte(enc.shared_buf, alt_len, BCF_BT_CHAR);
                enc.shared_buf.push(anchor.as_char() as u8);
                for b in inserted {
                    enc.shared_buf.push(b.as_char() as u8);
                }
            }
            Alleles::Deletion { anchor, .. } => {
                encode_type_byte(enc.shared_buf, 1, BCF_BT_CHAR);
                enc.shared_buf.push(anchor.as_char() as u8);
            }
            Alleles::Complex { alt_alleles, .. } => {
                for alt in alt_alleles {
                    encode_typed_string(enc.shared_buf, alt.as_bytes());
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use seqair_types::Base;

    /// Create encoder with in-memory buffers for testing (no BGZF).
    fn test_encoder<'a>(
        shared: &'a mut Vec<u8>,
        indiv: &'a mut Vec<u8>,
        bgzf_buf: &'a mut TestBgzf,
        info_tracker: &'a mut FieldTracker,
        fmt_tracker: &'a mut FieldTracker,
    ) -> BcfRecordEncoder<'a> {
        BcfRecordEncoder {
            shared_buf: shared,
            indiv_buf: indiv,
            bgzf: bgzf_buf,
            index: None,
            n_allele: 0,
            n_alt: 0,
            n_info: 0,
            n_fmt: 0,
            n_sample: 0,
            tid: 0,
            pos_0based: 0,
            rlen: 0,
            info_tracker,
            fmt_tracker,
        }
    }

    /// Mock BGZF for unit tests — just collects bytes.
    struct TestBgzf {
        data: Vec<u8>,
    }

    impl TestBgzf {
        fn new() -> Self {
            Self { data: Vec::new() }
        }
    }

    impl BgzfWrite for TestBgzf {
        fn virtual_offset(&self) -> VirtualOffset {
            VirtualOffset(self.data.len() as u64)
        }
        fn flush_if_needed(&mut self, _upcoming: usize) -> Result<(), crate::io::BgzfError> {
            Ok(())
        }
        fn write_all(&mut self, data: &[u8]) -> Result<(), crate::io::BgzfError> {
            self.data.extend_from_slice(data);
            Ok(())
        }
    }

    // r[verify bcf_encoder.begin_record]
    #[test]
    fn begin_record_writes_fixed_header_and_alleles() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut it = FieldTracker::default();
        let mut ft = FieldTracker::default();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf, &mut it, &mut ft);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let pos = Pos1::new(100).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), pos, Some(30.0)).unwrap();

        assert_eq!(enc.n_allele, 2);
        assert_eq!(enc.n_alt, 1);
        assert_eq!(enc.tid, 0);
        assert_eq!(enc.pos_0based, 99); // 0-based
        assert_eq!(enc.rlen, 1);
        assert!(shared.len() >= 24, "must have at least 24-byte header");
    }

    // ── FieldTracker unit tests ─────────────────────────────────────

    #[test]
    fn field_tracker_remove_duplicate_first_of_three() {
        // Simulate 3 BCF INFO fields: [AAAA][BBBB][CCCC]
        let mut buf = b"AAAABBBBCCCC".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0); // field 1 at offset 0
        t.push(2, 4); // field 2 at offset 4
        t.push(3, 8); // field 3 at offset 8

        assert!(t.remove_duplicate(&mut buf, 1, "F1"));
        assert_eq!(buf, b"BBBBCCCC");
        // Tracker should have 2 entries with adjusted offsets
        assert!(t.find(1).is_none());
        assert_eq!(t.find(2), Some(0));
        assert_eq!(t.find(3), Some(1));
    }

    #[test]
    fn field_tracker_remove_duplicate_middle_of_three() {
        let mut buf = b"AAAABBBBCCCC".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0);
        t.push(2, 4);
        t.push(3, 8);

        assert!(t.remove_duplicate(&mut buf, 2, "F2"));
        assert_eq!(buf, b"AAAACCCC");
        assert_eq!(t.find(1), Some(0));
        assert!(t.find(2).is_none());
        assert_eq!(t.find(3), Some(1));
    }

    #[test]
    fn field_tracker_remove_duplicate_last_of_three() {
        let mut buf = b"AAAABBBBCCCC".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0);
        t.push(2, 4);
        t.push(3, 8);

        assert!(t.remove_duplicate(&mut buf, 3, "F3"));
        assert_eq!(buf, b"AAAABBBB");
        assert_eq!(t.find(1), Some(0));
        assert_eq!(t.find(2), Some(1));
        assert!(t.find(3).is_none());
    }

    #[test]
    fn field_tracker_remove_duplicate_only_entry() {
        let mut buf = b"AAAA".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0);

        assert!(t.remove_duplicate(&mut buf, 1, "F1"));
        assert!(buf.is_empty());
        assert!(t.find(1).is_none());
    }

    #[test]
    fn field_tracker_no_duplicate_returns_false() {
        let mut buf = b"AAAA".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0);

        assert!(!t.remove_duplicate(&mut buf, 99, "missing"));
        assert_eq!(buf, b"AAAA"); // unchanged
    }

    #[test]
    fn field_tracker_vcf_remove_first_strips_separator() {
        // VCF: "DP=50;BQ=30" — DP at 0, ;BQ at 5
        let mut buf = b"DP=50;BQ=30".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0); // DP starts at 0
        t.push(2, 5); // ;BQ starts at 5
        let mut count = 2u16;

        assert!(t.remove_duplicate_vcf(&mut buf, &mut count, 1, "DP"));
        assert_eq!(buf, b"BQ=30");
        assert_eq!(count, 1);
        assert!(t.find(1).is_none());
        assert_eq!(t.find(2), Some(0));
    }

    #[test]
    fn field_tracker_vcf_remove_last() {
        let mut buf = b"DP=50;BQ=30".to_vec();
        let mut t = FieldTracker::default();
        t.push(1, 0);
        t.push(2, 5);
        let mut count = 2u16;

        assert!(t.remove_duplicate_vcf(&mut buf, &mut count, 2, "BQ"));
        assert_eq!(buf, b"DP=50");
        assert_eq!(count, 1);
    }

    #[test]
    fn field_tracker_clear_resets() {
        let mut t = FieldTracker::default();
        t.push(1, 0);
        t.push(2, 4);
        t.clear();
        assert!(t.find(1).is_none());
        assert!(t.find(2).is_none());
    }
}
