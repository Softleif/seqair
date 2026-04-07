//! Zero-allocation BCF direct encoder. Pre-resolved typed handles write values
//! straight into BCF buffers without constructing a `VcfRecord` intermediate.
//!
//! Usage:
//! ```ignore
//! let dp = writer.resolve_info_scalar::<i32>("DP")?;
//! let pass = writer.resolve_filter_pass();
//! let contig = writer.resolve_contig("chr1")?;
//! let mut enc = writer.record_encoder();
//! alleles.begin_record(&mut enc, contig, pos, qual)?;
//! pass.encode(&mut enc);
//! dp.encode(&mut enc, 50);
//! enc.emit()?;
//! ```

use super::bcf_encoding::*;
use super::error::{VcfError, VcfHeaderError};
use super::header::VcfHeader;
use super::index_builder::IndexBuilder;
use super::record::{Filters, Genotype, InfoValue, SampleValue, VcfRecord};
use crate::bam::bgzf::VirtualOffset;
use crate::bam::bgzf_writer::BgzfWriter;
use std::io::Write;
use std::marker::PhantomData;

// ── BcfValue trait ──────────────────────────────────────────────────────

// r[impl bcf_encoder.bcf_value]
/// Trait for types that can be directly encoded into BCF typed values.
pub trait BcfValue: Copy {
    /// The BCF type code used for scalar encoding.
    const TYPE_CODE: u8;
    /// BCF type code for a specific scalar value. Defaults to TYPE_CODE.
    /// Overridden by i32 to select the smallest fitting type per r[bcf_writer.smallest_int_type].
    fn scalar_type_code(self) -> u8 {
        Self::TYPE_CODE
    }
    /// Write this value's bytes into the buffer using the given type code.
    fn encode_bcf_as(self, buf: &mut Vec<u8>, type_code: u8);
    /// Write this value's bytes using the default type code.
    fn encode_bcf(self, buf: &mut Vec<u8>) {
        self.encode_bcf_as(buf, Self::TYPE_CODE)
    }
    /// Write the missing sentinel for this value type.
    fn encode_missing(buf: &mut Vec<u8>);
    /// Write the end-of-vector sentinel for this value type.
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

    /// Select smallest BCF int type for this value per r[bcf_writer.smallest_int_type].
    fn scalar_type_code(self) -> u8 {
        if (INT8_MIN..=INT8_MAX).contains(&self) {
            BCF_BT_INT8
        } else if (INT16_MIN..=INT16_MAX).contains(&self) {
            BCF_BT_INT16
        } else {
            BCF_BT_INT32
        }
    }

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

// ── Handle types ────────────────────────────────────────────────────────

// r[impl bcf_encoder.handles]

/// Pre-resolved contig (chromosome) handle. Carries the tid.
#[derive(Debug, Clone, Copy)]
pub struct ContigHandle(pub u32);

impl ContigHandle {
    pub fn tid(self) -> u32 {
        self.0
    }
}

/// Pre-resolved filter handle. PASS is always index 0.
#[derive(Debug, Clone, Copy)]
pub struct FilterHandle(pub u32);

impl FilterHandle {
    /// The PASS filter (always BCF dictionary index 0).
    pub const PASS: Self = Self(0);

    // r[impl bcf_encoder.handle_encode]
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>) {
        debug_assert!(self.0 <= i32::MAX as u32, "filter dict index overflow");
        encode_typed_int_vec(enc.shared_buf, &[self.0 as i32]);
    }
}

// r[impl bcf_encoder.handle_types]
// r[impl bcf_encoder.handle_value_type]

/// INFO field handle for `Number::Count(1)` — encodes exactly 1 value.
#[derive(Debug, Clone, Copy)]
pub struct ScalarInfoHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// INFO field handle for `Number::Count(0)` (Flag) — encodes no value.
#[derive(Debug, Clone, Copy)]
pub struct FlagInfoHandle {
    pub dict_idx: u32,
}

/// INFO field handle for `Number::A` — encodes `n_alt` values.
#[derive(Debug, Clone, Copy)]
pub struct PerAltInfoHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// INFO field handle for `Number::R` — encodes `n_allele` values.
#[derive(Debug, Clone, Copy)]
pub struct PerAlleleInfoHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// FORMAT field handle for `Number::Count(1)` — encodes 1 value per sample.
#[derive(Debug, Clone, Copy)]
pub struct ScalarFormatHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// FORMAT field handle for `Number::A` — encodes `n_alt` values per sample.
#[derive(Debug, Clone, Copy)]
pub struct PerAltFormatHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// FORMAT field handle for `Number::R` — encodes `n_allele` values per sample.
#[derive(Debug, Clone, Copy)]
pub struct PerAlleleFormatHandle<T: BcfValue> {
    pub dict_idx: u32,
    pub _marker: PhantomData<T>,
}

/// FORMAT field handle for GT (genotype) — special encoding.
#[derive(Debug, Clone, Copy)]
pub struct GtFormatHandle {
    pub dict_idx: u32,
}

// ── Handle encode implementations ──────────────────────────────────────

// r[impl bcf_encoder.handle_encode]

impl<T: BcfValue> ScalarInfoHandle<T> {
    // r[impl bcf_encoder.bcf_value_int] — uses scalar_type_code for smallest-int selection
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, value: T) {
        let type_code = value.scalar_type_code();
        encode_typed_int_key(enc.shared_buf, self.dict_idx);
        encode_type_byte(enc.shared_buf, 1, type_code);
        value.encode_bcf_as(enc.shared_buf, type_code);
        enc.n_info = enc.n_info.saturating_add(1);
    }
}

impl FlagInfoHandle {
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>) {
        encode_typed_int_key(enc.shared_buf, self.dict_idx);
        // Flag: type=0, count=0
        encode_type_byte(enc.shared_buf, 0, BCF_BT_NULL);
        enc.n_info = enc.n_info.saturating_add(1);
    }
}

impl<T: BcfValue> PerAltInfoHandle<T> {
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, values: &[T]) {
        debug_assert_eq!(
            values.len(),
            enc.n_alt as usize,
            "PerAlt field expects {} values (n_alt), got {}",
            enc.n_alt,
            values.len()
        );
        encode_typed_int_key(enc.shared_buf, self.dict_idx);
        encode_array_values(enc.shared_buf, values);
        enc.n_info = enc.n_info.saturating_add(1);
    }
}

impl<T: BcfValue> PerAlleleInfoHandle<T> {
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, values: &[T]) {
        debug_assert_eq!(
            values.len(),
            enc.n_allele as usize,
            "PerAllele field expects {} values (n_allele), got {}",
            enc.n_allele,
            values.len()
        );
        encode_typed_int_key(enc.shared_buf, self.dict_idx);
        encode_array_values(enc.shared_buf, values);
        enc.n_info = enc.n_info.saturating_add(1);
    }
}

impl<T: BcfValue> ScalarFormatHandle<T> {
    /// Encode a single-sample scalar FORMAT field.
    // r[impl bcf_encoder.bcf_value_int] — uses scalar_type_code for smallest-int selection
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, value: T) {
        let type_code = value.scalar_type_code();
        encode_typed_int_key(enc.indiv_buf, self.dict_idx);
        encode_type_byte(enc.indiv_buf, 1, type_code);
        value.encode_bcf_as(enc.indiv_buf, type_code);
        enc.n_fmt = enc.n_fmt.saturating_add(1);
    }
}

impl<T: BcfValue> PerAltFormatHandle<T> {
    /// Encode a single-sample per-alt FORMAT field.
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, values: &[T]) {
        debug_assert_eq!(values.len(), enc.n_alt as usize);
        encode_typed_int_key(enc.indiv_buf, self.dict_idx);
        encode_array_values(enc.indiv_buf, values);
        enc.n_fmt = enc.n_fmt.saturating_add(1);
    }
}

impl<T: BcfValue> PerAlleleFormatHandle<T> {
    /// Encode a single-sample per-allele FORMAT field.
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, values: &[T]) {
        debug_assert_eq!(values.len(), enc.n_allele as usize);
        encode_typed_int_key(enc.indiv_buf, self.dict_idx);
        encode_array_values(enc.indiv_buf, values);
        enc.n_fmt = enc.n_fmt.saturating_add(1);
    }
}

// r[impl bcf_encoder.handle_encode] (GT)
impl GtFormatHandle {
    /// Encode a single-sample GT FORMAT field.
    pub fn encode(&self, enc: &mut BcfRecordEncoder<'_>, gt: &Genotype) {
        encode_typed_int_key(enc.indiv_buf, self.dict_idx);

        let ploidy = gt.alleles.len();

        // Compute max encoded value to select smallest int type
        let max_allele: i32 =
            gt.alleles.iter().flatten().map(|idx| i32::from(*idx)).max().unwrap_or(0);
        let max_val = (max_allele.saturating_add(1)).saturating_mul(2).saturating_add(1);
        let typ = smallest_int_type(&[max_val]);
        encode_type_byte(enc.indiv_buf, ploidy, typ);

        for (i, allele_opt) in gt.alleles.iter().enumerate() {
            let encoded: i32 = match allele_opt {
                Some(idx) => {
                    let phased = if i == 0 {
                        false
                    } else {
                        gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                    };
                    (i32::from(*idx).saturating_add(1) << 1) | (phased as i32)
                }
                None => 0, // missing allele
            };
            encode_int_as(enc.indiv_buf, encoded, typ);
        }
        enc.n_fmt = enc.n_fmt.saturating_add(1);
    }
}

// ── BcfRecordEncoder ────────────────────────────────────────────────────

// r[impl bcf_encoder.encoder]
/// Zero-allocation BCF record encoder. Borrows the writer's buffers.
pub struct BcfRecordEncoder<'a> {
    pub(crate) shared_buf: &'a mut Vec<u8>,
    pub(crate) indiv_buf: &'a mut Vec<u8>,
    pub(crate) bgzf: &'a mut dyn BgzfWrite,
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
}

/// Trait to abstract over BgzfWriter<W> for different W types.
pub(crate) trait BgzfWrite {
    fn virtual_offset(&self) -> VirtualOffset;
    fn flush_if_needed(&mut self, upcoming: usize) -> Result<(), crate::bam::bgzf::BgzfError>;
    fn write_all(&mut self, data: &[u8]) -> Result<(), crate::bam::bgzf::BgzfError>;
}

impl<W: Write> BgzfWrite for BgzfWriter<W> {
    fn virtual_offset(&self) -> VirtualOffset {
        self.virtual_offset()
    }
    fn flush_if_needed(&mut self, upcoming: usize) -> Result<(), crate::bam::bgzf::BgzfError> {
        self.flush_if_needed(upcoming)
    }
    fn write_all(&mut self, data: &[u8]) -> Result<(), crate::bam::bgzf::BgzfError> {
        self.write_all(data)
    }
}

impl<'a> BcfRecordEncoder<'a> {
    /// Number of alleles (REF + ALTs).
    pub fn n_allele(&self) -> usize {
        self.n_allele as usize
    }

    /// Number of alt alleles.
    pub fn n_alt(&self) -> usize {
        self.n_alt as usize
    }

    /// Set number of samples for FORMAT encoding.
    pub fn begin_samples(&mut self, n_sample: u32) {
        self.n_sample = n_sample;
    }

    // r[impl bcf_encoder.emit]
    /// Patch the header, write the record to BGZF, push to index.
    pub fn emit(&mut self) -> Result<(), VcfError> {
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

// ── Alleles integration ─────────────────────────────────────────────────

use super::alleles::Alleles;
use seqair_types::{One, Pos};

impl Alleles {
    // r[impl bcf_encoder.begin_record]
    // r[impl bcf_encoder.checked_casts]
    /// Begin a BCF record: write the 24-byte fixed header, ID, and alleles.
    /// Sets n_allele/n_alt on the encoder for downstream field validation.
    pub fn begin_record(
        &self,
        enc: &mut BcfRecordEncoder<'_>,
        contig: ContigHandle,
        pos: Pos<One>,
        qual: Option<f32>,
    ) -> Result<(), VcfError> {
        enc.shared_buf.clear();
        enc.indiv_buf.clear();
        enc.n_info = 0;
        enc.n_fmt = 0;
        enc.n_sample = 0;

        enc.tid = i32::try_from(contig.0).map_err(|_| VcfError::ValueOverflow {
            field: "contig_tid",
            value: u64::from(contig.0),
            target_type: "i32",
        })?;
        enc.pos_0based =
            i32::try_from(pos.to_zero_based().get()).map_err(|_| VcfError::ValueOverflow {
                field: "pos",
                value: pos.to_zero_based().get() as u64,
                target_type: "i32",
            })?;
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

        let qual_bits = match qual {
            Some(q) => q.to_bits(),
            None => FLOAT_MISSING,
        };

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

// ── Encoding helpers ────────────────────────────────────────────────────

impl<'a> BcfRecordEncoder<'a> {
    /// Encode a complete `VcfRecord` into this encoder and emit it.
    pub fn write_vcf_record(
        &mut self,
        record: &VcfRecord,
        header: &VcfHeader,
    ) -> Result<(), VcfError> {
        let tid = header.contig_id(&record.contig)?;
        let tid_u32 = u32::try_from(tid).map_err(|_| VcfError::ValueOverflow {
            field: "contig_tid",
            value: tid as u64,
            target_type: "u32",
        })?;
        record.alleles.begin_record(self, ContigHandle(tid_u32), record.pos, record.qual)?;

        // Filter
        match &record.filters {
            Filters::Pass => {
                encode_type_byte(self.shared_buf, 1, BCF_BT_INT8);
                self.shared_buf.push(0);
            }
            Filters::Failed(ids) => {
                let mut indices = Vec::with_capacity(ids.len());
                for id in ids {
                    let idx = header.string_map().get(id).ok_or_else(|| {
                        VcfError::Header(VcfHeaderError::MissingFilter { id: id.clone() })
                    })?;
                    let idx_i32 = i32::try_from(idx).map_err(|_| VcfError::ValueOverflow {
                        field: "filter_dict_idx",
                        value: idx as u64,
                        target_type: "i32",
                    })?;
                    indices.push(idx_i32);
                }
                encode_typed_int_vec(self.shared_buf, &indices);
            }
            Filters::NotApplied => {
                encode_type_byte(self.shared_buf, 0, BCF_BT_NULL);
            }
        }

        // INFO fields
        for (key, value) in record.info.iter() {
            let idx = header
                .string_map()
                .get(key)
                .ok_or_else(|| VcfError::Header(VcfHeaderError::MissingInfo { id: key.clone() }))?;
            let dict_idx = i32::try_from(idx).map_err(|_| VcfError::ValueOverflow {
                field: "info_dict_idx",
                value: idx as u64,
                target_type: "i32",
            })?;
            encode_typed_int_vec(self.shared_buf, &[dict_idx]);
            encode_info_value(self.shared_buf, value);
            self.n_info = self.n_info.saturating_add(1);
        }

        // FORMAT / individual fields
        let n_sample = record.samples.values.len();
        if n_sample > 0 && !record.samples.format_keys.is_empty() {
            self.n_sample = u32::try_from(n_sample).map_err(|_| VcfError::ValueOverflow {
                field: "n_sample",
                value: n_sample as u64,
                target_type: "u32",
            })?;
            for (field_idx, key) in record.samples.format_keys.iter().enumerate() {
                let fmt_idx = header.string_map().get(key).ok_or_else(|| {
                    VcfError::Header(VcfHeaderError::MissingFormat { id: key.clone() })
                })?;
                let dict_idx = i32::try_from(fmt_idx).map_err(|_| VcfError::ValueOverflow {
                    field: "format_dict_idx",
                    value: fmt_idx as u64,
                    target_type: "i32",
                })?;
                encode_typed_int_vec(self.indiv_buf, &[dict_idx]);
                if key == "GT" {
                    encode_gt_field(self.indiv_buf, &record.samples.values, field_idx);
                } else {
                    encode_format_field(self.indiv_buf, &record.samples.values, field_idx);
                }
                self.n_fmt = self.n_fmt.saturating_add(1);
            }
        }

        self.emit()
    }
}

// ── Record-path encoding helpers ─────────────────────────────────────────

// r[impl bcf_writer.flag_encoding]
pub(crate) fn encode_info_value(buf: &mut Vec<u8>, value: &InfoValue) {
    match value {
        InfoValue::Integer(v) => encode_typed_int_vec(buf, &[*v]),
        InfoValue::Float(v) => {
            encode_type_byte(buf, 1, BCF_BT_FLOAT);
            buf.extend_from_slice(&v.to_le_bytes());
        }
        InfoValue::Flag => {
            encode_type_byte(buf, 0, BCF_BT_NULL);
        }
        InfoValue::String(s) => encode_typed_string(buf, s.as_bytes()),
        InfoValue::IntegerArray(arr) => {
            // r[impl bcf_writer.smallest_int_type]
            let typ = smallest_int_type_iter(arr.iter().filter_map(|v| *v));
            encode_type_byte(buf, arr.len(), typ);
            for v in arr {
                match v {
                    Some(n) => match typ {
                        BCF_BT_INT8 => buf.push(*n as u8),
                        BCF_BT_INT16 => buf.extend_from_slice(&(*n as i16).to_le_bytes()),
                        _ => buf.extend_from_slice(&n.to_le_bytes()),
                    },
                    // r[impl bcf_writer.missing_sentinels]
                    None => match typ {
                        BCF_BT_INT8 => buf.push(INT8_MISSING),
                        BCF_BT_INT16 => buf.extend_from_slice(&INT16_MISSING.to_le_bytes()),
                        _ => buf.extend_from_slice(&INT32_MISSING.to_le_bytes()),
                    },
                }
            }
        }
        InfoValue::FloatArray(arr) => {
            encode_type_byte(buf, arr.len(), BCF_BT_FLOAT);
            for v in arr {
                match v {
                    // r[impl bcf_writer.missing_sentinels]
                    Some(f) => buf.extend_from_slice(&f.to_le_bytes()),
                    None => buf.extend_from_slice(&FLOAT_MISSING.to_le_bytes()),
                }
            }
        }
        InfoValue::StringArray(arr) => {
            let joined: String =
                arr.iter().map(|v| v.as_deref().unwrap_or(".")).collect::<Vec<_>>().join(",");
            encode_typed_string(buf, joined.as_bytes());
        }
    }
}

// r[impl bcf_writer.gt_encoding]
pub(crate) fn encode_gt_field(
    buf: &mut Vec<u8>,
    samples: &[seqair_types::SmallVec<SampleValue, 6>],
    field_idx: usize,
) {
    let mut max_ploidy = 0usize;
    for sample in samples {
        if let Some(SampleValue::Genotype(gt)) = sample.get(field_idx) {
            max_ploidy = max_ploidy.max(gt.alleles.len());
        }
    }
    if max_ploidy == 0 {
        max_ploidy = 2;
    }

    let mut max_allele: i32 = 0;
    for sample in samples {
        if let Some(SampleValue::Genotype(gt)) = sample.get(field_idx) {
            for idx in gt.alleles.iter().flatten() {
                max_allele = max_allele.max(i32::from(*idx));
            }
        }
    }

    let max_val = (max_allele.saturating_add(1)).saturating_mul(2).saturating_add(1);
    let typ = smallest_int_type(&[max_val]);

    encode_type_byte(buf, max_ploidy, typ);

    for sample in samples {
        let gt = match sample.get(field_idx) {
            Some(SampleValue::Genotype(g)) => Some(g),
            _ => None,
        };

        for i in 0..max_ploidy {
            let encoded = if let Some(g) = gt {
                if let Some(allele_opt) = g.alleles.get(i) {
                    match allele_opt {
                        Some(idx) => {
                            let phased = if i == 0 {
                                false
                            } else {
                                g.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                            };
                            ((i32::from(*idx).saturating_add(1)) << 1) | (phased as i32)
                        }
                        None => 0,
                    }
                } else {
                    match typ {
                        BCF_BT_INT8 => {
                            buf.push(INT8_END_OF_VECTOR);
                            continue;
                        }
                        BCF_BT_INT16 => {
                            buf.extend_from_slice(&INT16_END_OF_VECTOR.to_le_bytes());
                            continue;
                        }
                        _ => {
                            buf.extend_from_slice(&INT32_END_OF_VECTOR.to_le_bytes());
                            continue;
                        }
                    }
                }
            } else {
                0
            };

            match typ {
                BCF_BT_INT8 => buf.push(encoded as u8),
                BCF_BT_INT16 => buf.extend_from_slice(&(encoded as i16).to_le_bytes()),
                _ => buf.extend_from_slice(&encoded.to_le_bytes()),
            }
        }
    }
}

pub(crate) fn encode_format_field(
    buf: &mut Vec<u8>,
    samples: &[seqair_types::SmallVec<SampleValue, 6>],
    field_idx: usize,
) {
    let first_val = samples
        .iter()
        .filter_map(|s| s.get(field_idx))
        .find(|v| !matches!(v, SampleValue::Missing));

    match first_val {
        Some(SampleValue::Integer(_)) | None => {
            // r[impl bcf_writer.smallest_int_type]
            let typ =
                smallest_int_type_iter(samples.iter().filter_map(|s| match s.get(field_idx) {
                    Some(SampleValue::Integer(v)) => Some(*v),
                    _ => None,
                }));
            encode_type_byte(buf, 1, typ);
            for sample in samples {
                match sample.get(field_idx) {
                    Some(SampleValue::Integer(v)) => match typ {
                        BCF_BT_INT8 => buf.push(*v as u8),
                        BCF_BT_INT16 => buf.extend_from_slice(&(*v as i16).to_le_bytes()),
                        _ => buf.extend_from_slice(&v.to_le_bytes()),
                    },
                    _ => match typ {
                        BCF_BT_INT8 => buf.push(INT8_MISSING),
                        BCF_BT_INT16 => buf.extend_from_slice(&INT16_MISSING.to_le_bytes()),
                        _ => buf.extend_from_slice(&INT32_MISSING.to_le_bytes()),
                    },
                }
            }
        }
        Some(SampleValue::Float(_)) => {
            encode_type_byte(buf, 1, BCF_BT_FLOAT);
            for sample in samples {
                match sample.get(field_idx) {
                    Some(SampleValue::Float(v)) => buf.extend_from_slice(&v.to_le_bytes()),
                    _ => buf.extend_from_slice(&FLOAT_MISSING.to_le_bytes()),
                }
            }
        }
        Some(SampleValue::String(s)) => {
            let max_len = samples
                .iter()
                .filter_map(|sv| match sv.get(field_idx) {
                    Some(SampleValue::String(s)) => Some(s.len()),
                    _ => None,
                })
                .max()
                .unwrap_or(s.len());

            encode_type_byte(buf, max_len, BCF_BT_CHAR);
            for sample in samples {
                match sample.get(field_idx) {
                    Some(SampleValue::String(sv)) => {
                        buf.extend_from_slice(sv.as_bytes());
                        for _ in sv.len()..max_len {
                            buf.push(0);
                        }
                    }
                    _ => {
                        buf.push(b'.');
                        for _ in 1..max_len {
                            buf.push(0);
                        }
                    }
                }
            }
        }
        Some(SampleValue::IntegerArray(arr)) => {
            let max_len = samples
                .iter()
                .filter_map(|sv| match sv.get(field_idx) {
                    Some(SampleValue::IntegerArray(a)) => Some(a.len()),
                    _ => None,
                })
                .max()
                .unwrap_or(arr.len());

            // r[impl bcf_writer.smallest_int_type]
            let typ = smallest_int_type_iter(
                samples
                    .iter()
                    .filter_map(|s| match s.get(field_idx) {
                        Some(SampleValue::IntegerArray(a)) => Some(a),
                        _ => None,
                    })
                    .flat_map(|a| a.iter().filter_map(|v| *v)),
            );
            encode_type_byte(buf, max_len, typ);
            for sample in samples {
                match sample.get(field_idx) {
                    Some(SampleValue::IntegerArray(a)) => {
                        for v in a.iter().take(max_len) {
                            encode_int_value_or_missing(buf, *v, typ);
                        }
                        // r[impl bcf_writer.end_of_vector]
                        for _ in a.len()..max_len {
                            encode_int_eov(buf, typ);
                        }
                    }
                    _ => {
                        encode_int_missing(buf, typ);
                        for _ in 1..max_len {
                            encode_int_eov(buf, typ);
                        }
                    }
                }
            }
        }
        Some(SampleValue::FloatArray(arr)) => {
            let max_len = samples
                .iter()
                .filter_map(|sv| match sv.get(field_idx) {
                    Some(SampleValue::FloatArray(a)) => Some(a.len()),
                    _ => None,
                })
                .max()
                .unwrap_or(arr.len());

            encode_type_byte(buf, max_len, BCF_BT_FLOAT);
            for sample in samples {
                match sample.get(field_idx) {
                    Some(SampleValue::FloatArray(a)) => {
                        for v in a.iter().take(max_len) {
                            match v {
                                Some(f) => buf.extend_from_slice(&f.to_le_bytes()),
                                None => buf.extend_from_slice(&FLOAT_MISSING.to_le_bytes()),
                            }
                        }
                        for _ in a.len()..max_len {
                            buf.extend_from_slice(&FLOAT_END_OF_VECTOR.to_le_bytes());
                        }
                    }
                    _ => {
                        buf.extend_from_slice(&FLOAT_MISSING.to_le_bytes());
                        for _ in 1..max_len {
                            buf.extend_from_slice(&FLOAT_END_OF_VECTOR.to_le_bytes());
                        }
                    }
                }
            }
        }
        _ => {
            encode_type_byte(buf, 0, BCF_BT_NULL);
        }
    }
}

/// Encode an array of BcfValue items. For i32, scans all values to select
/// the smallest BCF int type that fits ALL values per r[bcf_writer.smallest_int_type].
fn encode_array_values<T: BcfValue>(buf: &mut Vec<u8>, values: &[T]) {
    if values.is_empty() {
        encode_type_byte(buf, 0, T::TYPE_CODE);
        return;
    }
    // Find the widest type needed across all values.
    let mut type_code: u8 = 0;
    for v in values {
        let vtc = v.scalar_type_code();
        if vtc > type_code {
            type_code = vtc;
        }
    }
    if type_code == 0 {
        type_code = T::TYPE_CODE;
    }
    encode_type_byte(buf, values.len(), type_code);
    for v in values {
        v.encode_bcf_as(buf, type_code);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::record::Genotype;
    use seqair_types::Base;

    /// Create encoder with in-memory buffers for testing (no BGZF).
    fn test_encoder<'a>(
        shared: &'a mut Vec<u8>,
        indiv: &'a mut Vec<u8>,
        bgzf_buf: &'a mut TestBgzf,
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
        fn flush_if_needed(&mut self, _upcoming: usize) -> Result<(), crate::bam::bgzf::BgzfError> {
            Ok(())
        }
        fn write_all(&mut self, data: &[u8]) -> Result<(), crate::bam::bgzf::BgzfError> {
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
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let pos = Pos::<One>::new(100).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), pos, Some(30.0)).unwrap();

        assert_eq!(enc.n_allele, 2);
        assert_eq!(enc.n_alt, 1);
        assert_eq!(enc.tid, 0);
        assert_eq!(enc.pos_0based, 99); // 0-based
        assert_eq!(enc.rlen, 1);
        assert!(shared.len() >= 24, "must have at least 24-byte header");
    }

    // r[verify bcf_encoder.handle_encode]
    #[test]
    fn scalar_info_encode() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();

        let dp = ScalarInfoHandle::<i32> { dict_idx: 1, _marker: PhantomData };
        dp.encode(&mut enc, 50);
        assert_eq!(enc.n_info, 1);
    }

    // r[verify bcf_encoder.handle_encode]
    #[test]
    fn flag_info_encode() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::reference(Base::A);
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();

        let db = FlagInfoHandle { dict_idx: 4 };
        db.encode(&mut enc);
        assert_eq!(enc.n_info, 1);
    }

    // r[verify bcf_encoder.handle_encode]
    #[test]
    fn per_allele_info_encode() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();

        let ad = PerAlleleInfoHandle::<i32> { dict_idx: 3, _marker: PhantomData };
        ad.encode(&mut enc, &[30, 20]); // n_allele = 2
        assert_eq!(enc.n_info, 1);
    }

    // r[verify bcf_encoder.handle_encode]
    #[test]
    fn gt_format_encode() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();
        enc.begin_samples(1);

        let gt_handle = GtFormatHandle { dict_idx: 0 };
        let gt = Genotype::unphased(0, 1);
        gt_handle.encode(&mut enc, &gt);
        assert_eq!(enc.n_fmt, 1);
        // GT bytes in indiv_buf: key + type_byte + 2 allele bytes
        // Allele 0 unphased: (0+1)<<1|0 = 2
        // Allele 1 unphased: (1+1)<<1|0 = 4
        assert!(indiv.len() >= 4);
    }

    // r[verify bcf_encoder.emit]
    #[test]
    fn emit_patches_header_and_writes() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles
            .begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(100).unwrap(), Some(30.0))
            .unwrap();
        FilterHandle::PASS.encode(&mut enc);

        let dp = ScalarInfoHandle::<i32> { dict_idx: 1, _marker: PhantomData };
        dp.encode(&mut enc, 50);

        enc.begin_samples(1);
        let gt_handle = GtFormatHandle { dict_idx: 0 };
        gt_handle.encode(&mut enc, &Genotype::unphased(0, 1));

        enc.emit().unwrap();

        // Verify header was patched: n_info=1, n_allele=2, n_fmt=1, n_sample=1
        let patched_info_allele =
            u32::from_le_bytes([shared[16], shared[17], shared[18], shared[19]]);
        assert_eq!(patched_info_allele >> 16, 2); // n_allele
        assert_eq!(patched_info_allele & 0xFFFF, 1); // n_info

        let patched_fmt_sample =
            u32::from_le_bytes([shared[20], shared[21], shared[22], shared[23]]);
        assert_eq!(patched_fmt_sample >> 24, 1); // n_fmt
        assert_eq!(patched_fmt_sample & 0xFFFFFF, 1); // n_sample

        // Verify data was written to bgzf
        assert!(!bgzf.data.is_empty());
    }

    // r[verify bcf_encoder.handle_encode]
    #[test]
    fn filter_handle_encodes_large_dict_index() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::reference(Base::A);
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();

        // Dict index 200 exceeds int8 range [-120, 127] — must use int16
        let filter = FilterHandle(200);
        filter.encode(&mut enc);

        // The encoded filter should contain 200 as int16, not truncated to u8
        // Find the filter bytes after the 24-byte header + ID + alleles
        // The filter is encoded as a typed int vec: type_byte + value
        let filter_region = &shared[24..]; // after fixed header
        // Skip ID (type_byte + '.') and REF allele (type_byte + 'A')
        // ID: 1 byte type + 1 byte '.' = 2 bytes
        // REF: 1 byte type + 1 byte 'A' = 2 bytes
        // Filter starts at offset 4 from end of header
        let filter_start = 4; // 2 (ID) + 2 (REF allele)
        let type_byte = filter_region[filter_start];
        let type_code = type_byte & 0x0F;
        assert_eq!(type_code, BCF_BT_INT16, "filter dict index 200 must use int16, not int8");
        // Read the int16 value
        let val =
            i16::from_le_bytes([filter_region[filter_start + 1], filter_region[filter_start + 2]]);
        assert_eq!(val, 200, "filter dict index must be 200");
    }

    // r[verify bcf_encoder.info_counting]
    #[test]
    fn info_counting_increments() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None).unwrap();

        assert_eq!(enc.n_info, 0);
        ScalarInfoHandle::<i32> { dict_idx: 1, _marker: PhantomData }.encode(&mut enc, 10);
        assert_eq!(enc.n_info, 1);
        FlagInfoHandle { dict_idx: 4 }.encode(&mut enc);
        assert_eq!(enc.n_info, 2);
        ScalarInfoHandle::<f32> { dict_idx: 2, _marker: PhantomData }.encode(&mut enc, 0.5);
        assert_eq!(enc.n_info, 3);
    }
}
