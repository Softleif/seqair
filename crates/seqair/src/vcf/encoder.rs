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

use super::error::VcfError;
use super::index_builder::IndexBuilder;
use super::record::Genotype;
use crate::bam::bgzf::VirtualOffset;
use crate::bam::bgzf_writer::BgzfWriter;
use std::io::Write;
use std::marker::PhantomData;

// Re-use BCF constants from bcf_writer
const BCF_BT_NULL: u8 = 0;
const BCF_BT_INT8: u8 = 1;
const BCF_BT_INT16: u8 = 2;
const BCF_BT_INT32: u8 = 3;
const BCF_BT_FLOAT: u8 = 5;
const BCF_BT_CHAR: u8 = 7;

// Sentinel values — some are used by BcfValue trait impls, others reserved for multi-sample
#[allow(dead_code)]
const INT8_MISSING: u8 = 0x80;
#[allow(dead_code)]
const INT16_MISSING: u16 = 0x8000;
const INT32_MISSING: u32 = 0x80000000;
const FLOAT_MISSING: u32 = 0x7F800001;

#[allow(dead_code)]
const INT8_END_OF_VECTOR: u8 = 0x81;
#[allow(dead_code)]
const INT16_END_OF_VECTOR: u16 = 0x8001;
#[allow(dead_code)]
const INT32_END_OF_VECTOR: u32 = 0x80000001;
const FLOAT_END_OF_VECTOR: u32 = 0x7F800002;

// Int ranges for BCF type selection — MIN values reserved for multi-sample smallest_int_type
#[allow(dead_code)]
const INT8_MIN: i32 = -120;
const INT8_MAX: i32 = 127;
#[allow(dead_code)]
const INT16_MIN: i32 = -32760;
const INT16_MAX: i32 = 32767;

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
        // Single int8 vector containing the filter index
        encode_type_byte(enc.shared_buf, 1, BCF_BT_INT8);
        enc.shared_buf.push(self.0 as u8);
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
        // GT values fit int8 for allele indices < 63
        encode_type_byte(enc.indiv_buf, ploidy, BCF_BT_INT8);

        for (i, allele_opt) in gt.alleles.iter().enumerate() {
            let encoded: u8 = match allele_opt {
                Some(idx) => {
                    let phased = if i == 0 {
                        false
                    } else {
                        gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                    };
                    let val = (i32::from(*idx).saturating_add(1) << 1) | (phased as i32);
                    val as u8
                }
                None => 0, // missing allele
            };
            enc.indiv_buf.push(encoded);
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

        let l_shared = self.shared_buf.len() as u32;
        let l_indiv = self.indiv_buf.len() as u32;
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
            index
                .push(self.tid, beg, end, self.bgzf.virtual_offset())
                .map_err(|e| VcfError::Io(std::io::Error::other(e.to_string())))?;
        }

        Ok(())
    }
}

// ── Alleles integration ─────────────────────────────────────────────────

use super::alleles::Alleles;
use seqair_types::{One, Pos};

impl Alleles {
    // r[impl bcf_encoder.begin_record]
    /// Begin a BCF record: write the 24-byte fixed header, ID, and alleles.
    /// Sets n_allele/n_alt on the encoder for downstream field validation.
    pub fn begin_record(
        &self,
        enc: &mut BcfRecordEncoder<'_>,
        contig: ContigHandle,
        pos: Pos<One>,
        qual: Option<f32>,
    ) {
        enc.shared_buf.clear();
        enc.indiv_buf.clear();
        enc.n_info = 0;
        enc.n_fmt = 0;
        enc.n_sample = 0;

        enc.tid = contig.0 as i32;
        enc.pos_0based = pos.to_zero_based().get() as i32;
        enc.rlen = self.rlen() as i32;
        enc.n_allele = self.n_allele() as u16;
        enc.n_alt = self.n_allele().saturating_sub(1) as u16;

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

        // REF allele (zero-alloc)
        let ref_start = enc.shared_buf.len();
        // Reserve space for type byte, fill in after
        enc.shared_buf.push(0); // placeholder type byte
        self.write_ref_into(enc.shared_buf);
        let ref_len = enc.shared_buf.len().saturating_sub(ref_start).saturating_sub(1);
        // Patch the type byte
        if let Some(b) = enc.shared_buf.get_mut(ref_start) {
            *b = ((ref_len.min(14) as u8) << 4) | BCF_BT_CHAR;
        }

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
    }
}

// ── Encoding helpers ────────────────────────────────────────────────────

fn encode_type_byte(buf: &mut Vec<u8>, count: usize, type_code: u8) {
    if count < 15 {
        buf.push(((count as u8) << 4) | type_code);
    } else {
        buf.push((15 << 4) | type_code);
        if count <= INT8_MAX as usize {
            buf.push((1 << 4) | BCF_BT_INT8);
            buf.push(count as u8);
        } else if count <= INT16_MAX as usize {
            buf.push((1 << 4) | BCF_BT_INT16);
            buf.extend_from_slice(&(count as u16).to_le_bytes());
        } else {
            buf.push((1 << 4) | BCF_BT_INT32);
            buf.extend_from_slice(&(count as u32).to_le_bytes());
        }
    }
}

/// Encode an array of BcfValue items. For i32, scans all values to select
/// the smallest BCF int type per r[bcf_writer.smallest_int_type].
/// For f32, uses BCF_BT_FLOAT directly (no type selection needed).
/// Encode an array of BcfValue items. For i32, scans all values to select
/// the smallest BCF int type that fits ALL values per r[bcf_writer.smallest_int_type].
fn encode_array_values<T: BcfValue>(buf: &mut Vec<u8>, values: &[T]) {
    if values.is_empty() {
        encode_type_byte(buf, 0, T::TYPE_CODE);
        return;
    }
    // Find the widest type needed across all values.
    // Start from 0 (narrowest) and widen as needed.
    let mut type_code: u8 = 0;
    for v in values {
        let vtc = v.scalar_type_code();
        if vtc > type_code {
            type_code = vtc;
        }
    }
    // If no values contributed (shouldn't happen since !is_empty), fall back to TYPE_CODE
    if type_code == 0 {
        type_code = T::TYPE_CODE;
    }
    encode_type_byte(buf, values.len(), type_code);
    for v in values {
        v.encode_bcf_as(buf, type_code);
    }
}

fn encode_typed_string(buf: &mut Vec<u8>, s: &[u8]) {
    encode_type_byte(buf, s.len(), BCF_BT_CHAR);
    buf.extend_from_slice(s);
}

/// Encode a dictionary index as a typed int (used for INFO/FORMAT/FILTER keys).
fn encode_typed_int_key(buf: &mut Vec<u8>, dict_idx: u32) {
    if dict_idx <= INT8_MAX as u32 {
        encode_type_byte(buf, 1, BCF_BT_INT8);
        buf.push(dict_idx as u8);
    } else if dict_idx <= INT16_MAX as u32 {
        encode_type_byte(buf, 1, BCF_BT_INT16);
        buf.extend_from_slice(&(dict_idx as u16).to_le_bytes());
    } else {
        encode_type_byte(buf, 1, BCF_BT_INT32);
        buf.extend_from_slice(&dict_idx.to_le_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::header::{
        ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader,
    };
    use crate::vcf::record::Genotype;
    use seqair_types::{Base, SmolStr};
    use std::sync::Arc;

    fn test_header() -> Arc<VcfHeader> {
        Arc::new(
            VcfHeader::builder()
                .add_contig("chr1", ContigDef { length: Some(1000) })
                .unwrap()
                .add_info(
                    "DP",
                    InfoDef {
                        number: Number::Count(1),
                        typ: ValueType::Integer,
                        description: SmolStr::from("Depth"),
                    },
                )
                .unwrap()
                .add_info(
                    "AF",
                    InfoDef {
                        number: Number::AlternateBases,
                        typ: ValueType::Float,
                        description: SmolStr::from("Allele Frequency"),
                    },
                )
                .unwrap()
                .add_info(
                    "AD",
                    InfoDef {
                        number: Number::ReferenceAlternateBases,
                        typ: ValueType::Integer,
                        description: SmolStr::from("Allele Depth"),
                    },
                )
                .unwrap()
                .add_info(
                    "DB",
                    InfoDef {
                        number: Number::Count(0),
                        typ: ValueType::Flag,
                        description: SmolStr::from("dbSNP"),
                    },
                )
                .unwrap()
                .add_format(
                    "GT",
                    FormatDef {
                        number: Number::Count(1),
                        typ: ValueType::String,
                        description: SmolStr::from("Genotype"),
                    },
                )
                .unwrap()
                .add_format(
                    "DP",
                    FormatDef {
                        number: Number::Count(1),
                        typ: ValueType::Integer,
                        description: SmolStr::from("Read Depth"),
                    },
                )
                .unwrap()
                .add_sample("S1")
                .unwrap()
                .build()
                .unwrap(),
        )
    }

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
        alleles.begin_record(&mut enc, ContigHandle(0), pos, Some(30.0));

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
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None);

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
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None);

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
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None);

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
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None);
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
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(100).unwrap(), Some(30.0));
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

    // r[verify bcf_encoder.info_counting]
    #[test]
    fn info_counting_increments() {
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut bgzf = TestBgzf::new();
        let mut enc = test_encoder(&mut shared, &mut indiv, &mut bgzf);

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        alleles.begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(1).unwrap(), None);

        assert_eq!(enc.n_info, 0);
        ScalarInfoHandle::<i32> { dict_idx: 1, _marker: PhantomData }.encode(&mut enc, 10);
        assert_eq!(enc.n_info, 1);
        FlagInfoHandle { dict_idx: 4 }.encode(&mut enc);
        assert_eq!(enc.n_info, 2);
        ScalarInfoHandle::<f32> { dict_idx: 2, _marker: PhantomData }.encode(&mut enc, 0.5);
        assert_eq!(enc.n_info, 3);
    }
}
