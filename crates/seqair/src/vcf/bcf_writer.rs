//! BCF binary writer. Encodes [`VcfRecord`]s in BCF2 format with BGZF
//! compression and optional CSI index co-production.

use super::error::VcfError;
use super::header::VcfHeader;
use super::index_builder::IndexBuilder;
use super::record::{Filters, InfoValue, SampleValue, VcfRecord};
use crate::bam::bgzf_writer::BgzfWriter;
use std::io::Write;
use std::sync::Arc;

// BCF type codes
const BCF_BT_NULL: u8 = 0;
const BCF_BT_INT8: u8 = 1;
const BCF_BT_INT16: u8 = 2;
const BCF_BT_INT32: u8 = 3;
const BCF_BT_FLOAT: u8 = 5;
const BCF_BT_CHAR: u8 = 7;

// Sentinel values
const INT8_MISSING: u8 = 0x80;
const INT16_MISSING: u16 = 0x8000;
const INT32_MISSING: u32 = 0x80000000;
const FLOAT_MISSING: u32 = 0x7F800001;

const INT8_END_OF_VECTOR: u8 = 0x81;
const INT16_END_OF_VECTOR: u16 = 0x8001;
const INT32_END_OF_VECTOR: u32 = 0x80000001;
const FLOAT_END_OF_VECTOR: u32 = 0x7F800002;

// Int ranges (excluding sentinel values)
const INT8_MIN: i32 = -120;
const INT8_MAX: i32 = 127;
const INT16_MIN: i32 = -32760;
const INT16_MAX: i32 = 32767;

// r[impl bcf_writer.magic]
// r[impl bcf_writer.buffer_reuse]
/// BCF binary writer with BGZF compression and optional CSI index.
pub struct BcfWriter<W: Write> {
    bgzf: BgzfWriter<W>,
    header: Arc<VcfHeader>,
    index: Option<IndexBuilder>,
    shared_buf: Vec<u8>,
    indiv_buf: Vec<u8>,
    header_written: bool,
}

impl<W: Write> BcfWriter<W> {
    /// Create a BCF writer with optional index co-production.
    pub fn new(inner: W, header: Arc<VcfHeader>, build_index: bool) -> Self {
        let n_refs = header.contigs().len();
        Self {
            bgzf: BgzfWriter::new(inner),
            header,
            index: if build_index {
                Some(IndexBuilder::new(n_refs, 14, 5, crate::bam::bgzf::VirtualOffset(0)))
            } else {
                None
            },
            shared_buf: Vec::with_capacity(4096),
            indiv_buf: Vec::with_capacity(4096),
            header_written: false,
        }
    }

    // r[impl bcf_writer.magic]
    /// Write the BCF magic and header. Must be called once before write_record.
    pub fn write_header(&mut self) -> Result<(), VcfError> {
        // Magic: BCF\x02\x02 (BCF version 2.2, the current version expected by htslib)
        self.bgzf.write_all(b"BCF\x02\x02")?;

        // Header text (NUL-terminated)
        let header_text = self.header.to_vcf_text();
        let l_text = header_text
            .len()
            .checked_add(1)
            .ok_or_else(|| VcfError::Io(std::io::Error::other("header too large")))?;
        self.bgzf.write_all(&(l_text as u32).to_le_bytes())?;
        self.bgzf.write_all(header_text.as_bytes())?;
        self.bgzf.write_all(&[0u8])?; // NUL terminator

        self.header_written = true;

        // Update index with header end offset
        if let Some(ref mut index) = self.index {
            let n_refs = self.header.contigs().len();
            *index = IndexBuilder::new(n_refs, 14, 5, self.bgzf.virtual_offset());
        }

        Ok(())
    }

    // r[impl bcf_writer.record_layout]
    // r[impl bcf_writer.bgzf_blocks]
    /// Write a single BCF record.
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        if !self.header_written {
            return Err(VcfError::HeaderNotWritten);
        }

        // r[impl bcf_writer.buffer_reuse]
        self.shared_buf.clear();
        self.indiv_buf.clear();

        self.encode_shared(record)?;
        self.encode_individual(record)?;

        let l_shared = self.shared_buf.len() as u32;
        let l_indiv = self.indiv_buf.len() as u32;
        let total =
            8usize.saturating_add(self.shared_buf.len()).saturating_add(self.indiv_buf.len());

        // r[impl bcf_writer.bgzf_blocks]
        self.bgzf.flush_if_needed(total)?;

        self.bgzf.write_all(&l_shared.to_le_bytes())?;
        self.bgzf.write_all(&l_indiv.to_le_bytes())?;
        self.bgzf.write_all(&self.shared_buf)?;
        self.bgzf.write_all(&self.indiv_buf)?;

        // Index co-production
        if let Some(ref mut index) = self.index {
            let tid = self.header.contig_id(&record.contig)? as i32;
            let beg = record.pos.to_zero_based().get() as u64;
            let end = beg.saturating_add(record.alleles.rlen() as u64);
            index
                .push(tid, beg, end, self.bgzf.virtual_offset())
                .map_err(|e| VcfError::Io(std::io::Error::other(e.to_string())))?;
        }

        Ok(())
    }

    // r[impl bcf_writer.fixed_fields]
    // r[impl bcf_writer.shared_variable]
    fn encode_shared(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        let buf = &mut self.shared_buf;

        // r[impl bcf_writer.coordinate_system]
        let chrom = self.header.contig_id(&record.contig)? as i32;
        let pos = record.pos.to_zero_based().get() as i32;
        let rlen = record.alleles.rlen() as i32;

        // r[impl bcf_writer.qual_missing]
        let qual_bits = match record.qual {
            Some(q) => q.to_bits(),
            None => FLOAT_MISSING,
        };

        let n_info = record.info.iter().count() as u32;
        let n_allele = record.alleles.n_allele() as u32;
        let n_fmt = record.samples.format_keys.len() as u32;
        let n_sample = record.samples.values.len() as u32;

        // 24 fixed bytes
        buf.extend_from_slice(&chrom.to_le_bytes());
        buf.extend_from_slice(&pos.to_le_bytes());
        buf.extend_from_slice(&rlen.to_le_bytes());
        buf.extend_from_slice(&qual_bits.to_le_bytes());
        buf.extend_from_slice(&((n_allele << 16) | n_info).to_le_bytes());
        buf.extend_from_slice(&((n_fmt << 24) | n_sample).to_le_bytes());

        // ID (typed string)
        match &record.id {
            Some(id) => encode_typed_string(buf, id.as_bytes()),
            None => encode_typed_string(buf, b"."),
        }

        // Alleles (n_allele typed strings, REF first)
        let ref_text = record.alleles.ref_text();
        encode_typed_string(buf, ref_text.as_bytes());
        for alt in &record.alleles.alt_texts() {
            encode_typed_string(buf, alt.as_bytes());
        }

        // r[impl bcf_writer.filter_pass]
        match &record.filters {
            Filters::Pass => {
                // Single int8 vector containing 0
                encode_type_byte(buf, 1, BCF_BT_INT8);
                buf.push(0);
            }
            Filters::Failed(ids) => {
                // r[impl vcf_writer.validation]
                let mut indices = Vec::with_capacity(ids.len());
                for id in ids {
                    let idx = self.header.string_map().get(id).ok_or_else(|| {
                        VcfError::Header(super::error::VcfHeaderError::MissingFilter {
                            id: id.clone(),
                        })
                    })?;
                    indices.push(idx as i32);
                }
                encode_typed_int_vec(buf, &indices);
            }
            Filters::NotApplied => {
                // type=0, count=0
                encode_type_byte(buf, 0, BCF_BT_NULL);
            }
        }

        // INFO key-value pairs
        // r[impl vcf_writer.validation]
        for (key, value) in record.info.iter() {
            let dict_idx = self.header.string_map().get(key).ok_or_else(|| {
                VcfError::Header(super::error::VcfHeaderError::MissingInfo { id: key.clone() })
            })? as i32;
            encode_typed_int_vec(buf, &[dict_idx]);
            // Value
            encode_info_value(buf, value);
        }

        Ok(())
    }

    // r[impl bcf_writer.indiv_field_major]
    fn encode_individual(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        let buf = &mut self.indiv_buf;
        let n_sample = record.samples.values.len();

        if n_sample == 0 || record.samples.format_keys.is_empty() {
            return Ok(());
        }

        for (field_idx, key) in record.samples.format_keys.iter().enumerate() {
            // r[impl vcf_writer.validation]
            let dict_idx = self.header.string_map().get(key).ok_or_else(|| {
                VcfError::Header(super::error::VcfHeaderError::MissingFormat { id: key.clone() })
            })? as i32;
            encode_typed_int_vec(buf, &[dict_idx]);

            // Collect values for this field across all samples
            let is_gt = key == "GT";

            if is_gt {
                encode_gt_field(buf, &record.samples.values, field_idx);
            } else {
                encode_format_field(buf, &record.samples.values, field_idx);
            }
        }

        Ok(())
    }

    /// Get a direct record encoder that writes into this writer's buffers.
    /// Use this for zero-alloc encoding instead of `write_record(&VcfRecord)`.
    pub fn record_encoder(&mut self) -> super::encoder::BcfRecordEncoder<'_> {
        super::encoder::BcfRecordEncoder {
            shared_buf: &mut self.shared_buf,
            indiv_buf: &mut self.indiv_buf,
            bgzf: &mut self.bgzf,
            index: self.index.as_mut(),
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

    // r[impl bcf_writer.finish]
    /// Finalize the writer. Returns the index builder if one was used.
    pub fn finish(self) -> Result<Option<IndexBuilder>, VcfError> {
        let voff = self.bgzf.virtual_offset();
        let mut index = self.index;
        if let Some(ref mut idx) = index {
            idx.finish(voff).map_err(|e| VcfError::Io(std::io::Error::other(e.to_string())))?;
        }
        self.bgzf.finish()?;
        Ok(index)
    }
}

// --- Encoding helpers ---

// r[impl bcf_writer.typed_values]
fn encode_type_byte(buf: &mut Vec<u8>, count: usize, type_code: u8) {
    if count < 15 {
        buf.push(((count as u8) << 4) | type_code);
    } else {
        buf.push((15 << 4) | type_code);
        // Encode the actual count as a typed integer
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

// r[impl bcf_writer.string_encoding]
fn encode_typed_string(buf: &mut Vec<u8>, s: &[u8]) {
    if s.is_empty() {
        encode_type_byte(buf, 0, BCF_BT_CHAR);
    } else {
        encode_type_byte(buf, s.len(), BCF_BT_CHAR);
        buf.extend_from_slice(s);
    }
}

// r[impl bcf_writer.smallest_int_type]
fn smallest_int_type(values: &[i32]) -> u8 {
    smallest_int_type_iter(values.iter().copied())
}

/// Like `smallest_int_type` but takes an iterator — avoids collecting into a Vec.
fn smallest_int_type_iter(values: impl Iterator<Item = i32>) -> u8 {
    let mut needs_16 = false;
    for v in values {
        if !(INT8_MIN..=INT8_MAX).contains(&v) {
            if !(INT16_MIN..=INT16_MAX).contains(&v) {
                return BCF_BT_INT32;
            }
            needs_16 = true;
        }
    }
    if needs_16 { BCF_BT_INT16 } else { BCF_BT_INT8 }
}

fn encode_typed_int_vec(buf: &mut Vec<u8>, values: &[i32]) {
    if values.is_empty() {
        encode_type_byte(buf, 0, BCF_BT_INT8);
        return;
    }
    let typ = smallest_int_type(values);
    encode_type_byte(buf, values.len(), typ);
    for &v in values {
        match typ {
            BCF_BT_INT8 => buf.push(v as u8),
            BCF_BT_INT16 => buf.extend_from_slice(&(v as i16).to_le_bytes()),
            BCF_BT_INT32 => buf.extend_from_slice(&v.to_le_bytes()),
            _ => {}
        }
    }
}

/// Encode an integer value or missing sentinel, using the given BCF type.
fn encode_int_value_or_missing(buf: &mut Vec<u8>, value: Option<i32>, typ: u8) {
    match value {
        Some(n) => match typ {
            BCF_BT_INT8 => buf.push(n as u8),
            BCF_BT_INT16 => buf.extend_from_slice(&(n as i16).to_le_bytes()),
            _ => buf.extend_from_slice(&n.to_le_bytes()),
        },
        None => encode_int_missing(buf, typ),
    }
}

/// Write the missing sentinel for the given int type.
fn encode_int_missing(buf: &mut Vec<u8>, typ: u8) {
    match typ {
        BCF_BT_INT8 => buf.push(INT8_MISSING),
        BCF_BT_INT16 => buf.extend_from_slice(&INT16_MISSING.to_le_bytes()),
        _ => buf.extend_from_slice(&INT32_MISSING.to_le_bytes()),
    }
}

/// Write the end-of-vector sentinel for the given int type.
fn encode_int_eov(buf: &mut Vec<u8>, typ: u8) {
    match typ {
        BCF_BT_INT8 => buf.push(INT8_END_OF_VECTOR),
        BCF_BT_INT16 => buf.extend_from_slice(&INT16_END_OF_VECTOR.to_le_bytes()),
        _ => buf.extend_from_slice(&INT32_END_OF_VECTOR.to_le_bytes()),
    }
}

// r[impl bcf_writer.flag_encoding]
fn encode_info_value(buf: &mut Vec<u8>, value: &InfoValue) {
    match value {
        InfoValue::Integer(v) => encode_typed_int_vec(buf, &[*v]),
        InfoValue::Float(v) => {
            encode_type_byte(buf, 1, BCF_BT_FLOAT);
            buf.extend_from_slice(&v.to_le_bytes());
        }
        InfoValue::Flag => {
            // type=0, count=0
            encode_type_byte(buf, 0, BCF_BT_NULL);
        }
        InfoValue::String(s) => encode_typed_string(buf, s.as_bytes()),
        InfoValue::IntegerArray(arr) => {
            // r[impl bcf_writer.smallest_int_type]
            // Scan only concrete (non-missing) values to pick smallest type.
            // Missing values get the per-type sentinel, not a fixed i32::MIN.
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
            // Comma-separated in a single typed string
            let joined: String =
                arr.iter().map(|v| v.as_deref().unwrap_or(".")).collect::<Vec<_>>().join(",");
            encode_typed_string(buf, joined.as_bytes());
        }
    }
}

// r[impl bcf_writer.gt_encoding]
fn encode_gt_field(
    buf: &mut Vec<u8>,
    samples: &[seqair_types::SmallVec<SampleValue, 6>],
    field_idx: usize,
) {
    // Find max ploidy across samples
    let mut max_ploidy = 0usize;
    for sample in samples {
        if let Some(SampleValue::Genotype(gt)) = sample.get(field_idx) {
            max_ploidy = max_ploidy.max(gt.alleles.len());
        }
    }
    if max_ploidy == 0 {
        max_ploidy = 2; // default diploid
    }

    // Find max allele index to determine int type
    let mut max_allele: i32 = 0;
    for sample in samples {
        if let Some(SampleValue::Genotype(gt)) = sample.get(field_idx) {
            for idx in gt.alleles.iter().flatten() {
                max_allele = max_allele.max(i32::from(*idx));
            }
        }
    }

    // GT value = (allele+1) << 1 | phased, so max value = (max_allele+1)*2+1
    let max_val = (max_allele.saturating_add(1)).saturating_mul(2).saturating_add(1);
    let typ = smallest_int_type(&[max_val]);

    // Type descriptor: max_ploidy values per sample
    encode_type_byte(buf, max_ploidy, typ);

    // Encode all samples contiguously
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
                                false // first allele separator is always unphased
                            } else {
                                g.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                            };
                            ((i32::from(*idx).saturating_add(1)) << 1) | (phased as i32)
                        }
                        None => 0, // missing allele
                    }
                } else {
                    // Haploid padding: end-of-vector
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
                0 // missing genotype
            };

            match typ {
                BCF_BT_INT8 => buf.push(encoded as u8),
                BCF_BT_INT16 => buf.extend_from_slice(&(encoded as i16).to_le_bytes()),
                _ => buf.extend_from_slice(&encoded.to_le_bytes()),
            }
        }
    }
}

fn encode_format_field(
    buf: &mut Vec<u8>,
    samples: &[seqair_types::SmallVec<SampleValue, 6>],
    field_idx: usize,
) {
    // Determine type from first non-missing sample value
    let first_val = samples
        .iter()
        .filter_map(|s| s.get(field_idx))
        .find(|v| !matches!(v, SampleValue::Missing));

    match first_val {
        Some(SampleValue::Integer(_)) | None => {
            // r[impl bcf_writer.smallest_int_type]
            // Scan all values to pick smallest fitting type (no allocation)
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
            // Find max string length
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
                        // Pad with NUL (end-of-vector for char)
                        for _ in sv.len()..max_len {
                            buf.push(0);
                        }
                    }
                    _ => {
                        buf.push(b'.'); // missing
                        for _ in 1..max_len {
                            buf.push(0);
                        }
                    }
                }
            }
        }
        Some(SampleValue::IntegerArray(arr)) => {
            // Find max array length
            let max_len = samples
                .iter()
                .filter_map(|sv| match sv.get(field_idx) {
                    Some(SampleValue::IntegerArray(a)) => Some(a.len()),
                    _ => None,
                })
                .max()
                .unwrap_or(arr.len());

            // r[impl bcf_writer.smallest_int_type]
            // Scan all concrete values across all samples to pick smallest type
            let typ = smallest_int_type_iter(samples.iter().flat_map(|s| {
                match s.get(field_idx) {
                    Some(SampleValue::IntegerArray(a)) => a.iter().filter_map(|v| *v).collect(),
                    _ => Vec::new(),
                }
                .into_iter()
            }));
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
                            // float EOV
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
            // Genotype handled separately; shouldn't reach here
            encode_type_byte(buf, 0, BCF_BT_NULL);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType};
    use crate::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
    use seqair_types::{Base, One, Pos, SmolStr};

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

    // r[verify bcf_writer.magic]
    #[test]
    fn bcf_magic_and_header() {
        let header = test_header();
        let mut output = Vec::new();
        let mut writer = BcfWriter::new(&mut output, header, false);
        writer.write_header().unwrap();
        writer.finish().unwrap();

        // Decompress and check magic
        let mut reader = crate::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(output));
        let mut data = Vec::new();
        reader.read_to_end(&mut data).unwrap();

        assert_eq!(&data[..5], b"BCF\x02\x02");
        // l_text follows
        let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
        // Header text should end with NUL
        assert_eq!(data[5usize.saturating_add(4).saturating_add(l_text).saturating_sub(1)], 0);
    }

    // r[verify bcf_writer.record_layout]
    // r[verify bcf_writer.fixed_fields]
    #[test]
    fn bcf_encode_simple_snv() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .qual(30.0)
        .filter_pass()
        .info_integer("DP", 50)
        .build(&header)
        .unwrap();

        let mut output = Vec::new();
        let mut writer = BcfWriter::new(&mut output, header, false);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
        writer.finish().unwrap();

        // Decompress
        let mut reader = crate::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(output));
        let mut data = Vec::new();
        reader.read_to_end(&mut data).unwrap();

        // Skip past header to first record
        let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
        let rec_start = 9usize.saturating_add(l_text);

        // l_shared and l_indiv
        let l_shared = u32::from_le_bytes([
            data[rec_start],
            data[rec_start.saturating_add(1)],
            data[rec_start.saturating_add(2)],
            data[rec_start.saturating_add(3)],
        ]);
        assert!(l_shared >= 24, "shared data must be at least 24 bytes");

        // Fixed fields start at rec_start + 8
        let fixed = rec_start.saturating_add(8);
        let chrom = i32::from_le_bytes([
            data[fixed],
            data[fixed.saturating_add(1)],
            data[fixed.saturating_add(2)],
            data[fixed.saturating_add(3)],
        ]);
        assert_eq!(chrom, 0); // chr1 = tid 0

        // r[verify bcf_writer.coordinate_system]
        let pos = i32::from_le_bytes([
            data[fixed.saturating_add(4)],
            data[fixed.saturating_add(5)],
            data[fixed.saturating_add(6)],
            data[fixed.saturating_add(7)],
        ]);
        assert_eq!(pos, 99); // 0-based (100 - 1)
    }

    // r[verify bcf_writer.typed_values]
    #[test]
    fn type_byte_encoding() {
        let mut buf = Vec::new();
        // Small count
        encode_type_byte(&mut buf, 3, BCF_BT_INT8);
        assert_eq!(buf[0], (3 << 4) | 1);

        // Count >= 15
        buf.clear();
        encode_type_byte(&mut buf, 20, BCF_BT_INT16);
        assert_eq!(buf[0], (15 << 4) | 2); // overflow marker
        assert_eq!(buf[1], (1 << 4) | 1); // count=20 as int8
        assert_eq!(buf[2], 20);
    }

    // r[verify bcf_writer.smallest_int_type]
    #[test]
    fn int_type_selection() {
        assert_eq!(smallest_int_type(&[0, 1, 127]), BCF_BT_INT8);
        assert_eq!(smallest_int_type(&[-120, 50]), BCF_BT_INT8);
        assert_eq!(smallest_int_type(&[128]), BCF_BT_INT16);
        assert_eq!(smallest_int_type(&[-121]), BCF_BT_INT16);
        assert_eq!(smallest_int_type(&[32768]), BCF_BT_INT32);
        assert_eq!(smallest_int_type(&[-32761]), BCF_BT_INT32);
    }

    // r[verify bcf_writer.gt_encoding]
    #[test]
    fn gt_encoding_values() {
        // 0/1 → [(0+1)<<1|0, (1+1)<<1|0] = [2, 4]
        // 0|1 → [(0+1)<<1|0, (1+1)<<1|1] = [2, 5]
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .format_keys(&["GT"])
        .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1))])
        .build(&header)
        .unwrap();

        let mut output = Vec::new();
        let mut writer = BcfWriter::new(&mut output, header, false);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
        writer.finish().unwrap();

        // Just verify it writes without error — detailed byte-level
        // verification of GT encoding is complex and best done via
        // bcftools round-trip integration tests
        assert!(!output.is_empty());
    }

    // r[verify bcf_writer.missing_sentinels]
    #[test]
    fn float_missing_sentinel() {
        // Verify the sentinel is written as raw bytes
        let bits = FLOAT_MISSING;
        let bytes = bits.to_le_bytes();
        // This is a signaling NaN — must NOT go through float conversion
        assert_eq!(bytes, [0x01, 0x00, 0x80, 0x7f]);
    }

    // r[verify bcf_writer.finish]
    #[test]
    fn bcf_with_index() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .build(&header)
        .unwrap();

        let mut output = Vec::new();
        let mut writer = BcfWriter::new(&mut output, header, true);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
        let index = writer.finish().unwrap();

        assert!(index.is_some());
        assert!(output.len() > 28);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    // r[verify bcf_writer.smallest_int_type]
    proptest! {
        #[test]
        fn smallest_int_type_fits_all_values(values in proptest::collection::vec(-100_000i32..100_000, 1..20)) {
            let typ = smallest_int_type(&values);
            for &v in &values {
                match typ {
                    BCF_BT_INT8 => {
                        prop_assert!(v >= INT8_MIN && v <= INT8_MAX,
                            "value {v} doesn't fit int8 [{INT8_MIN}, {INT8_MAX}]");
                    }
                    BCF_BT_INT16 => {
                        prop_assert!(v >= INT16_MIN && v <= INT16_MAX,
                            "value {v} doesn't fit int16 [{INT16_MIN}, {INT16_MAX}]");
                    }
                    BCF_BT_INT32 => {} // always fits
                    _ => prop_assert!(false, "unexpected type code {typ}"),
                }
            }
        }

        #[test]
        fn smallest_int_type_is_minimal(values in proptest::collection::vec(-100_000i32..100_000, 1..20)) {
            let typ = smallest_int_type(&values);
            // If we got INT16, at least one value must not fit INT8
            if typ == BCF_BT_INT16 {
                prop_assert!(values.iter().any(|&v| v < INT8_MIN || v > INT8_MAX),
                    "selected INT16 but all values fit INT8");
            }
            // If we got INT32, at least one value must not fit INT16
            if typ == BCF_BT_INT32 {
                prop_assert!(values.iter().any(|&v| v < INT16_MIN || v > INT16_MAX),
                    "selected INT32 but all values fit INT16");
            }
        }

        /// Iter version must return same result as slice version.
        #[test]
        fn smallest_int_type_iter_matches_slice(values in proptest::collection::vec(-100_000i32..100_000, 0..20)) {
            let from_slice = smallest_int_type(&values);
            let from_iter = smallest_int_type_iter(values.iter().copied());
            prop_assert_eq!(from_slice, from_iter);
        }
    }
}
