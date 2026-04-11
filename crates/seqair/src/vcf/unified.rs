//! Unified VCF/BCF writer with typestate enforcement.

use super::OutputFormat;
use super::alleles::Alleles;
use super::bcf_encoding::*;
use super::encoder::{BcfRecordEncoder, BcfValue, BgzfWrite, ContigHandle};
use super::error::VcfError;
use super::header::VcfHeader;
use super::index_builder::IndexBuilder;
use super::record::Genotype;
use super::record_encoder::{ContigId, FieldId, FilterId, FormatEncoder, InfoEncoder};
use super::writer::{percent_encode_into, write_float_g};
use crate::bam::bgzf_writer::BgzfWriter;
use seqair_types::{One, Pos, SmallVec, SmolStr};
use std::io::Write;
use std::marker::PhantomData;

// ── Writer state types ─────────────────────────────────────────────────

/// Writer state: header has not been written yet.
pub struct Unstarted;

/// Writer state: header has been written, ready to encode records.
pub struct Ready;

// ── Record encoder state types ─────────────────────────────────────────

// r[impl record_encoder.typestate_states]
/// Record state: fixed fields written, awaiting filter.
pub struct Begun;

/// Record state: filter written, INFO fields may be encoded.
pub struct Filtered;

/// Record state: sample count declared, FORMAT fields may be encoded.
pub struct WithSamples;

// ── Writer ─────────────────────────────────────────────────────────────

// r[impl record_encoder.writer]
// r[impl record_encoder.writer_typestate]
/// Unified VCF/BCF writer. Output format is selected at construction time.
///
/// See the [module documentation](self) for complete usage examples.
pub struct Writer<W: Write, S = Unstarted> {
    inner: WriterInner<W>,
    _state: PhantomData<S>,
}

enum WriterInner<W: Write> {
    Vcf {
        output: W,
        buf: Vec<u8>,
        fmt_keys: Vec<SmolStr>,
        sample_bufs: Vec<Vec<u8>>,
        n_samples: u32,
    },
    VcfGz {
        bgzf: BgzfWriter<W>,
        index: Option<IndexBuilder>,
        buf: Vec<u8>,
        fmt_keys: Vec<SmolStr>,
        sample_bufs: Vec<Vec<u8>>,
        n_samples: u32,
    },
    Bcf {
        bgzf: BgzfWriter<W>,
        index: Option<IndexBuilder>,
        shared_buf: Vec<u8>,
        indiv_buf: Vec<u8>,
        n_samples: u32,
    },
}

// r[impl record_encoder.writer_new]
impl<W: Write> Writer<W> {
    /// Create a writer for the given output format. No header is needed at this point.
    pub fn new(inner: W, format: OutputFormat) -> Self {
        let inner = match format {
            OutputFormat::Vcf => WriterInner::Vcf {
                output: inner,
                buf: Vec::with_capacity(4096),
                fmt_keys: Vec::with_capacity(8),
                sample_bufs: Vec::new(),
                n_samples: 0,
            },
            OutputFormat::VcfGz => WriterInner::VcfGz {
                bgzf: BgzfWriter::new(inner),
                index: None,
                buf: Vec::with_capacity(4096),
                fmt_keys: Vec::with_capacity(8),
                sample_bufs: Vec::new(),
                n_samples: 0,
            },
            OutputFormat::Bcf => WriterInner::Bcf {
                bgzf: BgzfWriter::new(inner),
                index: None,
                shared_buf: Vec::with_capacity(4096),
                indiv_buf: Vec::with_capacity(4096),
                n_samples: 0,
            },
        };
        Writer { inner, _state: PhantomData }
    }

    // r[impl record_encoder.write_header]
    /// Write the file header. Consumes the `Unstarted` writer and returns `Ready`.
    pub fn write_header(mut self, header: &VcfHeader) -> Result<Writer<W, Ready>, VcfError> {
        let header_text = header.to_vcf_text();
        let n_refs = header.contigs().len();
        #[expect(
            clippy::cast_possible_truncation,
            reason = "sample count bounded by VcfHeader builder; realistic files have <100 samples"
        )]
        let header_n_samples = header.samples().len() as u32;

        match &mut self.inner {
            WriterInner::Vcf { output, n_samples, .. } => {
                output.write_all(header_text.as_bytes()).map_err(VcfError::Io)?;
                *n_samples = header_n_samples;
            }
            WriterInner::VcfGz { bgzf, index, n_samples, .. } => {
                bgzf.write_all(header_text.as_bytes())?;
                *index = Some(IndexBuilder::tbi(n_refs, bgzf.virtual_offset()));
                *n_samples = header_n_samples;
            }
            WriterInner::Bcf { bgzf, index, n_samples, .. } => {
                // r[impl bcf_writer.magic]
                bgzf.write_all(b"BCF\x02\x02")?;
                let l_text = header_text.len().checked_add(1).ok_or(VcfError::HeaderTooLarge)?;
                let l_text_u32 = u32::try_from(l_text).map_err(|_| VcfError::HeaderTooLarge)?;
                bgzf.write_all(&l_text_u32.to_le_bytes())?;
                bgzf.write_all(header_text.as_bytes())?;
                bgzf.write_all(&[0u8])?;
                *index = Some(IndexBuilder::new(n_refs, 14, 5, bgzf.virtual_offset()));
                *n_samples = header_n_samples;
            }
        }

        Ok(Writer { inner: self.inner, _state: PhantomData })
    }
}

impl<W: Write> Writer<W, Ready> {
    // r[impl record_encoder.begin_record]
    /// Begin encoding a new record. Returns an encoder in the `Begun` state.
    pub fn begin_record(
        &mut self,
        contig: &ContigId,
        pos: Pos<One>,
        alleles: &Alleles,
        qual: Option<f32>,
    ) -> Result<RecordEncoder<'_, Begun>, VcfError> {
        let inner = match &mut self.inner {
            WriterInner::Bcf { bgzf, index, shared_buf, indiv_buf, n_samples } => {
                // r[impl bcf_writer.buffer_reuse]
                let mut enc = BcfRecordEncoder {
                    shared_buf,
                    indiv_buf,
                    bgzf,
                    index: index.as_mut(),
                    n_allele: 0,
                    n_alt: 0,
                    n_info: 0,
                    n_fmt: 0,
                    n_sample: *n_samples,
                    tid: 0,
                    pos_0based: 0,
                    rlen: 0,
                };
                alleles.begin_record(&mut enc, ContigHandle(contig.tid()), pos, qual)?;
                EncoderInner::Bcf(enc)
            }
            WriterInner::Vcf { output, buf, fmt_keys, sample_bufs, n_samples } => {
                let output = VcfOutput::Plain(output);
                EncoderInner::Vcf(begin_vcf_record(
                    buf,
                    fmt_keys,
                    sample_bufs,
                    output,
                    contig,
                    pos,
                    alleles,
                    qual,
                    *n_samples,
                )?)
            }
            WriterInner::VcfGz { bgzf, index, buf, fmt_keys, sample_bufs, n_samples } => {
                let output = VcfOutput::Bgzf { bgzf, index: index.as_mut() };
                EncoderInner::Vcf(begin_vcf_record(
                    buf,
                    fmt_keys,
                    sample_bufs,
                    output,
                    contig,
                    pos,
                    alleles,
                    qual,
                    *n_samples,
                )?)
            }
        };

        Ok(RecordEncoder { inner, _state: PhantomData })
    }

    // r[impl record_encoder.finish]
    /// Finalize the writer. Returns the index builder if one was co-produced.
    pub fn finish(self) -> Result<Option<IndexBuilder>, VcfError> {
        match self.inner {
            WriterInner::Vcf { mut output, .. } => {
                output.flush().map_err(VcfError::Io)?;
                Ok(None)
            }
            WriterInner::VcfGz { bgzf, mut index, .. } => {
                // r[impl vcf_writer.finish]
                let voff = bgzf.virtual_offset();
                if let Some(ref mut idx) = index {
                    idx.finish(voff)?;
                }
                bgzf.finish()?;
                Ok(index)
            }
            WriterInner::Bcf { bgzf, mut index, .. } => {
                // r[impl bcf_writer.finish]
                let voff = bgzf.virtual_offset();
                if let Some(ref mut idx) = index {
                    idx.finish(voff)?;
                }
                bgzf.finish()?;
                Ok(index)
            }
        }
    }
}

// ── VCF text output abstraction ────────────────────────────────────────

enum VcfOutput<'a> {
    Plain(&'a mut dyn Write),
    Bgzf { bgzf: &'a mut dyn BgzfWrite, index: Option<&'a mut IndexBuilder> },
}

impl VcfOutput<'_> {
    fn write_line(&mut self, buf: &[u8]) -> Result<(), VcfError> {
        match self {
            VcfOutput::Plain(w) => w.write_all(buf).map_err(VcfError::Io),
            VcfOutput::Bgzf { bgzf, .. } => {
                bgzf.flush_if_needed(buf.len())?;
                bgzf.write_all(buf)?;
                Ok(())
            }
        }
    }

    fn push_index(&mut self, tid: i32, beg: u64, end: u64) -> Result<(), VcfError> {
        match self {
            VcfOutput::Plain(_) => Ok(()),
            VcfOutput::Bgzf { bgzf, index } => {
                if let Some(idx) = index {
                    idx.push(tid, beg, end, bgzf.virtual_offset())?;
                }
                Ok(())
            }
        }
    }
}

// ── VCF text encoder fields ────────────────────────────────────────────

struct VcfEncoderFields<'a> {
    buf: &'a mut Vec<u8>,
    fmt_keys: &'a mut Vec<SmolStr>,
    sample_bufs: &'a mut Vec<Vec<u8>>,
    output: VcfOutput<'a>,
    n_allele: u16,
    n_alt: u16,
    info_count: u16,
    tid: i32,
    pos_0based: u32,
    rlen: u32,
    n_samples: u32,
    filter_written: bool,
}

// r[impl record_encoder.buffer_reuse]
// r[impl vcf_writer.buffer_reuse]
#[expect(clippy::too_many_arguments, reason = "internal helper, borrows writer fields")]
fn begin_vcf_record<'a>(
    buf: &'a mut Vec<u8>,
    fmt_keys: &'a mut Vec<SmolStr>,
    sample_bufs: &'a mut Vec<Vec<u8>>,
    output: VcfOutput<'a>,
    contig: &ContigId,
    pos: Pos<One>,
    alleles: &Alleles,
    qual: Option<f32>,
    n_samples: u32,
) -> Result<VcfEncoderFields<'a>, VcfError> {
    buf.clear();
    fmt_keys.clear();
    // Retain existing Vec allocations for buffer reuse across records (r[record_encoder.buffer_reuse]).
    // Only clear contents; begin_samples will resize if needed.
    for sb in sample_bufs.iter_mut() {
        sb.clear();
    }

    let n_allele = u16::try_from(alleles.n_allele()).map_err(|_| VcfError::ValueOverflow {
        field: "n_allele",
        value: alleles.n_allele() as u64,
        target_type: "u16",
    })?;
    let n_alt = n_allele.saturating_sub(1);
    let tid = i32::try_from(contig.tid()).map_err(|_| VcfError::ValueOverflow {
        field: "contig_tid",
        value: u64::from(contig.tid()),
        target_type: "i32",
    })?;
    let pos_0based = pos.to_zero_based().get();
    let rlen = u32::try_from(alleles.rlen()).map_err(|_| VcfError::ValueOverflow {
        field: "rlen",
        value: alleles.rlen() as u64,
        target_type: "u32",
    })?;

    // r[impl vcf_writer.tab_delimited]
    buf.extend_from_slice(contig.name().as_bytes());
    buf.push(b'\t');
    // r[impl vcf_writer.integer_format]
    let mut itoa_buf = itoa::Buffer::new();
    buf.extend_from_slice(itoa_buf.format(pos.get()).as_bytes());
    buf.push(b'\t');
    buf.push(b'.');
    buf.push(b'\t');
    alleles.write_ref_into(buf);
    buf.push(b'\t');
    let n_alts = alleles.write_alts_into(buf);
    if n_alts == 0 {
        buf.push(b'.');
    }
    buf.push(b'\t');
    match qual {
        Some(q) => {
            // r[impl vcf_writer.float_precision]
            write_float_g(buf, q).map_err(|source| VcfError::FailedToWriteFormattedString {
                field: SmolStr::from("QUAL"),
                source,
            })?;
        }
        // r[impl vcf_writer.missing_dot]
        None => buf.push(b'.'),
    }
    buf.push(b'\t');

    Ok(VcfEncoderFields {
        buf,
        fmt_keys,
        sample_bufs,
        output,
        n_allele,
        n_alt,
        info_count: 0,
        tid,
        pos_0based,
        rlen,
        n_samples,
        filter_written: false,
    })
}

// ── RecordEncoder ──────────────────────────────────────────────────────

// r[impl record_encoder.typestate]
// r[impl record_encoder.typestate_must_use]
#[must_use = "record is silently discarded if emit() is not called"]
pub struct RecordEncoder<'a, S> {
    // r[impl record_encoder.typestate_w_erased]
    inner: EncoderInner<'a>,
    _state: PhantomData<S>,
}

// r[impl record_encoder.typestate_w_erased]
enum EncoderInner<'a> {
    Bcf(BcfRecordEncoder<'a>),
    Vcf(VcfEncoderFields<'a>),
}

// ── Begun → Filtered ───────────────────────────────────────────────────

// r[impl record_encoder.typestate_transitions]
// r[impl record_encoder.filters]
impl<'a> RecordEncoder<'a, Begun> {
    pub fn filter_pass(mut self) -> RecordEncoder<'a, Filtered> {
        match &mut self.inner {
            // r[impl bcf_writer.filter_pass]
            EncoderInner::Bcf(enc) => encode_typed_int_vec(enc.shared_buf, &[0i32]),
            EncoderInner::Vcf(vcf) => {
                vcf.filter_written = true;
                vcf.buf.extend_from_slice(b"PASS");
                vcf.buf.push(b'\t');
            }
        }
        RecordEncoder { inner: self.inner, _state: PhantomData }
    }

    pub fn filter_fail(mut self, filters: &[&FilterId]) -> RecordEncoder<'a, Filtered> {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let indices: SmallVec<i32, 3> = filters
                    .iter()
                    .map(|f| {
                        i32::try_from(f.dict_idx())
                            .expect("filter dict_idx fits i32 (validated at registration)")
                    })
                    .collect();
                encode_typed_int_vec(enc.shared_buf, &indices);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.filter_written = true;
                for (i, f) in filters.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b';');
                    }
                    vcf.buf.extend_from_slice(f.name().as_bytes());
                }
                vcf.buf.push(b'\t');
            }
        }
        RecordEncoder { inner: self.inner, _state: PhantomData }
    }

    pub fn no_filter(mut self) -> RecordEncoder<'a, Filtered> {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => encode_type_byte(enc.shared_buf, 0, BCF_BT_NULL),
            EncoderInner::Vcf(vcf) => {
                vcf.filter_written = true;
                vcf.buf.push(b'.');
                vcf.buf.push(b'\t');
            }
        }
        RecordEncoder { inner: self.inner, _state: PhantomData }
    }
}

// ── Filtered: InfoEncoder ──────────────────────────────────────────────

// r[impl record_encoder.info_encoder]
// r[impl record_encoder.bcf_encoding]
// r[impl record_encoder.vcf_encoding]
impl InfoEncoder for RecordEncoder<'_, Filtered> {
    fn info_int(&mut self, id: &FieldId, value: i32) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1); // saturates at u16::MAX; real VCFs have <100 INFO fields
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                vcf.buf.extend_from_slice(b.format(value).as_bytes());
            }
        }
    }
    fn info_float(&mut self, id: &FieldId, value: f32) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.float_precision]
                write_float_g(vcf.buf, value)
                    .expect("f32 with 6 significant digits never exceeds 32 chars");
            }
        }
    }
    fn info_ints(&mut self, id: &FieldId, values: &[i32]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b',');
                    }
                    vcf.buf.extend_from_slice(b.format(*v).as_bytes());
                }
            }
        }
    }
    fn info_floats(&mut self, id: &FieldId, values: &[f32]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b',');
                    }
                    // r[impl vcf_writer.float_precision]
                    write_float_g(vcf.buf, *v)
                        .expect("f32 with 6 significant digits never exceeds 32 chars");
                }
            }
        }
    }
    fn info_flag(&mut self, id: &FieldId) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.flag_encoding]
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 0, BCF_BT_NULL);
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
            }
        }
    }
    fn info_string(&mut self, id: &FieldId, value: &str) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_typed_string(enc.shared_buf, value.as_bytes());
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.percent_encoding]
                percent_encode_into(vcf.buf, value.as_bytes());
            }
        }
    }
    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                // r[impl bcf_writer.missing_sentinels]
                let typ = smallest_int_type_iter(values.iter().filter_map(|v| *v));
                encode_type_byte(enc.shared_buf, values.len(), typ);
                for &v in values {
                    encode_int_value_or_missing(enc.shared_buf, v, typ);
                }
                // r[impl bcf_encoder.info_counting]
                enc.n_info = enc.n_info.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                vcf_info_separator(vcf);
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                let mut b = itoa::Buffer::new();
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b',');
                    }
                    match v {
                        // r[impl vcf_writer.integer_format]
                        Some(n) => vcf.buf.extend_from_slice(b.format(*n).as_bytes()),
                        // r[impl vcf_writer.missing_dot]
                        None => vcf.buf.push(b'.'),
                    }
                }
            }
        }
    }
    // r[impl record_encoder.info_state_queries]
    fn n_allele(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_allele as usize,
            EncoderInner::Vcf(v) => v.n_allele as usize,
        }
    }
    fn n_alt(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_alt as usize,
            EncoderInner::Vcf(v) => v.n_alt as usize,
        }
    }
}

// r[impl record_encoder.typestate_transitions]
impl<'a> RecordEncoder<'a, Filtered> {
    /// Transition to the `WithSamples` state for encoding FORMAT fields.
    ///
    /// The sample count is derived from the header (set during `write_header`).
    pub fn begin_samples(mut self) -> RecordEncoder<'a, WithSamples> {
        match &mut self.inner {
            EncoderInner::Bcf(_enc) => {
                // n_sample already set from header in begin_record
            }
            EncoderInner::Vcf(vcf) => {
                let n = vcf.n_samples as usize;
                // Grow sample_bufs if needed, reuse existing Vec allocations
                vcf.sample_bufs.resize_with(n, Vec::new);
                // Truncate if header has fewer samples than a previous use
                vcf.sample_bufs.truncate(n);
                for buf in vcf.sample_bufs.iter_mut() {
                    buf.clear();
                }
            }
        }
        RecordEncoder { inner: self.inner, _state: PhantomData }
    }
    // r[impl record_encoder.emit]
    // r[impl record_encoder.emit_no_samples]
    pub fn emit(self) -> Result<(), VcfError> {
        emit_inner(self.inner)
    }
}

// ── WithSamples: FormatEncoder ─────────────────────────────────────────

// r[impl record_encoder.format_encoder]
#[allow(clippy::indexing_slicing, reason = "sample_bufs length validated by debug_assert")]
impl FormatEncoder for RecordEncoder<'_, WithSamples> {
    fn format_gt(&mut self, id: &FieldId, gts: &[Genotype]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.gt_encoding]
                // r[impl bcf_writer.indiv_field_major]
                // r[impl bcf_encoder.format_field_major]
                debug_assert_eq!(gts.len(), enc.n_sample as usize);
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                // Ploidy is taken from the first sample. Mixed ploidy (e.g., haploid +
                // diploid on chrX) is not yet supported — BCF allows per-sample padding
                // with end-of-vector sentinels but we don't encode that. For now, callers
                // must pass uniform ploidy; mismatched samples will produce incorrect BCF.
                let ploidy = gts.first().map_or(0, |g| g.alleles.len());
                let max_allele: i32 = gts
                    .iter()
                    .flat_map(|g| g.alleles.iter())
                    .flatten()
                    .map(|idx| i32::from(*idx))
                    .max()
                    .unwrap_or(0);
                let max_val = (max_allele.saturating_add(1)).saturating_mul(2).saturating_add(1);
                let typ = smallest_int_type(&[max_val]);
                encode_type_byte(enc.indiv_buf, ploidy, typ);
                for gt in gts {
                    for (i, allele_opt) in gt.alleles.iter().enumerate() {
                        let encoded: i32 = match allele_opt {
                            Some(idx) => {
                                let phased = if i == 0 {
                                    false
                                } else {
                                    gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                                };
                                (i32::from(*idx).saturating_add(1) << 1) | i32::from(phased)
                            }
                            None => 0,
                        };
                        encode_int_as(enc.indiv_buf, encoded, typ);
                    }
                }
                enc.n_fmt = enc.n_fmt.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                // r[impl vcf_writer.genotype_serialization]
                debug_assert_eq!(gts.len(), vcf.n_samples as usize);
                vcf_begin_format_field(vcf, id);
                let mut b = itoa::Buffer::new();
                for (si, gt) in gts.iter().enumerate() {
                    for (i, allele) in gt.alleles.iter().enumerate() {
                        if i > 0 {
                            let phased =
                                gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false);
                            vcf.sample_bufs[si].push(if phased { b'|' } else { b'/' });
                        }
                        match allele {
                            Some(idx) => {
                                vcf.sample_bufs[si].extend_from_slice(b.format(*idx).as_bytes())
                            }
                            None => vcf.sample_bufs[si].push(b'.'),
                        }
                    }
                }
            }
        }
    }
    fn format_int(&mut self, id: &FieldId, values: &[i32]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.indiv_field_major]
                debug_assert_eq!(values.len(), enc.n_sample as usize);
                let tc = smallest_int_type(values);
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    v.encode_bcf_as(enc.indiv_buf, tc);
                }
                enc.n_fmt = enc.n_fmt.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                debug_assert_eq!(values.len(), vcf.n_samples as usize);
                vcf_begin_format_field(vcf, id);
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                for (i, &v) in values.iter().enumerate() {
                    vcf.sample_bufs[i].extend_from_slice(b.format(v).as_bytes());
                }
            }
        }
    }
    fn format_float(&mut self, id: &FieldId, values: &[f32]) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.indiv_field_major]
                debug_assert_eq!(values.len(), enc.n_sample as usize);
                // f32::scalar_type_code() always returns BCF_BT_FLOAT regardless of value,
                // so checking only the first element is harmless (unlike format_int which
                // must scan all values to find the smallest fitting integer type).
                let tc = values.first().map_or(BCF_BT_FLOAT, |v| v.scalar_type_code());
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    v.encode_bcf_as(enc.indiv_buf, tc);
                }
                enc.n_fmt = enc.n_fmt.saturating_add(1);
            }
            EncoderInner::Vcf(vcf) => {
                debug_assert_eq!(values.len(), vcf.n_samples as usize);
                vcf_begin_format_field(vcf, id);
                // r[impl vcf_writer.float_precision]
                for (i, &v) in values.iter().enumerate() {
                    write_float_g(&mut vcf.sample_bufs[i], v)
                        .expect("f32 with 6 significant digits never exceeds 32 chars");
                }
            }
        }
    }
    // r[impl record_encoder.format_state_queries]
    fn n_allele(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_allele as usize,
            EncoderInner::Vcf(v) => v.n_allele as usize,
        }
    }
    fn n_alt(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_alt as usize,
            EncoderInner::Vcf(v) => v.n_alt as usize,
        }
    }
    fn n_samples(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_sample as usize,
            EncoderInner::Vcf(v) => v.n_samples as usize,
        }
    }
}

impl<'a> RecordEncoder<'a, WithSamples> {
    // r[impl record_encoder.emit]
    pub fn emit(self) -> Result<(), VcfError> {
        emit_inner(self.inner)
    }
}

// ── Emit ───────────────────────────────────────────────────────────────

fn emit_inner(mut inner: EncoderInner<'_>) -> Result<(), VcfError> {
    match &mut inner {
        EncoderInner::Bcf(enc) => enc.emit_inner(),
        EncoderInner::Vcf(vcf) => vcf_emit(vcf),
    }
}

fn vcf_emit(vcf: &mut VcfEncoderFields<'_>) -> Result<(), VcfError> {
    if !vcf.fmt_keys.is_empty() && vcf.n_samples == 0 {
        return Err(VcfError::FormatDataWithoutSamples);
    }
    if !vcf.filter_written {
        vcf.buf.extend_from_slice(b".\t");
    }
    // r[impl vcf_writer.missing_dot]
    if vcf.info_count == 0 {
        vcf.buf.push(b'.');
    }
    // r[impl vcf_writer.format_serialization]
    if !vcf.fmt_keys.is_empty() {
        vcf.buf.push(b'\t');
        for (i, key) in vcf.fmt_keys.iter().enumerate() {
            if i > 0 {
                vcf.buf.push(b':');
            }
            vcf.buf.extend_from_slice(key.as_bytes());
        }
        for sample_buf in vcf.sample_bufs.iter() {
            vcf.buf.push(b'\t');
            vcf.buf.extend_from_slice(sample_buf);
        }
    }
    vcf.buf.push(b'\n');
    let beg = u64::from(vcf.pos_0based);
    let end = beg.saturating_add(u64::from(vcf.rlen));
    vcf.output.write_line(vcf.buf)?;
    vcf.output.push_index(vcf.tid, beg, end)?;
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────

// r[impl vcf_writer.info_serialization]
fn vcf_info_separator(vcf: &mut VcfEncoderFields<'_>) {
    if !vcf.filter_written {
        vcf.filter_written = true;
        vcf.buf.extend_from_slice(b".\t");
    }
    if vcf.info_count > 0 {
        vcf.buf.push(b';');
    }
    vcf.info_count = vcf.info_count.saturating_add(1);
}

// r[impl vcf_writer.format_serialization]
fn vcf_begin_format_field(vcf: &mut VcfEncoderFields<'_>, id: &FieldId) {
    if !vcf.fmt_keys.is_empty() {
        for buf in vcf.sample_bufs.iter_mut() {
            buf.push(b':');
        }
    }
    vcf.fmt_keys.push(id.name().into());
}

fn encode_array_values<T: BcfValue>(buf: &mut Vec<u8>, values: &[T]) {
    let mut type_code = 0u8;
    for v in values {
        let tc = v.scalar_type_code();
        if tc > type_code {
            type_code = tc;
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
