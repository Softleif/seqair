//! Unified VCF/BCF writer with typestate enforcement.

use super::OutputFormat;
use super::alleles::Alleles;
use super::bcf_encoding::*;
use super::encoder::{BcfRecordEncoder, BcfValue, BgzfWrite, ContigHandle, FieldTracker};
use super::error::VcfError;
use super::header::VcfHeader;
use super::record::Genotype;
use super::record_encoder::{ContigId, FieldId, FilterId, FormatEncoder, InfoEncoder};
use super::writer::{percent_encode_into, write_float_g};
use crate::io::BgzfWriter;
use crate::io::IndexBuilder;
use seqair_types::{Pos1, SmallVec, SmolStr};
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
///
/// `begin_record` is only available after the header has been written:
///
/// ```compile_fail
/// use seqair::vcf::{OutputFormat, Writer};
///
/// let mut buf = Vec::new();
/// let mut writer = Writer::new(&mut buf, OutputFormat::Vcf);
/// // ERROR: begin_record requires Writer<_, Ready>, not Unstarted
/// writer.begin_record(todo!(), todo!(), todo!(), None).unwrap();
/// ```
///
/// `write_header` cannot be called twice:
///
/// ```compile_fail
/// use seqair::vcf::{OutputFormat, VcfHeader, Writer};
/// use std::sync::Arc;
///
/// let header = Arc::new(VcfHeader::builder().build().unwrap());
/// let mut buf = Vec::new();
/// let writer = Writer::new(&mut buf, OutputFormat::Vcf);
/// let writer = writer.write_header(&header).unwrap();
/// // ERROR: write_header is on Unstarted, not Ready
/// writer.write_header(&header).unwrap();
/// ```
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
        info_tracker: FieldTracker,
        fmt_tracker: FieldTracker,
    },
    VcfGz {
        bgzf: BgzfWriter<W>,
        index: Option<IndexBuilder>,
        buf: Vec<u8>,
        fmt_keys: Vec<SmolStr>,
        sample_bufs: Vec<Vec<u8>>,
        n_samples: u32,
        info_tracker: FieldTracker,
        fmt_tracker: FieldTracker,
    },
    Bcf {
        bgzf: BgzfWriter<W>,
        index: Option<IndexBuilder>,
        shared_buf: Vec<u8>,
        indiv_buf: Vec<u8>,
        n_samples: u32,
        info_tracker: FieldTracker,
        fmt_tracker: FieldTracker,
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
                info_tracker: FieldTracker::default(),
                fmt_tracker: FieldTracker::default(),
            },
            OutputFormat::VcfGz => WriterInner::VcfGz {
                bgzf: BgzfWriter::new(inner),
                index: None,
                buf: Vec::with_capacity(4096),
                fmt_keys: Vec::with_capacity(8),
                sample_bufs: Vec::new(),
                n_samples: 0,
                info_tracker: FieldTracker::default(),
                fmt_tracker: FieldTracker::default(),
            },
            OutputFormat::Bcf => WriterInner::Bcf {
                bgzf: BgzfWriter::new(inner),
                index: None,
                shared_buf: Vec::with_capacity(4096),
                indiv_buf: Vec::with_capacity(4096),
                n_samples: 0,
                info_tracker: FieldTracker::default(),
                fmt_tracker: FieldTracker::default(),
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
        pos: Pos1,
        alleles: &Alleles,
        qual: Option<f32>,
    ) -> Result<RecordEncoder<'_, Begun>, VcfError> {
        let inner = match &mut self.inner {
            WriterInner::Bcf {
                bgzf,
                index,
                shared_buf,
                indiv_buf,
                n_samples,
                info_tracker,
                fmt_tracker,
            } => {
                info_tracker.clear();
                fmt_tracker.clear();
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
                    info_tracker,
                    fmt_tracker,
                };
                alleles.begin_record(&mut enc, ContigHandle(contig.tid()), pos, qual)?;
                EncoderInner::Bcf(enc)
            }
            WriterInner::Vcf {
                output,
                buf,
                fmt_keys,
                sample_bufs,
                n_samples,
                info_tracker,
                fmt_tracker,
            } => {
                info_tracker.clear();
                fmt_tracker.clear();
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
                    info_tracker,
                )?)
            }
            WriterInner::VcfGz {
                bgzf,
                index,
                buf,
                fmt_keys,
                sample_bufs,
                n_samples,
                info_tracker,
                fmt_tracker,
            } => {
                info_tracker.clear();
                fmt_tracker.clear();
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
                    info_tracker,
                )?)
            }
        };

        Ok(RecordEncoder { inner, _state: PhantomData })
    }

    // r[impl record_encoder.finish]
    /// Finalize the writer, returning the inner writer and optional index builder.
    ///
    /// For VCF.gz and BCF, all buffered data is flushed and the BGZF EOF marker
    /// is written before this returns — the caller does not need to flush `W`.
    /// The index builder is `Some` for `VcfGz` (TBI) and `Bcf` (CSI), `None`
    /// for plain VCF. This signature matches [`crate::bam::BamWriter::finish`].
    pub fn finish(self) -> Result<(W, Option<IndexBuilder>), VcfError> {
        match self.inner {
            WriterInner::Vcf { mut output, .. } => {
                output.flush().map_err(VcfError::Io)?;
                Ok((output, None))
            }
            WriterInner::VcfGz { bgzf, mut index, .. } => {
                // r[impl vcf_writer.finish]
                let voff = bgzf.virtual_offset();
                if let Some(ref mut idx) = index {
                    idx.finish(voff)?;
                }
                let inner = bgzf.finish()?;
                Ok((inner, index))
            }
            WriterInner::Bcf { bgzf, mut index, .. } => {
                // r[impl bcf_writer.finish]
                let voff = bgzf.virtual_offset();
                if let Some(ref mut idx) = index {
                    idx.finish(voff)?;
                }
                let inner = bgzf.finish()?;
                Ok((inner, index))
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
    info_tracker: &'a mut FieldTracker,
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
    pos: Pos1,
    alleles: &Alleles,
    qual: Option<f32>,
    n_samples: u32,
    info_tracker: &'a mut FieldTracker,
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
    let pos_0based = *pos.to_zero_based();
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
    buf.extend_from_slice(itoa_buf.format(*pos).as_bytes());
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
            write_float_g(buf, q)
                .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
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
        info_tracker,
    })
}

// ── RecordEncoder ──────────────────────────────────────────────────────

// r[impl record_encoder.typestate]
// r[impl record_encoder.typestate_must_use]
/// Record encoder with typestate enforcement of the
/// `Begun` → `Filtered` → `WithSamples` → `emit()` chain.
///
/// INFO fields cannot be encoded before the filter is set:
///
/// ```compile_fail
/// use seqair::vcf::{Begun, RecordEncoder};
/// use seqair::vcf::record_encoder::InfoEncoder;
///
/// fn check(enc: &mut RecordEncoder<'_, Begun>) {
///     // ERROR: InfoEncoder is not implemented for RecordEncoder<Begun>
///     enc.info_int(todo!(), 42);
/// }
/// ```
///
/// FORMAT fields cannot be encoded before `begin_samples`:
///
/// ```compile_fail
/// use seqair::vcf::{Filtered, RecordEncoder};
/// use seqair::vcf::record_encoder::FormatEncoder;
///
/// fn check(enc: &mut RecordEncoder<'_, Filtered>) {
///     // ERROR: FormatEncoder is not implemented for RecordEncoder<Filtered>
///     enc.format_int(todo!(), &[42]).unwrap();
/// }
/// ```
///
/// INFO fields cannot be encoded after `begin_samples`:
///
/// ```compile_fail
/// use seqair::vcf::{WithSamples, RecordEncoder};
/// use seqair::vcf::record_encoder::InfoEncoder;
///
/// fn check(enc: &mut RecordEncoder<'_, WithSamples>) {
///     // ERROR: InfoEncoder is not implemented for RecordEncoder<WithSamples>
///     enc.info_int(todo!(), 42);
/// }
/// ```
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

// ── Field dedup helpers ───────────────────────────────────────────────

// r[impl record_encoder.info_dedup]
impl VcfEncoderFields<'_> {
    /// Handle duplicate removal and separator, returning the buffer offset for
    /// the caller to start writing key=value bytes.
    fn prepare_info_field(&mut self, id: &FieldId) {
        self.info_tracker.remove_duplicate_vcf(
            self.buf,
            &mut self.info_count,
            id.dict_idx(),
            id.name(),
        );
        if !self.filter_written {
            self.filter_written = true;
            self.buf.extend_from_slice(b".\t");
        }
        let start = self.buf.len();
        if self.info_count > 0 {
            self.buf.push(b';');
        }
        self.info_count = self.info_count.saturating_add(1);
        self.info_tracker.push(id.dict_idx(), start);
    }
}

// r[impl record_encoder.info_dedup]
impl BcfRecordEncoder<'_> {
    /// Handle duplicate removal, track offset. Returns `true` if replacing
    /// (caller should skip incrementing `n_info`).
    fn prepare_info_field(&mut self, id: &FieldId) -> bool {
        let replacing =
            self.info_tracker.remove_duplicate(self.shared_buf, id.dict_idx(), id.name());
        self.info_tracker.push(id.dict_idx(), self.shared_buf.len());
        replacing
    }

    // r[impl record_encoder.format_dedup]
    fn prepare_format_field(&mut self, id: &FieldId) -> bool {
        let replacing = self.fmt_tracker.remove_duplicate(self.indiv_buf, id.dict_idx(), id.name());
        self.fmt_tracker.push(id.dict_idx(), self.indiv_buf.len());
        replacing
    }
}

// r[impl record_encoder.format_dedup]
impl VcfEncoderFields<'_> {
    /// Handle FORMAT duplicate removal. Removes the previous key from `fmt_keys`
    /// and the corresponding colon-delimited data from each sample buffer.
    fn prepare_format_field(&mut self, id: &FieldId) {
        if let Some(idx) = self.fmt_keys.iter().position(|k| k.as_str() == id.name()) {
            tracing::warn!(
                field = id.name(),
                "FORMAT field encoded twice; overwriting previous value"
            );
            let total = self.fmt_keys.len();
            self.fmt_keys.remove(idx);
            for buf in self.sample_bufs.iter_mut() {
                remove_colon_delimited_field(buf, idx, total);
            }
        }
        // Normal begin_format_field logic
        if !self.fmt_keys.is_empty() {
            for buf in self.sample_bufs.iter_mut() {
                buf.push(b':');
            }
        }
        self.fmt_keys.push(id.name().into());
    }
}

/// Remove the `field_idx`-th colon-delimited field from `buf`.
///
/// `total_fields` is the field count _before_ removal. For a buffer like
/// `"GT:DP:GQ"` with `total_fields=3`, removing field 1 yields `"GT:GQ"`.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "field_idx < total_fields guaranteed by caller; \
              +1/-1 are safe because non-first fields always have a preceding ':'"
)]
fn remove_colon_delimited_field(buf: &mut Vec<u8>, field_idx: usize, total_fields: usize) {
    if total_fields == 0 {
        return;
    }
    // Walk to the target field by counting ':' separators.
    let mut current = 0;
    let mut start = 0;
    for i in 0..buf.len() {
        if buf[i] == b':' {
            if current == field_idx {
                // Field is at [start..i].
                drain_delimited(buf, start, i, field_idx, total_fields);
                return;
            }
            current += 1;
            start = i + 1;
        }
    }
    // Target is the last field: [start..buf.len()].
    if current == field_idx {
        drain_delimited(buf, start, buf.len(), field_idx, total_fields);
    }
}

/// Drain bytes for one colon-delimited field including exactly one separator.
#[allow(clippy::arithmetic_side_effects, reason = "see remove_colon_delimited_field")]
fn drain_delimited(
    buf: &mut Vec<u8>,
    field_start: usize,
    field_end: usize,
    field_idx: usize,
    total_fields: usize,
) {
    let range = if total_fields == 1 {
        field_start..field_end
    } else if field_idx == 0 {
        // First of multiple: remove field + following ':'
        field_start..field_end + 1
    } else {
        // Non-first: remove preceding ':' + field
        field_start - 1..field_end
    };
    buf.drain(range);
}

// ── Filtered: InfoEncoder ──────────────────────────────────────────────

// r[impl record_encoder.info_encoder]
// r[impl record_encoder.bcf_encoding]
// r[impl record_encoder.vcf_encoding]
impl InfoEncoder for RecordEncoder<'_, Filtered> {
    fn info_int(&mut self, id: &FieldId, value: i32) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let replacing = enc.prepare_info_field(id);
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
                let replacing = enc.prepare_info_field(id);
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
                let replacing = enc.prepare_info_field(id);
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
                let replacing = enc.prepare_info_field(id);
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
                let replacing = enc.prepare_info_field(id);
                // r[impl bcf_writer.flag_encoding]
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 0, BCF_BT_NULL);
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
                vcf.buf.extend_from_slice(id.name().as_bytes());
            }
        }
    }
    fn info_string(&mut self, id: &FieldId, value: &str) {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let replacing = enc.prepare_info_field(id);
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_typed_string(enc.shared_buf, value.as_bytes());
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
                let replacing = enc.prepare_info_field(id);
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                // r[impl bcf_writer.missing_sentinels]
                let typ = smallest_int_type_iter(values.iter().filter_map(|v| *v));
                encode_type_byte(enc.shared_buf, values.len(), typ);
                for &v in values {
                    encode_int_value_or_missing(enc.shared_buf, v, typ);
                }
                // r[impl bcf_encoder.info_counting]
                if !replacing {
                    enc.n_info = enc.n_info.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_info_field(id);
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
    fn format_gt(&mut self, id: &FieldId, gts: &[Genotype]) -> Result<(), VcfError> {
        // Ploidy is taken from the first sample. Mixed ploidy (e.g., haploid +
        // diploid on chrX) is not yet supported — BCF allows per-sample padding
        // with end-of-vector sentinels but we don't encode that.
        let ploidy = gts.first().map_or(0, |g| g.alleles.len());
        if let Some((i, g)) = gts.iter().enumerate().find(|(_, g)| g.alleles.len() != ploidy) {
            return Err(VcfError::MixedPloidy {
                first_ploidy: ploidy,
                mismatch_index: i,
                mismatch_ploidy: g.alleles.len(),
            });
        }

        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let replacing = enc.prepare_format_field(id);
                // r[impl bcf_writer.gt_encoding]
                // r[impl bcf_writer.indiv_field_major]
                // r[impl bcf_encoder.format_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
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
                if !replacing {
                    enc.n_fmt = enc.n_fmt.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                // r[impl vcf_writer.genotype_serialization]
                vcf.prepare_format_field(id);
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
        Ok(())
    }
    fn format_int(&mut self, id: &FieldId, values: &[i32]) -> Result<(), VcfError> {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let replacing = enc.prepare_format_field(id);
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.indiv_field_major]
                let tc = smallest_int_type(values);
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    v.encode_bcf_as(enc.indiv_buf, tc);
                }
                if !replacing {
                    enc.n_fmt = enc.n_fmt.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_format_field(id);
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                for (i, &v) in values.iter().enumerate() {
                    vcf.sample_bufs[i].extend_from_slice(b.format(v).as_bytes());
                }
            }
        }
        Ok(())
    }
    fn format_float(&mut self, id: &FieldId, values: &[f32]) -> Result<(), VcfError> {
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let replacing = enc.prepare_format_field(id);
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.indiv_field_major]
                // f32::scalar_type_code() always returns BCF_BT_FLOAT regardless of value,
                // so checking only the first element is harmless (unlike format_int which
                // must scan all values to find the smallest fitting integer type).
                let tc = values.first().map_or(BCF_BT_FLOAT, |v| v.scalar_type_code());
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    v.encode_bcf_as(enc.indiv_buf, tc);
                }
                if !replacing {
                    enc.n_fmt = enc.n_fmt.saturating_add(1);
                }
            }
            EncoderInner::Vcf(vcf) => {
                vcf.prepare_format_field(id);
                // r[impl vcf_writer.float_precision]
                for (i, &v) in values.iter().enumerate() {
                    write_float_g(&mut vcf.sample_bufs[i], v)
                        .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
                }
            }
        }
        Ok(())
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
