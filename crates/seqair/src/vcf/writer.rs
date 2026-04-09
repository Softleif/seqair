//! VCF text writer. Serializes [`VcfRecord`]s as tab-delimited VCF lines,
//! optionally BGZF-compressed with TBI index co-production.

use super::error::VcfError;
use super::header::VcfHeader;
use super::index_builder::IndexBuilder;
use super::record::{Filters, InfoValue, SampleValue, VcfRecord};
use crate::bam::bgzf_writer::BgzfWriter;
use seqair_types::SmolStr;
use seqair_types::smol_str::format_smolstr;
use std::io::Write;
use std::sync::Arc;

/// Output mode for the VCF writer.
enum Output<W: Write> {
    Plain(W),
    Bgzf { bgzf: BgzfWriter<W>, index: IndexBuilder },
}

// r[impl vcf_writer.header_first]
// r[impl vcf_writer.tab_delimited]
// r[impl vcf_writer.buffer_reuse]
/// VCF text writer with optional BGZF compression and TBI index co-production.
pub struct VcfWriter<W: Write> {
    output: Output<W>,
    header: Arc<VcfHeader>,
    buf: Vec<u8>,
    header_written: bool,
    /// FORMAT key accumulator (reused across records by record_encoder).
    fmt_keys: Vec<SmolStr>,
    /// FORMAT value buffer (reused across records by record_encoder).
    fmt_values: Vec<u8>,
}

impl<W: Write> VcfWriter<W> {
    /// Create a plain (uncompressed) VCF writer.
    pub fn new(inner: W, header: Arc<VcfHeader>) -> Self {
        Self {
            output: Output::Plain(inner),
            header,
            buf: Vec::with_capacity(4096),
            header_written: false,
            fmt_keys: Vec::with_capacity(8),
            fmt_values: Vec::with_capacity(256),
        }
    }

    /// Create a BGZF-compressed VCF writer with TBI index co-production.
    pub fn bgzf(inner: W, header: Arc<VcfHeader>) -> Self {
        let n_refs = header.contigs().len();
        let bgzf = BgzfWriter::new(inner);
        Self {
            output: Output::Bgzf {
                bgzf,
                index: IndexBuilder::tbi(n_refs, crate::bam::bgzf::VirtualOffset(0)),
            },
            header,
            buf: Vec::with_capacity(4096),
            header_written: false,
            fmt_keys: Vec::with_capacity(8),
            fmt_values: Vec::with_capacity(256),
        }
    }

    // r[impl vcf_writer.header_first]
    /// Write the VCF header. Must be called exactly once before any `write_record`.
    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; VcfWrite trait delegates to it for dyn dispatch"
    )]
    pub fn write_header(&mut self) -> Result<(), VcfError> {
        let text = self.header.to_vcf_text();
        match &mut self.output {
            Output::Plain(w) => w.write_all(text.as_bytes()).map_err(VcfError::Io)?,
            Output::Bgzf { bgzf, .. } => bgzf.write_all(text.as_bytes())?,
        }
        self.header_written = true;

        // Update index with header end offset for BGZF mode
        if let Output::Bgzf { bgzf, index } = &mut self.output {
            *index = IndexBuilder::tbi(self.header.contigs().len(), bgzf.virtual_offset());
        }
        Ok(())
    }

    // r[impl vcf_writer.tab_delimited]
    // r[impl vcf_writer.buffer_reuse]
    // r[impl vcf_writer.validation]
    /// Write a single VCF record.
    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; VcfWrite trait delegates to it for dyn dispatch"
    )]
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        if !self.header_written {
            return Err(VcfError::HeaderNotWritten);
        }

        // r[impl vcf_writer.buffer_reuse]
        self.buf.clear();

        // CHROM
        self.buf.extend_from_slice(record.contig.as_bytes());
        self.buf.push(b'\t');

        // POS (1-based)
        let mut itoa_buf = itoa::Buffer::new();
        self.buf.extend_from_slice(itoa_buf.format(record.pos.get()).as_bytes());
        self.buf.push(b'\t');

        // ID
        match &record.id {
            Some(id) => self.buf.extend_from_slice(id.as_bytes()),
            None => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // REF — zero-alloc, writes directly into buffer
        record.alleles.write_ref_into(&mut self.buf);
        self.buf.push(b'\t');

        // ALT — zero-alloc, writes directly into buffer
        let n_alts = record.alleles.write_alts_into(&mut self.buf);
        if n_alts == 0 {
            self.buf.push(b'.');
        }
        self.buf.push(b'\t');

        // QUAL
        // r[impl vcf_writer.float_precision]
        match record.qual {
            Some(q) => write_float_g(&mut self.buf, q).map_err(|source| {
                VcfError::FailedToWriteFormattedString {
                    field: SmolStr::new_static("qual"),
                    source,
                }
            })?,
            None => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // FILTER
        // r[impl vcf_writer.missing_dot]
        match &record.filters {
            Filters::Pass => self.buf.extend_from_slice(b"PASS"),
            Filters::Failed(ids) => {
                for (i, id) in ids.iter().enumerate() {
                    if i > 0 {
                        self.buf.push(b';');
                    }
                    self.buf.extend_from_slice(id.as_bytes());
                }
            }
            Filters::NotApplied => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        // INFO
        // r[impl vcf_writer.info_serialization]
        if record.info.is_empty() {
            self.buf.push(b'.');
        } else {
            for (i, (key, value)) in record.info.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b';');
                }
                self.buf.extend_from_slice(key.as_bytes());
                match value {
                    InfoValue::Flag => {} // key only, no =
                    _ => {
                        self.buf.push(b'=');
                        write_info_value(&mut self.buf, value).map_err(|source| {
                            VcfError::FailedToWriteFormattedString {
                                field: format_smolstr!("INFO:{key}"),
                                source,
                            }
                        })?;
                    }
                }
            }
        }

        // FORMAT + samples
        // r[impl vcf_writer.format_serialization]
        if !record.samples.format_keys.is_empty() && !record.samples.values.is_empty() {
            self.buf.push(b'\t');
            for (i, key) in record.samples.format_keys.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b':');
                }
                self.buf.extend_from_slice(key.as_bytes());
            }

            for sample_values in &record.samples.values {
                self.buf.push(b'\t');
                for (i, val) in sample_values.iter().enumerate() {
                    if i > 0 {
                        self.buf.push(b':');
                    }
                    write_sample_value(&mut self.buf, val).map_err(|source| {
                        VcfError::FailedToWriteFormattedString {
                            field: format_smolstr!(
                                "FORMAT:{}",
                                record
                                    .samples
                                    .format_keys
                                    .get(i)
                                    .unwrap_or(&SmolStr::new_static("?"))
                            ),
                            source,
                        }
                    })?;
                }
            }
        }

        self.buf.push(b'\n');

        // Write the serialized record to the output
        match &mut self.output {
            Output::Plain(w) => {
                w.write_all(&self.buf).map_err(VcfError::Io)?;
            }
            Output::Bgzf { bgzf, index } => {
                bgzf.flush_if_needed(self.buf.len())?;
                bgzf.write_all(&self.buf)?;

                let tid_usize = self.header.contig_id(&record.contig)?;
                let tid = i32::try_from(tid_usize).map_err(|_| VcfError::ValueOverflow {
                    field: "contig_tid",
                    value: tid_usize as u64,
                    target_type: "i32",
                })?;
                let beg = u64::from(record.pos.to_zero_based().get());
                let end = beg.saturating_add(record.alleles.rlen() as u64);
                index.push(tid, beg, end, bgzf.virtual_offset())?;
            }
        }

        Ok(())
    }

    // r[impl vcf_writer.finish]
    /// Flush and finalize the writer.
    /// For BGZF: writes EOF marker. Returns the index builder if present.
    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; VcfWrite trait delegates to it for dyn dispatch"
    )]
    pub fn finish(self) -> Result<Option<IndexBuilder>, VcfError> {
        match self.output {
            Output::Plain(mut w) => {
                w.flush().map_err(VcfError::Io)?;
                Ok(None)
            }
            Output::Bgzf { bgzf, mut index } => {
                let voff = bgzf.virtual_offset();
                index.finish(voff)?;
                bgzf.finish()?;
                Ok(Some(index))
            }
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum WriteError {
    #[error("failed to write formatted string")]
    WriteFormattedString { source: std::io::Error },
    #[error("Failed to strip trailing zeros from formatted float")]
    FailedToStripTrailingZeros,
    #[error("formatted float exceeds 32 bytes")]
    FormattedFloatLongerThan32Chars,
}

type FmtResult = Result<(), WriteError>;

// r[impl vcf_writer.integer_format]
fn write_info_value(buf: &mut Vec<u8>, value: &InfoValue) -> FmtResult {
    let mut itoa_buf = itoa::Buffer::new();
    match value {
        InfoValue::Integer(v) => buf.extend_from_slice(itoa_buf.format(*v).as_bytes()),
        InfoValue::Float(v) => write_float_g(buf, *v)?,
        InfoValue::Flag => {} // handled at call site
        InfoValue::String(s) => percent_encode_into(buf, s.as_bytes()),
        InfoValue::IntegerArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(n) => buf.extend_from_slice(itoa_buf.format(*n).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        InfoValue::FloatArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(f) => write_float_g(buf, *f)?,
                    None => buf.push(b'.'),
                }
            }
        }
        InfoValue::StringArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(s) => percent_encode_into(buf, s.as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
    }
    Ok(())
}

// r[impl vcf_writer.float_precision]
/// Format a float like C's `%g` with 6 significant digits — no trailing zeros,
/// no trailing decimal point. This matches htslib/bcftools VCF text output.
///
/// Examples: 35.89775 → "35.8978", 60.0 → "60", 0.777778 → "0.777778"
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "magnitude is log10 of f64 (range -308..=308), precision is a literal 6, cursor position is bounded by the 32-byte tmp buffer; all casts are safe"
)]
fn write_float_g(buf: &mut Vec<u8>, v: f32) -> Result<(), WriteError> {
    // %g with precision P means P significant digits.
    // In fixed notation: decimal_places = P - floor(log10(|v|)) - 1
    // For very small or very large values, %g switches to scientific notation,
    // but VCF floats are typically in a reasonable range.
    let v = f64::from(v);
    let precision = 6usize; // %g default

    if v == 0.0 {
        buf.push(b'0');
        return Ok(());
    }

    let magnitude = v.abs().log10().floor() as i32;
    let decimal_places = (precision as i32).saturating_sub(magnitude).saturating_sub(1);
    let decimal_places = decimal_places.max(0) as usize;

    let mut tmp = [0u8; 32];
    let len = {
        use std::io::Write as _;
        let mut cursor = std::io::Cursor::new(&mut tmp[..]);
        write!(cursor, "{v:.decimal_places$}")
            .map_err(|source| WriteError::WriteFormattedString { source })?;
        cursor.position() as usize
    };
    let s = &tmp.get(..len).ok_or(WriteError::FormattedFloatLongerThan32Chars)?;
    // Strip trailing zeros after decimal point, then trailing decimal point
    if let Some(dot) = s.iter().position(|&b| b == b'.') {
        let mut end = len;
        while end > dot.saturating_add(1) && s.get(end.saturating_sub(1)) == Some(&b'0') {
            end = end.saturating_sub(1);
        }
        if end == dot.saturating_add(1) {
            end = dot;
        }
        buf.extend_from_slice(s.get(..end).ok_or(WriteError::FailedToStripTrailingZeros)?);
    } else {
        buf.extend_from_slice(s);
    }
    Ok(())
}

// r[impl vcf_writer.percent_encoding]
/// Percent-encode special characters in VCF field values.
fn percent_encode_into(buf: &mut Vec<u8>, data: &[u8]) {
    for &b in data {
        match b {
            b':' => buf.extend_from_slice(b"%3A"),
            b';' => buf.extend_from_slice(b"%3B"),
            b'=' => buf.extend_from_slice(b"%3D"),
            b'%' => buf.extend_from_slice(b"%25"),
            b',' => buf.extend_from_slice(b"%2C"),
            b'\t' => buf.extend_from_slice(b"%09"),
            b'\n' => buf.extend_from_slice(b"%0A"),
            b'\r' => buf.extend_from_slice(b"%0D"),
            _ => buf.push(b),
        }
    }
}

// r[impl vcf_writer.genotype_serialization]
fn write_sample_value(buf: &mut Vec<u8>, value: &SampleValue) -> FmtResult {
    let mut itoa_buf = itoa::Buffer::new();
    match value {
        SampleValue::Missing => buf.push(b'.'),
        SampleValue::Integer(v) => buf.extend_from_slice(itoa_buf.format(*v).as_bytes()),
        SampleValue::Float(v) => write_float_g(buf, *v)?,
        SampleValue::String(s) => percent_encode_into(buf, s.as_bytes()),
        SampleValue::Genotype(gt) => {
            for (i, allele) in gt.alleles.iter().enumerate() {
                if i > 0 {
                    // Separator: phased[i-1] determines | or /
                    let phased = gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false);
                    buf.push(if phased { b'|' } else { b'/' });
                }
                match allele {
                    Some(idx) => buf.extend_from_slice(itoa_buf.format(*idx).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        SampleValue::IntegerArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(n) => buf.extend_from_slice(itoa_buf.format(*n).as_bytes()),
                    None => buf.push(b'.'),
                }
            }
        }
        SampleValue::FloatArray(arr) => {
            for (i, v) in arr.iter().enumerate() {
                if i > 0 {
                    buf.push(b',');
                }
                match v {
                    Some(f) => write_float_g(buf, *f)?,
                    None => buf.push(b'.'),
                }
            }
        }
    }
    Ok(())
}

// ── VcfOutput trait ────────────────────────────────────────────────────

// r[impl record_encoder.vcf_text_output]
/// Abstracts over plain and BGZF VCF output for the record encoder.
pub(crate) trait VcfOutput {
    fn write_line(&mut self, buf: &[u8]) -> Result<(), VcfError>;
    fn push_index(&mut self, tid: i32, beg: u64, end: u64) -> Result<(), VcfError>;
}

struct PlainOutput<'a, W: Write>(&'a mut W);

impl<W: Write> VcfOutput for PlainOutput<'_, W> {
    fn write_line(&mut self, buf: &[u8]) -> Result<(), VcfError> {
        self.0.write_all(buf).map_err(VcfError::Io)
    }
    fn push_index(&mut self, _tid: i32, _beg: u64, _end: u64) -> Result<(), VcfError> {
        Ok(())
    }
}

struct BgzfOutput<'a, W: Write> {
    bgzf: &'a mut BgzfWriter<W>,
    index: &'a mut IndexBuilder,
}

impl<W: Write> VcfOutput for BgzfOutput<'_, W> {
    fn write_line(&mut self, buf: &[u8]) -> Result<(), VcfError> {
        self.bgzf.flush_if_needed(buf.len())?;
        self.bgzf.write_all(buf)?;
        Ok(())
    }
    fn push_index(&mut self, tid: i32, beg: u64, end: u64) -> Result<(), VcfError> {
        self.index.push(tid, beg, end, self.bgzf.virtual_offset())?;
        Ok(())
    }
}

// ── VcfRecordEncoder ───────────────────────────────────────────────────

use super::record_encoder::{ContigId, FieldId, FilterId, RecordEncoder};

// r[impl record_encoder.vcf_text_encoder]
// r[impl record_encoder.vcf_text_buffer_reuse]
/// VCF text record encoder. Borrows reusable buffers from the [`VcfWriter`]
/// to avoid per-record allocation after warmup.
///
/// # Calling contract
///
/// Methods must be called in order: `begin` → `filter_pass`/`filter_fail` →
/// `info_*` → `begin_samples` → `format_*` → `emit`. Debug builds enforce this
/// with assertions. For compile-time enforcement consider a typestate wrapper.
pub struct VcfRecordEncoder<'a> {
    /// Main line buffer — contains CHROM through INFO.
    buf: &'a mut Vec<u8>,
    /// FORMAT key accumulator (reused across records).
    fmt_keys: &'a mut Vec<SmolStr>,
    /// Formatted sample values buffer (reused across records).
    fmt_values: &'a mut Vec<u8>,
    /// Output sink (plain or BGZF).
    output: Box<dyn VcfOutput + 'a>,
    /// Record state.
    n_allele: u16,
    n_alt: u16,
    info_count: u16,
    tid: i32,
    pos_0based: u32,
    rlen: u32,
    /// Number of samples declared via `begin_samples`. 0 means not yet called.
    n_samples: u32,
    /// Whether `begin()` has been called for the current record.
    #[cfg(debug_assertions)]
    record_begun: bool,
}

// r[impl record_encoder.vcf_text_begin]
// r[impl record_encoder.vcf_text_info]
// r[impl record_encoder.vcf_text_format_accumulation]
// r[impl record_encoder.vcf_text_emit]
impl RecordEncoder for VcfRecordEncoder<'_> {
    fn begin(
        &mut self,
        contig: &ContigId,
        pos: seqair_types::Pos<seqair_types::One>,
        alleles: &super::alleles::Alleles,
        qual: Option<f32>,
    ) -> Result<(), VcfError> {
        self.buf.clear();
        self.fmt_keys.clear();
        self.fmt_values.clear();
        self.info_count = 0;
        self.n_samples = 0;
        #[cfg(debug_assertions)]
        {
            self.record_begun = true;
        }

        self.n_allele = u16::try_from(alleles.n_allele()).map_err(|_| VcfError::ValueOverflow {
            field: "n_allele",
            value: alleles.n_allele() as u64,
            target_type: "u16",
        })?;
        self.n_alt = self.n_allele.saturating_sub(1);
        self.tid = i32::try_from(contig.tid()).map_err(|_| VcfError::ValueOverflow {
            field: "contig_tid",
            value: u64::from(contig.tid()),
            target_type: "i32",
        })?;
        self.pos_0based = pos.to_zero_based().get();
        self.rlen = u32::try_from(alleles.rlen()).map_err(|_| VcfError::ValueOverflow {
            field: "rlen",
            value: alleles.rlen() as u64,
            target_type: "u32",
        })?;

        // CHROM
        self.buf.extend_from_slice(contig.name().as_bytes());
        self.buf.push(b'\t');

        // POS (1-based)
        let mut itoa_buf = itoa::Buffer::new();
        self.buf.extend_from_slice(itoa_buf.format(pos.get()).as_bytes());
        self.buf.push(b'\t');

        // ID
        self.buf.push(b'.');
        self.buf.push(b'\t');

        // REF
        alleles.write_ref_into(self.buf);
        self.buf.push(b'\t');

        // ALT
        let n_alts = alleles.write_alts_into(self.buf);
        if n_alts == 0 {
            self.buf.push(b'.');
        }
        self.buf.push(b'\t');

        // QUAL
        match qual {
            Some(q) => {
                write_float_g(self.buf, q).expect("writing to Vec<u8> is infallible");
            }
            None => self.buf.push(b'.'),
        }
        self.buf.push(b'\t');

        Ok(())
    }

    fn filter_pass(&mut self) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "filter_pass() called before begin()");
        self.buf.extend_from_slice(b"PASS");
        self.buf.push(b'\t');
    }

    fn filter_fail(&mut self, filters: &[&FilterId]) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "filter_fail() called before begin()");
        for (i, f) in filters.iter().enumerate() {
            if i > 0 {
                self.buf.push(b';');
            }
            self.buf.extend_from_slice(f.name().as_bytes());
        }
        self.buf.push(b'\t');
    }

    fn info_int(&mut self, id: &FieldId, value: i32) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_int() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        let mut itoa_buf = itoa::Buffer::new();
        self.buf.extend_from_slice(itoa_buf.format(value).as_bytes());
    }

    fn info_float(&mut self, id: &FieldId, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_float() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        write_float_g(self.buf, value).expect("writing to Vec<u8> is infallible");
    }

    fn info_ints(&mut self, id: &FieldId, values: &[i32]) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_ints() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        let mut itoa_buf = itoa::Buffer::new();
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                self.buf.push(b',');
            }
            self.buf.extend_from_slice(itoa_buf.format(*v).as_bytes());
        }
    }

    fn info_floats(&mut self, id: &FieldId, values: &[f32]) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_floats() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                self.buf.push(b',');
            }
            write_float_g(self.buf, *v).expect("writing to Vec<u8> is infallible");
        }
    }

    fn info_flag(&mut self, id: &FieldId) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_flag() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
    }

    fn info_string(&mut self, id: &FieldId, value: &str) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_string() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        percent_encode_into(self.buf, value.as_bytes());
    }

    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "info_int_opts() called before begin()");
        self.info_separator();
        self.buf.extend_from_slice(id.name().as_bytes());
        self.buf.push(b'=');
        let mut itoa_buf = itoa::Buffer::new();
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                self.buf.push(b',');
            }
            match v {
                Some(n) => self.buf.extend_from_slice(itoa_buf.format(*n).as_bytes()),
                None => self.buf.push(b'.'),
            }
        }
    }

    fn begin_samples(&mut self, n: u32) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "begin_samples() called before begin()");
        self.n_samples = n;
    }

    fn format_gt(&mut self, id: &FieldId, gt: &super::record::Genotype) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "format_gt() called before begin()");
        self.push_format_key(id);
        let mut itoa_buf = itoa::Buffer::new();
        for (i, allele) in gt.alleles.iter().enumerate() {
            if i > 0 {
                let phased = gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false);
                self.fmt_values.push(if phased { b'|' } else { b'/' });
            }
            match allele {
                Some(idx) => {
                    self.fmt_values.extend_from_slice(itoa_buf.format(*idx).as_bytes());
                }
                None => self.fmt_values.push(b'.'),
            }
        }
    }

    fn format_int(&mut self, id: &FieldId, value: i32) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "format_int() called before begin()");
        self.push_format_key(id);
        let mut itoa_buf = itoa::Buffer::new();
        self.fmt_values.extend_from_slice(itoa_buf.format(value).as_bytes());
    }

    fn format_float(&mut self, id: &FieldId, value: f32) {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "format_float() called before begin()");
        self.push_format_key(id);
        write_float_g(&mut self.fmt_values, value).expect("writing to Vec<u8> is infallible");
    }

    fn n_allele(&self) -> usize {
        self.n_allele as usize
    }

    fn n_alt(&self) -> usize {
        self.n_alt as usize
    }

    fn emit(&mut self) -> Result<(), VcfError> {
        #[cfg(debug_assertions)]
        debug_assert!(self.record_begun, "emit() called before begin()");

        if !self.fmt_keys.is_empty() && self.n_samples == 0 {
            return Err(VcfError::FormatDataWithoutSamples);
        }

        // If no INFO fields were written, write "."
        if self.info_count == 0 {
            self.buf.push(b'.');
        }

        // FORMAT + samples
        if !self.fmt_keys.is_empty() {
            // FORMAT column: key1:key2:key3
            self.buf.push(b'\t');
            for (i, key) in self.fmt_keys.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b':');
                }
                self.buf.extend_from_slice(key.as_bytes());
            }

            // Sample values column
            self.buf.push(b'\t');
            self.buf.extend_from_slice(&self.fmt_values);
        }

        self.buf.push(b'\n');

        // Write to output
        let beg = u64::from(self.pos_0based);
        let end = beg.saturating_add(u64::from(self.rlen));
        self.output.write_line(self.buf)?;
        self.output.push_index(self.tid, beg, end)?;

        Ok(())
    }
}

impl VcfRecordEncoder<'_> {
    /// Write INFO field separator (`;` between fields, nothing before the first).
    fn info_separator(&mut self) {
        if self.info_count > 0 {
            self.buf.push(b';');
        }
        self.info_count = self.info_count.saturating_add(1);
    }

    /// Push a FORMAT key and add `:` separator to values if not the first field.
    fn push_format_key(&mut self, id: &FieldId) {
        if !self.fmt_keys.is_empty() {
            self.fmt_values.push(b':');
        }
        self.fmt_keys.push(id.name.clone());
    }
}

// ── VcfWriter::record_encoder() ────────────────────────────────────────

impl<W: Write> VcfWriter<W> {
    /// Get a direct record encoder for zero-alloc VCF text encoding.
    pub fn record_encoder(&mut self) -> VcfRecordEncoder<'_> {
        let output: Box<dyn VcfOutput + '_> = match &mut self.output {
            Output::Plain(w) => Box::new(PlainOutput(w)),
            Output::Bgzf { bgzf, index } => Box::new(BgzfOutput { bgzf, index }),
        };

        VcfRecordEncoder {
            buf: &mut self.buf,
            fmt_keys: &mut self.fmt_keys,
            fmt_values: &mut self.fmt_values,
            output,
            n_allele: 0,
            n_alt: 0,
            info_count: 0,
            tid: 0,
            pos_0based: 0,
            rlen: 0,
            n_samples: 0,
            #[cfg(debug_assertions)]
            record_begun: false,
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

    /// Write records and return the output as a String.
    fn write_vcf(header: &Arc<VcfHeader>, records: &[VcfRecord]) -> String {
        let mut output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut output, header.clone());
            writer.write_header().unwrap();
            for record in records {
                writer.write_record(record).unwrap();
            }
            writer.finish().unwrap();
        }
        String::from_utf8(output).unwrap()
    }

    // r[verify vcf_writer.header_first]
    #[test]
    fn write_header_before_records() {
        let header = test_header();
        let text = write_vcf(&header, &[]);
        assert!(text.starts_with("##fileformat=VCFv4.3\n"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    }

    /// Get the last data line fields from VCF text output.
    fn last_line_fields(text: &str) -> Vec<String> {
        text.lines().last().unwrap().split('\t').map(String::from).collect()
    }

    // r[verify vcf_writer.tab_delimited]
    // r[verify vcf_writer.missing_dot]
    #[test]
    fn write_simple_snv() {
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

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[0], "chr1");
        assert_eq!(f[1], "100");
        assert_eq!(f[2], ".");
        assert_eq!(f[3], "A");
        assert_eq!(f[4], "T");
        assert_eq!(f[5], "30"); // %g format: 30.0 → "30"
        assert_eq!(f[6], "PASS");
        assert_eq!(f[7], "DP=50");
    }

    // r[verify vcf_writer.info_serialization]
    #[test]
    fn write_info_flag() {
        let header = test_header();
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_flag("DB")
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[7], "DB");
    }

    // r[verify vcf_writer.genotype_serialization]
    #[test]
    fn write_genotype() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .format_keys(&["GT", "DP"])
        .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1)), SampleValue::Integer(30)])
        .build(&header)
        .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[8], "GT:DP");
        assert_eq!(f[9], "0/1:30");
    }

    // r[verify vcf_writer.genotype_serialization]
    #[test]
    fn write_phased_genotype() {
        let header = test_header();
        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .filter_pass()
        .format_keys(&["GT"])
        .add_sample(vec![SampleValue::Genotype(Genotype::phased_diploid(0, 1))])
        .build(&header)
        .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[9], "0|1");
    }

    // r[verify vcf_writer.missing_dot]
    #[test]
    fn write_reference_site() {
        let header = test_header();
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::G))
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[2], ".");
        assert_eq!(f[4], ".");
        assert_eq!(f[5], ".");
        assert_eq!(f[6], ".");
        assert_eq!(f[7], ".");
    }

    // r[verify vcf_writer.output_formats]
    #[test]
    fn bgzf_writer_produces_valid_output() {
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
        let mut writer = VcfWriter::bgzf(&mut output, header);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
        let index = writer.finish().unwrap();

        // Should have produced BGZF output
        assert!(output.len() > 28);
        assert_eq!(output[0], 0x1f);
        assert_eq!(output[1], 0x8b);

        // Should have an index
        assert!(index.is_some());
    }

    // r[verify vcf_writer.percent_encoding]
    #[test]
    fn percent_encoding_special_chars() {
        let mut buf = Vec::new();
        percent_encode_into(&mut buf, b"hello:world;key=val%done,x\ty\nz\rend");
        let result = String::from_utf8(buf).unwrap();
        assert_eq!(result, "hello%3Aworld%3Bkey%3Dval%25done%2Cx%09y%0Az%0Dend");
    }

    // r[verify vcf_writer.percent_encoding]
    #[test]
    fn info_string_value_is_percent_encoded() {
        let header = Arc::new(
            VcfHeader::builder()
                .add_contig("chr1", ContigDef { length: Some(1000) })
                .unwrap()
                .add_info(
                    "ANN",
                    InfoDef {
                        number: Number::Count(1),
                        typ: ValueType::String,
                        description: SmolStr::from("Annotation"),
                    },
                )
                .unwrap()
                .build()
                .unwrap(),
        );
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_string("ANN", "val:with;special=chars")
                .build(&header)
                .unwrap();

        let text = write_vcf(&header, &[record]);
        let f = last_line_fields(&text);
        assert_eq!(f[7], "ANN=val%3Awith%3Bspecial%3Dchars");
    }
}
