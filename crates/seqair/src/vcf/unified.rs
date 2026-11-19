//! Unified VCF/BCF writer with typestate enforcement.

use super::alleles::Alleles;
use super::bcf_encoding::*;
use super::encoder::{BcfRecordEncoder, BcfValue, BgzfWrite, ContigHandle, FieldTracker};
use super::error::{RecordIdProblem, VcfError};
use super::field_set::FieldSet;
use super::header::{HeaderId, VcfHeader};
use super::reader::rlen::{RlenColumns, RlenKeys, rlen};
use super::record::Genotype;
use super::record_encoder::{
    ContigId, FieldId, FilterId, FormatEncoder, FormatKey, InfoEncoder, InfoKey,
};
use super::writer::{needs_percent, percent_encode_into, push_percent, write_float_g};
use super::{IndexKind, OutputFormat};
use crate::bam::bgzf_writer::BgzfSink;
use crate::io::IndexBuilder;
use crate::io::{BgzfWriter, CompressionLevel};
use seqair_types::{Pos1, SmallVec, SmolStr};
use std::io::Write;
use std::marker::PhantomData;

// ── Writer state types ─────────────────────────────────────────────────

/// Writer state: header has not been written yet.
pub struct Unstarted;

/// Writer state: header has been written, ready to encode records.
pub struct Ready;

// ── Record encoder state types ─────────────────────────────────────────

// r[impl record_encoder.typestate_states+2]
/// Record state: fixed fields written, awaiting filter.
pub struct Begun;

/// Record state: filter written, INFO fields may be encoded.
pub struct Filtered;

/// Record state: sample count declared, FORMAT fields may be encoded.
pub struct WithSamples;

// r[impl vcf_passthrough.edit+2]
// r[impl vcf_passthrough.edit.alleles]
/// Record state: a copied record ([`Writer::copy_record`]), whose INFO and
/// FORMAT fields may be encoded, replaced or removed in any order.
pub struct Edit;

// r[verify vcf_passthrough.edit.alleles]
/// A copied record open for edits: see [`Writer::copy_record`]. Its alleles
/// are the input's: every `Number=A`, `R` and `G` field and every GT names
/// them, so nothing changes REF or ALT.
///
/// ```compile_fail,E0599
/// use seqair::vcf::{Alleles, RecordEdit};
///
/// fn realign(edit: &mut RecordEdit<'_>, alleles: &Alleles) {
///     // ERROR: a copy's alleles cannot be replaced
///     edit.set_alleles(alleles);
/// }
/// ```
pub type RecordEdit<'a> = RecordEncoder<'a, Edit>;

mod phase {
    pub trait Sealed {}
    impl Sealed for super::Filtered {}
    impl Sealed for super::WithSamples {}
    impl Sealed for super::Edit {}
}

/// The record states that take INFO fields: [`Filtered`] and [`Edit`].
pub trait InfoPhase: phase::Sealed {}
impl InfoPhase for Filtered {}
impl InfoPhase for Edit {}

/// The record states that take FORMAT fields: [`WithSamples`] and [`Edit`].
pub trait FormatPhase: phase::Sealed {}
impl FormatPhase for WithSamples {}
impl FormatPhase for Edit {}

// ── Writer ─────────────────────────────────────────────────────────────

// r[impl record_encoder.writer]
// r[impl record_encoder.writer_typestate]
/// Unified VCF/BCF writer. Output format is selected at construction time.
///
/// See the [module documentation](self) for complete usage examples.
///
/// The legal order — and the guard for the two rejections below. Every name
/// they mention appears here in a use that has to compile, so renaming one
/// fails *this* example instead of quietly making them free:
///
/// ```
/// use seqair::vcf::{Alleles, ContigDef, OutputFormat, VcfHeader, Writer};
/// use seqair_types::{Base, Pos1};
/// use std::sync::Arc;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut builder = VcfHeader::builder();
/// let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000)))?;
/// let header = Arc::new(builder.build()?);
///
/// let mut buf = Vec::new();
/// let writer = Writer::new(&mut buf, OutputFormat::Vcf);
/// let mut writer = writer.write_header(&header)?;
///
/// let alleles = Alleles::snv(Base::A, Base::T)?;
/// let pos = Pos1::new(100).expect("100 is a valid 1-based position");
/// let enc = writer.begin_record(&chr1, pos, &alleles, None)?;
/// enc.filter_pass().begin_samples().emit()?;
/// # Ok(())
/// # }
/// ```
///
/// `begin_record` is only available after the header has been written:
///
/// ```compile_fail,E0599
/// use seqair::vcf::{Alleles, ContigDef, OutputFormat, VcfHeader, Writer};
/// use seqair_types::{Base, Pos1};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut builder = VcfHeader::builder();
/// let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000)))?;
/// let alleles = Alleles::snv(Base::A, Base::T)?;
/// let pos = Pos1::new(100).expect("100 is a valid 1-based position");
///
/// let mut buf = Vec::new();
/// let mut writer = Writer::new(&mut buf, OutputFormat::Vcf);
/// // ERROR: begin_record requires Writer<_, Ready>, not Unstarted
/// writer.begin_record(&chr1, pos, &alleles, None)?;
/// # Ok(())
/// # }
/// ```
///
/// `write_header` cannot be called twice:
///
/// ```compile_fail,E0599
/// use seqair::vcf::{OutputFormat, VcfHeader, Writer};
/// use std::sync::Arc;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let header = Arc::new(VcfHeader::builder().build()?);
/// let mut buf = Vec::new();
/// let writer = Writer::new(&mut buf, OutputFormat::Vcf);
/// let writer = writer.write_header(&header)?;
/// // ERROR: write_header is on Unstarted, not Ready
/// writer.write_header(&header)?;
/// # Ok(())
/// # }
/// ```
pub struct Writer<W: Write, S = Unstarted> {
    inner: WriterInner<W>,
    /// Contig names in rid order, captured from the header. Empty for plain
    /// VCF; used to label the coordinate index returned by [`Writer::finish`].
    contig_names: Vec<SmolStr>,
    /// The written header's id, set by `write_header`: `begin_record` refuses
    /// a contig resolved from another header.
    pub(crate) header: Option<HeaderId>,
    /// BGZF settings for `VcfGz` and `Bcf`, applied by `write_header`.
    compression: Compression,
    /// The index `VcfGz` co-produces; `Bcf` takes only CSI.
    index_kind: IndexKind,
    /// Per-record state outside the format's encoder, reused across records.
    pub(crate) record: RecordState,
    _state: PhantomData<S>,
}

/// What a record carries besides its bytes, for checks and for its length.
#[derive(Default)]
pub(crate) struct RecordState {
    rlen: RlenState,
    /// Each sample's GT ploidy once GT is encoded, for `Number=G` counts
    /// (`r[record_encoder.count_check+3]`); empty before GT.
    ploidies: Vec<usize>,
    /// Buffers a passthrough copy gathers values in, reused across records.
    pub(crate) copy: super::passthrough::CopyScratch,
    /// The written header and the field map `copy_record` built from it.
    pub(crate) implicit: super::passthrough::ImplicitMap,
}

impl RecordState {
    fn begin(&mut self, alleles: &Alleles) {
        self.ploidies.clear();
        self.rlen.begin(alleles);
    }
}

// r[impl record_encoder.rlen]
/// The record's INFO END/SVLEN and FORMAT LEN as VCF text, the only fields
/// besides REF that change its length, kept for [`rlen`] at `emit`.
#[derive(Default)]
struct RlenState {
    keys: RlenKeys,
    /// Dictionary indices of the three, where the header declares them
    /// `Type=Integer` (`None` otherwise: they cannot change the length).
    end: Option<u32>,
    svlen: Option<u32>,
    len: Option<u32>,
    /// The record's REF and ALT columns, written only when one of the three is declared.
    reference: Vec<u8>,
    alt: Vec<u8>,
    /// `END=…`, `SVLEN=…` and `LEN\t…\t…` as encoded; empty if not encoded.
    end_text: Vec<u8>,
    svlen_text: Vec<u8>,
    len_text: Vec<u8>,
    info: Vec<u8>,
}

impl RlenState {
    fn new(header: &VcfHeader) -> Self {
        let keys = RlenKeys::from_header(header);
        let idx = |declared: bool, name: &str| {
            let i = header.string_map().get(name).filter(|_| declared)?;
            u32::try_from(i).ok()
        };
        Self {
            end: idx(keys.end, "END"),
            svlen: idx(keys.svlen, "SVLEN"),
            len: idx(keys.len, "LEN"),
            keys,
            ..Self::default()
        }
    }

    fn declared(&self) -> bool {
        self.end.is_some() || self.svlen.is_some() || self.len.is_some()
    }

    fn begin(&mut self, alleles: &Alleles) {
        self.end_text.clear();
        self.svlen_text.clear();
        self.len_text.clear();
        if self.declared() {
            self.reference.clear();
            alleles.write_ref_into(&mut self.reference);
            self.alt.clear();
            if alleles.write_alts_into(&mut self.alt) == 0 {
                self.alt.push(b'.');
            }
        }
    }

    /// The buffer an INFO field's values go into, `END=`/`SVLEN=` written,
    /// if it is one the length depends on.
    fn info_target(&mut self, dict_idx: u32) -> Option<&mut Vec<u8>> {
        let (buf, key): (&mut Vec<u8>, &[u8]) = if self.end == Some(dict_idx) {
            (&mut self.end_text, b"END=")
        } else if self.svlen == Some(dict_idx) {
            (&mut self.svlen_text, b"SVLEN=")
        } else {
            return None;
        };
        buf.clear();
        buf.extend_from_slice(key);
        Some(buf)
    }

    // r[impl vcf_passthrough.edit.remove]
    /// Forget INFO field `dict_idx` if the length depends on it.
    fn forget_info(&mut self, dict_idx: u32) {
        if self.end == Some(dict_idx) {
            self.end_text.clear();
        }
        if self.svlen == Some(dict_idx) {
            self.svlen_text.clear();
        }
    }

    // r[impl vcf_passthrough.edit.remove]
    /// Forget FORMAT field `dict_idx` if it is LEN.
    fn forget_format(&mut self, dict_idx: u32) {
        if self.len == Some(dict_idx) {
            self.len_text.clear();
        }
    }

    /// The buffer FORMAT LEN's samples go into, each after a TAB.
    fn format_target(&mut self, dict_idx: u32) -> Option<&mut Vec<u8>> {
        if self.len != Some(dict_idx) {
            return None;
        }
        self.len_text.clear();
        self.len_text.extend_from_slice(b"LEN");
        Some(&mut self.len_text)
    }

    /// The record's reference length: `ref_len` (REF's) unless it encoded one
    /// of the three, else what a reader computes from the record's text.
    fn compute(&mut self, pos_0based: i64, ref_len: i64) -> i64 {
        if self.end_text.is_empty() && self.svlen_text.is_empty() && self.len_text.is_empty() {
            return ref_len;
        }
        self.info.clear();
        self.info.extend_from_slice(&self.end_text);
        if !self.end_text.is_empty() && !self.svlen_text.is_empty() {
            self.info.push(b';');
        }
        self.info.extend_from_slice(&self.svlen_text);
        if self.info.is_empty() {
            self.info.push(b'.');
        }
        let cols = RlenColumns {
            reference: &self.reference,
            alt: &self.alt,
            info: &self.info,
            format: &self.len_text,
        };
        rlen(&cols, self.keys, pos_0based)
    }
}

/// Integers as VCF text writes them, comma-separated, `.` for a missing one.
fn push_ints(buf: &mut Vec<u8>, values: impl IntoIterator<Item = Option<i32>>) {
    let mut itoa = itoa::Buffer::new();
    for (i, v) in values.into_iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        match v {
            Some(v) => buf.extend_from_slice(itoa.format(v).as_bytes()),
            None => buf.push(b'.'),
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Compression {
    level: CompressionLevel,
    threads: usize,
}

// r[impl vcf_writer.buffered_output]
/// Capacity of the plain-VCF output buffer.
const PLAIN_OUTPUT_BUFFER: usize = 128 * 1024;

// r[impl vcf_writer.buffered_output]
/// Plain VCF's buffered sink. Dropped without [`into_inner`](Self::into_inner),
/// it flushes and logs a failure — `BufWriter`'s own drop would discard the
/// error, and with it up to a buffer's worth of records, without a word.
struct PlainOutput<W: Write>(Option<std::io::BufWriter<W>>);

impl<W: Write> PlainOutput<W> {
    fn new(inner: W) -> Self {
        Self(Some(std::io::BufWriter::with_capacity(PLAIN_OUTPUT_BUFFER, inner)))
    }

    /// Flush the buffer and then the sink, and hand the sink back.
    fn into_inner(mut self) -> Result<W, VcfError> {
        let Some(buffered) = self.0.take() else {
            return Err(VcfError::Io(std::io::ErrorKind::NotConnected.into()));
        };
        let mut inner = buffered.into_inner().map_err(|e| VcfError::Io(e.into_error()))?;
        inner.flush().map_err(VcfError::Io)?;
        Ok(inner)
    }
}

// The `None` arms are unreachable: only `into_inner` empties the option, and it
// consumes `self`.
impl<W: Write> Write for PlainOutput<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.0 {
            Some(w) => w.write(buf),
            None => Err(std::io::ErrorKind::NotConnected.into()),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match &mut self.0 {
            Some(w) => w.write_all(buf),
            None => Err(std::io::ErrorKind::NotConnected.into()),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match &mut self.0 {
            Some(w) => w.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write> Drop for PlainOutput<W> {
    fn drop(&mut self) {
        if let Some(w) = &mut self.0
            && let Err(e) = w.flush()
        {
            tracing::warn!("vcf::Writer dropped without finish(): flushing VCF output failed: {e}");
        }
    }
}

enum WriterInner<W: Write> {
    Vcf {
        output: PlainOutput<W>,
        buf: Vec<u8>,
        fmt_keys: Vec<SmolStr>,
        sample_bufs: Vec<Vec<u8>>,
        n_samples: u32,
        info_tracker: FieldTracker,
        fmt_seen: FieldSet,
    },
    VcfGz {
        bgzf: BgzfSink<W>,
        index: Option<IndexBuilder>,
        buf: Vec<u8>,
        fmt_keys: Vec<SmolStr>,
        sample_bufs: Vec<Vec<u8>>,
        n_samples: u32,
        info_tracker: FieldTracker,
        fmt_seen: FieldSet,
    },
    Bcf {
        bgzf: BgzfSink<W>,
        index: Option<IndexBuilder>,
        shared_buf: Vec<u8>,
        indiv_buf: Vec<u8>,
        n_samples: u32,
        info_tracker: FieldTracker,
        fmt_tracker: FieldTracker,
    },
}

// r[impl record_encoder.writer_new+2]
impl<W: Write> Writer<W> {
    /// Create a writer for the given output format. No header is needed at this point.
    pub fn new(inner: W, format: OutputFormat) -> Self {
        let inner = match format {
            OutputFormat::Vcf => WriterInner::Vcf {
                output: PlainOutput::new(inner),
                buf: Vec::with_capacity(4096),
                fmt_keys: Vec::with_capacity(8),
                sample_bufs: Vec::new(),
                n_samples: 0,
                info_tracker: FieldTracker::default(),
                fmt_seen: FieldSet::default(),
            },
            OutputFormat::VcfGz => WriterInner::VcfGz {
                bgzf: BgzfSink::Serial(BgzfWriter::new(inner)),
                index: None,
                buf: Vec::with_capacity(4096),
                fmt_keys: Vec::with_capacity(8),
                sample_bufs: Vec::new(),
                n_samples: 0,
                info_tracker: FieldTracker::default(),
                fmt_seen: FieldSet::default(),
            },
            OutputFormat::Bcf => WriterInner::Bcf {
                bgzf: BgzfSink::Serial(BgzfWriter::new(inner)),
                index: None,
                shared_buf: Vec::with_capacity(4096),
                indiv_buf: Vec::with_capacity(4096),
                n_samples: 0,
                info_tracker: FieldTracker::default(),
                fmt_tracker: FieldTracker::default(),
            },
        };
        Writer {
            inner,
            contig_names: Vec::new(),
            header: None,
            compression: Compression::default(),
            index_kind: IndexKind::default(),
            record: RecordState::default(),
            _state: PhantomData,
        }
    }

    // r[impl record_encoder.index_kind]
    /// The coordinate index a `VcfGz` writer co-produces: [`IndexKind::Csi`]
    /// (the default) or [`IndexKind::Tbi`]. BCF is indexed by CSI only:
    /// [`write_header`](Self::write_header) refuses `Tbi` for it with
    /// [`VcfError::TbiForBcf`]. Plain VCF is not indexed and ignores this.
    ///
    /// A TBI addresses positions below 2^29 only; a record that ends past it
    /// fails its `emit` with [`IndexError::BeyondIndexRange`](crate::io::IndexError::BeyondIndexRange).
    #[must_use]
    pub fn index_kind(mut self, kind: IndexKind) -> Self {
        self.index_kind = kind;
        self
    }

    // r[impl record_encoder.compression_level]
    /// BGZF compression level for `VcfGz` and `Bcf`. Defaults to
    /// [`CompressionLevel::Level6`], htslib's default. Plain VCF is not
    /// compressed and ignores this.
    #[must_use]
    pub fn compression_level(mut self, level: CompressionLevel) -> Self {
        self.compression.level = level;
        self
    }

    // r[impl record_encoder.compression_threads]
    /// Compress BGZF blocks (`VcfGz`, `Bcf`) on `threads` worker threads.
    /// Defaults to `0`: compress on the calling thread. Plain VCF is not
    /// compressed and ignores this.
    ///
    /// The calling thread keeps encoding records and writes the compressed
    /// blocks in order, so the output — and the coordinate index — is
    /// byte-identical to the single-threaded writer's. The threads start
    /// with [`write_header`](Self::write_header), which reports a failure to
    /// start them.
    #[must_use]
    pub fn compression_threads(mut self, threads: usize) -> Self {
        self.compression.threads = threads;
        self
    }

    /// Rebuild the BGZF sink, still unwritten, with the configured settings.
    fn apply_compression(self) -> Result<Self, VcfError> {
        if self.compression == Compression::default() {
            return Ok(self);
        }
        let Writer { inner, contig_names, header, compression, index_kind, record, _state } = self;
        let Compression { level, threads } = compression;
        let inner = match inner {
            WriterInner::Vcf { .. } => inner,
            WriterInner::VcfGz {
                bgzf,
                index,
                buf,
                fmt_keys,
                sample_bufs,
                n_samples,
                info_tracker,
                fmt_seen,
            } => WriterInner::VcfGz {
                bgzf: bgzf.with_compression(level, threads)?,
                index,
                buf,
                fmt_keys,
                sample_bufs,
                n_samples,
                info_tracker,
                fmt_seen,
            },
            WriterInner::Bcf {
                bgzf,
                index,
                shared_buf,
                indiv_buf,
                n_samples,
                info_tracker,
                fmt_tracker,
            } => WriterInner::Bcf {
                bgzf: bgzf.with_compression(level, threads)?,
                index,
                shared_buf,
                indiv_buf,
                n_samples,
                info_tracker,
                fmt_tracker,
            },
        };
        Ok(Writer { inner, contig_names, header, compression, index_kind, record, _state })
    }

    // r[impl record_encoder.write_header+2]
    /// Write the file header. Consumes the `Unstarted` writer and returns `Ready`.
    pub fn write_header(self, header: &VcfHeader) -> Result<Writer<W, Ready>, VcfError> {
        // r[impl record_encoder.index_kind]
        if self.index_kind == IndexKind::Tbi && matches!(self.inner, WriterInner::Bcf { .. }) {
            return Err(VcfError::TbiForBcf);
        }
        let mut writer = self.apply_compression()?;
        // r[impl vcf_header.serialization+2]
        let header_text = match writer.inner {
            WriterInner::Bcf { .. } => header.to_bcf_text(),
            _ => header.to_vcf_text(),
        };
        // Indexed by rid: a parsed header's contig dictionary may have holes.
        let n_refs = header.contig_map().len();
        // Only compressed output is indexed; plain VCF needs no contig names.
        let contig_names: Vec<SmolStr> = match writer.inner {
            WriterInner::Vcf { .. } => Vec::new(),
            // A hole names no contig and no record can use its rid; the index
            // needs a name per rid, so it gets an empty one.
            _ => (0..n_refs)
                .map(|rid| header.contig_map().name(rid).map(SmolStr::from).unwrap_or_default())
                .collect(),
        };
        #[expect(
            clippy::cast_possible_truncation,
            reason = "sample count bounded by VcfHeader builder; realistic files have <100 samples"
        )]
        let header_n_samples = header.samples().len() as u32;
        // The longest reference decides how deep the CSI has to bin. A header
        // with contigs but no declared lengths is `None`, which takes htslib's
        // assumption rather than a guess of ours; an empty header indexes
        // nothing, so its depth is immaterial and stays at the BAI default.
        let max_ref_len: Option<u64> = if header.contig_defs().is_empty() {
            Some(0)
        } else {
            header.contig_defs().values().filter_map(|contig| contig.length).max()
        };

        match &mut writer.inner {
            WriterInner::Vcf { output, n_samples, .. } => {
                output.write_all(header_text.as_bytes()).map_err(VcfError::Io)?;
                *n_samples = header_n_samples;
            }
            WriterInner::VcfGz { bgzf, index, n_samples, .. } => {
                bgzf.write_all(header_text.as_bytes())?;
                *index = Some(match writer.index_kind {
                    IndexKind::Csi => {
                        IndexBuilder::csi(n_refs, 14, max_ref_len, bgzf.index_offset())
                    }
                    IndexKind::Tbi => IndexBuilder::tbi(n_refs, bgzf.index_offset()),
                });
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
                *index = Some(IndexBuilder::csi(n_refs, 14, max_ref_len, bgzf.index_offset()));
                *n_samples = header_n_samples;
            }
        }

        Ok(Writer {
            inner: writer.inner,
            contig_names,
            header: Some(header.id()),
            compression: writer.compression,
            index_kind: writer.index_kind,
            record: RecordState {
                rlen: RlenState::new(header),
                implicit: super::passthrough::ImplicitMap::new(header),
                ..RecordState::default()
            },
            _state: PhantomData,
        })
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
        // r[impl vcf_values.header_brand+2]
        if self.header != Some(contig.header) {
            return Err(VcfError::ForeignKey { name: contig.name.clone() });
        }
        // r[impl record_encoder.reserved_rejection+2]
        // A QUAL whose bits are the float MISSING marker would read back as an
        // absent QUAL — indistinguishable from `None`.
        reject_reserved_floats("QUAL", qual)?;
        self.record.begin(alleles);
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
                    alleles_start: 0,
                    editing: false,
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
                fmt_seen,
            } => {
                info_tracker.clear();
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
                    fmt_seen,
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
                fmt_seen,
            } => {
                info_tracker.clear();
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
                    fmt_seen,
                )?)
            }
        };

        Ok(RecordEncoder {
            inner,
            header: contig.header,
            record: &mut self.record,
            _state: PhantomData,
        })
    }

    // r[impl record_encoder.finish+2]
    /// Finalize the writer, returning the inner writer and an optional
    /// [`CoordinateIndex`].
    ///
    /// For VCF.gz and BCF, all buffered data is flushed and the BGZF EOF marker
    /// is written before this returns — the caller does not need to flush `W`.
    /// The index is `Some` for `VcfGz`/`Bcf` and `None` for plain VCF; it is
    /// self-describing, so callers write it with [`CoordinateIndex::write`]
    /// without choosing a serializer or re-supplying contig names.
    pub fn finish(self) -> Result<(W, Option<CoordinateIndex>), VcfError> {
        let contig_names = self.contig_names;
        let kind = self.index_kind;
        match self.inner {
            WriterInner::Vcf { output, .. } => {
                // r[impl vcf_writer.buffered_output]
                Ok((output.into_inner()?, None))
            }
            WriterInner::VcfGz { mut bgzf, mut index, .. } => {
                // r[impl vcf_writer.finish+2]
                // r[impl vcf_writer.index_amend_last]
                // End the last block first, so the index ends at the next
                // block's start, as htslib's `sam_idx_save` does.
                if let Some(ref mut idx) = index {
                    bgzf.flush_if_needed(usize::MAX)?;
                    let voff = bgzf.index_offset();
                    idx.amend_last(voff);
                    idx.finish(voff)?;
                }
                let inner = bgzf.finish(index.as_mut())?;
                let index = index.map(|builder| CoordinateIndex {
                    builder,
                    format: OutputFormat::VcfGz,
                    kind,
                    contig_names,
                });
                Ok((inner, index))
            }
            WriterInner::Bcf { mut bgzf, mut index, .. } => {
                // r[impl bcf_writer.finish]
                // r[impl vcf_writer.index_amend_last]
                // End the last block first, so the index ends at the next
                // block's start, as htslib's `sam_idx_save` does.
                if let Some(ref mut idx) = index {
                    bgzf.flush_if_needed(usize::MAX)?;
                    let voff = bgzf.index_offset();
                    idx.amend_last(voff);
                    idx.finish(voff)?;
                }
                let inner = bgzf.finish(index.as_mut())?;
                let index = index.map(|builder| CoordinateIndex {
                    builder,
                    format: OutputFormat::Bcf,
                    kind: IndexKind::Csi,
                    contig_names,
                });
                Ok((inner, index))
            }
        }
    }
}

/// A finalized coordinate index for a VCF.gz or BCF file, ready to serialize.
///
/// Returned by [`Writer::finish`] for compressed output. It bundles the
/// [`IndexBuilder`], the [`OutputFormat`], and the contig names in rid order, so
/// callers serialize with one [`write`](Self::write) call instead of selecting
/// the right `IndexBuilder` method and re-supplying names and ref counts.
///
/// Convention: the index lives next to the data file with [`suffix`](Self::suffix)
/// appended (e.g. `calls.bcf` → `calls.bcf.csi`, `calls.vcf.gz` → `calls.vcf.gz.tbi`
/// for a writer set to [`IndexKind::Tbi`]).
///
/// ```
/// use seqair_types::{Base, Pos1};
/// use seqair::vcf::{Alleles, ContigDef, CoordinateIndex, OutputFormat, VcfHeader, Writer};
/// use std::sync::Arc;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut builder = VcfHeader::builder();
/// let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1000)))?;
/// let header = Arc::new(builder.build()?);
///
/// let mut bcf = Vec::new();
/// let mut writer = Writer::new(&mut bcf, OutputFormat::Bcf).write_header(&header)?;
/// let alleles = Alleles::snv(Base::A, Base::T)?;
/// writer
///     .begin_record(&chr1, Pos1::new(100).unwrap(), &alleles, Some(30.0))?
///     .filter_pass()
///     .emit()?;
///
/// // Compressed output yields an index from finish(); plain VCF yields `None`.
/// let (_, index) = writer.finish()?;
/// let index = index.expect("BCF output is indexed");
///
/// // Write it next to the data file, e.g. `calls.bcf.{suffix}`:
/// assert_eq!(index.suffix(), "csi");
/// let mut csi = Vec::new();
/// index.write(&mut csi)?;
/// assert!(!csi.is_empty());
/// # Ok(()) }
/// ```
pub struct CoordinateIndex {
    builder: IndexBuilder,
    format: OutputFormat,
    kind: IndexKind,
    contig_names: Vec<SmolStr>,
}

impl CoordinateIndex {
    /// Which index this is: what the writer was set to for VCF.gz, CSI for BCF.
    pub fn kind(&self) -> IndexKind {
        self.kind
    }

    /// The index file's suffix, without the dot: it lives at
    /// `<data path>.<suffix>` ([`IndexKind::suffix`] of [`kind`](Self::kind)).
    pub fn suffix(&self) -> &'static str {
        self.kind.suffix()
    }

    // r[impl record_encoder.finish+2]
    /// Serialize the index (BGZF-compressed) to `writer`, choosing the correct
    /// on-disk format: TBI, or CSI with tabix metadata, for VCF.gz; plain CSI
    /// for BCF.
    pub fn write<W: Write>(&self, writer: W) -> Result<(), VcfError> {
        match (self.format, self.kind) {
            (OutputFormat::VcfGz, IndexKind::Tbi) => {
                self.builder.write_tbi(writer, &self.contig_names)?;
            }
            (OutputFormat::VcfGz, IndexKind::Csi) => {
                self.builder.write_csi_tabix(writer, &self.contig_names)?;
            }
            (OutputFormat::Bcf, _) => {
                self.builder.write_csi(writer, self.contig_names.len(), &[])?;
            }
            // finish() never builds a CoordinateIndex for plain VCF.
            (OutputFormat::Vcf, _) => {}
        }
        Ok(())
    }

    /// Serialize the index to `path`, creating or truncating the file.
    pub fn write_to_path(&self, path: impl AsRef<std::path::Path>) -> Result<(), VcfError> {
        let file = std::fs::File::create(path).map_err(VcfError::Io)?;
        self.write(std::io::BufWriter::new(file))
    }

    /// Borrow the underlying [`IndexBuilder`] for direct serialization (e.g. to
    /// compare formats in tests).
    pub fn builder(&self) -> &IndexBuilder {
        &self.builder
    }

    /// Consume the index, returning the underlying [`IndexBuilder`].
    pub fn into_builder(self) -> IndexBuilder {
        self.builder
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
            VcfOutput::Bgzf { bgzf, index } => {
                bgzf.flush_if_needed(buf.len())?;
                // r[impl vcf_writer.index_amend_last]
                // A new block may have started: name the previous record's
                // end by it, as htslib's `vcf_write` does.
                if let Some(idx) = index {
                    idx.amend_last(bgzf.index_offset());
                }
                bgzf.write_all(buf)?;
                Ok(())
            }
        }
    }

    fn push_index(&mut self, tid: u32, beg: u64, end: u64) -> Result<(), VcfError> {
        match self {
            VcfOutput::Plain(_) => Ok(()),
            VcfOutput::Bgzf { bgzf, index } => {
                if let Some(idx) = index {
                    idx.push(tid, beg, end, bgzf.index_offset())?;
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
    tid: u32,
    pos_0based: u32,
    rlen: u32,
    /// The ID column's bytes in `buf`, which `RecordEncoder::ids` replaces.
    id_range: (usize, usize),
    n_samples: u32,
    filter_written: bool,
    /// A passthrough edit (`r[vcf_passthrough.edit+2]`): replacing a field
    /// is the point, so a duplicate is not warned of.
    editing: bool,
    info_tracker: &'a mut FieldTracker,
    /// The `dict_idx` of every FORMAT field in `fmt_keys`, so the common
    /// first-write path needs no name comparison.
    fmt_seen: &'a mut FieldSet,
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
    fmt_seen: &'a mut FieldSet,
) -> Result<VcfEncoderFields<'a>, VcfError> {
    buf.clear();
    fmt_keys.clear();
    fmt_seen.clear();
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
    let pos_0based = pos.to_zero_based().as_u32();
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
    buf.extend_from_slice(itoa_buf.format(pos.as_u32()).as_bytes());
    buf.push(b'\t');
    let id_start = buf.len();
    buf.push(b'.');
    let id_range = (id_start, buf.len());
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
        tid: contig.tid(),
        pos_0based,
        rlen,
        id_range,
        n_samples,
        filter_written: false,
        editing: false,
        info_tracker,
        fmt_seen,
    })
}

// ── RecordEncoder ──────────────────────────────────────────────────────

// r[impl record_encoder.typestate]
// r[impl record_encoder.typestate_must_use]
/// Record encoder with typestate enforcement of the
/// `Begun` → `Filtered` → `WithSamples` → `emit()` chain.
///
/// The three rejections below are written against the handles a caller
/// actually holds — `InfoKey::encode` and `FormatKey::encode` — rather than
/// the `InfoEncoder`/`FormatEncoder` methods, whose [`FieldId`] argument has
/// no public constructor. This example is their guard: it uses every name
/// they use, legally, so a rename breaks it loudly instead of leaving three
/// tests that pass because their `use` no longer resolves.
///
/// [`FieldId`]: super::record_encoder::FieldId
///
/// ```
/// use seqair::vcf::{Begun, Filtered, FormatInt, InfoInt, RecordEncoder, VcfError, WithSamples};
///
/// fn filter_comes_first<'a>(enc: RecordEncoder<'a, Begun>) -> RecordEncoder<'a, Filtered> {
///     enc.filter_pass()
/// }
///
/// fn info_belongs_to_filtered(enc: &mut RecordEncoder<'_, Filtered>, dp: &InfoInt) {
///     dp.encode(enc, 42);
/// }
///
/// fn samples_come_next<'a>(enc: RecordEncoder<'a, Filtered>) -> RecordEncoder<'a, WithSamples> {
///     enc.begin_samples()
/// }
///
/// fn format_belongs_to_with_samples(
///     enc: &mut RecordEncoder<'_, WithSamples>,
///     dp: &FormatInt,
/// ) -> Result<(), VcfError> {
///     dp.encode(enc, &[42])
/// }
/// ```
///
/// INFO fields cannot be encoded before the filter is set:
///
/// ```compile_fail,E0277
/// use seqair::vcf::{Begun, InfoInt, RecordEncoder};
///
/// fn check(enc: &mut RecordEncoder<'_, Begun>, dp: &InfoInt) {
///     // ERROR: InfoEncoder is not implemented for RecordEncoder<Begun>
///     dp.encode(enc, 42);
/// }
/// ```
///
/// FORMAT fields cannot be encoded before `begin_samples`:
///
/// ```compile_fail,E0277
/// use seqair::vcf::{Filtered, FormatInt, RecordEncoder, VcfError};
///
/// fn check(enc: &mut RecordEncoder<'_, Filtered>, dp: &FormatInt) -> Result<(), VcfError> {
///     // ERROR: FormatEncoder is not implemented for RecordEncoder<Filtered>
///     dp.encode(enc, &[42])
/// }
/// ```
///
/// INFO fields cannot be encoded after `begin_samples`:
///
/// ```compile_fail,E0277
/// use seqair::vcf::{InfoInt, RecordEncoder, WithSamples};
///
/// fn check(enc: &mut RecordEncoder<'_, WithSamples>, dp: &InfoInt) {
///     // ERROR: InfoEncoder is not implemented for RecordEncoder<WithSamples>
///     dp.encode(enc, 42);
/// }
/// ```
#[must_use = "record is silently discarded if emit() is not called"]
pub struct RecordEncoder<'a, S> {
    // r[impl record_encoder.typestate_w_erased]
    inner: EncoderInner<'a>,
    /// The written header's id, which every key must carry
    /// (`r[record_encoder.foreign_keys]`).
    header: HeaderId,
    /// GT ploidies and what the record's length depends on.
    pub(crate) record: &'a mut RecordState,
    _state: PhantomData<S>,
}

impl<'a, S> RecordEncoder<'a, S> {
    // r[impl vcf_passthrough.edit+2]
    /// This record as an edit of a copy, whose replacements are not warned of.
    pub(crate) fn into_edit(mut self) -> RecordEncoder<'a, Edit> {
        match &mut self.inner {
            EncoderInner::Bcf(e) => e.editing = true,
            EncoderInner::Vcf(v) => v.editing = true,
        }
        self.into_state()
    }

    pub(crate) fn into_state<T>(self) -> RecordEncoder<'a, T> {
        RecordEncoder {
            inner: self.inner,
            header: self.header,
            record: self.record,
            _state: PhantomData,
        }
    }

    // r[impl record_encoder.foreign_keys]
    fn check_key(&self, id: &FieldId) -> Result<(), VcfError> {
        if id.header == self.header {
            return Ok(());
        }
        Err(VcfError::ForeignKey { name: id.name.clone() })
    }

    fn n_allele_count(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => usize::from(e.n_allele),
            EncoderInner::Vcf(v) => usize::from(v.n_allele),
        }
    }

    // r[impl record_encoder.count_check+3]
    /// Check every non-empty sample's count; an empty sample is missing (`.`).
    fn check_sample_counts<T>(&self, id: &FieldId, per_sample: &[&[T]]) -> Result<(), VcfError> {
        let n_allele = self.n_allele_count();
        per_sample.iter().enumerate().filter(|(_, s)| !s.is_empty()).try_for_each(|(i, s)| {
            let ploidy = self.record.ploidies.get(i).copied().filter(|&p| p > 0);
            id.check_count(n_allele, ploidy, s.len())
        })
    }
}

// r[impl record_encoder.typestate_w_erased]
enum EncoderInner<'a> {
    Bcf(BcfRecordEncoder<'a>),
    Vcf(VcfEncoderFields<'a>),
}

/// BCF's `n_fmt` is a `u8`; htslib's `MAX_N_FMT` holds VCF text to the same.
const MAX_FORMAT_FIELDS: u8 = u8::MAX;

impl EncoderInner<'_> {
    // r[impl record_encoder.format_count_limit]
    /// Reject a new FORMAT field once the record holds [`MAX_FORMAT_FIELDS`],
    /// before any of its bytes are written. A duplicate always fits.
    fn check_format_room(&self, id: &FieldId) -> Result<(), VcfError> {
        let (count, present) = match self {
            Self::Bcf(enc) => (usize::from(enc.n_fmt), enc.fmt_tracker.contains(id.dict_idx())),
            Self::Vcf(vcf) => (vcf.fmt_keys.len(), vcf.fmt_seen.contains(id.dict_idx())),
        };
        if count < usize::from(MAX_FORMAT_FIELDS) || present {
            return Ok(());
        }
        Err(VcfError::ValueOverflow {
            field: "n_fmt",
            value: u64::from(MAX_FORMAT_FIELDS).saturating_add(1),
            target_type: "u8",
        })
    }
}

// ── Begun → Filtered ───────────────────────────────────────────────────

// r[impl record_encoder.typestate_transitions+4]
// r[impl record_encoder.filters+3]
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
        self.into_state()
    }

    /// Mark the record as failing the given filters. Accepts any iterator of
    /// `&FilterId` (e.g. a slice's `.iter().copied()`, or a `map` over a
    /// caller's lookup table) so callers need not materialize a `&[&FilterId]`.
    ///
    /// # Errors
    /// [`VcfError::EmptyFilterList`] for no filters, [`VcfError::ForeignKey`]
    /// for a filter of another header.
    pub fn filter_fail<'f>(
        mut self,
        filters: impl IntoIterator<Item = &'f FilterId>,
    ) -> Result<RecordEncoder<'a, Filtered>, VcfError> {
        let filters: SmallVec<&FilterId, 3> = filters.into_iter().collect();
        if filters.is_empty() {
            return Err(VcfError::EmptyFilterList);
        }
        // r[impl record_encoder.foreign_keys]
        // PASS (no header) is index 0 in every header.
        if let Some(f) = filters.iter().find(|f| f.header.is_some_and(|h| h != self.header)) {
            return Err(VcfError::ForeignKey { name: f.name.clone() });
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let indices = filters
                    .iter()
                    .map(|f| {
                        i32::try_from(f.dict_idx()).map_err(|_| VcfError::ValueOverflow {
                            field: "FILTER",
                            value: u64::from(f.dict_idx()),
                            target_type: "i32",
                        })
                    })
                    .collect::<Result<SmallVec<i32, 3>, _>>()?;
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
        Ok(self.into_state())
    }

    // r[impl record_encoder.record_id]
    /// Set the record's ID column: `ids` written `;`-separated, or `.` if there
    /// are none. Calling it again replaces the column.
    ///
    /// # Errors
    /// [`VcfError::InvalidRecordId`] for an identifier that is empty, `.`, or
    /// holds whitespace or `;`, and [`VcfError::DuplicateRecordId`] for one
    /// given twice, before anything is written; the record is not emitted.
    pub fn ids<'i>(mut self, ids: impl IntoIterator<Item = &'i str>) -> Result<Self, VcfError> {
        let ids: SmallVec<&str, 4> = ids.into_iter().collect();
        for (i, &id) in ids.iter().enumerate() {
            let problem = if id.is_empty() {
                Some(RecordIdProblem::Empty)
            } else if id == "." {
                Some(RecordIdProblem::Missing)
            } else if id.chars().any(char::is_whitespace) {
                Some(RecordIdProblem::Whitespace)
            } else if id.contains(';') {
                Some(RecordIdProblem::Semicolon)
            } else {
                None
            };
            if let Some(problem) = problem {
                return Err(VcfError::InvalidRecordId { id: SmolStr::from(id), problem });
            }
            if ids.get(..i).is_some_and(|before| before.contains(&id)) {
                return Err(VcfError::DuplicateRecordId { id: SmolStr::from(id) });
            }
        }
        // Written at the buffer's end, then moved over the old column: no allocation.
        let write_column = |buf: &mut Vec<u8>| {
            for (i, id) in ids.iter().enumerate() {
                if i > 0 {
                    buf.push(b';');
                }
                buf.extend_from_slice(id.as_bytes());
            }
            if ids.is_empty() {
                buf.push(b'.');
            }
        };
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let from = enc.shared_buf.len();
                // r[impl bcf_writer.shared_variable+2]
                let len = ids.iter().map(|id| id.len().saturating_add(1)).sum::<usize>();
                encode_type_byte(enc.shared_buf, len.saturating_sub(1), BCF_BT_CHAR);
                if !ids.is_empty() {
                    write_column(enc.shared_buf);
                }
                let new_len = splice_appended(enc.shared_buf, 24..enc.alleles_start, from);
                enc.alleles_start = 24usize.saturating_add(new_len);
            }
            EncoderInner::Vcf(vcf) => {
                let from = vcf.buf.len();
                write_column(vcf.buf);
                let (start, end) = vcf.id_range;
                let new_len = splice_appended(vcf.buf, start..end, from);
                vcf.id_range = (start, start.saturating_add(new_len));
            }
        }
        Ok(self)
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
        self.into_state()
    }
}

/// Move the bytes appended to `buf` at `from..` over `range`, which they
/// replace, keeping everything between in place; returns their length.
fn splice_appended(buf: &mut Vec<u8>, range: std::ops::Range<usize>, from: usize) -> usize {
    debug_assert!(range.start <= range.end && range.end <= from && from <= buf.len());
    let new_len = buf.len().saturating_sub(from);
    if let Some(tail) = buf.get_mut(range.start..) {
        tail.rotate_right(new_len);
    }
    let old = range.start.saturating_add(new_len)..range.end.saturating_add(new_len);
    buf.drain(old);
    new_len
}

// ── Field dedup helpers ───────────────────────────────────────────────

// r[impl record_encoder.info_dedup+2]
impl VcfEncoderFields<'_> {
    /// Open an INFO field: ensure the FILTER column is closed, then emit the
    /// `;` separator (when this is not the first field) and return the byte
    /// offset of the field unit. The caller writes `key=value` after this, then
    /// calls [`commit_info_field`](Self::commit_info_field).
    fn begin_info_field(&mut self) -> usize {
        if !self.filter_written {
            self.filter_written = true;
            self.buf.extend_from_slice(b".\t");
        }
        let start = self.buf.len();
        if self.info_count > 0 {
            self.buf.push(b';');
        }
        start
    }

    /// Close the INFO field opened at `start`. A duplicate is spliced over its
    /// previous occurrence in place (dropping the just-written leading `;` when
    /// the overwritten field is first); a new field bumps `info_count`.
    fn commit_info_field(&mut self, id: &FieldId, start: usize) {
        if self.info_tracker.commit(self.buf, id.dict_idx(), start, true) {
            self.info_count = self.info_count.saturating_add(1);
        } else {
            warn_twice(self.editing, id);
        }
    }
}

// r[impl record_encoder.info_dedup+2]
impl BcfRecordEncoder<'_> {
    /// Commit an INFO field whose bytes were just appended at `[start..]` in
    /// `shared_buf`. On a duplicate, splices the new bytes over the previous
    /// occurrence in place; otherwise bumps `n_info`.
    fn commit_info_field(&mut self, id: &FieldId, start: usize) {
        if self.info_tracker.commit(self.shared_buf, id.dict_idx(), start, false) {
            self.n_info = self.n_info.saturating_add(1);
        } else {
            warn_twice(self.editing, id);
        }
    }

    // r[impl record_encoder.format_dedup+2]
    /// FORMAT counterpart of [`commit_info_field`](Self::commit_info_field),
    /// operating on `indiv_buf` and `n_fmt`.
    fn commit_format_field(&mut self, id: &FieldId, start: usize) {
        if self.fmt_tracker.commit(self.indiv_buf, id.dict_idx(), start, false) {
            self.n_fmt = self.n_fmt.saturating_add(1);
        } else {
            warn_twice(self.editing, id);
        }
    }
}

// r[impl vcf_passthrough.edit+2]
/// A field encoded twice into a new record is likely a bug; in an edit of a
/// copy, replacing it is what the caller asked for.
fn warn_twice(editing: bool, id: &FieldId) {
    if !editing {
        tracing::warn!(field = id.name(), "field encoded twice; overwriting previous value");
    }
}

// r[impl record_encoder.format_dedup+2]
impl VcfEncoderFields<'_> {
    /// Open a FORMAT field: append the `:` separator to every sample buffer when
    /// this is not the first field. The caller then writes each sample's value,
    /// followed by [`commit_format_field`](Self::commit_format_field).
    fn begin_format_field(&mut self) {
        if !self.fmt_keys.is_empty() {
            for buf in self.sample_bufs.iter_mut() {
                buf.push(b':');
            }
        }
    }

    /// Close a FORMAT field. A new key is appended to `fmt_keys`. A duplicate
    /// keeps its original column position: the freshly-appended value (now the
    /// last colon-field in each sample buffer) is moved over the previous
    /// occurrence in place, preserving key order to match the BCF/htslib path.
    fn commit_format_field(&mut self, id: &FieldId) {
        if self.fmt_seen.insert(id.dict_idx()) {
            // A clone shares or copies the header's name; no re-validation.
            self.fmt_keys.push(id.name.clone());
            return;
        }
        match self.fmt_keys.iter().position(|k| k.as_str() == id.name()) {
            Some(idx) => {
                warn_twice(self.editing, id);
                // `+1`: the appended occupies a transient extra column we now fold away.
                let total = self.fmt_keys.len().saturating_add(1);
                for buf in self.sample_bufs.iter_mut() {
                    move_last_colon_field_to(buf, idx, total);
                }
            }
            None => self.fmt_keys.push(id.name.clone()),
        }
    }
}

/// Move the last colon-delimited field of `buf` over the `target_idx`-th field,
/// in place, keeping every other field's position. `total_fields` counts the
/// fields currently present (including the appended last one).
///
/// Example: `buf = "0/1:45:99"`, `target_idx = 1`, `total_fields = 3` → the last
/// field `99` replaces field 1, yielding `"0/1:99"`.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "target_idx < total_fields-1 so the target lies strictly before the final ':'; \
              last_sep+1 is in bounds because a duplicate implies >= 2 colon-fields"
)]
fn move_last_colon_field_to(buf: &mut Vec<u8>, target_idx: usize, total_fields: usize) {
    if total_fields <= 1 || target_idx + 1 >= total_fields {
        return;
    }
    let Some(last_sep) = buf.iter().rposition(|&b| b == b':') else {
        return;
    };
    // Capture the appended value (everything after the final ':'), then drop it.
    let value: SmallVec<u8, 16> = buf[last_sep + 1..].iter().copied().collect();
    buf.truncate(last_sep);
    // Locate the target field's value range among the remaining fields and splice.
    let (ts, te) = nth_colon_field(buf, target_idx);
    buf.splice(ts..te, value.iter().copied());
}

/// Byte range of the `n`-th colon-delimited field's value (separators excluded).
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "i < buf.len(); start/idx only advance past a found ':' so stay in bounds"
)]
fn nth_colon_field(buf: &[u8], n: usize) -> (usize, usize) {
    let mut idx = 0;
    let mut start = 0;
    for (i, &byte) in buf.iter().enumerate() {
        if byte == b':' {
            if idx == n {
                return (start, i);
            }
            idx += 1;
            start = i + 1;
        }
    }
    (start, buf.len())
}

// ── Filtered: InfoEncoder ──────────────────────────────────────────────

// r[impl record_encoder.info_encoder+2]
// r[impl record_encoder.bcf_encoding]
// r[impl record_encoder.vcf_encoding]
impl<S: InfoPhase> InfoEncoder for RecordEncoder<'_, S> {
    // r[impl record_encoder.reserved_rejection+2]
    fn info_int(&mut self, id: &FieldId, value: i32) -> Result<(), VcfError> {
        self.check_key(id)?;
        reject_reserved_ints(id.name(), [value])?;
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.info_target(id.dict_idx()) {
            push_ints(buf, [Some(value)]);
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                vcf.buf.extend_from_slice(b.format(value).as_bytes());
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn info_float(&mut self, id: &FieldId, value: f32) -> Result<(), VcfError> {
        self.check_key(id)?;
        reject_reserved_floats(id.name(), [value])?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                // r[impl bcf_writer.smallest_int_type]
                let tc = value.scalar_type_code();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 1, tc);
                value.encode_bcf_as(enc.shared_buf, tc);
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                // r[impl vcf_writer.float_precision]
                write_float_g(vcf.buf, value)
                    .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn info_ints(&mut self, id: &FieldId, values: &[i32]) -> Result<(), VcfError> {
        self.check_key(id)?;
        id.check_count(self.n_allele_count(), None, values.len())?;
        reject_reserved_ints(id.name(), values.iter().copied())?;
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.info_target(id.dict_idx()) {
            push_ints(buf, values.iter().copied().map(Some));
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
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
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn info_floats(&mut self, id: &FieldId, values: &[f32]) -> Result<(), VcfError> {
        self.check_key(id)?;
        id.check_count(self.n_allele_count(), None, values.len())?;
        reject_reserved_floats(id.name(), values.iter().copied())?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_array_values(enc.shared_buf, values);
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b',');
                    }
                    // r[impl vcf_writer.float_precision]
                    write_float_g(vcf.buf, *v)
                        .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
                }
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    fn info_flag(&mut self, id: &FieldId) -> Result<(), VcfError> {
        self.check_key(id)?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                // r[impl bcf_writer.flag_encoding]
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, 0, BCF_BT_NULL);
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    fn info_string(&mut self, id: &FieldId, value: &str) -> Result<(), VcfError> {
        self.info_str(id, value.as_bytes(), Escape::Percent)
    }
    // r[impl record_encoder.string_list]
    fn info_strings(&mut self, id: &FieldId, values: &[&str]) -> Result<(), VcfError> {
        self.check_key(id)?;
        id.check_count(self.n_allele_count(), None, values.len())?;
        let mut joined = Vec::new();
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                joined.push(b',');
            }
            percent_encode_into(&mut joined, v.as_bytes());
        }
        self.info_str(id, if values.is_empty() { b"." } else { &joined }, Escape::Encoded)
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]) -> Result<(), VcfError> {
        self.check_key(id)?;
        id.check_count(self.n_allele_count(), None, values.len())?;
        reject_reserved_ints(id.name(), values.iter().filter_map(|v| *v))?;
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.info_target(id.dict_idx()) {
            push_ints(buf, values.iter().copied());
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                // r[impl bcf_writer.missing_sentinels]
                let typ = smallest_int_type_iter(values.iter().filter_map(|v| *v));
                encode_type_byte(enc.shared_buf, values.len(), typ);
                for &v in values {
                    encode_int_value_or_missing(enc.shared_buf, v, typ);
                }
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
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
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
    // r[impl vcf_passthrough.opt_arrays]
    // r[impl record_encoder.reserved_rejection+2]
    fn info_float_opts(&mut self, id: &FieldId, values: &[Option<f32>]) -> Result<(), VcfError> {
        self.check_key(id)?;
        id.check_count(self.n_allele_count(), None, values.len())?;
        reject_reserved_floats(id.name(), values.iter().filter_map(|v| *v))?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_type_byte(enc.shared_buf, values.len(), BCF_BT_FLOAT);
                for v in values {
                    // r[impl bcf_writer.missing_sentinels]
                    match v {
                        Some(v) => v.encode_bcf_as(enc.shared_buf, BCF_BT_FLOAT),
                        None => f32::encode_missing(enc.shared_buf),
                    }
                }
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        vcf.buf.push(b',');
                    }
                    match v {
                        // r[impl vcf_writer.float_precision]
                        Some(v) => write_float_g(vcf.buf, *v)
                            .map_err(|source| VcfError::FailedToWriteFormattedString { source })?,
                        // r[impl vcf_writer.missing_dot]
                        None => vcf.buf.push(b'.'),
                    }
                }
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
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

// r[impl record_encoder.typestate_transitions+4]
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
        self.into_state()
    }
    // r[impl record_encoder.emit]
    // r[impl record_encoder.emit_no_samples]
    pub fn emit(self) -> Result<(), VcfError> {
        emit_inner(self.inner, self.record)
    }
}

/// How a string value's bytes are written into a VCF text line.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Escape {
    /// Every character VCF gives a meaning (`r[vcf_writer.percent_encoding]`):
    /// the colon is the FORMAT subfield delimiter, so a raw `:` (or `;`, TAB,
    /// …) would corrupt the line and disagree with the length-prefixed BCF
    /// encoding.
    Percent,
    // r[impl vcf_passthrough.copy.values]
    /// A copied value, which the reader returns as stored: only TAB, LF, CR
    /// and the section's separator (`;` for INFO, `:` for FORMAT), which
    /// would break the line.
    Structure(u8),
    /// Already percent-encoded (a list joined from encoded elements).
    Encoded,
}

impl Escape {
    // r[impl bcf_writer.string_escape]
    /// The bytes BCF stores: what VCF text would hold, so that a BCF and a
    /// VCF of the same record read back the same text. A copied value is
    /// already as stored.
    fn stored(self, bytes: &[u8]) -> std::borrow::Cow<'_, [u8]> {
        match self {
            Self::Percent if bytes.iter().any(|&b| needs_percent(b)) => {
                let mut buf = Vec::with_capacity(bytes.len().saturating_add(8));
                percent_encode_into(&mut buf, bytes);
                std::borrow::Cow::Owned(buf)
            }
            _ => std::borrow::Cow::Borrowed(bytes),
        }
    }

    fn write(self, buf: &mut Vec<u8>, bytes: &[u8]) {
        match self {
            Self::Encoded => buf.extend_from_slice(bytes),
            // r[impl vcf_writer.percent_encoding]
            Self::Percent => percent_encode_into(buf, bytes),
            Self::Structure(sep) => {
                for &b in bytes {
                    if matches!(b, b'\t' | b'\n' | b'\r') || b == sep {
                        push_percent(buf, b);
                    } else {
                        buf.push(b);
                    }
                }
            }
        }
    }
}

impl<S: InfoPhase> RecordEncoder<'_, S> {
    pub(crate) fn info_str(
        &mut self,
        id: &FieldId,
        value: &[u8],
        escape: Escape,
    ) -> Result<(), VcfError> {
        self.check_key(id)?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.shared_buf.len();
                encode_typed_int_key(enc.shared_buf, id.dict_idx());
                encode_typed_string(enc.shared_buf, &escape.stored(value));
                // r[impl bcf_encoder.info_counting]
                enc.commit_info_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                let start = vcf.begin_info_field();
                vcf.buf.extend_from_slice(id.name().as_bytes());
                vcf.buf.push(b'=');
                escape.write(vcf.buf, value);
                vcf.commit_info_field(id, start);
            }
        }
        Ok(())
    }
}

#[allow(clippy::indexing_slicing, reason = "i < values.len() == n_samples() == sample_bufs.len()")]
impl<S: FormatPhase> RecordEncoder<'_, S> {
    fn n_samples_inner(&self) -> usize {
        match &self.inner {
            EncoderInner::Bcf(e) => e.n_sample as usize,
            EncoderInner::Vcf(v) => v.n_samples as usize,
        }
    }

    pub(crate) fn format_str(
        &mut self,
        id: &FieldId,
        values: &[&str],
        escape: Escape,
    ) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(values.len(), self.n_samples_inner(), "one string per sample");
        // A missing (empty) sample is stored as the literal `.` char (htslib's BCF
        // string-missing convention), so the column width must leave room for it.
        // Note: a *non-empty* value of exactly "." is indistinguishable from missing
        // on read-back — that ambiguity is inherent to BCF/VCF, not specific to us.
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let stored: SmallVec<std::borrow::Cow<'_, [u8]>, 8> =
                    values.iter().map(|s| escape.stored(s.as_bytes())).collect();
                let width = stored
                    .iter()
                    .map(|s| if s.is_empty() { 1 } else { s.len() })
                    .max()
                    .unwrap_or(0);
                let start = enc.indiv_buf.len();
                // BCF FORMAT char vectors are fixed-width: every sample stores `width`
                // bytes, shorter strings NUL-padded (htslib bcf_enc_vchar).
                // r[impl bcf_writer.string_encoding]
                // r[impl bcf_writer.indiv_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, width, BCF_BT_CHAR);
                for s in &stored {
                    // Match the VCF text path's `.`: an empty string is stored as the
                    // literal `.` char, which bcftools/noodles render as missing.
                    let bytes: &[u8] = if s.is_empty() { b"." } else { s };
                    enc.indiv_buf.extend_from_slice(bytes);
                    for _ in bytes.len()..width {
                        enc.indiv_buf.push(0);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                for (i, s) in values.iter().enumerate() {
                    if s.is_empty() {
                        // r[impl vcf_writer.missing_dot]
                        vcf.sample_bufs[i].push(b'.');
                    } else {
                        escape.write(&mut vcf.sample_bufs[i], s.as_bytes());
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
}

// ── WithSamples: FormatEncoder ─────────────────────────────────────────

// r[impl record_encoder.format_encoder+2]
// `sample_bufs[i]` is indexed with `i` in `0..values.len()`; each method asserts
// `values.len() == n_samples()` up front, and `sample_bufs` is sized to
// `n_samples()` in `begin_samples()`. The typed `FormatKey::encode` wrappers
// return `SampleCountMismatch` for the same condition; the `debug_assert`s catch
// direct (untyped) trait calls that bypass that check.
#[allow(clippy::indexing_slicing, reason = "i < values.len() == n_samples() == sample_bufs.len()")]
impl<S: FormatPhase> FormatEncoder for RecordEncoder<'_, S> {
    // r[impl record_encoder.mixed_ploidy]
    fn format_gt(&mut self, id: &FieldId, gts: &[Genotype]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(gts.len(), self.n_samples(), "one GT per sample");

        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
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
                // The column is as wide as the largest ploidy; a genotype with
                // no alleles still takes one slot, its missing call.
                let width = gts.iter().map(|g| g.alleles.len()).max().unwrap_or(0).max(1);
                encode_type_byte(enc.indiv_buf, width, typ);
                for gt in gts {
                    let ploidy = gt.alleles.len();
                    for (i, allele_opt) in gt.alleles.iter().enumerate() {
                        // r[impl vcf_record.gt_first_phase+3]
                        // The first allele carries a phase bit too, and from
                        // VCF 4.4 it is read rather than ignored: the spec
                        // makes the leading indicator "`/` if any phasing
                        // indicators are `/` and `|` otherwise". A haploid call
                        // has no separator, so it is phased.
                        let phased = if i == 0 {
                            gt.phased.iter().all(|p| *p)
                        } else {
                            gt.phased.get(i.saturating_sub(1)).copied().unwrap_or(false)
                        };
                        let encoded: i32 = match allele_opt {
                            Some(idx) => {
                                (i32::from(*idx).saturating_add(1) << 1) | i32::from(phased)
                            }
                            // A missing allele keeps its phase bit, as htslib
                            // encodes `.|.`; a missing haploid call is `0x00`.
                            None => i32::from(phased && ploidy > 1),
                        };
                        encode_int_as(enc.indiv_buf, encoded, typ);
                    }
                    if ploidy == 0 {
                        encode_int_as(enc.indiv_buf, 0, typ);
                    }
                    // r[impl bcf_writer.end_of_vector]
                    for _ in ploidy.max(1)..width {
                        encode_int_eov(enc.indiv_buf, typ);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                // r[impl vcf_writer.genotype_serialization+2]
                vcf.begin_format_field();
                let mut b = itoa::Buffer::new();
                for (si, gt) in gts.iter().enumerate() {
                    if gt.alleles.is_empty() {
                        vcf.sample_bufs[si].push(b'.');
                    }
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
                vcf.commit_format_field(id);
            }
        }
        self.record.ploidies.clear();
        self.record.ploidies.extend(gts.iter().map(|g| g.alleles.len()));
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn format_int(&mut self, id: &FieldId, values: &[i32]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(values.len(), self.n_samples(), "one value per sample");
        reject_reserved_ints(id.name(), values.iter().copied())?;
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.format_target(id.dict_idx()) {
            for &v in values {
                buf.push(b'\t');
                push_ints(buf, [Some(v)]);
            }
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.indiv_field_major]
                let tc = smallest_int_type(values);
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    v.encode_bcf_as(enc.indiv_buf, tc);
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                // r[impl vcf_writer.integer_format]
                let mut b = itoa::Buffer::new();
                for (i, &v) in values.iter().enumerate() {
                    vcf.sample_bufs[i].extend_from_slice(b.format(v).as_bytes());
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn format_float(&mut self, id: &FieldId, values: &[f32]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(values.len(), self.n_samples(), "one value per sample");
        reject_reserved_floats(id.name(), values.iter().copied())?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
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
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                // r[impl vcf_writer.float_precision]
                for (i, &v) in values.iter().enumerate() {
                    write_float_g(&mut vcf.sample_bufs[i], v)
                        .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.format_scalar_opt]
    // r[impl record_encoder.reserved_rejection+2]
    fn format_int_opt(&mut self, id: &FieldId, values: &[Option<i32>]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(values.len(), self.n_samples(), "one value per sample");
        reject_reserved_ints(id.name(), values.iter().filter_map(|v| *v))?;
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.format_target(id.dict_idx()) {
            for &v in values {
                buf.push(b'\t');
                push_ints(buf, [v]);
            }
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.missing_sentinels]
                // Every sample missing leaves the iterator empty: int8, as htslib picks.
                let tc = smallest_int_type_iter(values.iter().filter_map(|v| *v));
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, tc);
                for &v in values {
                    encode_int_value_or_missing(enc.indiv_buf, v, tc);
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                let mut b = itoa::Buffer::new();
                for (i, v) in values.iter().enumerate() {
                    match v {
                        // r[impl vcf_writer.integer_format]
                        Some(v) => vcf.sample_bufs[i].extend_from_slice(b.format(*v).as_bytes()),
                        // r[impl vcf_writer.missing_dot]
                        None => vcf.sample_bufs[i].push(b'.'),
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.format_scalar_opt]
    // r[impl record_encoder.reserved_rejection+2]
    fn format_float_opt(&mut self, id: &FieldId, values: &[Option<f32>]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(values.len(), self.n_samples(), "one value per sample");
        reject_reserved_floats(id.name(), values.iter().filter_map(|v| *v))?;
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, 1, BCF_BT_FLOAT);
                for &v in values {
                    // r[impl bcf_writer.missing_sentinels]
                    match v {
                        Some(v) => v.encode_bcf_as(enc.indiv_buf, BCF_BT_FLOAT),
                        None => f32::encode_missing(enc.indiv_buf),
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                for (i, v) in values.iter().enumerate() {
                    match v {
                        // r[impl vcf_writer.float_precision]
                        Some(v) => write_float_g(&mut vcf.sample_bufs[i], *v)
                            .map_err(|source| VcfError::FailedToWriteFormattedString { source })?,
                        // r[impl vcf_writer.missing_dot]
                        None => vcf.sample_bufs[i].push(b'.'),
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn format_ints(&mut self, id: &FieldId, per_sample: &[&[i32]]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(per_sample.len(), self.n_samples(), "one slice per sample");
        self.check_sample_counts(id, per_sample)?;
        reject_reserved_ints(id.name(), per_sample.iter().flat_map(|s| s.iter().copied()))?;
        let width = per_sample.iter().map(|s| s.len()).max().unwrap_or(0);
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.format_target(id.dict_idx()) {
            for s in per_sample {
                buf.push(b'\t');
                push_ints(buf, s.iter().copied().map(Some));
            }
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // BCF stores one integer type for the whole FORMAT column, so the
                // type must fit every value across every sample, not just one sample.
                // r[impl bcf_writer.smallest_int_type]
                let typ = smallest_int_type_iter(per_sample.iter().flat_map(|s| s.iter().copied()));
                // r[impl bcf_writer.indiv_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, width, typ);
                for s in per_sample {
                    if s.is_empty() {
                        // Match the VCF text path's `.`: a sample with no values is
                        // the missing sentinel followed by end-of-vector padding,
                        // which bcftools/noodles render as missing (not an empty cell).
                        if width > 0 {
                            encode_int_missing(enc.indiv_buf, typ);
                            for _ in 1..width {
                                encode_int_eov(enc.indiv_buf, typ);
                            }
                        }
                        continue;
                    }
                    for &v in *s {
                        encode_int_as(enc.indiv_buf, v, typ);
                    }
                    // r[impl bcf_writer.end_of_vector]
                    for _ in s.len()..width {
                        encode_int_eov(enc.indiv_buf, typ);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                let mut b = itoa::Buffer::new();
                for (i, s) in per_sample.iter().enumerate() {
                    if s.is_empty() {
                        // r[impl vcf_writer.missing_dot]
                        vcf.sample_bufs[i].push(b'.');
                        continue;
                    }
                    for (j, &v) in s.iter().enumerate() {
                        if j > 0 {
                            vcf.sample_bufs[i].push(b',');
                        }
                        // r[impl vcf_writer.integer_format]
                        vcf.sample_bufs[i].extend_from_slice(b.format(v).as_bytes());
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn format_floats(&mut self, id: &FieldId, per_sample: &[&[f32]]) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(per_sample.len(), self.n_samples(), "one slice per sample");
        self.check_sample_counts(id, per_sample)?;
        reject_reserved_floats(id.name(), per_sample.iter().flat_map(|s| s.iter().copied()))?;
        let width = per_sample.iter().map(|s| s.len()).max().unwrap_or(0);
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // r[impl bcf_writer.indiv_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, width, BCF_BT_FLOAT);
                for s in per_sample {
                    if s.is_empty() {
                        // Match the VCF text path's `.`: a sample with no values is
                        // the missing sentinel followed by end-of-vector padding.
                        if width > 0 {
                            f32::encode_missing(enc.indiv_buf);
                            for _ in 1..width {
                                f32::encode_end_of_vector(enc.indiv_buf);
                            }
                        }
                        continue;
                    }
                    for &v in *s {
                        v.encode_bcf_as(enc.indiv_buf, BCF_BT_FLOAT);
                    }
                    // r[impl bcf_writer.end_of_vector]
                    for _ in s.len()..width {
                        f32::encode_end_of_vector(enc.indiv_buf);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                for (i, s) in per_sample.iter().enumerate() {
                    if s.is_empty() {
                        // r[impl vcf_writer.missing_dot]
                        vcf.sample_bufs[i].push(b'.');
                        continue;
                    }
                    for (j, &v) in s.iter().enumerate() {
                        if j > 0 {
                            vcf.sample_bufs[i].push(b',');
                        }
                        // r[impl vcf_writer.float_precision]
                        write_float_g(&mut vcf.sample_bufs[i], v)
                            .map_err(|source| VcfError::FailedToWriteFormattedString { source })?;
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl record_encoder.reserved_rejection+2]
    fn format_float_opts(
        &mut self,
        id: &FieldId,
        per_sample: &[&[Option<f32>]],
    ) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(per_sample.len(), self.n_samples(), "one slice per sample");
        self.check_sample_counts(id, per_sample)?;
        reject_reserved_floats(
            id.name(),
            per_sample.iter().flat_map(|s| s.iter().filter_map(|v| *v)),
        )?;
        let width = per_sample.iter().map(|s| s.len()).max().unwrap_or(0);
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // r[impl bcf_writer.indiv_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, width, BCF_BT_FLOAT);
                for s in per_sample {
                    if s.is_empty() && width > 0 {
                        f32::encode_missing(enc.indiv_buf);
                        for _ in 1..width {
                            f32::encode_end_of_vector(enc.indiv_buf);
                        }
                        continue;
                    }
                    for v in *s {
                        // r[impl bcf_writer.missing_sentinels]
                        match v {
                            Some(v) => v.encode_bcf_as(enc.indiv_buf, BCF_BT_FLOAT),
                            None => f32::encode_missing(enc.indiv_buf),
                        }
                    }
                    // r[impl bcf_writer.end_of_vector]
                    for _ in s.len()..width {
                        f32::encode_end_of_vector(enc.indiv_buf);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                for (i, s) in per_sample.iter().enumerate() {
                    if s.is_empty() {
                        // r[impl vcf_writer.missing_dot]
                        vcf.sample_bufs[i].push(b'.');
                        continue;
                    }
                    for (j, v) in s.iter().enumerate() {
                        if j > 0 {
                            vcf.sample_bufs[i].push(b',');
                        }
                        match v {
                            // r[impl vcf_writer.float_precision]
                            Some(v) => {
                                write_float_g(&mut vcf.sample_bufs[i], *v).map_err(|source| {
                                    VcfError::FailedToWriteFormattedString { source }
                                })?
                            }
                            // r[impl vcf_writer.missing_dot]
                            None => vcf.sample_bufs[i].push(b'.'),
                        }
                    }
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    // r[impl vcf_passthrough.opt_arrays]
    // r[impl record_encoder.reserved_rejection+2]
    fn format_int_opts(
        &mut self,
        id: &FieldId,
        per_sample: &[&[Option<i32>]],
    ) -> Result<(), VcfError> {
        self.check_key(id)?;
        self.inner.check_format_room(id)?;
        debug_assert_eq!(per_sample.len(), self.n_samples(), "one slice per sample");
        self.check_sample_counts(id, per_sample)?;
        let concrete = || per_sample.iter().flat_map(|s| s.iter().filter_map(|v| *v));
        reject_reserved_ints(id.name(), concrete())?;
        let width = per_sample.iter().map(|s| s.len()).max().unwrap_or(0);
        // r[impl record_encoder.rlen]
        if let Some(buf) = self.record.rlen.format_target(id.dict_idx()) {
            for s in per_sample {
                buf.push(b'\t');
                push_ints(buf, s.iter().copied());
            }
        }
        match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let start = enc.indiv_buf.len();
                // r[impl bcf_writer.smallest_int_type]
                // r[impl bcf_writer.missing_sentinels]
                let typ = smallest_int_type_iter(concrete());
                // r[impl bcf_writer.indiv_field_major]
                encode_typed_int_key(enc.indiv_buf, id.dict_idx());
                encode_type_byte(enc.indiv_buf, width, typ);
                for s in per_sample {
                    if s.is_empty() {
                        if width > 0 {
                            encode_int_missing(enc.indiv_buf, typ);
                            for _ in 1..width {
                                encode_int_eov(enc.indiv_buf, typ);
                            }
                        }
                        continue;
                    }
                    for &v in *s {
                        encode_int_value_or_missing(enc.indiv_buf, v, typ);
                    }
                    // r[impl bcf_writer.end_of_vector]
                    for _ in s.len()..width {
                        encode_int_eov(enc.indiv_buf, typ);
                    }
                }
                enc.commit_format_field(id, start);
            }
            EncoderInner::Vcf(vcf) => {
                vcf.begin_format_field();
                for (i, s) in per_sample.iter().enumerate() {
                    if s.is_empty() {
                        // r[impl vcf_writer.missing_dot]
                        vcf.sample_bufs[i].push(b'.');
                        continue;
                    }
                    push_ints(&mut vcf.sample_bufs[i], s.iter().copied());
                }
                vcf.commit_format_field(id);
            }
        }
        Ok(())
    }
    fn format_string(&mut self, id: &FieldId, values: &[&str]) -> Result<(), VcfError> {
        self.format_str(id, values, Escape::Percent)
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
        emit_inner(self.inner, self.record)
    }
}

// ── Edit: removal and emit ─────────────────────────────────────────────

impl RecordEncoder<'_, Edit> {
    // r[impl vcf_passthrough.edit.remove]
    /// Remove INFO field `key` from the record, keeping the others in order;
    /// whether it was there.
    ///
    /// # Errors
    /// [`VcfError::ForeignKey`] for a key of another header.
    pub fn remove_info<V>(&mut self, key: &InfoKey<V>) -> Result<bool, VcfError> {
        let id = &key.0;
        self.check_key(id)?;
        self.record.rlen.forget_info(id.dict_idx());
        Ok(match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let removed = enc.info_tracker.remove(enc.shared_buf, id.dict_idx(), false);
                enc.n_info = enc.n_info.saturating_sub(u16::from(removed));
                removed
            }
            EncoderInner::Vcf(vcf) => {
                let removed = vcf.info_tracker.remove(vcf.buf, id.dict_idx(), true);
                vcf.info_count = vcf.info_count.saturating_sub(u16::from(removed));
                removed
            }
        })
    }

    // r[impl vcf_passthrough.edit.remove]
    /// Remove FORMAT field `key` from the record, for every sample, keeping
    /// the others in order; whether it was there.
    ///
    /// # Errors
    /// [`VcfError::ForeignKey`] for a key of another header.
    pub fn remove_format<V>(&mut self, key: &FormatKey<V>) -> Result<bool, VcfError> {
        let id = &key.0;
        self.check_key(id)?;
        self.record.rlen.forget_format(id.dict_idx());
        let removed = match &mut self.inner {
            EncoderInner::Bcf(enc) => {
                let removed = enc.fmt_tracker.remove(enc.indiv_buf, id.dict_idx(), false);
                enc.n_fmt = enc.n_fmt.saturating_sub(u8::from(removed));
                removed
            }
            EncoderInner::Vcf(vcf) => {
                let position = vcf.fmt_keys.iter().position(|k| k.as_str() == id.name());
                match position.filter(|_| vcf.fmt_seen.remove(id.dict_idx())) {
                    Some(idx) => {
                        let total = vcf.fmt_keys.len();
                        vcf.fmt_keys.remove(idx);
                        for buf in vcf.sample_bufs.iter_mut() {
                            remove_colon_field(buf, idx, total);
                        }
                        true
                    }
                    None => false,
                }
            }
        };
        if removed && id.name() == "GT" {
            self.record.ploidies.clear();
        }
        Ok(removed)
    }

    // r[impl vcf_passthrough.edit.emit]
    /// Write the record (see [`RecordEncoder::<WithSamples>::emit`]).
    pub fn emit(self) -> Result<(), VcfError> {
        emit_inner(self.inner, self.record)
    }
}

/// Remove the `idx`-th of `total` colon-delimited fields of `buf` with one
/// of its separators: the one before it, or for the first field the one after.
fn remove_colon_field(buf: &mut Vec<u8>, idx: usize, total: usize) {
    if total <= 1 {
        buf.clear();
        return;
    }
    let (start, end) = nth_colon_field(buf, idx);
    let range = if idx == 0 {
        start..end.saturating_add(1).min(buf.len())
    } else {
        start.saturating_sub(1)..end
    };
    buf.drain(range);
}

// ── Emit ───────────────────────────────────────────────────────────────

// r[impl record_encoder.rlen]
fn emit_inner(mut inner: EncoderInner<'_>, record: &mut RecordState) -> Result<(), VcfError> {
    match &mut inner {
        EncoderInner::Bcf(enc) => {
            let rlen = record.rlen.compute(i64::from(enc.pos_0based), i64::from(enc.rlen));
            enc.set_rlen(i32::try_from(rlen).map_err(|_| overflow("rlen", rlen, "i32"))?);
            enc.emit_inner()
        }
        EncoderInner::Vcf(vcf) => {
            let rlen = record.rlen.compute(i64::from(vcf.pos_0based), i64::from(vcf.rlen));
            vcf.rlen = u32::try_from(rlen).map_err(|_| overflow("rlen", rlen, "u32"))?;
            vcf_emit(vcf)
        }
    }
}

fn overflow(field: &'static str, value: i64, target_type: &'static str) -> VcfError {
    VcfError::ValueOverflow { field, value: value.unsigned_abs(), target_type }
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
    // r[impl vcf_writer.format_serialization+2]
    if vcf.fmt_keys.is_empty() && vcf.n_samples > 0 {
        vcf.buf.extend_from_slice(b"\t.");
        for _ in 0..vcf.n_samples {
            vcf.buf.extend_from_slice(b"\t.");
        }
    }
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
