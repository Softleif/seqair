//! [`VcfFile`]: a VCF/BCF file opened once — format, header, index and the
//! tables every reader of it shares (`r[vcf_reader.file+2]`).

use super::error::{VcfFormat, VcfReadError};
use super::fields::Fields;
use super::keys::FieldTable;
use super::query::Partition;
use super::region_index::VcfIndex;
use super::rlen::RlenKeys;
use super::stream::ReaderBuilder;
use crate::bam::bgzf::{BgzfStream, VirtualOffset, find_bsize};
use crate::vcf::header::{MAX_HEADER_TEXT, VcfHeader, VcfVersion};
use crate::vcf::record_encoder::{ContigId, FilterId};
use seqair_types::SmolStr;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, RwLock};
use tracing::warn;

const BCF_MAGIC: &[u8; 5] = b"BCF\x02\x02";
const BCF_21_MAGIC: &[u8; 5] = b"BCF\x02\x01";
const VCF_MAGIC: &[u8] = b"##fileformat=";
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Bytes each read from a file asks for: BGZF blocks are at most 64 KiB, so
/// this keeps the block parser's small reads off the system call path.
pub(crate) const FILE_BUFFER: usize = 128 * 1024;

// r[impl vcf_reader.file+2]
// r[impl vcf_reader.file.shared+2]
/// A VCF or BCF file: its header, index and contig table, parsed once and
/// shared by every reader opened from it.
///
/// `Send + Sync + Clone`; a clone shares the parsed state and re-reads
/// nothing, so a rayon task can open its own reader from a clone.
#[derive(Clone)]
pub struct VcfFile(pub(crate) Arc<FileShared>);

impl std::fmt::Debug for VcfFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// Where a file's bytes come from.
#[derive(Debug)]
pub(crate) enum Source {
    Path(PathBuf),
    Bytes(Arc<[u8]>),
    /// A [`VcfStream`](super::VcfStream)'s input, read once by its reader.
    Stream,
}

/// A reader's input: a file handle, a view of in-memory bytes, or a stream
/// that cannot seek.
pub(crate) enum Input {
    File(BufReader<File>),
    Bytes(Cursor<Arc<[u8]>>),
    Stream(BufReader<Box<dyn Read + Send>>),
}

impl std::fmt::Debug for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::File(file) => f.debug_tuple("File").field(file).finish(),
            Self::Bytes(_) => f.write_str("Bytes"),
            Self::Stream(_) => f.write_str("Stream"),
        }
    }
}

impl Read for Input {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::File(f) => f.read(buf),
            Self::Bytes(b) => b.read(buf),
            Self::Stream(s) => s.read(buf),
        }
    }
}

impl BufRead for Input {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            Self::File(f) => f.fill_buf(),
            Self::Bytes(b) => b.fill_buf(),
            Self::Stream(s) => s.fill_buf(),
        }
    }

    fn consume(&mut self, n: usize) {
        match self {
            Self::File(f) => f.consume(n),
            Self::Bytes(b) => b.consume(n),
            Self::Stream(s) => s.consume(n),
        }
    }
}

impl Seek for Input {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::File(f) => f.seek(pos),
            Self::Bytes(b) => b.seek(pos),
            Self::Stream(_) => Err(std::io::ErrorKind::Unsupported.into()),
        }
    }
}

/// Where a file's records start.
#[derive(Debug, Clone, Copy)]
pub(crate) enum DataStart {
    /// At this virtual offset.
    At(VirtualOffset),
    /// This many (decompressed) bytes from the start: in a plain file, or
    /// at an offset the stream could not name.
    After(u64),
}

/// What reading a header finds: the format, the header, where records start,
/// and (VCF text) how many lines the header took.
pub(crate) type HeaderStart = (VcfFormat, VcfHeader, DataStart, u64);

/// Everything a [`VcfFile`]'s readers share.
pub(crate) struct FileShared {
    source: Source,
    pub(crate) format: VcfFormat,
    pub(crate) header: VcfHeader,
    pub(crate) data_start: DataStart,
    pub(crate) index: Option<VcfIndex>,
    pub(crate) contigs: ContigTable,
    /// VCF text: the header's lines, so a record's line number is known.
    pub(crate) header_lines: u64,
    /// FILTER ids by dictionary index; `None` where the index names no FILTER.
    pub(crate) filters: Vec<Option<FilterId>>,
    /// The INFO/FORMAT keys records can name, by id and (VCF text) by name.
    pub(crate) fields: FieldTable,
    /// Whether a GT's first phase bit is read (VCF 4.4 and later).
    pub(crate) first_phase_is_read: bool,
    /// VCF text: which of END, SVLEN and LEN set the reference length.
    pub(crate) rlen_keys: RlenKeys,
    warned_duplicate: AtomicBool,
    warned_undeclared: AtomicBool,
}

// r[impl vcf_reader.debug]
/// The file, not its tables: every record component holds one.
impl std::fmt::Debug for FileShared {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut d = f.debug_struct("VcfFile");
        match &self.source {
            Source::Path(path) => d.field("path", path),
            Source::Bytes(bytes) => d.field("bytes", &bytes.len()),
            Source::Stream => d.field("stream", &true),
        };
        d.field("format", &self.format).finish_non_exhaustive()
    }
}

impl Source {
    fn open(&self) -> Result<Input, VcfReadError> {
        match self {
            Source::Path(path) => File::open(path)
                .map(|f| Input::File(BufReader::with_capacity(FILE_BUFFER, f)))
                .map_err(|source| VcfReadError::Open { path: path.clone(), source }),
            Source::Bytes(bytes) => Ok(Input::Bytes(Cursor::new(Arc::clone(bytes)))),
            // Its one reader holds the input; nothing reopens it but a
            // query, which fails first for want of an index.
            Source::Stream => Err(VcfReadError::NoIndex),
        }
    }
}

impl FileShared {
    /// A fresh handle on the file's bytes, at their start.
    pub(crate) fn open_input(&self) -> Result<Input, VcfReadError> {
        self.source.open()
    }

    /// The file's path, unless it was read from memory.
    pub(crate) fn path(&self) -> Option<&Path> {
        match &self.source {
            Source::Path(path) => Some(path),
            Source::Bytes(_) | Source::Stream => None,
        }
    }

    /// Whether records are VCF text (plain or BGZF).
    pub(crate) fn is_text(&self) -> bool {
        self.format != VcfFormat::Bcf
    }

    // r[impl vcf_reader.undeclared.filter]
    /// Log, once per file, that a record names a FILTER the header does not
    /// declare.
    pub(crate) fn warn_undeclared_filter(&self, name: &[u8]) {
        if !self.warned_undeclared.swap(true, Ordering::Relaxed) {
            warn!(
                filter = std::str::from_utf8(name).ok(),
                path = ?self.path(),
                "a record names a FILTER the header does not declare; the record fails it \
                 (further undeclared FILTERs in this file are not reported)"
            );
        }
    }

    // r[impl vcf_values.duplicates+2]
    /// Log, once per file, that a record holds `name` twice.
    pub(crate) fn warn_duplicate(&self, name: &str) {
        if !self.warned_duplicate.swap(true, Ordering::Relaxed) {
            warn!(
                field = name,
                path = ?self.path(),
                "a record holds the same field twice; the first occurrence is read \
                 (further duplicates in this file are not reported)"
            );
        }
    }
}

impl VcfFile {
    /// Open the VCF or BCF file at `path`: detect its format from its bytes,
    /// parse its header, and parse the index next to it if there is one.
    ///
    /// # Errors
    /// The file cannot be read, is not a supported format, its header is
    /// malformed, or an index exists and does not parse.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, VcfReadError> {
        let path = path.as_ref();
        let source = Source::Path(path.to_path_buf());
        let start = read_header(source.open()?)?;
        let index = match crate::utils::find_index(path, index_extensions(start.0)) {
            Some(index) => Some(
                std::fs::read(&index)
                    .map_err(|source| VcfReadError::Open { path: index, source })?,
            ),
            None => None,
        };
        Self::new(source, start, index.as_deref())
    }

    /// A file held in memory, with its index's bytes (TBI or CSI) if it has
    /// one. A plain-text VCF has no index: one given for it is not read.
    ///
    /// # Errors
    /// As [`Self::open`].
    pub fn from_bytes(
        data: impl Into<Arc<[u8]>>,
        index: Option<&[u8]>,
    ) -> Result<Self, VcfReadError> {
        let source = Source::Bytes(data.into());
        let start = read_header(source.open()?)?;
        Self::new(source, start, index)
    }

    pub(crate) fn new(
        source: Source,
        (format, header, data_start, header_lines): HeaderStart,
        index: Option<&[u8]>,
    ) -> Result<Self, VcfReadError> {
        let mut contigs = ContigTable::new(&header);
        let index = match index {
            Some(bytes) if format != VcfFormat::Vcf => {
                Some(VcfIndex::parse(bytes, &header, &mut contigs)?)
            }
            _ => None,
        };
        let mut fields = FieldTable::new(&header, format != VcfFormat::Bcf);
        if let Source::Path(path) = &source {
            fields.path = Some(path.clone());
        }
        let shared = FileShared {
            rlen_keys: RlenKeys::new(&fields),
            fields,
            contigs,
            filters: filter_table(&header),
            first_phase_is_read: header.version() >= VcfVersion::V4_4,
            header_lines,
            warned_duplicate: AtomicBool::new(false),
            warned_undeclared: AtomicBool::new(false),
            format,
            source,
            header,
            data_start,
            index,
        };
        Ok(Self(Arc::new(shared)))
    }

    // r[impl vcf_reader.passthrough_header]
    /// The file's header, as [`Writer::write_header`](crate::vcf::Writer::write_header) takes it.
    pub fn header(&self) -> &VcfHeader {
        &self.0.header
    }

    /// The file's format, as detected from its content.
    pub fn format(&self) -> VcfFormat {
        self.0.format
    }

    /// Whether an index was found (or given) and parsed.
    pub fn has_index(&self) -> bool {
        self.0.index.is_some()
    }

    // r[impl vcf_reader.partitions+2]
    /// Split the file into about `n` [`Partition`]s of roughly equal
    /// compressed bytes, in file order: whole regions of the indexed
    /// contigs, at least one per contig. Read each with
    /// [`ReaderBuilder::partition`]; every record is in exactly one.
    ///
    /// # Errors
    /// [`VcfReadError::NoIndex`] without an index.
    pub fn partitions(&self, n: NonZeroUsize) -> Result<Vec<Partition>, VcfReadError> {
        super::query::partitions(&self.0, n)
    }

    // r[impl vcf_reader.contig_lookup]
    /// The contig `name`, if the header or the index declares it or a
    /// record named it: whether a query of `name` resolves.
    pub fn contig(&self, name: &str) -> Option<ContigId> {
        self.0.contigs.find(&self.0.header, name)
    }

    // r[impl vcf_reader.reader]
    /// A reader of records into `F`, to configure and [`open`](ReaderBuilder::open).
    pub fn reader<F: Fields>(&self) -> ReaderBuilder<'_, F> {
        ReaderBuilder::new(self)
    }
}

/// The index names of `r[index.locate]` a format's index may have: a
/// plain-text VCF has none.
fn index_extensions(format: VcfFormat) -> &'static [&'static str] {
    match format {
        VcfFormat::Vcf => &[],
        VcfFormat::VcfGz => &["csi", "tbi"],
        VcfFormat::Bcf => &["csi"],
    }
}

/// FILTER ids by dictionary index.
fn filter_table(header: &VcfHeader) -> Vec<Option<FilterId>> {
    let mut table = vec![None; header.string_map().len()];
    for name in header.filters().keys() {
        if let (Some(slot), Ok(id)) =
            (header.string_map().get(name).and_then(|i| table.get_mut(i)), header.filter(name))
        {
            *slot = Some(id);
        }
    }
    table
}

/// How the first bytes of a file are framed.
pub(crate) enum Container {
    Bgzf,
    Gzip,
    Plain,
}

// r[impl vcf_reader.formats]
/// The container `prefix` (the file's first bytes) starts with.
pub(crate) fn container(prefix: &[u8]) -> Container {
    if !prefix.starts_with(&GZIP_MAGIC) {
        return Container::Plain;
    }
    // A gzip member with FEXTRA whose extra field holds BGZF's BC subfield.
    let flags = prefix.get(3).copied().unwrap_or(0);
    let xlen = prefix
        .get(10..)
        .and_then(<[u8]>::first_chunk::<2>)
        .map_or(0, |b| usize::from(u16::from_le_bytes(*b)));
    let extra = prefix.get(12..12usize.saturating_add(xlen));
    match extra {
        Some(extra) if flags & 4 != 0 && find_bsize(extra).is_some() => Container::Bgzf,
        _ => Container::Gzip,
    }
}

/// The format a decompressed (or plain) start of a file names.
pub(crate) fn format_of(start: &[u8], compressed: bool) -> Result<VcfFormat, VcfReadError> {
    if start.starts_with(BCF_MAGIC) {
        return if compressed { Ok(VcfFormat::Bcf) } else { Err(VcfReadError::UncompressedBcf) };
    }
    if start.starts_with(BCF_21_MAGIC) {
        return Err(VcfReadError::Bcf21);
    }
    if start.starts_with(VCF_MAGIC) {
        return Ok(if compressed { VcfFormat::VcfGz } else { VcfFormat::Vcf });
    }
    Err(VcfReadError::UnknownFormat)
}

/// Detect the format, read the header, and find where records start and
/// (VCF text) how many lines the header took.
fn read_header(mut input: Input) -> Result<HeaderStart, VcfReadError> {
    let mut prefix = Vec::with_capacity(64);
    (&mut input).take(64).read_to_end(&mut prefix).map_err(read_failed)?;
    input.seek(SeekFrom::Start(0)).map_err(read_failed)?;
    match container(&prefix) {
        Container::Gzip => Err(VcfReadError::NotBgzf),
        Container::Plain => {
            let format = format_of(&prefix, false)?;
            let (header, len, lines) = read_text_header(&mut input)?;
            Ok((format, header, DataStart::After(len), lines))
        }
        Container::Bgzf => {
            let mut stream = BgzfStream::new(input);
            let format = format_of(stream.fill()?, true)?;
            let (header, len, lines) = if format == VcfFormat::Bcf {
                let (header, len) = read_bcf_header(&mut stream)?;
                (header, len, 0)
            } else {
                read_text_header(&mut stream)?
            };
            // r[impl vcf_reader.file.header_only]
            // No fill: where the header ends in its block, or at a block's end
            // (the next block's start), is known without inflating more.
            let start = stream.virtual_offset().map_or(DataStart::After(len), DataStart::At);
            Ok((format, header, start, lines))
        }
    }
}

/// A BCF header after its magic, and the bytes it took.
pub(crate) fn read_bcf_header(
    stream: &mut BgzfStream<Input>,
) -> Result<(VcfHeader, u64), VcfReadError> {
    stream.advance(BCF_MAGIC.len());
    let mut l_text = [0u8; 4];
    stream.read_exact_into(&mut l_text)?;
    let l_text = u32::from_le_bytes(l_text) as usize;
    // r[impl vcf_reader.limits]
    if l_text > MAX_HEADER_TEXT {
        return Err(
            crate::vcf::VcfHeaderError::TextTooLarge { len: l_text, max: MAX_HEADER_TEXT }.into()
        );
    }
    let mut text = vec![0u8; l_text];
    stream.read_exact_into(&mut text)?;
    let text =
        std::str::from_utf8(&text).map_err(|source| VcfReadError::HeaderNotUtf8 { source })?;
    let header = VcfHeader::parse(text)?;
    let len = (BCF_MAGIC.len() as u64).saturating_add(4).saturating_add(l_text as u64);
    Ok((header, len))
}

/// A VCF text header: the `##` lines through the `#CHROM` line, and the
/// bytes and lines they took. A line that is not a header line ends it, and
/// [`VcfHeader::parse`] says what is missing.
pub(crate) fn read_text_header(
    input: &mut impl BufRead,
) -> Result<(VcfHeader, u64, u64), VcfReadError> {
    let mut text = Vec::new();
    let mut lines = 0u64;
    loop {
        let start = text.len();
        // r[impl vcf_reader.limits]
        let room = (MAX_HEADER_TEXT.saturating_sub(start) as u64).saturating_add(1);
        let peek = input.fill_buf().map_err(read_failed)?;
        if peek.first().is_some_and(|&b| b != b'#') {
            break;
        }
        let n = (&mut *input).take(room).read_until(b'\n', &mut text).map_err(read_failed)?;
        lines = lines.saturating_add(u64::from(n > 0));
        if text.len() > MAX_HEADER_TEXT {
            let len = text.len();
            return Err(
                crate::vcf::VcfHeaderError::TextTooLarge { len, max: MAX_HEADER_TEXT }.into()
            );
        }
        if n == 0 || text.get(start..).is_some_and(|line| line.starts_with(b"#CHROM")) {
            break;
        }
    }
    let len = text.len() as u64;
    let text =
        std::str::from_utf8(&text).map_err(|source| VcfReadError::HeaderNotUtf8 { source })?;
    Ok((VcfHeader::parse(text)?, len, lines))
}

pub(crate) fn read_failed(source: std::io::Error) -> VcfReadError {
    VcfReadError::Bgzf(crate::bam::bgzf::BgzfError::ReadFailed { source })
}

// r[impl vcf_reader.contigs]
/// A file's contigs: the header's, by contig dictionary index, then (for VCF
/// text) contigs records name that neither the header nor the index
/// declares, appended in order of first use and shared by every reader.
#[derive(Debug, Default)]
pub(crate) struct ContigTable {
    by_rid: Vec<Option<ContigId>>,
    discovered: RwLock<Discovered>,
}

#[derive(Debug, Default)]
struct Discovered {
    contigs: Vec<ContigId>,
    warned: bool,
}

impl ContigTable {
    fn new(header: &VcfHeader) -> Self {
        let mut by_rid = vec![None; header.contig_map().len()];
        for (rid, name) in header.contig_map().iter() {
            if let (Some(slot), Ok(id)) = (by_rid.get_mut(rid), header.contig(name)) {
                *slot = Some(id);
            }
        }
        Self { by_rid, discovered: RwLock::default() }
    }

    /// Add `name`, which an index lists, if the header does not declare
    /// it (`r[vcf_reader.contigs]`: such contigs follow the header's); its tid.
    pub(crate) fn add_indexed(&mut self, header: &VcfHeader, name: &str) -> u32 {
        if let Some(id) = self.find(header, name) {
            return id.tid;
        }
        let discovered =
            self.discovered.get_mut().unwrap_or_else(std::sync::PoisonError::into_inner);
        let tid = self.by_rid.len().saturating_add(discovered.contigs.len());
        let tid = u32::try_from(tid).unwrap_or(u32::MAX);
        discovered.contigs.push(ContigId { tid, name: SmolStr::from(name), header: header.id() });
        tid
    }

    /// The contig named `name`, if the header, the index or a record named it.
    pub(crate) fn find(&self, header: &VcfHeader, name: &str) -> Option<ContigId> {
        if let Some(id) =
            header.contig_map().get(name).and_then(|rid| self.by_rid.get(rid)?.clone())
        {
            return Some(id);
        }
        let read = self.discovered.read().unwrap_or_else(std::sync::PoisonError::into_inner);
        read.contigs.iter().find(|c| c.name() == name).cloned()
    }

    /// The contig with tid `tid`.
    pub(crate) fn get(&self, tid: u32) -> Option<ContigId> {
        if let Some(id) = self.by_rid.iter().flatten().find(|c| c.tid == tid) {
            return Some(id.clone());
        }
        let read = self.discovered.read().unwrap_or_else(std::sync::PoisonError::into_inner);
        read.contigs.iter().find(|c| c.tid == tid).cloned()
    }

    /// The contig a BCF record's `CHROM` index names.
    pub(crate) fn by_rid(&self, rid: i32) -> Option<&ContigId> {
        self.by_rid.get(usize::try_from(rid).ok()?)?.as_ref()
    }

    /// The name of contig `tid`, header-declared or discovered.
    pub(crate) fn name(&self, tid: u32) -> Option<SmolStr> {
        let tid = tid as usize;
        if let Some(id) = self.by_rid.get(tid) {
            return id.as_ref().map(|id| id.name.clone());
        }
        let read = self.discovered.read().unwrap_or_else(std::sync::PoisonError::into_inner);
        read.contigs.get(tid.checked_sub(self.by_rid.len())?).map(|id| id.name.clone())
    }

    /// The contig a VCF text record's `CHROM` names, added (with one `warn!`
    /// per file) if neither the header nor the index declares it.
    pub(crate) fn by_name(&self, header: &VcfHeader, name: &str) -> ContigId {
        if let Some(id) =
            header.contig_map().get(name).and_then(|rid| self.by_rid.get(rid)?.clone())
        {
            return id;
        }
        let read = self.discovered.read().unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(id) = read.contigs.iter().find(|c| c.name() == name) {
            return id.clone();
        }
        drop(read);
        let mut write = self.discovered.write().unwrap_or_else(std::sync::PoisonError::into_inner);
        // Another reader may have added it meanwhile.
        if let Some(id) = write.contigs.iter().find(|c| c.name() == name) {
            return id.clone();
        }
        if !std::mem::replace(&mut write.warned, true) {
            warn!(
                contig = name,
                "a record names a contig the header does not declare; added (further ones are not reported)"
            );
        }
        let tid = self.by_rid.len().saturating_add(write.contigs.len());
        let id = ContigId {
            tid: u32::try_from(tid).unwrap_or(u32::MAX),
            name: SmolStr::from(name),
            header: header.id(),
        };
        write.contigs.push(id.clone());
        id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_reader.contigs]
    #[test]
    fn undeclared_contigs_are_appended_once() {
        let mut builder = VcfHeader::builder();
        builder.register_contig("chr1", crate::vcf::ContigDef::new(None)).unwrap();
        let header = builder.build().unwrap();
        let table = ContigTable::new(&header);
        assert_eq!(table.by_rid(0).unwrap().name(), "chr1");
        assert!(table.by_rid(1).is_none() && table.by_rid(-1).is_none());
        assert_eq!(table.by_name(&header, "chr1").tid(), 0);
        let un = table.by_name(&header, "chrUn");
        assert_eq!((un.tid(), un.name().as_str()), (1, "chrUn"));
        assert_eq!(table.by_name(&header, "chrUn"), un);
        assert_eq!(table.by_name(&header, "chrM").tid(), 2);
    }

    // r[verify vcf_reader.formats]
    #[test]
    fn formats_are_detected_from_content() {
        assert!(matches!(container(b"BCF\x02\x02"), Container::Plain));
        assert!(matches!(format_of(b"BCF\x02\x02", false), Err(VcfReadError::UncompressedBcf)));
        assert!(matches!(format_of(b"BCF\x02\x01", true), Err(VcfReadError::Bcf21)));
        assert!(matches!(format_of(b"##fileformat=VCFv4.2", false), Ok(VcfFormat::Vcf)));
        assert!(matches!(format_of(b"##fileformat=VCFv4.2", true), Ok(VcfFormat::VcfGz)));
        assert!(matches!(format_of(b"BCF\x02\x02", true), Ok(VcfFormat::Bcf)));
        assert!(matches!(format_of(b"hello", true), Err(VcfReadError::UnknownFormat)));
        // A gzip member without BGZF's extra field.
        assert!(matches!(container(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3]), Container::Gzip));
    }
}
