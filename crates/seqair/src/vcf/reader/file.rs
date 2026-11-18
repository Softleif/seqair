//! [`VcfFile`]: a VCF/BCF file opened once — format, header, index and the
//! tables every reader of it shares (`r[vcf_reader.file]`).

use super::error::{VcfFormat, VcfReadError};
use super::fields::Fields;
use super::stream::ReaderBuilder;
use crate::bam::bgzf::{BgzfStream, VirtualOffset, find_bsize};
use crate::bam::csi_index::CsiIndex;
use crate::vcf::field_set::FieldSet;
use crate::vcf::header::{MAX_HEADER_TEXT, VcfHeader, VcfVersion};
use crate::vcf::record_encoder::{ContigId, FilterId};
use seqair_types::SmolStr;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
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
const FILE_BUFFER: usize = 128 * 1024;

// r[impl vcf_reader.file]
// r[impl vcf_reader.file.shared]
/// A VCF or BCF file: its header, index and contig table, parsed once and
/// shared by every reader opened from it.
///
/// `Send + Sync + Clone`; a clone shares the parsed state and re-reads
/// nothing, so a rayon task can open its own reader from a clone.
#[derive(Debug, Clone)]
pub struct VcfFile(pub(crate) Arc<FileShared>);

/// Where a file's bytes come from.
#[derive(Debug)]
enum Source {
    Path(PathBuf),
    Bytes(Arc<[u8]>),
}

/// A reader's input: a file handle or a view of in-memory bytes.
#[derive(Debug)]
pub(crate) enum Input {
    File(BufReader<File>),
    Bytes(Cursor<Arc<[u8]>>),
}

impl Read for Input {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::File(f) => f.read(buf),
            Self::Bytes(b) => b.read(buf),
        }
    }
}

impl Seek for Input {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::File(f) => f.seek(pos),
            Self::Bytes(b) => b.seek(pos),
        }
    }
}

/// Where a file's records start.
#[derive(Debug, Clone, Copy)]
pub(crate) enum DataStart {
    /// At this virtual offset.
    At(VirtualOffset),
    /// This many decompressed bytes from the start (an offset the stream
    /// could not name).
    After(u64),
}

/// Everything a [`VcfFile`]'s readers share.
#[derive(Debug)]
pub(crate) struct FileShared {
    source: Source,
    pub(crate) header: VcfHeader,
    pub(crate) data_start: DataStart,
    index: Option<CsiIndex>,
    pub(crate) contigs: ContigTable,
    /// FILTER ids by dictionary index; `None` where the index names no FILTER.
    pub(crate) filters: Vec<Option<FilterId>>,
    /// Dictionary indices the header declares as INFO / FORMAT fields.
    pub(crate) info_ids: FieldSet,
    pub(crate) format_ids: FieldSet,
    /// Whether a GT's first phase bit is read (VCF 4.4 and later).
    pub(crate) first_phase_is_read: bool,
    warned_duplicate: AtomicBool,
}

impl Source {
    fn open(&self) -> Result<Input, VcfReadError> {
        match self {
            Source::Path(path) => File::open(path)
                .map(|f| Input::File(BufReader::with_capacity(FILE_BUFFER, f)))
                .map_err(|source| VcfReadError::Open { path: path.clone(), source }),
            Source::Bytes(bytes) => Ok(Input::Bytes(Cursor::new(Arc::clone(bytes)))),
        }
    }
}

impl FileShared {
    /// A fresh handle on the file's bytes, at their start.
    pub(crate) fn open_input(&self) -> Result<Input, VcfReadError> {
        self.source.open()
    }

    // r[impl vcf_values.duplicates]
    /// Log, once per file, that a record holds `name` twice.
    pub(crate) fn warn_duplicate(&self, name: &str) {
        if !self.warned_duplicate.swap(true, Ordering::Relaxed) {
            warn!(
                field = name,
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
        let index = match crate::utils::find_index(path, &["csi"]) {
            Some(index) => Some(CsiIndex::from_path(&index)?),
            None => None,
        };
        Self::new(Source::Path(path.to_path_buf()), index)
    }

    /// A file held in memory, with its index's bytes if it has one.
    ///
    /// # Errors
    /// As [`Self::open`].
    pub fn from_bytes(
        data: impl Into<Arc<[u8]>>,
        index: Option<&[u8]>,
    ) -> Result<Self, VcfReadError> {
        let index = index.map(CsiIndex::from_bytes).transpose()?;
        Self::new(Source::Bytes(data.into()), index)
    }

    fn new(source: Source, index: Option<CsiIndex>) -> Result<Self, VcfReadError> {
        let (header, data_start) = read_header(source.open()?)?;
        let shared = FileShared {
            info_ids: declared(&header, header.infos().keys()),
            format_ids: declared(&header, header.formats().keys()),
            contigs: ContigTable::new(&header),
            filters: filter_table(&header),
            first_phase_is_read: header.version() >= VcfVersion::V4_4,
            warned_duplicate: AtomicBool::new(false),
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

    /// Whether an index was found (or given) and parsed.
    pub fn has_index(&self) -> bool {
        self.0.index.is_some()
    }

    // r[impl vcf_reader.reader]
    /// A reader of records into `F`, to configure and [`open`](ReaderBuilder::open).
    pub fn reader<F: Fields>(&self) -> ReaderBuilder<'_, F> {
        ReaderBuilder::new(self)
    }
}

/// The dictionary indices of `names`.
fn declared<'h>(header: &VcfHeader, names: impl Iterator<Item = &'h SmolStr>) -> FieldSet {
    let mut set = FieldSet::default();
    for id in names.filter_map(|n| header.string_map().get(n)) {
        if let Ok(id) = u32::try_from(id) {
            set.insert(id);
        }
    }
    set
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
enum Container {
    Bgzf,
    Gzip,
    Plain,
}

// r[impl vcf_reader.formats]
/// The container `prefix` (the file's first bytes) starts with.
fn container(prefix: &[u8]) -> Container {
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
fn format_of(start: &[u8], compressed: bool) -> Result<VcfFormat, VcfReadError> {
    if start.starts_with(BCF_MAGIC) {
        return if compressed { Ok(VcfFormat::Bcf) } else { Err(VcfReadError::UncompressedBcf) };
    }
    if start.starts_with(BCF_21_MAGIC) {
        return Err(VcfReadError::Bcf21);
    }
    if start.starts_with(VCF_MAGIC) {
        let format = if compressed { VcfFormat::VcfGz } else { VcfFormat::Vcf };
        return Err(VcfReadError::Unsupported { format });
    }
    Err(VcfReadError::UnknownFormat)
}

/// Detect the format, read the header, and find where records start.
fn read_header(mut input: Input) -> Result<(VcfHeader, DataStart), VcfReadError> {
    let mut prefix = Vec::with_capacity(64);
    (&mut input).take(64).read_to_end(&mut prefix).map_err(read_failed)?;
    input.seek(SeekFrom::Start(0)).map_err(read_failed)?;
    match container(&prefix) {
        Container::Gzip => return Err(VcfReadError::NotBgzf),
        Container::Plain => {
            format_of(&prefix, false)?;
            return Err(VcfReadError::UnknownFormat);
        }
        Container::Bgzf => {}
    }
    let mut stream = BgzfStream::new(input);
    format_of(stream.fill()?, true)?;
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
    stream.fill()?;
    let header_len = (BCF_MAGIC.len() as u64).saturating_add(4).saturating_add(l_text as u64);
    let start = stream.virtual_offset().map_or(DataStart::After(header_len), DataStart::At);
    Ok((header, start))
}

fn read_failed(source: std::io::Error) -> VcfReadError {
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

    /// The contig a BCF record's `CHROM` index names.
    pub(crate) fn by_rid(&self, rid: i32) -> Option<&ContigId> {
        self.by_rid.get(usize::try_from(rid).ok()?)?.as_ref()
    }

    /// The contig a VCF text record's `CHROM` names, added (with one `warn!`
    /// per file) if neither the header nor the index declares it.
    #[cfg_attr(not(test), expect(dead_code, reason = "VCF text records (next phase)"))]
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
        assert_eq!((un.tid(), un.name()), (1, "chrUn"));
        assert_eq!(table.by_name(&header, "chrUn"), un);
        assert_eq!(table.by_name(&header, "chrM").tid(), 2);
    }

    // r[verify vcf_reader.formats]
    #[test]
    fn formats_are_detected_from_content() {
        assert!(matches!(container(b"BCF\x02\x02"), Container::Plain));
        assert!(matches!(format_of(b"BCF\x02\x02", false), Err(VcfReadError::UncompressedBcf)));
        assert!(matches!(format_of(b"BCF\x02\x01", true), Err(VcfReadError::Bcf21)));
        assert!(matches!(
            format_of(b"##fileformat=VCFv4.2", false),
            Err(VcfReadError::Unsupported { format: VcfFormat::Vcf })
        ));
        assert!(matches!(format_of(b"BCF\x02\x02", true), Ok(VcfFormat::Bcf)));
        assert!(matches!(format_of(b"hello", true), Err(VcfReadError::UnknownFormat)));
        // A gzip member without BGZF's extra field.
        assert!(matches!(container(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3]), Container::Gzip));
    }
}
