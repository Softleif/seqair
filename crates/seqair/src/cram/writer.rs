//! CRAM writer: CRAM 3.1 (or 3.0) against a reference FASTA, next to
//! [`BamWriter`](crate::bam::BamWriter) and taking the same inputs.
//!
//! The writer is the inverse of seqair's CRAM reader: each read is stored as
//! its differences from the reference, every data series and tag in a block
//! of its own compressed by whichever codec is smallest, one slice per
//! container, mates linked within a slice. Input must be coordinate-sorted.
//! A record CRAM cannot store as it is gets a typed [`CramWriteError`] and
//! leaves the writer usable; what reads back differently (`=`/`X` CIGAR
//! operations as `M`, `MD`/`NM` regenerated, …) is listed in the
//! `r[cram_writer.fidelity]` rule of the spec.
//!
//! Records written and read back:
//!
//! ```
//! use seqair::bam::{BamHeader, OwnedBamRecord};
//! use seqair::bam::cigar::{CigarOp, CigarOpType};
//! use seqair::cram::stream::CramStreamReader;
//! use seqair::cram::writer::CramWriterBuilder;
//! use seqair::fasta::IndexedFastaReader;
//! use seqair_types::{IupacBase, Pos0};
//! # let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
//! #     .join("../../tests/data/unplaced_multiref.fa");
//!
//! // The FASTA's chr1 starts ACGTACGTAC.
//! let header = BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:200\n@SQ\tSN:chr2\tLN:200\n")?;
//! let chr1 = header.tid("chr1");
//! let record = OwnedBamRecord::builder(chr1, Pos0::new(2), b"read1".to_vec())
//!     .cigar(vec![CigarOp::new(CigarOpType::Match, 6)])
//!     .seq(IupacBase::from_ascii_slice(b"GTTCGT")) // one substitution
//!     .build()?;
//!
//! let reference = IndexedFastaReader::open(&fasta)?;
//! let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
//!     .write_index(true)
//!     .build()?;
//! writer.write(&record)?;
//! let (cram, index) = writer.finish()?;
//! assert_eq!(index.map(|i| i.entries().len()), Some(1)); // one slice
//!
//! let mut reader = CramStreamReader::new(&cram[..])?;
//! let mut decoder = reader.decoder(&fasta)?;
//! let mut records = Vec::new();
//! while let Some(container) = reader.next_container()? {
//!     decoder.decode(&container, &mut records)?;
//! }
//! assert_eq!(records[0].qname(), b"read1");
//! assert_eq!(records[0].pos(), Pos0::new(2));
//! // The decode regenerates MD and NM from the reference.
//! assert_eq!(records[0].aux().as_bytes(), b"MDZ2A3\0NMC\x01");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! A coordinate-sorted BAM rewritten as CRAM, with its index, on four
//! compression threads (see also `examples/cram_write.rs`):
//!
//! ```no_run
//! use seqair::bam::{BamHeader, RawBamRecord};
//! use seqair::bam::stream::BamStreamReader;
//! use seqair::cram::writer::CramWriterBuilder;
//! use seqair::fasta::IndexedFastaReader;
//! use std::fs::File;
//! use std::path::Path;
//!
//! let mut bam = BamStreamReader::new(File::open("in.bam")?)?;
//! let header = BamHeader::from_template(bam.header());
//! let reference = IndexedFastaReader::open(Path::new("ref.fa"))?;
//! let mut writer = CramWriterBuilder::to_path(Path::new("out.cram"), &header, reference)
//!     .compression_threads(4)
//!     .write_index(true)
//!     .build()?;
//! let mut record = RawBamRecord::default();
//! while bam.read_record(&mut record)? {
//!     writer.write_raw(&record)?;
//! }
//! let (_file, index) = writer.finish()?;
//! if let Some(index) = index {
//!     index.write_crai(File::create("out.cram.crai")?)?;
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod arith;
mod container;
mod fqzcomp;
mod header;
mod mates;
mod metrics;
mod pool;
mod rans;
mod rans_nx16;
mod record;
mod tok3;

use self::container::{Assembled, SliceBuilder, SliceKind};
use self::header::CramHeader;
use self::metrics::Metrics;
use self::pool::{Done, Pool};
use self::record::{RecordView, Scratch};
use super::index::{CraiEntry, CramIndex};
use crate::bam::aux::AuxTag;
use crate::bam::header::BamHeader;
use crate::bam::owned_record::OwnedBamRecord;
use crate::bam::raw_record::RawBamRecord;
use crate::bam::record_idx::RecordIdx;
use crate::bam::record_store::RecordStore;
use crate::bam::writer::{ToPath, ToWriter};
use crate::fasta::{FastaError, IndexedFastaReader};
use crate::io::CompressionLevel;
use crate::reader::Tid;
use seqair_types::{BaseQuality, Pos0, SmolStr};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use thiserror::Error;

pub use self::record::Placement;

// r[impl cram_writer.validate_before_write]
/// Errors from CRAM writing. Errors about a record leave the writer usable
/// (the record is not written); errors writing to the sink poison it.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CramWriteError {
    #[error("I/O error")]
    Io {
        #[from]
        source: io::Error,
    },

    #[error("reference error")]
    Fasta {
        #[from]
        source: FastaError,
    },

    /// The writer was poisoned by an earlier error writing to its sink.
    #[error("writer poisoned after previous error")]
    Poisoned,

    #[error(
        "the header's @SQ line {index} names {sq_name:?}, but target {index} is {target:?}; \
         CRAM keeps the references only in the header text"
    )]
    HeaderTargetsMismatch { index: usize, sq_name: Option<SmolStr>, target: Option<SmolStr> },

    #[error("@SQ {contig} has LN:{header}, but the reference FASTA's is {fasta} bases long")]
    ReferenceLengthMismatch { contig: SmolStr, header: u64, fasta: u64 },

    #[error("@SQ {contig} has M5:{header}, but the reference FASTA's bases have MD5 {fasta}")]
    ReferenceMd5Mismatch { contig: SmolStr, header: SmolStr, fasta: SmolStr },

    // r[impl cram_writer.sort_order]
    #[error(
        "record at {record:?} comes after one at {previous:?}; CRAM input must be coordinate-sorted"
    )]
    Unsorted { previous: Placement, record: Placement },

    #[error(
        "record at {placement:?} has a reference without a position, or is mapped without a reference"
    )]
    InvalidPlacement { placement: Placement },

    #[error("record starts at {pos:?}, past the end of its reference ({length} bases)")]
    PositionPastReferenceEnd { pos: Pos0, length: u64 },

    #[error("record references target {tid}, but the header has only {n_targets} targets")]
    TidOutOfRange { tid: u32, n_targets: usize },

    // r[impl cram_writer.record.name]
    #[error("record has an empty read name")]
    EmptyReadName,

    /// `RN` is NUL-terminated.
    #[error("read name holds a NUL byte")]
    NulInReadName,

    #[error("tag name {tag} is not a letter and a letter or digit")]
    InvalidTagName { tag: AuxTag },

    #[error("CIGAR covers {cigar} query bases, but SEQ has {seq}")]
    CigarSeqLengthMismatch { cigar: u64, seq: usize },

    #[error("record has {qual} qualities for {seq} bases")]
    QualLengthMismatch { qual: usize, seq: usize },

    #[error("CIGAR operation code {code} is not one of MIDNSHP=X")]
    UnknownCigarOp { code: u8 },

    #[error("aux data is malformed")]
    InvalidAux,

    #[error("tag {tag} has type {ty:?}, which CRAM does not store")]
    UnsupportedTagType { tag: AuxTag, ty: char },

    // r[impl cram_writer.tags]
    #[error("{tag} holds a TAB, which SAM does not allow in a string")]
    TabInStringTag { tag: AuxTag },

    #[error("{tag} is reserved: htslib's CRAM decoder takes it for its own flags")]
    ReservedTag { tag: AuxTag },

    #[error("{what} is {value}, more than CRAM can store")]
    ValueTooLarge { what: &'static str, value: u64 },

    #[error("could not start a compression thread")]
    ThreadSpawn { source: io::Error },

    #[error("a compression thread stopped")]
    CompressionWorkerLost,

    /// [`CramWriter::write_store_record`] was given an index the store does
    /// not hold.
    #[error("no record {idx} in the store")]
    NoSuchRecord { idx: RecordIdx },
}

/// A record's qualities, `None` when it has none: empty, or `0xFF` first,
/// as [`RawBamRecord::qual`] reads them.
fn qualities(qual: &[BaseQuality]) -> Option<&[u8]> {
    let bytes = BaseQuality::slice_to_bytes(qual);
    bytes.first().filter(|&&q| q != 0xFF).map(|_| bytes)
}

/// Bytes written so far, which the index's container offsets count.
struct Counting<W> {
    inner: W,
    written: u64,
}

impl<W: Write> Counting<W> {
    fn write_all(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.inner.write_all(bytes)?;
        self.written = self.written.saturating_add(bytes.len() as u64);
        Ok(())
    }
}

// r[impl cram_writer.scope]
/// CRAM 3.0 writer. Build with [`CramWriterBuilder`], write coordinate-sorted
/// records, then [`finish`](Self::finish).
pub struct CramWriter<W: Write> {
    out: Counting<W>,
    header: CramHeader,
    fasta: IndexedFastaReader,
    /// The uppercased bases of the reference the current slice is on.
    reference: Option<(Tid, Vec<u8>)>,
    slice: SliceBuilder,
    last: Option<Placement>,
    record_counter: u64,
    /// Data containers so far.
    slices: u64,
    index: Option<Vec<CraiEntry>>,
    codecs: container::Codecs,
    /// Compression workers, with `compression_threads`.
    pool: Option<Pool>,
    /// Slices back from the pool, to fill again.
    spare: Vec<SliceBuilder>,
    records_per_slice: u32,
    warned_unmapped_alignment: bool,
    scratch: Scratch,
    poisoned: bool,
}

impl<W: Write> CramWriter<W> {
    fn new(
        inner: W,
        file_id: [u8; 20],
        header: &BamHeader,
        mut fasta: IndexedFastaReader,
        options: Options,
    ) -> Result<Self, CramWriteError> {
        let (header, hashed) = CramHeader::new(header, &mut fasta)?;
        let mut out = Counting { inner, written: 0 };
        let metrics = Arc::new(Metrics::new(options.exhaustive_codec_search));
        let mut codecs =
            container::Codecs::new(options.level, options.version, Arc::clone(&metrics));
        // r[impl cram_writer.file_definition]
        let mut def = Vec::with_capacity(26);
        def.extend_from_slice(b"CRAM\x03");
        def.push(match options.version {
            CramVersion::V3_0 => 0,
            CramVersion::V3_1(_) => 1,
        });
        def.extend_from_slice(&file_id);
        out.write_all(&def)?;
        out.write_all(&container::header_container(header.text(), &mut codecs)?)?;
        Ok(Self {
            out,
            header,
            fasta,
            reference: hashed,
            slice: SliceBuilder::default(),
            last: None,
            record_counter: 0,
            slices: 0,
            index: options.write_index.then(Vec::new),
            codecs,
            pool: (options.threads > 0)
                .then(|| Pool::new(options.threads, options.level, options.version, &metrics))
                .transpose()?,
            spare: Vec::new(),
            records_per_slice: options.records_per_slice.unwrap_or(match options.version {
                CramVersion::V3_1(CramProfile::Archive) => container::ARCHIVE_RECORDS,
                _ => container::MAX_RECORDS,
            }),
            warned_unmapped_alignment: false,
            scratch: Scratch::default(),
            poisoned: false,
        })
    }

    // r[impl cram_writer.inputs]
    /// Write one record.
    ///
    /// # Errors
    /// A record CRAM cannot store as it is, or one out of coordinate order,
    /// is refused and leaves the writer as it was. An error writing a full
    /// slice to the sink poisons the writer.
    pub fn write(&mut self, record: &OwnedBamRecord) -> Result<(), CramWriteError> {
        let view = RecordView {
            tid: record.tid,
            pos: record.pos,
            flags: record.flags,
            mapq: record.mapq,
            mate: record.mate_fields,
            qname: &record.qname,
            cigar: &record.cigar,
            seq: &record.seq,
            qual: qualities(&record.qual),
            aux: record.aux.as_bytes(),
        };
        self.write_view(&view)
    }

    /// Write a record as a streaming reader read it. A long CIGAR held in a
    /// `CG` tag is written as the record's CIGAR, without the tag.
    ///
    /// # Errors
    /// As [`Self::write`].
    pub fn write_raw(&mut self, record: &RawBamRecord) -> Result<(), CramWriteError> {
        let mut input = std::mem::take(&mut self.scratch.input);
        input.cigar.clear();
        record.cigar().extend_into(&mut input.cigar);
        input.seq.clear();
        record.iupac_bases_into(&mut input.seq);
        let view = RecordView {
            tid: record.tid(),
            pos: record.pos(),
            flags: record.flags(),
            mapq: record.mapq(),
            mate: record.mate_fields(),
            qname: record.qname(),
            cigar: &input.cigar,
            seq: &input.seq,
            qual: record.qual().map(BaseQuality::slice_to_bytes),
            aux: record.aux().as_bytes(),
        };
        let result = self.write_view(&view);
        self.scratch.input = input;
        result
    }

    /// Write a record of a [`RecordStore`] without copying it out first.
    ///
    /// # Errors
    /// [`CramWriteError::NoSuchRecord`] for an index the store lacks, and the
    /// errors of [`Self::write`].
    pub fn write_store_record(
        &mut self,
        store: &RecordStore,
        idx: RecordIdx,
    ) -> Result<(), CramWriteError> {
        let rec = store.record(idx).ok_or(CramWriteError::NoSuchRecord { idx })?;
        let view = RecordView {
            tid: rec.tid,
            pos: rec.tid.map(|_| rec.pos),
            flags: rec.flags,
            mapq: rec.mapq,
            mate: rec.mate_fields,
            qname: rec.qname(),
            cigar: rec.cigar(),
            seq: rec.seq(),
            qual: qualities(rec.qual()),
            aux: rec.aux(),
        };
        self.write_view(&view)
    }

    fn write_view(&mut self, view: &RecordView<'_>) -> Result<(), CramWriteError> {
        if self.poisoned {
            return Err(CramWriteError::Poisoned);
        }
        let checked = record::check(view, &self.header, self.last)?;
        let kind = match view.tid {
            Some(tid) => SliceKind::Mapped { tid, reference: self.header.reference(tid).is_some() },
            None => SliceKind::Unplaced,
        };
        // r[impl cram_writer.slicing+2]
        if self.slice.ends_before(kind, checked.read_len, self.records_per_slice) {
            self.flush_slice().inspect_err(|_| self.poisoned = true)?;
        }
        if let SliceKind::Mapped { tid, reference: true } = kind {
            self.load_reference(tid)?;
        }
        let reference = match kind {
            SliceKind::Mapped { reference: true, .. } => {
                self.reference.as_ref().map(|(_, bases)| bases.as_slice())
            }
            _ => None,
        };
        record::encode(
            &mut self.slice,
            kind,
            view,
            &checked,
            reference,
            &self.header,
            &mut self.scratch,
        );
        self.last = Some(checked.placement);
        if checked.drops_alignment && !self.warned_unmapped_alignment {
            self.warned_unmapped_alignment = true;
            tracing::warn!(
                qname = std::str::from_utf8(view.qname).unwrap_or("<not UTF-8>"),
                "unmapped record with a CIGAR or MAPQ; CRAM stores neither, so it reads back \
                 with MAPQ 0 and no CIGAR (warned once per writer)"
            );
        }
        Ok(())
    }

    /// Make `self.reference` the bases of `tid`, checking them against the
    /// header's `M5` the first time.
    fn load_reference(&mut self, tid: Tid) -> Result<(), CramWriteError> {
        if self.reference.as_ref().is_some_and(|(loaded, _)| *loaded == tid) {
            return Ok(());
        }
        let mut bases = self.reference.take().map(|(_, b)| b).unwrap_or_default();
        // Not a chr1-sized buffer for a small contig for the rest of the file.
        let len = usize::try_from(self.header.target_len(tid)).unwrap_or(usize::MAX);
        if bases.capacity() / 2 > len {
            bases = Vec::new();
        }
        let result = self.header.fetch(tid, &mut self.fasta, &mut bases);
        result?;
        self.reference = Some((tid, bases));
        Ok(())
    }

    /// Assemble the slice and write it, or hand it to the pool and write
    /// whatever containers have come back.
    fn flush_slice(&mut self) -> Result<(), CramWriteError> {
        let Some(kind) = self.slice.kind() else { return Ok(()) };
        let reference = match kind {
            SliceKind::Mapped { tid, reference: true } => self
                .reference
                .as_ref()
                .filter(|(loaded, _)| *loaded == tid)
                .map(|(_, bases)| bases.as_slice()),
            _ => None,
        };
        let length = match kind {
            SliceKind::Mapped { tid, .. } => self.header.target_len(tid),
            SliceKind::Unplaced => 0,
        };
        let frame = self.slice.frame(self.slices, self.record_counter, length, reference);
        self.slices = self.slices.saturating_add(1);
        self.record_counter = self.record_counter.saturating_add(u64::from(self.slice.n_records()));
        match &mut self.pool {
            None => {
                let assembled = self.slice.assemble(&frame, &mut self.codecs)?;
                self.emit(kind, &assembled)?;
                self.slice.clear();
            }
            Some(pool) => {
                let fresh = self.spare.pop().unwrap_or_default();
                pool.submit(std::mem::replace(&mut self.slice, fresh), frame)?;
                self.drain(false)?;
            }
        }
        Ok(())
    }

    /// Write the containers the pool has finished, in order; with `all`,
    /// every one submitted.
    fn drain(&mut self, all: bool) -> Result<(), CramWriteError> {
        while let Some(done) = self.pool.as_mut().map(|p| p.next(all)).transpose()?.flatten() {
            let Done { mut slice, result, .. } = done;
            self.emit(slice.kind().unwrap_or(SliceKind::Unplaced), &result?)?;
            slice.clear();
            self.spare.push(slice);
        }
        Ok(())
    }

    /// Write a container, and index its slice.
    fn emit(&mut self, kind: SliceKind, assembled: &Assembled) -> Result<(), CramWriteError> {
        let container_offset = self.out.written;
        self.out.write_all(&assembled.bytes)?;
        // r[impl cram_writer.crai]
        if let Some(index) = &mut self.index {
            index.push(CraiEntry {
                ref_id: match kind {
                    SliceKind::Mapped { tid, .. } => Some(tid.as_u32()),
                    SliceKind::Unplaced => None,
                },
                alignment_start: assembled.alignment_start,
                alignment_span: assembled.alignment_span,
                container_offset,
                slice_offset: assembled.slice_offset,
                slice_size: assembled.slice_size,
            });
        }
        Ok(())
    }

    // r[impl cram_writer.finish]
    // r[impl cram_writer.eof]
    /// Write the last slice and the EOF container, and return the inner
    /// writer and, with [`CramWriterBuilder::write_index`], the index of what
    /// was written; [`CramIndex::write_crai`] writes it.
    ///
    /// # Errors
    /// [`CramWriteError::Poisoned`] after an earlier sink error, and errors
    /// writing the last slice.
    pub fn finish(mut self) -> Result<(W, Option<CramIndex>), CramWriteError> {
        if self.poisoned {
            return Err(CramWriteError::Poisoned);
        }
        self.flush_slice()?;
        self.drain(true)?;
        self.out.write_all(&container::EOF_CONTAINER)?;
        self.out.inner.flush()?;
        Ok((self.out.inner, self.index.map(CramIndex::from_entries)))
    }
}

// r[impl cram_writer.builder]
/// Configuration for a [`CramWriter`]: [`Self::to_path`] or
/// [`Self::to_writer`], then [`Self::build`].
pub struct CramWriterBuilder<'a, T> {
    target: T,
    header: &'a BamHeader,
    reference: IndexedFastaReader,
    options: Options,
}

/// What a [`CramWriterBuilder`] sets.
#[derive(Debug, Clone, Copy)]
struct Options {
    write_index: bool,
    level: CompressionLevel,
    threads: usize,
    /// `None`: the profile's.
    records_per_slice: Option<u32>,
    version: CramVersion,
    exhaustive_codec_search: bool,
}

// r[impl cram_writer.version+2]
/// The CRAM version a writer writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CramVersion {
    /// gzip and rANS 4×8 only, for readers without CRAM 3.1 (htsjdk before
    /// 4.x, older tools).
    V3_0,
    /// htslib's default: rANS Nx16 and the name tokeniser as well, and
    /// whatever else the profile adds.
    V3_1(CramProfile),
}

impl Default for CramVersion {
    fn default() -> Self {
        Self::V3_1(CramProfile::default())
    }
}

// r[impl cram_writer.profile]
/// Which CRAM 3.1 codecs the writer tries; it keeps whichever stores a block
/// smallest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum CramProfile {
    /// htslib's default: rANS Nx16 and the name tokeniser.
    #[default]
    Normal,
    /// Also the adaptive arithmetic coder, and fqzcomp for qualities: a few
    /// percent smaller, slower to write and to read.
    Archive,
}

impl<'a> CramWriterBuilder<'a, ToPath<'a>> {
    /// Write the records of `header` to a file at `path`, against
    /// `reference`.
    pub fn to_path(path: &'a Path, header: &'a BamHeader, reference: IndexedFastaReader) -> Self {
        Self::with_target(ToPath { path }, header, reference)
    }

    /// Create the file and write the file definition and header container.
    ///
    /// # Errors
    /// A header whose references disagree with the FASTA, and I/O errors.
    pub fn build(self) -> Result<CramWriter<BufWriter<File>>, CramWriteError> {
        let mut file_id = [0u8; 20];
        let name = self.target.path.file_name().map(|n| n.as_encoded_bytes()).unwrap_or_default();
        for (dst, &src) in file_id.iter_mut().zip(name) {
            *dst = src;
        }
        let file = BufWriter::new(File::create(self.target.path)?);
        CramWriter::new(file, file_id, self.header, self.reference, self.options)
    }
}

impl<'a, W: Write> CramWriterBuilder<'a, ToWriter<W>> {
    /// Write the records of `header` to `writer`, against `reference`.
    pub fn to_writer(writer: W, header: &'a BamHeader, reference: IndexedFastaReader) -> Self {
        Self::with_target(ToWriter { writer }, header, reference)
    }

    /// Write the file definition and header container to the writer. An
    /// index's offsets count from the first byte written.
    ///
    /// # Errors
    /// As [`CramWriterBuilder::build`] for a path.
    pub fn build(self) -> Result<CramWriter<W>, CramWriteError> {
        CramWriter::new(self.target.writer, [0; 20], self.header, self.reference, self.options)
    }
}

impl<'a, T> CramWriterBuilder<'a, T> {
    fn with_target(target: T, header: &'a BamHeader, reference: IndexedFastaReader) -> Self {
        Self {
            target,
            header,
            reference,
            options: Options {
                write_index: false,
                level: CompressionLevel::default(),
                threads: 0,
                records_per_slice: None,
                version: CramVersion::default(),
                exhaustive_codec_search: false,
            },
        }
    }

    /// Build a [`CramIndex`] (`.crai`) of the slices written; `finish`
    /// returns it. Defaults to `false`.
    #[must_use]
    pub fn write_index(mut self, write_index: bool) -> Self {
        self.options.write_index = write_index;
        self
    }

    // r[impl cram_writer.compression+2]
    /// The gzip level blocks are tried at. Defaults to
    /// [`CompressionLevel::Level6`].
    #[must_use]
    pub fn compression_level(mut self, level: CompressionLevel) -> Self {
        self.options.level = level;
        self
    }

    // r[impl cram_writer.threads+2]
    /// Compress containers on `threads` worker threads while the calling
    /// thread encodes the next slice. Defaults to `0`: compress on the
    /// calling thread. The output is the same either way.
    #[must_use]
    pub fn compression_threads(mut self, threads: usize) -> Self {
        self.options.threads = threads;
        self
    }

    /// The CRAM version to write. Defaults to [`CramVersion::V3_1`] with
    /// [`CramProfile::Normal`], as
    /// htslib does.
    #[must_use]
    pub fn version(mut self, version: CramVersion) -> Self {
        self.options.version = version;
        self
    }

    /// End slices after `records` records rather than htslib's 10,000, for
    /// tests that want many slices from few records.
    #[doc(hidden)]
    #[must_use]
    pub fn records_per_slice(mut self, records: u32) -> Self {
        self.options.records_per_slice = Some(records.max(1));
        self
    }

    // r[impl cram_writer.compression.learning]
    /// Try every codec on every block rather than the ones that did best on
    /// the trial slices, for tests and measurements.
    #[doc(hidden)]
    #[must_use]
    pub fn exhaustive_codec_search(mut self, exhaustive: bool) -> Self {
        self.options.exhaustive_codec_search = exhaustive;
        self
    }
}

/// The archive codecs' encoders, for their fuzz targets.
#[cfg(feature = "fuzz")]
#[doc(hidden)]
pub mod fuzz {
    /// The arith encoder with the combination named by `flags`' order,
    /// STRIPE, RLE and PACK bits.
    pub fn arith_encode(data: &[u8], flags: u8) -> Option<Vec<u8>> {
        super::arith::encode(data, super::arith::ArithFlags::from_byte(flags))
    }

    /// The fqzcomp encoder with `(len, reverse, read2)` records and
    /// strategy `strategy % 4` (htslib's `FQZ`, `FQZ_b`, `FQZ_c`, `FQZ_d`).
    pub fn fqzcomp_encode(
        qualities: &[u8],
        records: &[(u32, bool, bool)],
        strategy: u8,
    ) -> Option<Vec<u8>> {
        use super::fqzcomp::{QualRecord, Strategy};
        let records: Vec<QualRecord> = records
            .iter()
            .map(|&(len, reverse, read2)| QualRecord { len, reverse, read2 })
            .collect();
        let strategy = Strategy::ALL.get(usize::from(strategy % 4)).copied()?;
        super::fqzcomp::encode(qualities, &records, strategy)
    }
}

#[cfg(test)]
mod tests;
