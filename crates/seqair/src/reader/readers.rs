use super::{ReaderError, indexed::IndexedReader};
use crate::{
    bam::{BamHeader, pileup::PileupEngine, record_store::RecordStore},
    fasta::{FastaError, IndexedFastaReader},
};
use seqair_types::Base;
use seqair_types::Pos0;
use std::path::Path;
use std::rc::Rc;
use tracing::instrument;

/// Alignment + reference reader bundle.
///
/// Bundles an [`IndexedReader`] (BAM/SAM/CRAM) with an [`IndexedFastaReader`]
/// so that CRAM has access to the reference it needs and all formats have
/// uniform open/fork/fetch semantics.
///
/// # Quick start: pileup at a genomic region
///
/// ```no_run
/// use seqair::reader::Readers;
/// use seqair::bam::pileup::PileupOp;
/// use seqair_types::Pos0;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Open BAM + FASTA — auto-detects BAM/SAM/CRAM
/// let mut readers = Readers::open(Path::new("sample.bam"), Path::new("reference.fa"))?;
///
/// let tid = readers.header().tid("chr19").expect("contig not found");
/// let start = Pos0::new(6_100_000).unwrap();
/// let end   = Pos0::new(6_200_000).unwrap();
/// let mut pileup = readers.pileup(tid, start, end)?;
///
/// for column in pileup.by_ref() {
///     let _pos      = column.pos();
///     let _ref_base = column.reference_base();
///     // depth() counts all alignments; match_depth() excludes deletions/ref-skips
///     let _depth = column.depth();
///
///     for aln in column.alignments() {
///         match &aln.op {
///             PileupOp::Match { base, qual, .. } => {
///                 let _ = (base, qual);
///             }
///             PileupOp::Insertion { base, qual, insert_len, .. } => {
///                 let _ = (base, qual, insert_len);
///             }
///             PileupOp::Deletion { del_len } => {
///                 let _ = del_len;
///             }
///             PileupOp::ComplexIndel { del_len, insert_len, .. } => {
///                 let _ = (del_len, insert_len);
///             }
///             PileupOp::RefSkip => {}
///         }
///     }
/// }
///
/// // Return the RecordStore to Readers for reuse on the next region
/// readers.recover_store(&mut pileup);
/// # Ok(())
/// # }
/// ```
///
/// # Multi-threaded pileup
///
/// [`fork`](Readers::fork) gives each thread a fresh file handle while sharing
/// the parsed index and header via `Arc` — no locking, no re-parsing.
///
/// ```no_run
/// # use seqair::reader::Readers;
/// # use seqair_types::Pos0;
/// # use std::path::Path;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let readers = Readers::open(Path::new("sample.bam"), Path::new("ref.fa"))?;
///
/// std::thread::scope(|s| {
///     for _ in 0..4 {
///         let mut forked = readers.fork().unwrap();
///         s.spawn(move || {
///             let pileup = forked.pileup(
///                 0,
///                 Pos0::new(0).unwrap(),
///                 Pos0::new(100_000).unwrap(),
///             ).unwrap();
///             for _col in pileup { /* process */ }
///         });
///     }
/// });
/// # Ok(())
/// # }
/// ```
// r[impl unified.readers_struct]
pub struct Readers {
    pub(crate) alignment: IndexedReader,
    pub(crate) fasta: IndexedFastaReader,
    pub(crate) store: RecordStore,
    pub(crate) fasta_buf: Vec<u8>,
}

impl std::fmt::Debug for Readers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Readers").field("alignment", &self.alignment).finish()
    }
}

impl Readers {
    /// Open an alignment file (BAM/SAM/CRAM) and a FASTA reference.
    ///
    /// Auto-detects the alignment format. For CRAM, the FASTA path is passed
    /// to the CRAM reader for sequence reconstruction.
    // r[impl unified.readers_open]
    #[instrument(level = "debug", fields(alignment = %alignment_path.display(), fasta = %fasta_path.display()))]
    pub fn open(alignment_path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        let fasta = IndexedFastaReader::open(fasta_path)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        let alignment = IndexedReader::open_with_fasta(alignment_path, fasta_path)?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    /// Fork both the alignment reader and the FASTA reader.
    // r[impl unified.readers_fork]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        let alignment = self.alignment.fork()?;
        let fasta = self.fasta.fork().map_err(|source| ReaderError::FastaFork { source })?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    // r[impl unified.readers_accessors]
    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    /// Configure the internal [`RecordStore`] used by `pileup()`.
    pub fn configure_store(&mut self, f: impl FnOnce(&mut RecordStore)) {
        f(&mut self.store);
    }

    /// Fetch records for a region and return a [`PileupEngine`] ready for iteration.
    ///
    /// Uses an internal [`RecordStore`] whose capacity is retained across calls.
    /// After iterating the engine, call [`recover_store`](Self::recover_store) to
    /// return the store for reuse. If not called, the next `pileup()` call
    /// allocates a fresh store (small perf hit, not a correctness issue).
    pub fn pileup(
        &mut self,
        tid: impl ResolveTid, // TODO: trait with `.resolve_tid(&self, &bam_header) -> Result<Tid>` impl'd on &str and u32, and where `struct Tid(u32)` represents a resolved and valid tid that fetch_into accepts
        start: Pos0,
        end: Pos0,
    ) -> Result<PileupEngine, ReaderError> {
        self.alignment.fetch_into(tid, start, end, &mut self.store)?;

        let store = std::mem::take(&mut self.store);
        Ok(PileupEngine::new(store, start, end))
    }

    // r[impl pileup.extras.recover_store]
    /// Recover the [`RecordStore`] from a consumed [`PileupEngine`] for reuse.
    ///
    /// Call this after iteration is complete. The store retains its allocated
    /// capacity, avoiding ~39 MB of re-allocation on the next `pileup()` call.
    /// Accepts `PileupEngine<U>` for any `U` — extras are stripped during recovery.
    pub fn recover_store<U>(&mut self, engine: &mut PileupEngine<U>) {
        if let Some(store) = engine.take_store() {
            self.store = store.strip_extras();
        }
    }

    pub fn fasta(&self) -> &IndexedFastaReader {
        &self.fasta
    }

    pub fn fasta_mut(&mut self) -> &mut IndexedFastaReader {
        &mut self.fasta
    }

    // r[impl fasta.fetch.buffer_reuse]
    /// Fetch a reference sequence region and return it as `Rc<[Base]>`.
    ///
    /// Uses an internal buffer whose capacity is retained across calls,
    /// avoiding per-segment allocation. The returned `Rc<[Base]>` owns its
    /// own copy — the internal buffer is reused on the next call.
    pub fn fetch_base_seq(
        &mut self,
        name: &str,
        start: Pos0,
        stop: Pos0,
    ) -> Result<Rc<[Base]>, FastaError> {
        self.fasta.fetch_seq_into(name, start, stop, &mut self.fasta_buf)?;
        // Take the buffer so from_ascii_vec can reinterpret it in-place (safe),
        // then create the Rc (which copies). The buffer capacity is lost but
        // fasta_buf re-grows on the next call.
        let buf = std::mem::take(&mut self.fasta_buf);
        let bases = Base::from_ascii_vec(buf);
        Ok(Rc::from(bases))
    }

    pub fn alignment(&self) -> &IndexedReader {
        &self.alignment
    }

    pub fn alignment_mut(&mut self) -> &mut IndexedReader {
        &mut self.alignment
    }
}
