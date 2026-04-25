use super::{
    ReaderError,
    indexed::IndexedReader,
    resolve::{ResolveTid, Tid},
};
use crate::{
    bam::{
        BamHeader,
        pileup::{PileupEngine, RefSeq},
        record_store::{CustomizeRecordStore, RecordStore},
    },
    fasta::{FastaError, IndexedFastaReader},
};
use seqair_types::{Base, Pos0, RegionString, SmolStr};
use std::path::Path;
use std::rc::Rc;
use tracing::instrument;

/// Alignment + reference reader bundle.
///
/// Bundles an [`IndexedReader`] (BAM/SAM/CRAM) with an [`IndexedFastaReader`]
/// so that CRAM has access to the reference it needs and all formats have
/// uniform open/fork/fetch semantics.
///
/// The optional type parameter `E` attaches a [`CustomizeRecordStore`]
/// value — user code defines per-record filtering (`keep_record`) and
/// per-record data computation (`compute`), and accesses the computed data
/// during pileup iteration via
/// [`AlignmentView::extra`](crate::bam::pileup::AlignmentView::extra). When
/// `E = ()` (the default), no records are filtered and no extras are
/// computed (zero-cost).
///
/// # Quick start: pileup at a genomic region
///
/// ```no_run
/// use seqair::reader::Readers;
/// use seqair::bam::pileup::PileupOp;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Open BAM + FASTA — auto-detects BAM/SAM/CRAM
/// let mut readers = Readers::open(Path::new("sample.bam"), Path::new("reference.fa"))?;
///
/// // Resolve "chr19:6_100_000-6_200_000" against the header.
/// let region = "chr19:6100000-6200000".parse()?;
/// let (tid, start, end) = readers.resolve_region(&region)?;
///
/// // `pileup()` fetches records AND the reference sequence for the region.
/// let mut pileup = readers.pileup(tid, start, end)?;
///
/// while let Some(column) = pileup.pileups() {
///     let _pos      = column.pos();
///     let _ref_base = column.reference_base();
///     let _depth    = column.depth();
///
///     for aln in column.alignments() {
///         match &aln.op {
///             PileupOp::Match { base, qual, .. } => { let _ = (base, qual); }
///             PileupOp::Insertion { base, qual, insert_len, .. } => {
///                 let _ = (base, qual, insert_len);
///             }
///             PileupOp::Deletion { del_len } => { let _ = del_len; }
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
///             let mut pileup = forked.pileup(
///                 0u32,
///                 Pos0::new(0).unwrap(),
///                 Pos0::new(100_000).unwrap(),
///             ).unwrap();
///             while let Some(_col) = pileup.pileups() { /* process */ }
///         });
///     }
/// });
/// # Ok(())
/// # }
/// ```
// r[impl unified.readers_struct]
pub struct Readers<E: CustomizeRecordStore = ()> {
    pub(crate) alignment: IndexedReader,
    pub(crate) fasta: IndexedFastaReader,
    pub(crate) store: RecordStore<()>,
    pub(crate) fasta_buf: Vec<u8>,
    pub(crate) customize: E,
}

impl<E: CustomizeRecordStore> std::fmt::Debug for Readers<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Readers").field("alignment", &self.alignment).finish()
    }
}

impl Readers<()> {
    /// Open an alignment file (BAM/SAM/CRAM) and a FASTA reference.
    ///
    /// Auto-detects the alignment format. For CRAM, the FASTA path is passed
    /// to the CRAM reader for sequence reconstruction.
    // r[impl unified.readers_open]
    #[instrument(level = "debug", fields(alignment = %alignment_path.display(), fasta = %fasta_path.display()))]
    pub fn open(alignment_path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        Self::open_customized(alignment_path, fasta_path, ())
    }
}

impl<E: CustomizeRecordStore> Readers<E> {
    // r[impl unified.readers_open_customized]
    /// Open like [`open`](Readers::open) but attach a [`CustomizeRecordStore`]
    /// value so per-record filtering and extras computation runs every time
    /// records are loaded.
    pub fn open_customized(
        alignment_path: &Path,
        fasta_path: &Path,
        customize: E,
    ) -> Result<Self, ReaderError> {
        let fasta = IndexedFastaReader::open(fasta_path)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        let alignment = IndexedReader::open_with_fasta(alignment_path, fasta_path)?;
        Ok(Readers {
            alignment,
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
            customize,
        })
    }

    /// Fork both the alignment reader and the FASTA reader.
    ///
    /// The customize value is cloned so each fork has its own copy.
    // r[impl unified.readers_fork]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        let alignment = self.alignment.fork()?;
        let fasta = self.fasta.fork().map_err(|source| ReaderError::FastaFork { source })?;
        Ok(Readers {
            alignment,
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
            customize: self.customize.clone(),
        })
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
        store: &mut RecordStore<()>,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    /// Access the customize value, e.g. to inspect any internal counters it carries.
    pub fn customize(&self) -> &E {
        &self.customize
    }

    /// Mutable access to the customize value, e.g. to reset state between regions.
    pub fn customize_mut(&mut self) -> &mut E {
        &mut self.customize
    }

    // r[impl unified.readers_resolve_region]
    /// Resolve a [`RegionString`] into `(tid, start, end)` for [`pileup`](Self::pileup).
    ///
    /// Missing start defaults to position 1 (1-based) → `Pos0(0)`. Missing end
    /// defaults to the contig's last base (inclusive) — so a bare `chr1`
    /// resolves to `[0, contig_len - 1]`. Returns [`ReaderError::EmptyContig`]
    /// if the contig has zero length.
    pub fn resolve_region(&self, region: &RegionString) -> Result<(Tid, Pos0, Pos0), ReaderError> {
        let tid = region.chromosome.as_str().resolve_tid(self.header())?;
        let start_pos0 = match region.start {
            Some(p1) => p1.to_zero_based(),
            None => Pos0::ZERO,
        };
        let end_pos0 = match region.end {
            Some(p1) => p1.to_zero_based(),
            None => {
                // Default end: last base of the contig (inclusive).
                let len = self.header().target_len(tid.as_u32()).unwrap_or(0);
                if len == 0 {
                    return Err(ReaderError::EmptyContig { name: region.chromosome.clone() });
                }
                let last = len.checked_sub(1).expect("len > 0 checked above");
                let last_u32 = u32::try_from(last)
                    .map_err(|_| ReaderError::RegionEndTooLarge { end: last })?;
                Pos0::new(last_u32).ok_or(ReaderError::RegionEndTooLarge { end: last })?
            }
        };
        Ok((tid, start_pos0, end_pos0))
    }

    // r[impl unified.readers_pileup]
    /// Fetch records for a region and return a [`PileupEngine`] ready for iteration.
    ///
    /// Accepts any [`ResolveTid`] value (a `u32`, a contig name `&str`, or a
    /// pre-resolved [`Tid`]). Always fetches the reference sequence for
    /// `[start, end]` and attaches it to the engine, so [`PileupColumn::reference_base`]
    /// returns real bases instead of `Unknown`.
    ///
    /// When `E != ()`, the extras provider is run once per record after the
    /// records are loaded, before the engine is constructed.
    ///
    /// Uses an internal [`RecordStore`] whose capacity is retained across calls.
    /// After iterating the engine, call [`recover_store`](Self::recover_store) to
    /// return the store for reuse. If not called, the next `pileup()` call
    /// allocates a fresh store (small perf hit, not a correctness issue).
    ///
    /// [`PileupColumn::reference_base`]: crate::bam::pileup::PileupColumn::reference_base
    pub fn pileup(
        &mut self,
        tid: impl ResolveTid,
        start: Pos0,
        end: Pos0,
    ) -> Result<PileupEngine<E::Extra>, ReaderError> {
        let tid = tid.resolve_tid(self.header())?;

        // Split borrows: fetch writes into `store` while customize is consulted
        // for keep_record. Three disjoint &mut borrows of self.
        let alignment = &mut self.alignment;
        let store = &mut self.store;
        let customize = &mut self.customize;
        alignment.fetch_into_customized(tid.as_u32(), start, end, store, customize)?;

        let contig = self
            .header()
            .target_name(tid.as_u32())
            .ok_or_else(|| super::resolve::TidError::TidOutOfRange {
                tid: tid.as_u32(),
                n_targets: u32::try_from(self.header().target_count()).unwrap_or(u32::MAX),
            })?
            .to_owned();
        // Fetch `[start, end]` (inclusive). FASTA APIs expect half-open [start, stop),
        // so stop = end + 1; clamp at i32::MAX.
        let stop =
            end.checked_add_offset(seqair_types::Offset::new(1)).unwrap_or_else(Pos0::max_value);
        let contig_name: SmolStr = contig.into();
        self.fasta.fetch_seq_into(&contig_name, start, stop, &mut self.fasta_buf).map_err(
            |source| ReaderError::FastaFetch {
                contig: contig_name.clone(),
                start: u32::try_from(start.as_i64().max(0)).unwrap_or(0),
                end: u32::try_from(end.as_i64().max(0)).unwrap_or(0),
                source,
            },
        )?;
        let fasta_buf = std::mem::take(&mut self.fasta_buf);
        let bases = Base::from_ascii_vec(fasta_buf);
        let ref_seq = RefSeq::new(Rc::from(bases), start);

        let base_store = std::mem::take(&mut self.store);
        let typed_store = base_store.apply_customize(&mut self.customize);

        let mut engine = PileupEngine::new(typed_store, start, end);
        engine.set_reference_seq(ref_seq);
        Ok(engine)
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
