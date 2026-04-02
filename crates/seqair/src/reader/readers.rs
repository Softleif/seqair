use super::{ReaderError, indexed::IndexedReader};
use crate::{
    bam::{BamHeader, pileup::PileupEngine, record_store::RecordStore},
    fasta::{FastaError, IndexedFastaReader},
};
use seqair_types::Base;
use seqair_types::Pos;
use seqair_types::Zero;
use std::path::Path;
use std::rc::Rc;
use tracing::instrument;

/// Alignment + reference reader bundle.
///
/// Bundles an [`IndexedReader`] (BAM/SAM/CRAM) with an [`IndexedFastaReader`]
/// so that CRAM has access to the reference it needs and all formats have
/// uniform open/fork/fetch semantics.
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
    #[instrument(level = "debug", fields(alignment = %alignment_path.display(), fasta = %fasta_path.display()), err)]
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
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    /// Fetch records for a region and return a [`PileupEngine`] ready for iteration.
    ///
    /// Uses an internal [`RecordStore`] whose capacity is retained across calls.
    /// After iterating the engine, call [`recover_store`](Self::recover_store) to
    /// return the store for reuse. If not called, the next `pileup()` call
    /// allocates a fresh store (small perf hit, not a correctness issue).
    pub fn pileup(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
    ) -> Result<PileupEngine, ReaderError> {
        self.alignment.fetch_into(tid, start, end, &mut self.store)?;

        let store = std::mem::take(&mut self.store);
        Ok(PileupEngine::new(store, start, end))
    }

    /// Recover the [`RecordStore`] from a consumed [`PileupEngine`] for reuse.
    ///
    /// Call this after iteration is complete. The store retains its allocated
    /// capacity, avoiding ~39 MB of re-allocation on the next `pileup()` call.
    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
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
        start: Pos<Zero>,
        stop: Pos<Zero>,
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
