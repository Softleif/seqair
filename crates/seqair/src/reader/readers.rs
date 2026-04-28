use super::{
    ReaderError,
    indexed::IndexedReader,
    segment::{IntoSegmentTarget, Segment, SegmentOptions, Segments},
};
use crate::{
    bam::{
        BamHeader,
        pileup::{PileupEngine, PileupGuard, RefSeq},
        record_store::{CustomizeRecordStore, RecordStore},
    },
    fasta::{FastaError, IndexedFastaReader},
};
use seqair_types::{Base, Pos0};
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
/// `pileup()` only takes a [`Segment`], which you obtain from
/// [`segments`](Self::segments). Pick a `max_len` that bounds the per-tile
/// memory footprint, iterate the plan, and pile up each segment.
///
/// ```no_run
/// use seqair::reader::{Readers, SegmentOptions};
/// use seqair::bam::pileup::PileupOp;
/// use seqair_types::RegionString;
/// use std::num::NonZeroU32;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Open BAM + FASTA — auto-detects BAM/SAM/CRAM
/// let mut readers = Readers::open(Path::new("sample.bam"), Path::new("reference.fa"))?;
///
/// // 100 kb tiles, no overlap.
/// let opts = SegmentOptions::new(NonZeroU32::new(100_000).unwrap());
///
/// // Plan a pileup over a parsed region. The iterator borrows `&self` only.
/// let region: RegionString = "chr19:6100000-6200000".parse()?;
/// let plan: Vec<_> = readers.segments(&region, opts)?.collect();
///
/// for segment in &plan {
///     // `pileup()` fetches records AND the reference sequence for the segment.
///     // The returned guard derefs to PileupEngine and on drop returns the
///     // RecordStore to Readers (no manual recover step needed).
///     let mut pileup = readers.pileup(segment)?;
///
///     while let Some(column) = pileup.pileups() {
///         let _pos      = column.pos();
///         let _ref_base = column.reference_base();
///         let _depth    = column.depth();
///
///         for aln in column.alignments() {
///             match &aln.op {
///                 PileupOp::Match { base, qual, .. } => { let _ = (base, qual); }
///                 PileupOp::Insertion { base, qual, insert_len, .. } => {
///                     let _ = (base, qual, insert_len);
///                 }
///                 PileupOp::Deletion { del_len } => { let _ = del_len; }
///                 PileupOp::ComplexIndel { del_len, insert_len, .. } => {
///                     let _ = (del_len, insert_len);
///                 }
///                 PileupOp::RefSkip => {}
///             }
///         }
///     }
///     // `pileup` drops here, returning the RecordStore to `readers` for
///     // the next segment.
/// }
/// # Ok(())
/// # }
/// ```
///
/// # Multi-threaded pileup
///
/// [`fork`](Readers::fork) gives each thread a fresh file handle while sharing
/// the parsed index and header via `Arc` — no locking, no re-parsing. Build
/// the segment plan once on the orchestrator, then ship segments to workers.
///
/// ```no_run
/// # use seqair::reader::{Readers, SegmentOptions};
/// # use std::num::NonZeroU32;
/// # use std::path::Path;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let readers = Readers::open(Path::new("sample.bam"), Path::new("ref.fa"))?;
/// let opts = SegmentOptions::new(NonZeroU32::new(100_000).unwrap());
/// let segments: Vec<_> = readers.segments("chr1", opts)?.collect();
///
/// std::thread::scope(|s| {
///     for chunk in segments.chunks(segments.len().max(1).div_ceil(4)) {
///         let mut forked = readers.fork().unwrap();
///         let chunk = chunk.to_vec();
///         s.spawn(move || {
///             for seg in &chunk {
///                 let mut pileup = forked.pileup(seg).unwrap();
///                 while let Some(_col) = pileup.pileups() { /* process */ }
///                 // `pileup` drops here; store goes back into `forked`.
///             }
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
    pub(crate) store: RecordStore<E::Extra>,
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

    // r[impl unified.readers_segments]
    /// Plan a pileup pass: produce an iterator of [`Segment`]s covering
    /// `target` according to `opts`.
    ///
    /// `target` is anything that implements [`IntoSegmentTarget`] — the
    /// trait is **sealed**, so the closed list of accepted targets is
    /// exactly: a contig name (`&str` / `String` / `SmolStr`), a
    /// pre-resolved [`Tid`] or `u32`, a parsed [`RegionString`], an
    /// explicit `(resolver, start, end)` tuple, or `()` for a whole-genome
    /// scan. Each yielded `Segment`'s **core** is at most `opts.max_len()`
    /// bases long; the full `[start, end]` includes `opts.overlap()` bases
    /// of context on each side (clipped to the requested range at edges),
    /// so internal tiles can be up to `max_len + 2 * overlap` bases.
    ///
    /// The returned iterator borrows `&self` only — feed each segment to
    /// [`pileup`](Self::pileup) which needs `&mut self`.
    ///
    /// ```no_run
    /// # use seqair::reader::{Readers, SegmentOptions};
    /// # use std::num::NonZeroU32;
    /// # use std::path::Path;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut readers = Readers::open(Path::new("sample.bam"), Path::new("ref.fa"))?;
    /// let opts = SegmentOptions::new(NonZeroU32::new(100_000).unwrap()).with_overlap(200)?;
    /// let plan: Vec<_> = readers.segments("chr19", opts)?.collect();
    /// for seg in plan {
    ///     let mut p = readers.pileup(&seg)?;
    ///     while let Some(_col) = p.pileups() { /* ... */ }
    ///     // `p` drops here; the RecordStore goes back into `readers`.
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn segments(
        &self,
        target: impl IntoSegmentTarget,
        opts: SegmentOptions,
    ) -> Result<Segments, ReaderError> {
        let ranges = super::segment::resolve_target(target, self.header())?;
        Ok(Segments::new(ranges, opts))
    }

    // r[impl unified.readers_pileup]
    /// Fetch records and reference sequence for `segment`, returning a
    /// [`PileupGuard`] that derefs to a [`PileupEngine`] ready for iteration.
    ///
    /// `segment` must come from [`segments`](Self::segments) — there is no
    /// other way to construct one. The segment's `tid`, `start`, `end`, and
    /// `contig` name are used directly.
    ///
    /// **Header consistency.** [`Segment`] is `Send + Clone`, so it can be
    /// shipped between threads or across `fork()`s. To catch the foot-gun
    /// of feeding a segment built against one [`Readers`] into a different
    /// [`Readers`] whose header doesn't match, this method validates that
    /// `segment.contig()` resolves to the same `tid` *and* that the
    /// header's contig length matches the segment's `contig_last_pos`. If
    /// the contig name resolves to a different tid (or is absent), returns
    /// [`ReaderError::SegmentHeaderMismatch`]; if the contig length differs
    /// (e.g. the segment was built against a 200 Mbp panel but the current
    /// header has the same contig at 100 Mbp), returns
    /// [`ReaderError::SegmentContigLengthMismatch`] rather than letting the
    /// fetch silently return zero records past the shorter contig's end.
    ///
    /// **FASTA contig name.** The contig name carried by the segment is
    /// the BAM-side name; it's passed verbatim to the FASTA reader. If
    /// your BAM uses `chr1` and your FASTA uses `1` (or vice versa), the
    /// FASTA fetch will fail. Either pre-normalize names when building the
    /// reference, or wrap [`Readers`] with your own translation layer.
    ///
    /// When `E != ()`, the extras provider runs once per record at fetch
    /// time, before the engine is constructed.
    ///
    /// **Buffer reuse.** Uses an internal [`RecordStore`] whose capacity
    /// is retained across calls. The returned [`PileupGuard`] borrows
    /// `&mut self.store`; on drop (end of scope, `?`-propagated error,
    /// `break` out of the loop) the populated-then-cleared store is
    /// moved back, so subsequent `pileup()` calls keep the ~39 MB
    /// allocation. No explicit recover step is needed.
    ///
    /// [`PileupColumn::reference_base`]: crate::bam::pileup::PileupColumn::reference_base
    pub fn pileup(&mut self, segment: &Segment) -> Result<PileupGuard<'_, E::Extra>, ReaderError> {
        let tid = segment.tid();
        let start = segment.start();
        let end = segment.end();

        // Header-consistency check: the segment's contig name MUST resolve to
        // the same tid against this Readers' header AND the segment's
        // `contig_last_pos` MUST equal the header's last position for that
        // contig. This catches two foot-guns:
        //
        // * Tid mismatch: shipping a segment from one Readers into another
        //   whose contig order differs.
        // * Length mismatch: shipping a segment built against one reference
        //   panel (e.g. chr1 = 200 Mbp) into another panel (e.g. chr1 = 100
        //   Mbp at the same tid). Without this check, fetches past the
        //   shorter contig's end would silently return zero records.
        let header = self.alignment.header();
        match header.tid(segment.contig().as_str()) {
            Some(t) if t == tid.as_u32() => {}
            _ => {
                return Err(ReaderError::SegmentHeaderMismatch {
                    contig: segment.contig().clone(),
                    expected_tid: tid.as_u32(),
                });
            }
        }
        let header_last_pos = header.target_len(tid.as_u32()).and_then(|len| len.checked_sub(1));
        if header_last_pos != Some(segment.contig_last_pos().as_u64()) {
            return Err(ReaderError::SegmentContigLengthMismatch {
                contig: segment.contig().clone(),
                segment_last_pos: segment.contig_last_pos().as_u64(),
                header_last_pos,
            });
        }

        // Split borrows: fetch writes into `store` while customize is consulted
        // for keep_record. Three disjoint &mut borrows of self.
        let alignment = &mut self.alignment;
        let store = &mut self.store;
        let customize = &mut self.customize;
        alignment.fetch_into_customized(tid.as_u32(), start, end, store, customize)?;

        // Fetch `[start, end]` (inclusive). FASTA APIs expect half-open [start, stop).
        // Use the u64 path so `end == Pos0::max_value()` doesn't truncate the
        // last reference base — `stop = end + 1` is i32::MAX + 1, which doesn't
        // fit in a Pos0 but does fit comfortably in a u64.
        let contig_name = segment.contig();
        let stop_u64 = end.as_u64().saturating_add(1);
        self.fasta
            .fetch_seq_into_u64(contig_name, start.as_u64(), stop_u64, &mut self.fasta_buf)
            .map_err(|source| ReaderError::FastaFetch {
            contig: contig_name.clone(),
            start: start.as_u64(),
            end: end.as_u64(),
            source,
        })?;
        // Convert in-place and copy into the Rc<[Base]> while keeping
        // `fasta_buf` (and its capacity) for the next pileup call.
        let bases: &[Base] = Base::convert_ascii_in_place_as_slice(&mut self.fasta_buf);
        let ref_seq = RefSeq::new(Rc::from(bases), start);

        // Move the populated store into the engine. After this `mem::take`,
        // `self.store` holds a default (empty) RecordStore — the slot the
        // guard's Drop will overwrite with the recovered store on scope exit.
        let store = std::mem::take(&mut self.store);

        let mut engine = PileupEngine::new(store, start, end);
        engine.set_reference_seq(ref_seq);
        Ok(PileupGuard::new(engine, &mut self.store))
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
    /// Uses an internal buffer whose capacity is retained across calls. The
    /// returned `Rc<[Base]>` owns its own copy of the bases — the internal
    /// buffer is converted in place and the slice is copied into the `Rc`
    /// allocation, so we pay one allocation for the `Rc` and zero for the
    /// buffer on subsequent calls.
    pub fn fetch_base_seq(
        &mut self,
        name: &str,
        start: Pos0,
        stop: Pos0,
    ) -> Result<Rc<[Base]>, FastaError> {
        self.fasta.fetch_seq_into(name, start, stop, &mut self.fasta_buf)?;
        // Convert ASCII → Base in place; the &[Base] borrow keeps fasta_buf
        // alive (and its capacity).
        let bases: &[Base] = Base::convert_ascii_in_place_as_slice(&mut self.fasta_buf);
        Ok(Rc::from(bases))
    }

    pub fn alignment(&self) -> &IndexedReader {
        &self.alignment
    }

    pub fn alignment_mut(&mut self) -> &mut IndexedReader {
        &mut self.alignment
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod tests {
    use super::*;
    use crate::reader::segment::Segment;
    use seqair_types::{Pos0, SmolStr};

    fn test_bam_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
    }

    fn test_fasta_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
    }

    // r[verify fasta.fetch.buffer_reuse]
    /// `fetch_base_seq` must keep `fasta_buf`'s capacity across calls so the
    /// next fetch doesn't reallocate. Pre-fix this test would have failed:
    /// `mem::take` left an empty Vec with zero capacity.
    #[test]
    fn fetch_base_seq_retains_buffer_capacity() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let start = Pos0::new(6_100_000).unwrap();
        let stop = Pos0::new(6_101_000).unwrap();
        let _first = readers.fetch_base_seq("chr19", start, stop).unwrap();
        let cap_after_first = readers.fasta_buf.capacity();
        assert!(cap_after_first > 0, "buffer should have grown to hold the fetched region");

        // Second call: capacity must not drop. (Could grow if region is bigger,
        // but for the same region must stay equal.)
        let _second = readers.fetch_base_seq("chr19", start, stop).unwrap();
        assert_eq!(
            cap_after_first,
            readers.fasta_buf.capacity(),
            "second fetch_base_seq must reuse the existing buffer allocation"
        );
    }

    // r[verify pileup.extras.recover_store]
    /// The `PileupGuard` returned by `pileup()` recovers the
    /// [`RecordStore`] into `Readers` on drop, retaining capacity for the
    /// next call. The user is *not* expected to call any explicit recover
    /// step. The test verifies that:
    ///
    /// 1. After a normal pileup loop, `readers.store` has non-zero
    ///    capacity (the guard's Drop wrote back).
    /// 2. After breaking out of the pileup loop early — the foot-gun the
    ///    old explicit `recover_store` API failed on — the same capacity
    ///    is recovered.
    #[test]
    fn pileup_guard_recovers_store_on_drop_even_with_early_break() {
        use std::num::NonZeroU32;
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let opts = SegmentOptions::new(NonZeroU32::new(3_000).unwrap());
        // Known-populated region from `tests/segments.rs`.
        let segments: Vec<_> = readers
            .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
            .unwrap()
            .collect();
        assert!(!segments.is_empty(), "test BAM should yield at least one segment");

        // Full loop: drain the engine, then drop.
        {
            let mut p = readers.pileup(&segments[0]).unwrap();
            while p.pileups().is_some() {}
            // p drops at end of scope.
        }
        let cap_after_full = readers.store.records_capacity();
        assert!(cap_after_full > 0, "guard's Drop must recover the store with non-zero capacity");

        // Early-break loop: pull one column then break, do *not* call any
        // recover step. The guard's Drop must still write the store back.
        {
            let mut p = readers.pileup(&segments[0]).unwrap();
            // Pull one column then break — simulating the "user found what
            // they were looking for and broke out" pattern that the old
            // explicit `recover_store` API failed on. Suppress
            // `clippy::never_loop` because the loop *intentionally* runs
            // at most one iteration to exercise the early-drop path.
            #[allow(clippy::never_loop, reason = "intentional early break to exercise guard Drop")]
            while let Some(_col) = p.pileups() {
                break;
            }
            // p drops here.
        }
        let cap_after_early_break = readers.store.records_capacity();
        assert!(
            cap_after_early_break > 0,
            "guard's Drop must recover the store even when the loop is broken early"
        );
        // Capacity should be at least as large as after the full drain
        // (records are cleared, allocations retained).
        assert_eq!(
            cap_after_early_break, cap_after_full,
            "early-break recovery should retain the same capacity as full-drain recovery"
        );
    }

    // r[verify unified.readers_pileup]
    /// Pileup must reject a `Segment` whose contig name doesn't resolve to
    /// the same tid in this `Readers`' header. Catches the foot-gun where a
    /// `Segment` built against one `Readers` is fed to another.
    #[test]
    fn pileup_rejects_segment_with_mismatched_contig_name() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        // Resolve a real tid so the segment's `tid()` is in range, but use a
        // contig name that doesn't exist in the header.
        let chr19_tid = readers.header().tid("chr19").expect("test BAM has chr19");
        let real_tid =
            crate::reader::ResolveTid::resolve_tid(&chr19_tid, readers.header()).unwrap();
        let bogus_contig: SmolStr = "chr_does_not_exist".into();
        let last_pos =
            Pos0::new(u32::try_from(readers.header().target_len(chr19_tid).unwrap() - 1).unwrap())
                .unwrap();
        let segment = Segment::new(
            real_tid,
            bogus_contig.clone(),
            Pos0::new(0).unwrap(),
            Pos0::new(99).unwrap(),
            0,
            0,
            last_pos,
        )
        .unwrap();

        let err = readers.pileup(&segment).unwrap_err();
        match err {
            ReaderError::SegmentHeaderMismatch { contig, expected_tid } => {
                assert_eq!(contig, bogus_contig);
                assert_eq!(expected_tid, real_tid.as_u32());
            }
            other => panic!("expected SegmentHeaderMismatch, got {other:?}"),
        }
    }

    // r[verify unified.readers_pileup]
    /// The mismatch check also catches the case where contig name and tid
    /// resolve correctly but the segment's `contig_last_pos` does not match
    /// the current header's contig length. This is the foot-gun of building
    /// a `Segment` against one reference panel (e.g. `chr1` = 200 Mbp) and
    /// feeding it to another panel (e.g. `chr1` = 100 Mbp at the same tid).
    /// Without this check, a fetch past the shorter contig's end would
    /// silently return zero records.
    #[test]
    fn pileup_rejects_segment_with_mismatched_contig_length() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").expect("test BAM has chr19");
        let real_tid =
            crate::reader::ResolveTid::resolve_tid(&chr19_tid, readers.header()).unwrap();
        let actual_last_pos =
            Pos0::new(u32::try_from(readers.header().target_len(chr19_tid).unwrap() - 1).unwrap())
                .unwrap();
        // Pretend the segment was built against a shorter version of chr19.
        let bogus_last_pos = Pos0::new(u32::from(actual_last_pos) / 2).unwrap();
        let segment = Segment::new(
            real_tid,
            "chr19".into(),
            Pos0::new(0).unwrap(),
            Pos0::new(99).unwrap(),
            0,
            0,
            bogus_last_pos,
        )
        .unwrap();

        let err = readers.pileup(&segment).unwrap_err();
        match err {
            ReaderError::SegmentContigLengthMismatch {
                contig,
                segment_last_pos,
                header_last_pos,
            } => {
                assert_eq!(contig.as_str(), "chr19");
                assert_eq!(segment_last_pos, bogus_last_pos.as_u64());
                assert_eq!(header_last_pos, Some(actual_last_pos.as_u64()));
            }
            other => panic!("expected SegmentContigLengthMismatch, got {other:?}"),
        }
    }

    // r[verify unified.readers_pileup]
    /// The mismatch check also catches the case where the contig name *does*
    /// exist in the header but resolves to a different numeric tid than the
    /// segment's pre-resolved tid (e.g., a segment built against a different
    /// header where contig order differs).
    #[test]
    fn pileup_rejects_segment_with_mismatched_tid() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        // Find two distinct tids in the header — a segment built with tid
        // for chrA but contig name for chrB triggers the same error path.
        let chr19_tid = readers.header().tid("chr19").expect("test BAM has chr19");
        // Pick a different tid (any other contig); fall back to chr19 again
        // if there's only one contig (then the test trivially passes by name
        // mismatch, not tid).
        let other_tid: u32 = (0..u32::try_from(readers.header().target_count()).unwrap())
            .find(|t| *t != chr19_tid)
            .expect("test BAM has more than one contig");
        let other_name: SmolStr = readers.header().target_name(other_tid).unwrap().into();

        // Resolve a Tid newtype against the wrong numeric id.
        let chr19_real_tid =
            crate::reader::ResolveTid::resolve_tid(&chr19_tid, readers.header()).unwrap();
        let last_pos =
            Pos0::new(u32::try_from(readers.header().target_len(chr19_tid).unwrap() - 1).unwrap())
                .unwrap();
        let segment = Segment::new(
            chr19_real_tid,
            other_name,
            Pos0::new(0).unwrap(),
            Pos0::new(99).unwrap(),
            0,
            0,
            last_pos,
        )
        .unwrap();

        let err = readers.pileup(&segment).unwrap_err();
        assert!(
            matches!(err, ReaderError::SegmentHeaderMismatch { .. }),
            "expected SegmentHeaderMismatch, got {err:?}"
        );
    }
}
