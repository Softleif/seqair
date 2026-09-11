use super::{
    ReaderError,
    indexed::IndexedReader,
    segment::{IntoSegmentTarget, Segment, SegmentOptions, Segments},
};
use crate::{
    bam::{
        BamHeader,
        depth_cap::{DepthCap, DepthLimit},
        pileup::{PileupEngine, PileupGuard, RefSeq},
        record_store::{CustomizeRecordStore, Prepared, RecordStore},
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
/// value — user code defines per-record filtering (`filter`) and
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
///     let mut pileup = readers.pileup(segment, seqair::reader::DepthLimit::Unlimited).run()?;
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
///                 PileupOp::SoftClip { base, qual, .. } => { let _ = (base, qual); }
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
///                 let mut pileup = forked.pileup(seg, seqair::reader::DepthLimit::Unlimited).run().unwrap();
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
    /// Pooled pileup scratch buffers, reused across `pileup()` calls.
    pub(crate) pileup_scratch: crate::bam::pileup::PileupScratch,
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

    /// Assemble a reference-aware reader from an already-open alignment handle
    /// and FASTA reader.
    ///
    /// This is the composition primitive behind [`open`](Self::open): it does
    /// no I/O and performs no format detection. Use it when you opened an
    /// [`IndexedReader`] first (e.g. to stream records) and later decided you
    /// want pileup, or to supply a FASTA reader obtained from elsewhere. The
    /// ergonomic spelling is [`IndexedReader::into_pileup`].
    // r[impl unified.readers_from_parts]
    pub fn from_parts(alignment: IndexedReader, fasta: IndexedFastaReader) -> Self {
        Self::from_parts_customized(alignment, fasta, ())
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
        let alignment = IndexedReader::open_with_reference(alignment_path, fasta_path)?;
        Ok(Self::from_parts_customized(alignment, fasta, customize))
    }

    /// Assemble a reader from already-open parts and a [`CustomizeRecordStore`].
    ///
    /// The generic form of [`from_parts`](Readers::from_parts); the shared
    /// builder behind [`open_customized`](Self::open_customized) and
    /// [`fork`](Self::fork). Does no I/O.
    // r[impl unified.readers_from_parts]
    pub fn from_parts_customized(
        alignment: IndexedReader,
        fasta: IndexedFastaReader,
        customize: E,
    ) -> Self {
        Readers {
            alignment,
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
            customize,
            pileup_scratch: crate::bam::pileup::PileupScratch::default(),
        }
    }

    /// Fork both the alignment reader and the FASTA reader.
    ///
    /// The customize value is cloned so each fork has its own copy.
    // r[impl unified.readers_fork]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        let alignment = self.alignment.fork()?;
        let fasta = self.fasta.fork().map_err(|source| ReaderError::FastaFork { source })?;
        Ok(Self::from_parts_customized(alignment, fasta, self.customize.clone()))
    }

    // r[impl unified.readers_accessors+1]
    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore<E::Extra>,
    ) -> Result<usize, ReaderError> {
        self.alignment
            .fetch_into_customized(tid, start, end, store, &mut self.customize)
            .map(|c| c.kept)
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
    /// pre-resolved [`Tid`](crate::reader::Tid) or `u32`, a parsed
    /// [`RegionString`](seqair_types::RegionString), an
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
    ///     let mut p = readers.pileup(&seg, seqair::reader::DepthLimit::Unlimited).run()?;
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

        // With a byte budget set, tile positionally and then subdivide each
        // tile so its estimated compressed load stays within budget. This is
        // eager (it queries the index per tile), but the per-tile cost is an
        // in-memory BAI lookup and the result is what callers collect anyway.
        // Backends that can't estimate bytes (CRAM) report 0 here, so the
        // split is a no-op and tiles pass through unchanged.
        let Some(budget) = opts.max_bytes() else {
            return Ok(Segments::new(ranges, opts));
        };
        let budget = budget.get();
        let overlap = opts.overlap();
        let mut estimate = |tid: u32, start: Pos0, end: Pos0| {
            self.alignment.estimate_region_bytes(tid, start, end).unwrap_or(0)
        };
        let mut out = Vec::new();
        for seg in Segments::new(ranges, opts) {
            super::segment::split_segment_by_bytes(&seg, overlap, budget, &mut estimate, &mut out);
        }
        Ok(Segments::precomputed(out))
    }

    /// Estimate the **compressed** bytes a `[start, end]` region query would
    /// load. `None` for CRAM (its slice reader bounds memory differently). This
    /// is the same estimate [`segments`](Self::segments) budgets against.
    #[must_use]
    pub fn estimate_region_bytes(&self, tid: u32, start: Pos0, end: Pos0) -> Option<u64> {
        self.alignment.estimate_region_bytes(tid, start, end)
    }

    // r[impl unified.readers_pileup+1]
    // r[impl unified.pileup_plan]
    /// Plan a pileup over `segment`: configure it, then [`run`](Pileup::run).
    ///
    /// ```ignore
    /// readers.pileup(&segment, depth).run()?;                       // plain
    /// readers.pileup(&segment, depth).mutate(realign).run()?;       // + hook
    /// readers.pileup(&segment, depth).with_reference(r).run()?;     // + reference
    /// readers.pileup(&segment, depth).with_reference(r).mutate(realign).run()?;
    /// ```
    ///
    /// Nothing is read until `run`, and the options are independent — which is
    /// why this is a plan and not four methods: the last line above had no
    /// spelling when each combination needed its own entry point.
    ///
    /// `segment` must come from [`segments`](Self::segments) — there is no
    /// other way to construct one. The segment's `tid`, `start`, `end`, and
    /// `contig` name are used directly.
    ///
    /// **Header consistency.** [`Segment`] is `Send + Clone`, so it can be
    /// shipped between threads or across `fork()`s. To catch the foot-gun
    /// of feeding a segment built against one [`Readers`] into a different
    /// [`Readers`] whose header doesn't match, `run` validates that
    /// `segment.contig()` resolves to the same `tid` *and* that the
    /// header's contig length matches the segment's `contig_last_pos`. If
    /// the contig name resolves to a different tid (or is absent), it returns
    /// [`ReaderError::SegmentHeaderMismatch`]; if the contig length differs
    /// (e.g. the segment was built against a 200 Mbp panel but the current
    /// header has the same contig at 100 Mbp), it returns
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
    /// moved back, so subsequent pileups keep the ~39 MB allocation. No
    /// explicit recover step is needed.
    pub fn pileup<'a>(&'a mut self, segment: &'a Segment, depth: DepthLimit) -> Pileup<'a, E> {
        Pileup { readers: self, segment, depth, reference: None, mutate: None, cover_reads: false }
    }

    fn pileup_impl<F>(
        &mut self,
        segment: &Segment,
        depth: DepthLimit,
        mutate: Option<F>,
        supplied_ref: Option<RefSeq>,
        cover_reads: bool,
    ) -> Result<PileupGuard<'_, E::Extra>, ReaderError>
    where
        F: FnMut(&mut RecordStore<E::Extra>, &RefSeq),
    {
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

        match depth.per_column() {
            None => {
                // Split borrows: fetch writes into `store` while customize is
                // consulted for filter. Three disjoint &mut borrows of self.
                let alignment = &mut self.alignment;
                let store = &mut self.store;
                let customize = &mut self.customize;
                alignment.fetch_into_customized(tid.as_u32(), start, end, store, customize)?;
            }
            Some(cap) => {
                // Cap reads overlapping any column at load time so an
                // ultra-deep pileup can't blow up the `RecordStore`. The
                // customize is cloned for the fetch and restored afterwards so
                // any state it accumulates survives.
                let mut capped = DepthCap::with_cap(self.customize.clone(), Some(cap));
                self.alignment.fetch_into_customized(
                    tid.as_u32(),
                    start,
                    end,
                    &mut self.store,
                    &mut capped,
                )?;
                self.customize = capped.into_inner();
            }
        }

        // r[impl unified.pileup_reference_covers_reads]
        // The records are in hand, so the span they actually cover is known
        // exactly — no padding constant, and a read that stops inside the
        // segment widens nothing.
        let (ref_start, ref_end) = if cover_reads {
            span_covering_records(&self.store, start, end, segment.contig_last_pos())
        } else {
            (start, end)
        };

        let ref_seq = match supplied_ref {
            // r[impl unified.readers_pileup_supplied_reference+1]
            Some(ref_seq) => {
                // An uncovered position reads as `Base::Unknown`, which is
                // indistinguishable from a genuine `N` in the reference — every
                // comparison against it would silently become a mismatch.
                let have_start = ref_seq.start_pos().as_u64();
                let have_last = have_start.saturating_add(ref_seq.len() as u64).saturating_sub(1);
                if have_start > start.as_u64() || have_last < end.as_u64() {
                    return Err(ReaderError::SuppliedReferenceTooSmall {
                        contig: segment.contig().clone(),
                        ref_start: have_start,
                        ref_end: have_last,
                        segment_start: start.as_u64(),
                        segment_end: end.as_u64(),
                    });
                }
                // r[impl unified.pileup_reference_covers_reads]
                // Nothing can be widened here, so the option becomes the
                // requirement it always was: silently handing the hook a
                // reference that stops short is the failure it exists to close.
                if cover_reads && (have_start > ref_start.as_u64() || have_last < ref_end.as_u64())
                {
                    return Err(ReaderError::SuppliedReferenceMissesReads {
                        contig: segment.contig().clone(),
                        ref_start: have_start,
                        ref_end: have_last,
                        reads_start: ref_start.as_u64(),
                        reads_end: ref_end.as_u64(),
                    });
                }
                ref_seq
            }
            // Fetch `[ref_start, ref_end]` (inclusive). FASTA APIs expect half-open
            // [start, stop). Use the u64 path so `end == Pos0::max_value()` doesn't
            // truncate the last reference base — `stop = end + 1` is i32::MAX + 1,
            // which doesn't fit in a Pos0 but does fit comfortably in a u64.
            None => {
                let contig_name = segment.contig();
                let stop_u64 = ref_end.as_u64().saturating_add(1);
                self.fasta
                    .fetch_seq_into_u64(
                        contig_name,
                        ref_start.as_u64(),
                        stop_u64,
                        &mut self.fasta_buf,
                    )
                    .map_err(|source| ReaderError::FastaFetch {
                        contig: contig_name.clone(),
                        start: ref_start.as_u64(),
                        end: ref_end.as_u64(),
                        source,
                    })?;
                // Convert in-place and copy into the Rc<[Base]> while keeping
                // `fasta_buf` (and its capacity) for the next pileup call.
                let bases: &[Base] = Base::convert_ascii_in_place_as_slice(&mut self.fasta_buf);
                RefSeq::new(Rc::from(bases), ref_start)
            }
        };

        // r[impl unified.readers_pileup_store_mutation+4]
        // In-place realignment hook: rewrite (pos, CIGAR) on the fetched store.
        // Runs after the reference fetch so the hook can score against the
        // same bases the engine will report. Restoring position order is
        // `prepare_for_pileup`'s job below — the hook only has to mark the
        // store as reordered, which `set_alignment` does.
        if let Some(mut mutate) = mutate {
            mutate(&mut self.store, &ref_seq);
        }

        // Move the populated store into the engine. After this `mem::take`,
        // `self.store` holds a default (empty) RecordStore — the slot the
        // guard's Drop will overwrite with the recovered store on scope exit.
        let store = std::mem::take(&mut self.store);
        let scratch = std::mem::take(&mut self.pileup_scratch);

        // r[impl record_store.link_mates.invalidated]
        // r[impl record_store.link_mates.stats]
        // Sorts (the realignment hook above may have moved records) and links,
        // in that order — a mate index is a store index and cannot survive a
        // reorder. Producing the engine's input is the only way to do this, so
        // neither step can be skipped.
        let Prepared { input, stats: mate_links } = store.prepare_for_pileup();
        if !mate_links.is_clean() {
            // Repeated primary qnames, or the same record present twice.
            // Neither is valid input; say so once per fetch, not per column.
            tracing::debug!(
                contig = %segment.contig(),
                start = start.as_u64(),
                pairs = mate_links.pairs,
                ambiguous_qnames = mate_links.ambiguous_qnames,
                "qnames are not unique here; those reads are not mate-deduplicated"
            );
        }

        let mut engine = PileupEngine::with_scratch(input, start, end, scratch);
        engine.set_reference_seq(ref_seq);
        Ok(PileupGuard::new(engine, &mut self.store, Some(&mut self.pileup_scratch)))
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

// r[impl unified.pileup_reference_covers_reads]
/// The span a reference must cover for every record in `store` to have its
/// bases: the segment `[start, end]`, widened to the records' own extent and
/// clamped to the contig's last position.
///
/// Never narrower than the segment — the engine reports a column for every
/// position in it — and never past the contig, where the FASTA has nothing.
/// Unmapped records are included: they carry a position too, and excluding
/// them would make the span depend on a distinction the caller cannot see
/// while walking `store.records()`.
fn span_covering_records<U>(
    store: &RecordStore<U>,
    start: Pos0,
    end: Pos0,
    contig_last: Pos0,
) -> (Pos0, Pos0) {
    let mut first = start;
    let mut last = end;
    for rec in store.records() {
        first = first.min(rec.pos);
        last = last.max(rec.end_pos);
    }
    // `end <= contig_last` for any segment the header produced, so clamping
    // `last` can never pull it back inside the segment.
    (first, last.min(contig_last))
}

/// A planned pileup: what [`Readers::pileup`] returns, configured by
/// [`with_reference`](Self::with_reference) and [`mutate`](Self::mutate) and
/// executed by [`run`](Self::run).
///
/// Nothing is read while the plan is being built — it holds the borrow of the
/// [`Readers`] and the options, and `run` is where the fetch, the reference,
/// the hook and `prepare_for_pileup` happen, in that order.
///
/// `F` is the mutator's type, defaulted to a function pointer for the plans
/// that have none; [`mutate`](Self::mutate) replaces it with the closure's own
/// type, so a hook costs no indirection.
// r[impl unified.pileup_plan]
#[must_use = "a pileup plan does nothing until `run` is called"]
pub struct Pileup<
    'a,
    E: CustomizeRecordStore,
    F = fn(&mut RecordStore<<E as CustomizeRecordStore>::Extra>, &RefSeq),
> {
    readers: &'a mut Readers<E>,
    segment: &'a Segment,
    depth: DepthLimit,
    reference: Option<RefSeq>,
    mutate: Option<F>,
    cover_reads: bool,
}

impl<'a, E: CustomizeRecordStore, F> Pileup<'a, E, F>
where
    F: FnMut(&mut RecordStore<E::Extra>, &RefSeq),
{
    // r[impl unified.readers_pileup_supplied_reference+1]
    /// Drive the engine from a reference the caller already holds instead of
    /// reading it from the FASTA.
    ///
    /// A caller that fetched a region's bases for its own use — and then tiles
    /// that region into several segments — otherwise pays one FASTA seek and
    /// decode per segment for bases it is already holding. [`RefSeq`] resolves
    /// absolute positions, so one reference spanning the whole region serves
    /// every segment within it.
    ///
    /// `ref_seq` MUST cover the segment; a reference that falls short is
    /// rejected by [`run`](Self::run) with
    /// [`ReaderError::SuppliedReferenceTooSmall`] rather than silently reading
    /// the uncovered positions as [`Base::Unknown`].
    ///
    /// This is also the reference a [`mutate`](Self::mutate) hook receives.
    ///
    /// [`Base::Unknown`]: seqair_types::Base::Unknown
    pub fn with_reference(mut self, ref_seq: RefSeq) -> Self {
        self.reference = Some(ref_seq);
        self
    }

    // r[impl unified.pileup_reference_covers_reads]
    /// Make the reference cover every record that was fetched, not just the
    /// segment.
    ///
    /// A read at a tile's edge reaches past it, and its outer bases then have
    /// no reference: [`RefSeq::try_base_at`] answers `None` there, and
    /// [`RefSeq::base_at`] answers [`Base::Unknown`], which a comparison reads
    /// as a mismatch. On this crate's 80 bp fixture that is 27 % of the reads
    /// at a 500 bp tile and 3 % at 5 kb — small per tile, but it is the reads
    /// at the edges, and a hook that realigns or rescores whole reads needs
    /// their bases.
    ///
    /// The widening is exact rather than a padding constant: the records are
    /// already loaded when the reference is read, so the span is their own
    /// `min(pos)`..`max(end_pos)`, unioned with the segment and clamped to the
    /// contig. A tile whose reads all stop inside it fetches exactly what it
    /// would have without this.
    ///
    /// Columns are unaffected. [`RefSeq`] resolves absolute positions, so
    /// [`PileupColumn::reference_base`] reports the same base at every column
    /// of the segment either way; only what lies *outside* the segment changes
    /// from absent to present.
    ///
    /// With [`with_reference`](Self::with_reference) there is nothing to
    /// widen, so this becomes a requirement on the reference the caller
    /// supplied: one that stops short of the records is rejected with
    /// [`ReaderError::SuppliedReferenceMissesReads`] rather than quietly
    /// handing the hook a reference that does not reach.
    ///
    /// [`PileupColumn::reference_base`]: crate::bam::pileup::PileupColumn::reference_base
    /// [`Base::Unknown`]: seqair_types::Base::Unknown
    pub fn reference_covers_reads(mut self) -> Self {
        self.cover_reads = true;
        self
    }

    // r[impl unified.readers_pileup_store_mutation+4]
    /// Run `mutate` on the freshly fetched [`RecordStore`] **before** the
    /// pileup engine is built, then re-sort the store by position.
    ///
    /// This is the hook for in-place local realignment: the mutator may call
    /// [`RecordStore::set_alignment`] to rewrite read (pos, CIGAR) pairs (e.g.
    /// from a POA consensus), and the subsequent pileup sees the corrected
    /// alignments. Buffer reuse is preserved exactly as without a hook.
    ///
    /// The mutator also receives the reference: the very [`RefSeq`] the engine
    /// is built with — whether fetched from the FASTA or supplied by
    /// [`with_reference`](Self::with_reference) — so `ref_seq.base_at(pos)`
    /// inside the hook is what [`PileupColumn::reference_base`] later reports
    /// at that column. A hook that normalises or rescores alignments against
    /// the reference reads it from here and never loads its own; one that does
    /// not need it ignores the argument.
    ///
    /// The reference covers exactly the segment, not the reads. A record
    /// reaching past either end has bases the `RefSeq` does not hold; read
    /// those with [`RefSeq::try_base_at`] and treat `None` as unavailable —
    /// [`RefSeq::base_at`] returns [`Base::Unknown`] there, which is
    /// indistinguishable from a genuine `N`.
    ///
    /// The mutator must preserve each record's query length (`set_alignment`
    /// enforces this); positions may change freely since the store is
    /// re-sorted afterward. It runs only once the reference is in hand: when
    /// the FASTA fetch fails, the hook has not run and the error is all the
    /// caller gets.
    ///
    /// [`PileupColumn::reference_base`]: crate::bam::pileup::PileupColumn::reference_base
    /// [`Base::Unknown`]: seqair_types::Base::Unknown
    pub fn mutate<G>(self, mutate: G) -> Pileup<'a, E, G>
    where
        G: FnMut(&mut RecordStore<E::Extra>, &RefSeq),
    {
        Pileup {
            readers: self.readers,
            segment: self.segment,
            depth: self.depth,
            reference: self.reference,
            mutate: Some(mutate),
            cover_reads: self.cover_reads,
        }
    }

    // r[impl unified.readers_pileup+1]
    // r[impl unified.pileup_plan]
    /// Fetch the records and the reference, run the hook if there is one, and
    /// hand back a [`PileupGuard`] that derefs to a [`PileupEngine`] ready for
    /// iteration.
    ///
    /// See [`Readers::pileup`] for what is validated here and what is reused
    /// between calls.
    pub fn run(self) -> Result<PileupGuard<'a, E::Extra>, ReaderError> {
        let Self { readers, segment, depth, reference, mutate, cover_reads } = self;
        readers.pileup_impl(segment, depth, mutate, reference, cover_reads)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod tests {
    use super::*;
    use crate::bam::cigar::{CigarOp, CigarOpType};
    use crate::bam::record_store::RecordIdx;
    use crate::reader::segment::Segment;
    use seqair_types::{Pos0, SmolStr};

    fn test_bam_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
    }

    fn test_fasta_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
    }

    // r[verify unified.segment_byte_budget]
    /// `segments()` with a byte budget subdivides a region (spanning several
    /// index leaf bins) so each emitted segment's estimated load stays within
    /// budget, while the cores still tile the requested range exactly.
    #[test]
    fn segments_byte_budget_subdivides_region() {
        use std::num::{NonZeroU32, NonZeroU64};
        let readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        // chr19 reads span ~6.10–6.14 Mb (several 16 kb leaf bins).
        let start = Pos0::new(6_100_000).unwrap();
        let end = Pos0::new(6_145_000).unwrap();
        let tid = readers.header().tid("chr19").expect("chr19 in header");
        let total = readers.estimate_region_bytes(tid, start, end).expect("BAM estimate");
        assert!(total > 0, "test region should contain reads");

        // One positional tile over the whole region (no byte budget).
        let big = NonZeroU32::new(1_000_000).unwrap();
        let positional: Vec<_> = readers
            .segments(("chr19", start, end), SegmentOptions::new(big).without_byte_budget())
            .unwrap()
            .collect();
        assert_eq!(positional.len(), 1, "without a budget the region is one tile");

        // A budget below the whole-region size (but above one ~16 kb leaf bin)
        // must subdivide it into several smaller, in-budget segments.
        let budget = NonZeroU64::new((total / 2).max(1)).unwrap();
        let budgeted: Vec<_> = readers
            .segments(("chr19", start, end), SegmentOptions::new(big).with_max_bytes(budget))
            .unwrap()
            .collect();
        assert!(budgeted.len() > 1, "budget {budget} (of {total}) must split the region");

        // Cores tile [start, end] exactly.
        assert_eq!(*budgeted.first().unwrap().core_range().start(), start);
        assert_eq!(*budgeted.last().unwrap().core_range().end(), end);
        for w in budgeted.windows(2) {
            assert_eq!(
                w[0].core_range().end().as_u64() + 1,
                w[1].core_range().start().as_u64(),
                "cores must be contiguous"
            );
        }

        // The peak segment load is strictly below the whole-region load (the
        // point of subdivision), and no segment exceeds the budget unless a
        // single index leaf bin alone does (none here — budget > one bin).
        for s in &budgeted {
            let bytes = readers.estimate_region_bytes(tid, s.start(), s.end()).unwrap();
            assert!(bytes < total, "segment {bytes} B not below whole-region {total} B");
            assert!(
                bytes <= budget.get(),
                "segment {:?} = {bytes} B over budget {budget}",
                (s.start().as_u64(), s.end().as_u64())
            );
        }
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
            let mut p = readers.pileup(&segments[0], DepthLimit::Unlimited).run().unwrap();
            while p.pileups().is_some() {}
            // p drops at end of scope.
        }
        let cap_after_full = readers.store.records_capacity();
        assert!(cap_after_full > 0, "guard's Drop must recover the store with non-zero capacity");

        // Early-break loop: pull one column then break, do *not* call any
        // recover step. The guard's Drop must still write the store back.
        {
            let mut p = readers.pileup(&segments[0], DepthLimit::Unlimited).run().unwrap();
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

    /// A region known to contain reads in the test BAM.
    #[cfg(test)]
    fn realign_test_segment(readers: &Readers) -> Segment {
        use std::num::NonZeroU32;
        let opts = SegmentOptions::new(NonZeroU32::new(5_000).unwrap());
        readers
            .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
            .unwrap()
            .next()
            .expect("test BAM should yield a segment")
    }

    #[cfg(test)]
    fn pileup_profile(p: &mut PileupGuard<'_, ()>) -> Vec<(u64, usize)> {
        let mut out = Vec::new();
        while let Some(col) = p.pileups() {
            out.push((col.pos().as_u64(), col.depth()));
        }
        out
    }

    /// Fetch the reference for a whole span as a single `RefSeq`, the way a
    /// caller that already needs those bases would hold them.
    ///
    /// `end` is inclusive, matching `Segment::end()`; `fetch_base_seq` is
    /// half-open, hence the `+ 1`.
    #[cfg(test)]
    fn whole_span_reference(readers: &mut Readers, contig: &str, start: Pos0, end: Pos0) -> RefSeq {
        let stop = Pos0::try_from(end.as_u64().saturating_add(1)).unwrap();
        let bases = readers.fetch_base_seq(contig, start, stop).unwrap();
        RefSeq::new(bases, start)
    }

    // r[verify unified.readers_pileup_supplied_reference+1]
    /// Supplying the reference must not change a single column: the engine sees
    /// the same bases it would have fetched itself.
    #[test]
    fn pileup_with_reference_matches_pileup() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        let expected = {
            let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
            pileup_profile(&mut p)
        };
        let expected_bases = {
            let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
            let mut out = Vec::new();
            while let Some(col) = p.pileups() {
                out.push(col.reference_base());
            }
            out
        };

        let ref_seq = whole_span_reference(&mut readers, "chr19", segment.start(), segment.end());
        let mut p =
            readers.pileup(&segment, DepthLimit::Unlimited).with_reference(ref_seq).run().unwrap();
        let mut actual = Vec::new();
        let mut actual_bases = Vec::new();
        while let Some(col) = p.pileups() {
            actual.push((col.pos().as_u64(), col.depth()));
            actual_bases.push(col.reference_base());
        }

        assert_eq!(actual, expected, "columns must be identical to a self-fetched pileup");
        assert_eq!(actual_bases, expected_bases, "reference bases must be identical");
    }

    // r[verify unified.readers_pileup_supplied_reference+1]
    /// One reference spanning a whole region serves every segment tiled inside
    /// it — `RefSeq` resolves absolute positions, so no per-segment slicing.
    #[test]
    fn one_reference_serves_every_segment_in_the_span() {
        use std::num::NonZeroU32;
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let (start, end) = (Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap());

        // Tile the span finely enough that several segments come out of it.
        let opts = SegmentOptions::new(NonZeroU32::new(500).unwrap());
        let segments: Vec<Segment> =
            readers.segments(("chr19", start, end), opts).unwrap().collect();
        assert!(segments.len() > 1, "span should tile into several segments");

        let expected: Vec<Vec<(u64, usize)>> = segments
            .iter()
            .map(|seg| {
                let mut p = readers.pileup(seg, DepthLimit::Unlimited).run().unwrap();
                pileup_profile(&mut p)
            })
            .collect();

        let ref_seq = whole_span_reference(&mut readers, "chr19", start, end);
        for (seg, want) in segments.iter().zip(&expected) {
            let mut p = readers
                .pileup(seg, DepthLimit::Unlimited)
                .with_reference(ref_seq.clone())
                .run()
                .unwrap();
            assert_eq!(&pileup_profile(&mut p), want, "segment {:?} diverged", seg.start());
        }
    }

    // r[verify unified.pileup_plan]
    // r[verify unified.readers_pileup_supplied_reference+1]
    // r[verify unified.readers_pileup_store_mutation+4]
    /// The combination the plan exists for, and which no method could express:
    /// a caller that holds the region's bases *and* rewrites the store against
    /// them. The hook must see the supplied reference — not a second fetch —
    /// and the columns must be the ones the same hook produces without it.
    #[test]
    fn a_supplied_reference_and_a_hook_compose() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);
        let (start, end) = (segment.start(), segment.end());

        // The same realignment, run once with the FASTA fetch and once with a
        // reference the caller already holds.
        let shift = |store: &mut RecordStore, ref_seq: &RefSeq, seen: &mut Vec<Base>| {
            *seen = (start.as_u64()..=end.as_u64())
                .map(|p| ref_seq.base_at(Pos0::try_from(p).unwrap()))
                .collect();
            let mut plan: Vec<(RecordIdx, Pos0, Vec<CigarOp>)> = Vec::new();
            for idx in store.indices() {
                let rec = store.record(idx);
                let Ok(cigar) = rec.cigar(store) else { continue };
                let Some(first) = cigar.first() else { continue };
                if rec.flags.is_unmapped()
                    || first.op_type() != CigarOpType::Match
                    || first.len() < 2
                {
                    continue;
                }
                let mut new = vec![
                    CigarOp::new(CigarOpType::SoftClip, 1),
                    CigarOp::new(CigarOpType::Match, first.len() - 1),
                ];
                new.extend_from_slice(cigar.get(1..).unwrap_or_default());
                let Ok(new_pos) = Pos0::try_from(rec.pos.as_u64().saturating_add(1)) else {
                    continue;
                };
                plan.push((idx, new_pos, new));
            }
            for (idx, pos, cigar) in &plan {
                store.set_alignment(*idx, *pos, cigar).unwrap();
            }
            plan.len()
        };

        let (mut fetched_ref, mut fetched_moved) = (Vec::new(), 0);
        let from_fasta = {
            let mut p = readers
                .pileup(&segment, DepthLimit::Unlimited)
                .mutate(|store, ref_seq| fetched_moved = shift(store, ref_seq, &mut fetched_ref))
                .run()
                .unwrap();
            pileup_profile(&mut p)
        };
        assert!(fetched_moved > 0, "the fixture has reads with a leading M");

        let supplied = whole_span_reference(&mut readers, "chr19", start, end);
        let (mut hook_ref, mut hook_moved) = (Vec::new(), 0);
        let from_supplied = {
            let mut p = readers
                .pileup(&segment, DepthLimit::Unlimited)
                .with_reference(supplied)
                .mutate(|store, ref_seq| hook_moved = shift(store, ref_seq, &mut hook_ref))
                .run()
                .unwrap();
            pileup_profile(&mut p)
        };

        assert_eq!(hook_moved, fetched_moved, "the hook ran over the same records");
        assert_eq!(hook_ref, fetched_ref, "the hook read the supplied reference, base for base");
        assert_eq!(
            from_supplied, from_fasta,
            "the columns must not depend on where the bases came from"
        );
    }

    // r[verify unified.pileup_reference_covers_reads]
    /// The widened reference must reach every read's ends — and must be
    /// exactly as wide as the reads are, not a padding guess. The fixture's
    /// reads are 80 bp, so a small tile has reads hanging off both edges.
    #[test]
    fn reference_covers_reads_reaches_every_read() {
        use std::num::NonZeroU32;
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let (from, to) = (Pos0::new(6_103_500).unwrap(), Pos0::new(6_104_500).unwrap());
        let opts = SegmentOptions::new(NonZeroU32::new(500).unwrap());
        let segments: Vec<Segment> = readers.segments(("chr19", from, to), opts).unwrap().collect();

        let mut checked = 0usize;
        for segment in &segments {
            let (seg_start, seg_end) = (segment.start(), segment.end());
            let mut overhanging = 0usize;
            let mut want = None;
            let mut holds = None;
            let mut p = readers
                .pileup(segment, DepthLimit::Unlimited)
                .reference_covers_reads()
                .mutate(|store, ref_seq| {
                    // Every record's own ends must resolve, which is the whole
                    // point: without the option these are `None`.
                    for rec in store.records() {
                        assert!(
                            ref_seq.try_base_at(rec.pos).is_some(),
                            "read start {} has no reference",
                            rec.pos.as_u64()
                        );
                        assert!(
                            ref_seq.try_base_at(rec.end_pos).is_some(),
                            "read end {} has no reference",
                            rec.end_pos.as_u64()
                        );
                        if rec.pos < seg_start || rec.end_pos > seg_end {
                            overhanging += 1;
                        }
                    }
                    // Exact, not padded: the span is the records' own extent
                    // unioned with the segment.
                    let lo = store.records().map(|r| r.pos).min().unwrap_or(seg_start);
                    let hi = store.records().map(|r| r.end_pos).max().unwrap_or(seg_end);
                    want = Some((lo.min(seg_start).as_u64(), hi.max(seg_end).as_u64()));
                    let first = ref_seq.start_pos().as_u64();
                    holds = Some((first, first + ref_seq.len() as u64 - 1));
                })
                .run()
                .unwrap();
            let _ = pileup_profile(&mut p);
            drop(p);

            assert_eq!(holds, want, "segment at {}", seg_start.as_u64());
            checked += overhanging;
        }
        assert!(checked > 0, "500 bp tiles over 80 bp reads must produce overhanging reads");

        // The option is doing the work: without it, the same reads have ends
        // the hook cannot read.
        let mut unreachable_ends = 0usize;
        for segment in &segments {
            let mut p = readers
                .pileup(segment, DepthLimit::Unlimited)
                .mutate(|store, ref_seq| {
                    unreachable_ends += store
                        .records()
                        .filter(|rec| {
                            ref_seq.try_base_at(rec.pos).is_none()
                                || ref_seq.try_base_at(rec.end_pos).is_none()
                        })
                        .count();
                })
                .run()
                .unwrap();
            let _ = pileup_profile(&mut p);
        }
        assert_eq!(
            unreachable_ends, checked,
            "every overhanging read must be unreachable without the option"
        );
    }

    // r[verify unified.pileup_reference_covers_reads]
    /// The columns are the same either way — `RefSeq` resolves absolute
    /// positions, so only what lies outside the segment changes.
    #[test]
    fn reference_covers_reads_does_not_change_any_column() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        let bases = |p: &mut PileupGuard<'_, ()>| {
            let mut out = Vec::new();
            while let Some(col) = p.pileups() {
                out.push((col.pos().as_u64(), col.reference_base(), col.depth()));
            }
            out
        };
        let narrow = {
            let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
            bases(&mut p)
        };
        let wide = {
            let mut p = readers
                .pileup(&segment, DepthLimit::Unlimited)
                .reference_covers_reads()
                .run()
                .unwrap();
            bases(&mut p)
        };
        assert_eq!(wide, narrow);
    }

    // r[verify unified.pileup_reference_covers_reads]
    /// With a supplied reference nothing can be widened, so the option is a
    /// requirement: one that covers the segment but stops short of the reads
    /// is a typed error, not a silent `None` inside the hook.
    #[test]
    fn reference_covers_reads_rejects_a_supplied_reference_that_stops_at_the_segment() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);
        let exactly_the_segment =
            whole_span_reference(&mut readers, "chr19", segment.start(), segment.end());

        // It covers the segment, so the plain path accepts it...
        assert!(
            readers
                .pileup(&segment, DepthLimit::Unlimited)
                .with_reference(exactly_the_segment.clone())
                .run()
                .is_ok()
        );
        // ...and asking for the reads' coverage rejects the same reference.
        assert!(matches!(
            readers
                .pileup(&segment, DepthLimit::Unlimited)
                .with_reference(exactly_the_segment)
                .reference_covers_reads()
                .run(),
            Err(ReaderError::SuppliedReferenceMissesReads { .. })
        ));

        // A reference that does reach the reads is accepted.
        let wide = whole_span_reference(
            &mut readers,
            "chr19",
            Pos0::new(segment.start().as_u32().saturating_sub(1_000)).unwrap(),
            Pos0::try_from(segment.end().as_u64().saturating_add(1_000)).unwrap(),
        );
        assert!(
            readers
                .pileup(&segment, DepthLimit::Unlimited)
                .with_reference(wide)
                .reference_covers_reads()
                .run()
                .is_ok()
        );
    }

    // r[verify unified.pileup_plan]
    /// Building a plan reads nothing: the store is still empty, and the
    /// segment is not validated until `run`. A plan that is dropped must leave
    /// the `Readers` exactly as it found it.
    #[test]
    fn a_plan_reads_nothing_until_run() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        let mut hook_ran = false;
        let plan = readers.pileup(&segment, DepthLimit::Unlimited).mutate(|_, _| hook_ran = true);
        drop(plan);
        assert!(!hook_ran, "a dropped plan must not run its hook");
        assert!(readers.store.is_empty(), "a dropped plan must not have fetched");

        // And the Readers is still usable, with its buffers intact.
        let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
        assert!(!pileup_profile(&mut p).is_empty(), "the Readers still works after a dropped plan");
    }

    // r[verify unified.readers_pileup_supplied_reference+1]
    /// A reference that falls short of the segment must be a typed error, not a
    /// run of `Unknown` bases that silently mismatch everything.
    #[test]
    fn pileup_with_short_reference_is_rejected() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        // Starts where the segment starts but stops well before its end.
        let short_end = Pos0::try_from(segment.start().as_u64().saturating_add(10)).unwrap();
        let short = whole_span_reference(&mut readers, "chr19", segment.start(), short_end);
        assert!(matches!(
            readers.pileup(&segment, DepthLimit::Unlimited).with_reference(short).run(),
            Err(ReaderError::SuppliedReferenceTooSmall { .. })
        ));

        // Covers the segment's end but begins after its start.
        let late_start = Pos0::try_from(segment.start().as_u64().saturating_add(10)).unwrap();
        let late = whole_span_reference(&mut readers, "chr19", late_start, segment.end());
        assert!(matches!(
            readers.pileup(&segment, DepthLimit::Unlimited).with_reference(late).run(),
            Err(ReaderError::SuppliedReferenceTooSmall { .. })
        ));
    }

    // r[verify unified.readers_pileup_store_mutation+4]
    /// A no-op `pileup_with` must produce exactly the same columns as `pileup`:
    /// the hook is transparent when the mutator changes nothing, and it sees the
    /// fully fetched store.
    #[test]
    fn pileup_with_noop_matches_pileup() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        let baseline = {
            let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
            pileup_profile(&mut p)
        };
        assert!(!baseline.is_empty(), "region should pile up some columns");

        let mut seen_records = 0usize;
        let mutated = {
            let mut p = readers
                .pileup(&segment, DepthLimit::Unlimited)
                .mutate(|store, _| {
                    seen_records = store.len();
                })
                .run()
                .unwrap();
            pileup_profile(&mut p)
        };

        assert!(seen_records > 0, "mutator must observe the fetched records");
        assert_eq!(baseline, mutated, "a no-op mutator must not change the pileup");
    }

    // r[verify unified.readers_pileup_store_mutation+4]
    /// A mutator that soft-clips the first aligned base of every read whose
    /// CIGAR starts with `M` (shifting `pos` right by one) must be reflected in
    /// the pileup — the corrected alignments are what the engine sees.
    #[test]
    fn pileup_with_realignment_is_observed() {
        use crate::bam::CigarOp;
        use crate::bam::cigar::CigarOpType;
        use crate::bam::record_store::RecordIdx;

        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);

        let baseline = {
            let mut p = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap();
            pileup_profile(&mut p)
        };

        let mut changed = 0usize;
        let realigned = {
            let mut p = readers
                .pileup(&segment, DepthLimit::Unlimited)
                .mutate(|store, _| {
                    let mut plan: Vec<(RecordIdx, Pos0, Vec<CigarOp>)> = Vec::new();
                    for idx in store.indices() {
                        let rec = store.record(idx);
                        if rec.flags.is_unmapped() {
                            continue;
                        }
                        let Ok(cigar) = rec.cigar(store) else { continue };
                        let Some(first) = cigar.first() else { continue };
                        if first.op_type() != CigarOpType::Match || first.len() < 2 {
                            continue;
                        }
                        let mut new = Vec::with_capacity(cigar.len() + 1);
                        new.push(CigarOp::new(CigarOpType::SoftClip, 1));
                        new.push(CigarOp::new(CigarOpType::Match, first.len() - 1));
                        new.extend_from_slice(&cigar[1..]);
                        let Ok(new_pos) = Pos0::try_from(rec.pos.as_u64().saturating_add(1)) else {
                            continue;
                        };
                        plan.push((idx, new_pos, new));
                    }
                    for (idx, pos, cig) in &plan {
                        store.set_alignment(*idx, *pos, cig).unwrap();
                    }
                    changed = plan.len();
                })
                .run()
                .unwrap();
            pileup_profile(&mut p)
        };

        assert!(changed > 0, "fixture should have reads with a leading M op to realign");
        assert_ne!(baseline, realigned, "soft-clipping leading bases must change the pileup");
    }

    // r[verify unified.readers_pileup_store_mutation+4]
    /// The mutator's `RefSeq` is the engine's: the base it reads at every
    /// position of the segment is the base the column later reports there,
    /// and it covers the segment exactly — nothing before, nothing after.
    #[test]
    fn pileup_with_sees_the_engines_reference() {
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let segment = realign_test_segment(&readers);
        let (start, end) = (segment.start(), segment.end());

        let mut seen: Vec<(u64, Base)> = Vec::new();
        let mut seen_records = 0usize;
        let mut p = readers
            .pileup(&segment, DepthLimit::Unlimited)
            .mutate(|store, ref_seq| {
                seen_records = store.len();
                assert_eq!(ref_seq.start_pos(), start, "reference starts at the segment");
                assert_eq!(
                    ref_seq.len() as u64,
                    end.as_u64() - start.as_u64() + 1,
                    "reference ends at the segment"
                );
                if let Some(before) = start.as_u32().checked_sub(1).and_then(Pos0::new) {
                    assert_eq!(ref_seq.try_base_at(before), None, "nothing before the segment");
                }
                let after = Pos0::try_from(end.as_u64() + 1).unwrap();
                assert_eq!(ref_seq.try_base_at(after), None, "nothing after the segment");
                for pos in start.as_u64()..=end.as_u64() {
                    let pos = Pos0::try_from(pos).unwrap();
                    seen.push((pos.as_u64(), ref_seq.base_at(pos)));
                }
            })
            .run()
            .unwrap();

        let mut columns = 0usize;
        while let Some(col) = p.pileups() {
            let idx = usize::try_from(col.pos().as_u64() - start.as_u64()).unwrap();
            let (pos, base) = seen[idx];
            assert_eq!(pos, col.pos().as_u64());
            assert_eq!(base, col.reference_base(), "mutator and column disagree at {pos}");
            columns += 1;
        }
        drop(p);
        assert!(seen_records > 0, "mutator must observe the fetched records");
        assert!(columns > 0, "region should pile up some columns");
        assert!(
            seen.iter().any(|(_, b)| *b != Base::Unknown),
            "the fixture reference must not be all N"
        );
    }

    // r[verify pileup.extras.recover_store]
    /// Documented footgun: the guard derefs to `PileupEngine`, so
    /// `guard.reclaim_allocation()` is reachable. If the user calls it, the
    /// guard's `Drop` finds nothing to recover and the next `pileup()`
    /// re-allocates from scratch. This test pins that behavior so a
    /// future change to the guard's recovery logic doesn't silently
    /// alter the contract documented on `PileupGuard`.
    #[test]
    fn pileup_guard_reclaim_allocation_via_deref_disables_recovery() {
        use std::num::NonZeroU32;
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let opts = SegmentOptions::new(NonZeroU32::new(3_000).unwrap());
        let segments: Vec<_> = readers
            .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
            .unwrap()
            .collect();
        assert!(!segments.is_empty());

        // Misuse: drain the store through Deref before the guard drops.
        {
            let mut p = readers.pileup(&segments[0], DepthLimit::Unlimited).run().unwrap();
            while p.pileups().is_some() {}
            let drained =
                p.reclaim_allocation().expect("populated store available for the first take");
            // Hold the drained store alive past the guard's drop.
            assert!(drained.records_capacity() > 0);
            // p drops here. reclaim_allocation on the engine returns None (already
            // taken), so the slot is left as the empty Default that
            // `Readers::pileup` put there at construction time.
        }
        assert_eq!(
            readers.store.records_capacity(),
            0,
            "after the user drains via Deref, the slot is left as an empty Default — \
             the next pileup() call allocates fresh. This is the documented contract."
        );
    }

    // r[verify unified.readers_pileup+1]
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

        let err = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap_err();
        match err {
            ReaderError::SegmentHeaderMismatch { contig, expected_tid } => {
                assert_eq!(contig, bogus_contig);
                assert_eq!(expected_tid, real_tid.as_u32());
            }
            other => panic!("expected SegmentHeaderMismatch, got {other:?}"),
        }
    }

    // r[verify unified.readers_pileup+1]
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

        let err = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap_err();
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

    // r[verify unified.readers_pileup+1]
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

        let err = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap_err();
        assert!(
            matches!(err, ReaderError::SegmentHeaderMismatch { .. }),
            "expected SegmentHeaderMismatch, got {err:?}"
        );
    }

    // r[verify unified.open_with_reference]
    /// The record-iteration handle streams records from any format without a
    /// pileup engine: opened via `IndexedReader::open_with_reference`, a plain
    /// `fetch_into` into a caller-owned `RecordStore` populates records for a
    /// known-covered region.
    #[test]
    fn open_with_reference_streams_records_without_pileup() {
        let mut reader =
            IndexedReader::open_with_reference(test_bam_path(), test_fasta_path()).unwrap();
        let tid = reader.header().tid("chr19").expect("test BAM has chr19");
        let mut store = RecordStore::default();
        let n = reader
            .fetch_into(
                tid,
                Pos0::new(6_103_500).unwrap(),
                Pos0::new(6_106_500).unwrap(),
                &mut store,
            )
            .unwrap();
        assert!(n > 0, "known-covered region must yield records");
        assert_eq!(n, store.len(), "fetch_into return value must match records pushed");
    }

    // r[verify unified.readers_from_parts]
    /// `from_parts` assembles a `Readers` from already-open handles with no I/O
    /// and starts with fresh, empty record/scratch buffers — it MUST NOT carry
    /// over any buffers from the source readers.
    #[test]
    fn from_parts_starts_with_empty_buffers() {
        let alignment =
            IndexedReader::open_with_reference(test_bam_path(), test_fasta_path()).unwrap();
        let fasta = IndexedFastaReader::open(test_fasta_path()).unwrap();
        let readers = Readers::from_parts(alignment, fasta);
        assert!(readers.store.is_empty(), "assembled store must start empty");
        assert_eq!(readers.store.records_capacity(), 0, "no record buffer should be pre-allocated");
        assert_eq!(readers.fasta_buf.capacity(), 0, "no fasta buffer should be pre-allocated");
    }

    // r[verify unified.into_pileup]
    /// `IndexedReader::into_pileup` is an alternate construction path for a
    /// reference-aware `Readers`. Independent oracle: a pileup driven through
    /// the `open_with_reference(...).into_pileup(...)` path must produce
    /// exactly the same per-column `(pos, depth)` sequence as the same pileup
    /// driven through `Readers::open`.
    #[test]
    fn into_pileup_matches_readers_open() {
        use std::num::NonZeroU32;

        fn column_depths(readers: &mut Readers, seg: &Segment) -> Vec<(u64, usize)> {
            let mut out = Vec::new();
            let mut p = readers.pileup(seg, DepthLimit::Unlimited).run().unwrap();
            while let Some(col) = p.pileups() {
                out.push((col.pos().as_u64(), col.depth()));
            }
            out
        }

        let opts = SegmentOptions::new(NonZeroU32::new(3_000).unwrap());
        let region = ("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap());

        // Path A: the canonical Readers::open.
        let mut via_open = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let seg_a = via_open.segments(region, opts).unwrap().next().expect("at least one segment");
        let depths_open = column_depths(&mut via_open, &seg_a);

        // Path B: record-iteration handle upgraded with into_pileup.
        let alignment =
            IndexedReader::open_with_reference(test_bam_path(), test_fasta_path()).unwrap();
        let fasta = IndexedFastaReader::open(test_fasta_path()).unwrap();
        let mut via_into = alignment.into_pileup(fasta);
        let seg_b = via_into.segments(region, opts).unwrap().next().expect("at least one segment");
        let depths_into = column_depths(&mut via_into, &seg_b);

        assert!(!depths_open.is_empty(), "test region must produce pileup columns");
        assert_eq!(
            depths_open, depths_into,
            "into_pileup construction must yield identical pileup to Readers::open"
        );
    }
}
