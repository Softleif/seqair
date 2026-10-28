//! [`StreamPileup`]: pileup columns over a whole coordinate-sorted BAM or SAM
//! stream, read front to back without an index. See
//! `docs/spec/4-pileup-stream.md`.

use super::{AlignmentStream, AlignmentStreamError, DepthCap, DepthLimit, Tid};
use crate::bam::{
    BamHeader, DecodeError, RawBamRecord, RecordIdx,
    pileup::{PileupColumn, PileupEngine, PileupScratch, RefSeq},
    record_store::{CustomizeRecordStore, Prepared, RecordStore},
};
use crate::fasta::{FastaError, IndexedFastaReader};
use core::range::RangeInclusive;
use seqair_types::{Base, Offset, Pos0, SmolStr};
use std::{num::NonZeroU32, rc::Rc};

/// Why a [`StreamPileup`] stopped.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum StreamPileupError {
    #[error("reading the alignment stream failed")]
    Stream {
        #[from]
        source: AlignmentStreamError,
    },

    /// A mapped record precedes the last mapped record before it
    /// (`r[pileup_stream.sorted+1]`). A position of `None` is a mapped record
    /// on a reference with no position, which htslib orders as `-1`.
    #[error(
        "record {record} (reference {tid}, position {pos:?}) precedes the record before it \
         (reference {previous_tid}, position {previous_pos:?}): the input is not coordinate-sorted"
    )]
    Unsorted {
        record: u64,
        tid: Tid,
        pos: Option<Pos0>,
        previous_tid: Tid,
        previous_pos: Option<Pos0>,
    },

    /// The record store refused a record.
    #[error("record {record} could not be stored")]
    Record { record: u64, source: DecodeError },

    /// More reads start at one position than a store can hold (its `u32`
    /// offsets would overflow), so no window boundary can split them
    /// (`r[pileup_stream.budget]`). [`DepthLimit::PerColumn`] bounds it.
    #[error("too many reads start at reference {tid}, position {pos} for one record store")]
    PileTooDeep { tid: Tid, pos: Pos0 },

    #[error("reading the reference {contig}:{start}-{end} (0-based, inclusive) failed")]
    FastaFetch { contig: SmolStr, start: u64, end: u64, source: FastaError },

    /// An earlier call returned an error; the pileup cannot go on.
    #[error("the streaming pileup already failed; an earlier call returned the error")]
    Poisoned,
}

// r[impl pileup_stream.window_len]
/// The default window length (`r[pileup_stream.window_len]`).
pub const DEFAULT_WINDOW_LEN: NonZeroU32 = NonZeroU32::new(64 * 1024).expect("non-zero");

// r[impl pileup_stream.budget]
/// When a window closes before its length is reached: once its store holds
/// more than `records` records, or its longest slab (bytes, or CIGAR
/// operations) is longer than `slab_len` (`r[pileup_stream.budget]`).
///
/// A soft limit: a window never closes inside one position, so a pile of
/// reads starting at one position grows past it — up to the store's `u32`
/// limits, where the pileup fails with [`StreamPileupError::PileTooDeep`].
/// The default is half those limits, so the reads carried into the next
/// window always fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowBudget {
    pub records: NonZeroU32,
    pub slab_len: NonZeroU32,
}

impl Default for WindowBudget {
    fn default() -> Self {
        let half = NonZeroU32::new(1 << 31).expect("non-zero");
        Self { records: half, slab_len: half }
    }
}

/// A column and the reference it is on, as [`StreamPileup::pileups`] yields it.
pub type StreamColumn<'a, U> = (Tid, PileupColumn<'a, U>);

/// The window being filled: positions `start..=end` of `tid`.
#[derive(Debug, Clone, Copy)]
struct Window {
    tid: Tid,
    start: Pos0,
    end: Pos0,
}

/// The window whose columns are being yielded.
struct Yielding<U> {
    tid: Tid,
    engine: PileupEngine<U>,
}

// r[impl pileup_stream.columns]
// r[impl pileup_stream.memory]
/// Pileup columns over a whole coordinate-sorted alignment stream — `stdin`,
/// a pipe — with no index, as `samtools mpileup` reads unindexed input.
///
/// Yields, contig by contig, the columns [`Readers::pileup`] yields for a
/// segment spanning each whole contig, by running the same
/// [`PileupEngine`] over a sliding window of positions. Columns depend on
/// neither the window length nor the [`WindowBudget`]. Three things differ
/// from an indexed pileup, by design:
///
/// - **Store indices are per window.** `record_idx`, `mate_idx` and
///   `fragment_idx` index the current window's store. Mates that overlap
///   each other are always in the same window and always linked, so
///   [`PileupColumn::mate_of`] and `in_mate_overlap` are exactly the indexed
///   pileup's. Mates that do not overlap are linked only when they happen to
///   share a window, so for them `mate_idx` may be `None` (and
///   `fragment_idx` the read's own index) where a whole-contig pileup links
///   them; do not rely on either for such mates. `samtools mpileup` has the
///   same limit: it pairs a read with its mate only while both are in its
///   buffer.
/// - **The depth cap is a whole-contig one.** With
///   [`DepthLimit::PerColumn`] the [`DepthCap`] sees each contig's records
///   once, in order, and keeps exactly the reads an indexed pileup over one
///   segment per contig keeps. A *tiled* indexed pileup restarts the cap at
///   every segment, so near tile boundaries it can keep different reads.
/// - **The sort check is htslib's.** Mapped records on a reference that the
///   customize value keeps must be in `(tid, pos)` order, else
///   [`StreamPileupError::Unsorted`]; unplaced and unmapped records, and the
///   ones the customize value rejects, are not checked, as `bam_plp_push`
///   does not check them.
///
/// ```no_run
/// use seqair::reader::{AlignmentStream, DepthLimit, StreamPileup};
///
/// let stream = AlignmentStream::from_reader(std::io::stdin())?;
/// let mut pileup = StreamPileup::new(stream, (), DepthLimit::Unlimited);
/// while let Some((tid, column)) = pileup.pileups()? {
///     println!("{tid}\t{}\t{}", column.pos(), column.depth());
/// }
/// # Ok::<(), seqair::reader::StreamPileupError>(())
/// ```
///
/// [`Readers::pileup`]: crate::reader::Readers::pileup
pub struct StreamPileup<E: CustomizeRecordStore = ()> {
    stream: AlignmentStream,
    customize: DepthCap<E>,
    fasta: Option<IndexedFastaReader>,
    fasta_buf: Vec<u8>,
    window_len: NonZeroU32,
    budget: WindowBudget,
    record: RawBamRecord,
    /// Records read so far, for error messages (1-based numbering).
    records_read: u64,
    /// The last checked record's `(tid, pos)`, for the sort check; `None`
    /// sorts first, as htslib's `-1`.
    last: Option<(Tid, Option<Pos0>)>,
    filling: Option<Window>,
    /// The store of the window being filled.
    store: RecordStore<E::Extra>,
    /// A cleared store kept for the next window, so two stores' allocations
    /// are reused from window to window.
    spare: RecordStore<E::Extra>,
    yielding: Option<Yielding<E::Extra>>,
    scratch: PileupScratch,
    input_ended: bool,
    failed: bool,
}

impl<E: CustomizeRecordStore> std::fmt::Debug for StreamPileup<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StreamPileup")
            .field("stream", &self.stream)
            .field("window_len", &self.window_len)
            .field("records_read", &self.records_read)
            .field("filling", &self.filling)
            .finish_non_exhaustive()
    }
}

impl<E> StreamPileup<E>
where
    E: CustomizeRecordStore,
    E::Extra: Clone,
{
    // r[impl pileup_stream.depth_limit+1]
    /// A pileup over `stream`, pushing every placed record through
    /// `customize` once, and capping reads per column at load time as
    /// [`Readers::pileup`](crate::reader::Readers::pileup) does for `depth`.
    pub fn new(stream: AlignmentStream, customize: E, depth: DepthLimit) -> Self {
        Self {
            stream,
            customize: DepthCap::with_cap(customize, depth.per_column()),
            fasta: None,
            fasta_buf: Vec::new(),
            window_len: DEFAULT_WINDOW_LEN,
            budget: WindowBudget::default(),
            record: RawBamRecord::default(),
            records_read: 0,
            last: None,
            filling: None,
            store: RecordStore::new(),
            spare: RecordStore::new(),
            yielding: None,
            scratch: PileupScratch::default(),
            input_ended: false,
            failed: false,
        }
    }

    // r[impl pileup_stream.reference]
    /// Report each column's reference base from `fasta`, read by the
    /// contig's name. Without it every reference base is [`Base::Unknown`].
    #[must_use]
    pub fn with_reference(mut self, fasta: IndexedFastaReader) -> Self {
        self.fasta = Some(fasta);
        self
    }

    /// Set the window length (`r[pileup_stream.window_len]`). Columns do not
    /// depend on it; memory and per-window work do.
    #[must_use]
    pub fn with_window_len(mut self, window_len: NonZeroU32) -> Self {
        self.window_len = window_len;
        self
    }

    /// Set the window budget (`r[pileup_stream.budget]`). Columns do not
    /// depend on it.
    #[must_use]
    pub fn with_window_budget(mut self, budget: WindowBudget) -> Self {
        self.budget = budget;
        self
    }

    /// The stream's header.
    pub fn header(&self) -> &BamHeader {
        self.stream.header()
    }

    /// The customize value, for inspecting state it accumulated.
    pub fn customize(&self) -> &E {
        self.customize.inner()
    }

    // r[impl pileup_stream.lending]
    // r[impl pileup_stream.errors]
    /// The next column and the reference it is on, or `None` once the input
    /// is exhausted. The column borrows the pileup until the next call.
    ///
    /// # Errors
    /// Reading the stream, an unsorted record, a record the store refuses,
    /// a pile too deep for one store, or the reference. Any error ends the pileup: later calls return
    /// [`StreamPileupError::Poisoned`].
    pub fn pileups(&mut self) -> Result<Option<StreamColumn<'_, E::Extra>>, StreamPileupError> {
        if self.failed {
            return Err(StreamPileupError::Poisoned);
        }
        let step = loop {
            if let Some(yielding) = self.yielding.as_mut()
                && let Some(step) = yielding.engine.advance()
            {
                break step;
            }
            match self.next_window() {
                Ok(true) => {}
                Ok(false) => return Ok(None),
                Err(e) => {
                    self.failed = true;
                    return Err(e);
                }
            }
        };
        // `break step` above came from this very engine.
        let Some(yielding) = self.yielding.as_ref() else { return Ok(None) };
        Ok(Some((yielding.tid, yielding.engine.column(step))))
    }

    // r[impl pileup_stream.window+2]
    // r[impl pileup_stream.records+1]
    /// Recycle the exhausted engine, then read records until a window
    /// closes and start an engine over it. `false` once nothing is left.
    fn next_window(&mut self) -> Result<bool, StreamPileupError> {
        if let Some(Yielding { mut engine, .. }) = self.yielding.take() {
            self.scratch = engine.take_scratch();
            if let Some(store) = engine.reclaim_allocation() {
                self.spare = store;
            }
        }
        loop {
            if self.input_ended {
                let Some(window) = self.filling.take() else {
                    // Only placed-unmapped records, if anything: no column.
                    self.store.clear();
                    return Ok(false);
                };
                let store = std::mem::take(&mut self.store);
                let last = self.final_last(window.tid, &store);
                self.start_engine(window.tid, store, window.start, last)?;
                return Ok(true);
            }
            if !self.stream.read_record(&mut self.record)? {
                self.input_ended = true;
                continue;
            }
            self.records_read = self.records_read.saturating_add(1);
            let Some(tid) = self.record.tid() else { continue };
            let mapped = !self.record.flags().is_unmapped();
            let Some(pos) = self.record.pos() else {
                // Checked at -1 as htslib does, then skipped: a store holds
                // placed records only.
                if mapped {
                    self.check_order(tid, None)?;
                }
                continue;
            };
            // Unmapped and rejected records go into the filling store (or
            // nowhere) without moving a window boundary.
            if !self.push_record(tid, pos)? || !mapped {
                continue;
            }
            self.check_order(tid, Some(pos))?;

            // The record is the filling store's last.
            let window_end = pos
                .checked_add_offset(Offset::new(i64::from(self.window_len.get().saturating_sub(1))))
                .unwrap_or(Pos0::MAX);
            let Some(window) = self.filling else {
                self.filling = Some(Window { tid, start: pos, end: window_end });
                continue;
            };
            let same_contig = window.tid == tid;
            // r[impl pileup_stream.budget]
            // Only a record past the window's start can end it early: there is
            // no boundary inside one position.
            let end = if same_contig && pos <= window.end {
                if pos == window.start || !self.over_budget() {
                    continue;
                }
                pos.checked_sub_offset(Offset::new(1)).unwrap_or(Pos0::ZERO)
            } else {
                window.end
            };

            // The record closes `window`: every record that can reach a
            // position up to `end` (or any of a finished contig) is in
            // `self.store`. It moves to the next store, after the carried
            // reads.
            let mut next = std::mem::take(&mut self.spare);
            next.clear();
            let closer = self.store.len().saturating_sub(1);
            let start = if same_contig {
                self.carry(&mut next, end, closer, tid)?;
                if next.is_empty() {
                    pos
                } else {
                    end.checked_add_offset(Offset::new(1)).unwrap_or(Pos0::MAX)
                }
            } else {
                pos
            };
            self.copy_into(&mut next, closer, tid)?;
            self.store.pop_last();
            let last = if same_contig {
                end.min(self.contig_last(window.tid))
            } else {
                self.final_last(window.tid, &self.store)
            };
            let full = std::mem::replace(&mut self.store, next);
            self.filling = Some(Window { tid, start, end: window_end });
            self.start_engine(window.tid, full, window.start, last)?;
            return Ok(true);
        }
    }

    // r[impl pileup_stream.sorted+1]
    /// Fail when `(tid, pos)` precedes the last checked record, as
    /// `bam_plp_push` does; `None` is htslib's position `-1`.
    fn check_order(&mut self, tid: Tid, pos: Option<Pos0>) -> Result<(), StreamPileupError> {
        if let Some((previous_tid, previous_pos)) = self.last
            && (tid, pos) < (previous_tid, previous_pos)
        {
            return Err(StreamPileupError::Unsorted {
                record: self.records_read,
                tid,
                pos,
                previous_tid,
                previous_pos,
            });
        }
        self.last = Some((tid, pos));
        Ok(())
    }

    /// Whether the filling store has outgrown the budget.
    fn over_budget(&self) -> bool {
        self.store.len() > self.budget.records.get() as usize
            || self.store.longest_slab() > self.budget.slab_len.get() as usize
    }

    /// Push the record just read, at `tid:pos`, into the filling store
    /// through the customize value. `true` when it was kept.
    fn push_record(&mut self, tid: Tid, pos: Pos0) -> Result<bool, StreamPileupError> {
        let Self { store, record, stream, customize, records_read, .. } = self;
        match store.push_raw(record.as_payload(), stream.header(), customize) {
            Ok(kept) => Ok(kept.is_some()),
            Err(DecodeError::SlabOverflow) => Err(StreamPileupError::PileTooDeep { tid, pos }),
            Err(source) => Err(StreamPileupError::Record { record: *records_read, source }),
        }
    }

    // r[impl pileup_stream.carry]
    /// Copy the records of the filling store before `closer` that extend
    /// past `end` into `next`, in store order.
    fn carry(
        &self,
        next: &mut RecordStore<E::Extra>,
        end: Pos0,
        closer: usize,
        tid: Tid,
    ) -> Result<(), StreamPileupError> {
        for i in 0..closer {
            let crosses = RecordIdx::from_usize(i)
                .and_then(|idx| self.store.record(idx))
                .is_some_and(|rec| rec.end_pos > end);
            if crosses {
                self.copy_into(next, i, tid)?;
            }
        }
        Ok(())
    }

    /// Copy record `i` of the filling store into `next`.
    fn copy_into(
        &self,
        next: &mut RecordStore<E::Extra>,
        i: usize,
        tid: Tid,
    ) -> Result<(), StreamPileupError> {
        let Some(idx) = RecordIdx::from_usize(i) else { return Ok(()) };
        let pos = self.store.record(idx).map_or(Pos0::ZERO, |rec| rec.pos);
        next.copy_record_from(&self.store, idx).map_err(|source| match source {
            DecodeError::SlabOverflow => StreamPileupError::PileTooDeep { tid, pos },
            source => StreamPileupError::Record { record: self.records_read, source },
        })?;
        Ok(())
    }

    /// Where the last window of `tid` stops: the furthest end of its records,
    /// within the contig.
    fn final_last(&self, tid: Tid, store: &RecordStore<E::Extra>) -> Pos0 {
        let furthest = store
            .records()
            .filter(|rec| !rec.flags.is_unmapped())
            .map(|rec| rec.end_pos)
            .max()
            .unwrap_or(Pos0::ZERO);
        furthest.min(self.contig_last(tid))
    }

    /// The last position of `tid` per the header.
    fn contig_last(&self, tid: Tid) -> Pos0 {
        self.stream
            .header()
            .target_len(tid)
            .and_then(|len| len.checked_sub(1))
            .and_then(|last| u32::try_from(last).ok())
            .and_then(Pos0::new)
            .unwrap_or(Pos0::MAX)
    }

    // r[impl pileup_stream.mates+1]
    /// Link `store`'s mates and start an engine yielding `start..=last`.
    fn start_engine(
        &mut self,
        tid: Tid,
        store: RecordStore<E::Extra>,
        start: Pos0,
        last: Pos0,
    ) -> Result<(), StreamPileupError> {
        let reference = self.reference(tid, start, last, store.is_empty())?;
        let Prepared { input, stats } = store.prepare_for_pileup();
        if !stats.is_clean() {
            tracing::debug!(
                tid = tid.as_u32(),
                start = start.as_u64(),
                pairs = stats.pairs,
                ambiguous_qnames = stats.ambiguous_qnames,
                "qnames are not unique here; those reads are not mate-deduplicated"
            );
        }
        let span = RangeInclusive { start, last };
        let mut engine = PileupEngine::with_scratch(input, span, std::mem::take(&mut self.scratch));
        if let Some(reference) = reference {
            engine.set_reference_seq(reference);
        }
        self.yielding = Some(Yielding { tid, engine });
        Ok(())
    }

    // r[impl pileup_stream.reference]
    /// The reference over `start..=last`, when there is a FASTA and the
    /// window can yield a column at all.
    fn reference(
        &mut self,
        tid: Tid,
        start: Pos0,
        last: Pos0,
        empty: bool,
    ) -> Result<Option<RefSeq>, StreamPileupError> {
        let Some(fasta) = self.fasta.as_mut() else { return Ok(None) };
        if empty || start > last {
            return Ok(None);
        }
        let contig = self.stream.header().target_name(tid).unwrap_or_default();
        fasta.fetch_seq_into(contig, RangeInclusive { start, last }, &mut self.fasta_buf).map_err(
            |source| StreamPileupError::FastaFetch {
                contig: SmolStr::new(contig),
                start: start.as_u64(),
                end: last.as_u64(),
                source,
            },
        )?;
        let bases: &[Base] = Base::convert_ascii_in_place_as_slice(&mut self.fasta_buf);
        Ok(Some(RefSeq::new(Rc::from(bases), start)))
    }
}
