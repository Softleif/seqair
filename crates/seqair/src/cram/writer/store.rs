//! `write_store`: the slices a store's records complete, encoded from the
//! store and compressed whole on worker threads, while the calling thread
//! checks records and decides where slices end.

use super::container::{self, Assembled, Codecs, SliceBuilder, SliceKind};
use super::header::CramHeader;
use super::metrics::MAX_IN_FLIGHT;
use super::record::{self, Checked, RefBases, Scratch, SliceRef};
use super::{CramWriteError, CramWriter, embed, encode_held, store_view};
use crate::bam::record_idx::RecordIdx;
use crate::bam::record_store::RecordStore;
use crate::reader::Tid;
use std::collections::BTreeMap;
use std::io::Write;
use std::sync::{Arc, Mutex, PoisonError, mpsc};

/// The bases a job's records on one reference are stored against.
enum Bases {
    /// The writer's whole reference, shared: a single-reference slice.
    #[allow(clippy::rc_buffer, reason = "the writer's buffer, reused once no slice holds it")]
    Whole(Arc<RefBases>),
    /// A copy of the window a multi-reference slice's records on it cover,
    /// so a slice across many references holds none of them whole.
    Window { start: usize, bases: Vec<u8>, iupac: bool },
}

impl Bases {
    fn slice_ref(&self) -> SliceRef<'_> {
        match self {
            Self::Whole(bases) => SliceRef::fasta(bases),
            Self::Window { start, bases, iupac } => {
                SliceRef { bases, start: *start, embedded: false, iupac: *iupac }
            }
        }
    }
}

/// A slice the store's records complete, with what encoding it needs.
struct Job {
    kind: SliceKind,
    records: Vec<(RecordIdx, Checked)>,
    references: Vec<(Tid, Bases)>,
    n_bases: u64,
}

impl Job {
    fn new(kind: SliceKind) -> Self {
        Self { kind, records: Vec::new(), references: Vec::new(), n_bases: 0 }
    }

    fn n_records(&self) -> u32 {
        u32::try_from(self.records.len()).unwrap_or(u32::MAX)
    }

    fn reference(&self, tid: Option<Tid>) -> Option<SliceRef<'_>> {
        let tid = tid?;
        self.references.iter().find(|(t, _)| *t == tid).map(|(_, b)| b.slice_ref())
    }

    /// The runs of records on one reference, as a multi-reference slice
    /// counts them.
    fn reference_runs(&self, store: &RecordStore) -> usize {
        let mut runs = 0_usize;
        let mut last = None;
        for &(idx, _) in &self.records {
            let tid = store.record(idx).and_then(|r| r.tid);
            if runs == 0 || tid != last {
                runs = runs.saturating_add(1);
                last = tid;
            }
        }
        runs
    }

    /// Replace the whole references of a multi-reference slice by the
    /// windows its records cover.
    fn window(&mut self, store: &RecordStore) {
        if self.kind != SliceKind::Multi {
            return;
        }
        for (tid, bases) in &mut self.references {
            let Bases::Whole(whole) = bases else { continue };
            let mut extent: Option<(usize, usize)> = None;
            for (idx, checked) in &self.records {
                let Some(rec) = store.record(*idx).filter(|r| r.tid == Some(*tid)) else {
                    continue;
                };
                let start = rec.pos.as_usize();
                let end = usize::try_from(checked.end).unwrap_or(usize::MAX);
                let (s, e) = extent.get_or_insert((start, end));
                *s = (*s).min(start);
                *e = (*e).max(end);
            }
            let (start, end) = extent.unwrap_or((0, 0));
            let end = end.min(whole.bases.len());
            let start = start.min(end);
            let window = whole.bases.get(start..end).unwrap_or_default().to_vec();
            let iupac = whole.iupac || record::all_iupac(&window);
            *bases = Bases::Window { start, bases: window, iupac };
        }
    }
}

/// A job numbered for writing.
struct Numbered {
    seq: u64,
    slice_index: u64,
    record_counter: u64,
    job: Job,
}

/// Encode and assemble a job's slice into `slice`.
fn run(
    numbered: Numbered,
    store: &RecordStore,
    header: &CramHeader,
    codecs: &mut Codecs,
    scratch: &mut Scratch,
    slice: &mut SliceBuilder,
) -> Result<Assembled, CramWriteError> {
    // Clean whatever an earlier job left on an error or a panic.
    slice.clear();
    let Numbered { slice_index, record_counter, job, .. } = numbered;
    let kind = job.kind;
    for (idx, checked) in &job.records {
        let view = store_view(store, *idx).ok_or(CramWriteError::NoSuchRecord { idx: *idx })?;
        match kind {
            SliceKind::Mapped { reference: false, .. } => {
                slice.hold(kind, embed::Held::new(&view, checked.clone()));
            }
            SliceKind::Unplaced => {
                record::encode(slice, kind, &view, checked, None, header, scratch)
            }
            SliceKind::Mapped { .. } | SliceKind::Multi => {
                let reference = job.reference(view.tid);
                record::encode(slice, kind, &view, checked, reference, header, scratch);
            }
        }
    }
    let length = match kind {
        SliceKind::Mapped { tid, reference: false } => {
            encode_held(slice, tid, header, scratch);
            header.target_len(tid)
        }
        SliceKind::Mapped { tid, .. } => header.target_len(tid),
        SliceKind::Multi | SliceKind::Unplaced => 0,
    };
    let reference = match (kind, job.references.first()) {
        (SliceKind::Mapped { reference: true, .. }, Some((_, Bases::Whole(bases)))) => {
            Some(bases.bases.as_slice())
        }
        _ => None,
    };
    let frame = slice.frame(slice_index, record_counter, length, reference);
    slice.assemble(&frame, codecs)
}

/// Jobs handed out, and their containers written back in order.
struct Pipeline {
    jobs: Option<mpsc::SyncSender<Numbered>>,
    done: mpsc::Receiver<(u64, Result<Assembled, CramWriteError>)>,
    ready: BTreeMap<u64, Result<Assembled, CramWriteError>>,
    next_seq: u64,
    next_out: u64,
    limit: u64,
}

impl<W: Write> CramWriter<W> {
    // r[impl cram_writer.store]
    /// Write every record of `store`, in store order, as
    /// [`write_store_record`](Self::write_store_record) would one by one,
    /// to the same bytes. With
    /// [`compression_threads`](super::CramWriterBuilder::compression_threads),
    /// the slices the store's records complete are encoded from the store
    /// and compressed whole on that many threads, without copying a record
    /// into the writer; the slice its last records leave open takes them as
    /// `write_store_record` does, for later records to complete. A store of
    /// fewer than two slices' records, or a writer without threads, takes
    /// them one by one.
    ///
    /// # Errors
    /// The first record the writer refuses ends the call with its error:
    /// the records before it are written, it and the ones after are not.
    /// An error writing a slice poisons the writer.
    pub fn write_store(&mut self, store: &RecordStore) -> Result<(), CramWriteError> {
        if self.poisoned {
            return Err(CramWriteError::Poisoned);
        }
        let mut indices = store.indices();
        let few = store.len() / 2 < usize::try_from(self.records_per_slice).unwrap_or(usize::MAX);
        if self.threads == 0 || few {
            // The persistent pool overlaps with the caller better here.
            for idx in indices {
                let view = store_view(store, idx).ok_or(CramWriteError::NoSuchRecord { idx })?;
                self.write_view(&view)?;
            }
            return Ok(());
        }
        // Records that continue the slice already open go in it one by one,
        // up to the first that ends it.
        let mut first = None;
        for idx in indices.by_ref() {
            let view = store_view(store, idx).ok_or(CramWriteError::NoSuchRecord { idx })?;
            let checked = record::check(&view, &self.header, self.last)?;
            let kind = self.kind_of(view.tid);
            if self.slice.kind().is_none()
                || self.slice.ends_before(kind, checked.read_len, self.records_per_slice)
            {
                first = Some((idx, checked));
                break;
            }
            self.add(kind, &view, checked)?;
        }
        let Some(first) = first else { return Ok(()) };
        self.flush_slice().inspect_err(|_| self.poisoned = true)?;
        // The pool's containers precede the ones encoded here.
        self.drain(true).inspect_err(|_| self.poisoned = true)?;

        let threads = self.threads;
        let header = Arc::clone(&self.header);
        let (job_tx, job_rx) = mpsc::sync_channel::<Numbered>(threads);
        let (done_tx, done_rx) = mpsc::channel();
        let job_rx = Mutex::new(job_rx);
        let limit = u64::try_from(threads.saturating_mul(2)).unwrap_or(u64::MAX);
        let mut pipeline = Pipeline {
            jobs: Some(job_tx),
            done: done_rx,
            ready: BTreeMap::new(),
            next_seq: 0,
            next_out: 0,
            limit: limit.clamp(1, MAX_IN_FLIGHT),
        };
        std::thread::scope(|scope| {
            for i in 0..threads {
                let (job_rx, done_tx, header) = (&job_rx, done_tx.clone(), &*header);
                let mut codecs = self.codecs.fork();
                let spawned = std::thread::Builder::new()
                    .name(format!("seqair-cram-store-{i}"))
                    .spawn_scoped(scope, move || {
                        let (mut scratch, mut slice) =
                            (Scratch::default(), SliceBuilder::default());
                        loop {
                            // Hold the lock only to take a job.
                            let job = job_rx.lock().unwrap_or_else(PoisonError::into_inner).recv();
                            let Ok(job) = job else { return };
                            let seq = job.seq;
                            // A panic comes back as an error: the caller
                            // writes every container in order.
                            let result =
                                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                    run(job, store, header, &mut codecs, &mut scratch, &mut slice)
                                }))
                                .unwrap_or(Err(CramWriteError::CompressionWorkerLost));
                            if done_tx.send((seq, result)).is_err() {
                                return;
                            }
                        }
                    });
                if let Err(source) = spawned {
                    // Close the jobs, or the threads spawned wait forever.
                    pipeline.jobs = None;
                    return Err(CramWriteError::ThreadSpawn { source });
                }
            }
            drop(done_tx);
            let refused = self.decide(store, first, &mut indices, &mut pipeline);
            pipeline.jobs = None;
            let written = self.write_ready(&mut pipeline, true);
            written.and(refused)
        })
    }

    /// Decide where the rest of `store`'s slices end, from `first` on,
    /// handing each complete slice to `pipeline`; the slice left open takes
    /// its records into the writer.
    fn decide(
        &mut self,
        store: &RecordStore,
        first: (RecordIdx, Checked),
        indices: &mut impl Iterator<Item = RecordIdx>,
        pipeline: &mut Pipeline,
    ) -> Result<(), CramWriteError> {
        let mut open: Option<Job> = None;
        let mut next = Some(first);
        let mut refused = Ok(());
        loop {
            let record = match next.take() {
                Some((idx, checked)) => store_view(store, idx)
                    .map(|view| (idx, checked, view))
                    .ok_or(CramWriteError::NoSuchRecord { idx }),
                None => {
                    let Some(idx) = indices.next() else { break };
                    store_view(store, idx).ok_or(CramWriteError::NoSuchRecord { idx }).and_then(
                        |view| {
                            let checked = record::check(&view, &self.header, self.last)?;
                            Ok((idx, checked, view))
                        },
                    )
                }
            };
            let (idx, checked, view) = match record {
                Ok(record) => record,
                Err(err) => {
                    refused = Err(err);
                    break;
                }
            };
            let mut kind = self.kind_of(view.tid);
            // r[impl cram_writer.slicing+2]
            if let Some(mut job) = open.take_if(|job| {
                job.kind != kind
                    || container::is_full(
                        job.n_records(),
                        job.n_bases,
                        checked.read_len,
                        self.records_per_slice,
                    )
            }) {
                self.slice_ends(job.kind, job.n_records(), job.reference_runs(store));
                job.window(store);
                self.submit(job, pipeline).inspect_err(|_| self.poisoned = true)?;
                // The slice's end may have switched multi-reference slices on or off.
                kind = self.kind_of(view.tid);
            }
            let tid = match kind {
                SliceKind::Mapped { tid, reference: true } => Some(tid),
                SliceKind::Multi => view.tid,
                SliceKind::Mapped { reference: false, .. } | SliceKind::Unplaced => None,
            };
            let job = open.get_or_insert_with(|| Job::new(kind));
            if let Some(tid) = tid
                && !job.references.iter().any(|(t, _)| *t == tid)
            {
                // A multi-reference slice is done with the references
                // before this one: keep only what its records cover.
                job.window(store);
                if let Err(err) = self.load_reference(tid) {
                    refused = Err(err);
                    break;
                }
                if let Some((_, bases)) = &self.reference {
                    job.references.push((tid, Bases::Whole(Arc::clone(bases))));
                }
            }
            job.n_bases = job.n_bases.saturating_add(u64::from(checked.read_len.unsigned_abs()));
            self.last = Some(checked.placement);
            self.warn_once(checked.drops_alignment, &view);
            job.records.push((idx, checked));
        }
        // The open slice stays open in the writer, for later records,
        // encoded against the bases the job holds, so this cannot fail.
        if let Some(job) = open {
            for (idx, checked) in &job.records {
                let Some(view) = store_view(store, *idx) else { continue };
                match job.kind {
                    SliceKind::Mapped { reference: false, .. } => {
                        self.slice.hold(job.kind, embed::Held::new(&view, checked.clone()));
                    }
                    kind => {
                        let reference =
                            job.reference(view.tid).filter(|_| kind != SliceKind::Unplaced);
                        record::encode(
                            &mut self.slice,
                            kind,
                            &view,
                            checked,
                            reference,
                            &self.header,
                            &mut self.scratch,
                        );
                    }
                }
            }
        }
        refused
    }

    /// Number `job` as the next slice and hand it out, writing whatever has
    /// come back.
    fn submit(&mut self, job: Job, pipeline: &mut Pipeline) -> Result<(), CramWriteError> {
        let numbered = Numbered {
            seq: pipeline.next_seq,
            slice_index: self.slices,
            record_counter: self.record_counter,
            job,
        };
        pipeline.next_seq = pipeline.next_seq.saturating_add(1);
        self.slices = self.slices.saturating_add(1);
        let n = u64::from(numbered.job.n_records());
        self.record_counter = self.record_counter.saturating_add(n);
        let jobs = pipeline.jobs.as_ref().ok_or(CramWriteError::CompressionWorkerLost)?;
        jobs.send(numbered).map_err(|_| CramWriteError::CompressionWorkerLost)?;
        self.write_ready(pipeline, false)
    }

    /// Write the containers that are back, in order: with `all`, every one
    /// handed out; else waiting only while too many are out. A poisoned
    /// writer writes nothing more, but waits for every container out.
    fn write_ready(&mut self, pipeline: &mut Pipeline, all: bool) -> Result<(), CramWriteError> {
        let mut failed = Ok(());
        loop {
            while let Some(result) = pipeline.ready.remove(&pipeline.next_out) {
                pipeline.next_out = pipeline.next_out.saturating_add(1);
                if self.poisoned {
                    continue;
                }
                if let Err(err) = result.and_then(|assembled| self.emit(&assembled)) {
                    self.poisoned = true;
                    failed = Err(err);
                }
            }
            let out = pipeline.next_seq.saturating_sub(pipeline.next_out);
            let wait = ((all || self.poisoned) && out > 0) || out >= pipeline.limit;
            let (seq, result) = if wait {
                if let Ok(done) = pipeline.done.recv() {
                    done
                } else {
                    self.poisoned = true;
                    return failed.and(Err(CramWriteError::CompressionWorkerLost));
                }
            } else {
                match pipeline.done.try_recv() {
                    Ok(done) => done,
                    Err(_) => return failed,
                }
            };
            pipeline.ready.insert(seq, result);
        }
    }
}
