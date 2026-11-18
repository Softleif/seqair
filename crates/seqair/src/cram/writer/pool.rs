//! Containers assembled (compressed) on worker threads while the caller
//! encodes the next slice, and handed back in the order they were submitted.

use super::container::{Assembled, Codecs, Frame, SliceBuilder};
use super::metrics::{MAX_IN_FLIGHT, Metrics};
use super::{CramVersion, CramWriteError};
use crate::io::CompressionLevel;
use std::collections::VecDeque;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;

struct Job {
    index: u64,
    slice: SliceBuilder,
    frame: Frame,
}

/// An assembled container, with its slice to reuse.
pub(super) struct Done {
    index: u64,
    pub slice: SliceBuilder,
    pub result: Result<Assembled, CramWriteError>,
}

// r[impl cram_writer.threads+2]
pub(super) struct Pool {
    jobs: Option<mpsc::SyncSender<Job>>,
    done: mpsc::Receiver<Done>,
    workers: Vec<thread::JoinHandle<()>>,
    /// Containers that came back before an earlier one; slot `i` holds
    /// container `next_out + i`.
    ready: VecDeque<Option<Done>>,
    next_in: u64,
    next_out: u64,
    max_in_flight: u64,
}

impl Pool {
    pub(super) fn new(
        threads: usize,
        level: CompressionLevel,
        version: CramVersion,
        metrics: &Arc<Metrics>,
    ) -> Result<Self, CramWriteError> {
        let threads = threads.max(1);
        let (job_tx, job_rx) = mpsc::sync_channel::<Job>(threads);
        let (done_tx, done_rx) = mpsc::channel::<Done>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        let mut workers = Vec::with_capacity(threads);
        for i in 0..threads {
            let jobs = Arc::clone(&job_rx);
            let done = done_tx.clone();
            let metrics = Arc::clone(metrics);
            let worker = thread::Builder::new()
                .name(format!("seqair-cram-{i}"))
                .spawn(move || {
                    let mut codecs = Codecs::new(level, version, metrics);
                    // Hold the lock only to take a job, not while compressing.
                    while let Ok(Job { index, slice, frame }) =
                        jobs.lock().map_err(drop).and_then(|rx| rx.recv().map_err(drop))
                    {
                        // A panic must come back as an error: the caller
                        // waits for every container in order.
                        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            slice.assemble(&frame, &mut codecs)
                        }))
                        .unwrap_or(Err(CramWriteError::CompressionWorkerLost));
                        if done.send(Done { index, slice, result }).is_err() {
                            return;
                        }
                    }
                })
                .map_err(|source| CramWriteError::ThreadSpawn { source })?;
            workers.push(worker);
        }
        Ok(Self {
            jobs: Some(job_tx),
            done: done_rx,
            workers,
            ready: VecDeque::new(),
            next_in: 0,
            next_out: 0,
            // Each worker busy and one queued: no worker waits for the caller.
            // At most `MAX_IN_FLIGHT`, which the learned codec choices rely on.
            max_in_flight: u64::try_from(threads.saturating_mul(2))
                .unwrap_or(u64::MAX)
                .min(MAX_IN_FLIGHT),
        })
    }

    pub(super) fn submit(
        &mut self,
        slice: SliceBuilder,
        frame: Frame,
    ) -> Result<(), CramWriteError> {
        let jobs = self.jobs.as_ref().ok_or(CramWriteError::CompressionWorkerLost)?;
        jobs.send(Job { index: self.next_in, slice, frame })
            .map_err(|_| CramWriteError::CompressionWorkerLost)?;
        self.next_in = self.next_in.saturating_add(1);
        Ok(())
    }

    /// The next container in submission order: waiting for it when `all`, or
    /// when too many are outstanding, else only if it is already back.
    pub(super) fn next(&mut self, all: bool) -> Result<Option<Done>, CramWriteError> {
        loop {
            if let Some(Some(_)) = self.ready.front() {
                self.next_out = self.next_out.saturating_add(1);
                return Ok(self.ready.pop_front().flatten());
            }
            let outstanding = self.next_in.saturating_sub(self.next_out);
            let wait = if all { outstanding > 0 } else { outstanding >= self.max_in_flight };
            let done = if wait {
                self.done.recv().map_err(|_| CramWriteError::CompressionWorkerLost)?
            } else {
                match self.done.try_recv() {
                    Ok(done) => done,
                    Err(mpsc::TryRecvError::Empty) => return Ok(None),
                    Err(mpsc::TryRecvError::Disconnected) => {
                        return Err(CramWriteError::CompressionWorkerLost);
                    }
                }
            };
            let slot = usize::try_from(done.index.saturating_sub(self.next_out))
                .map_err(|_| CramWriteError::CompressionWorkerLost)?;
            if self.ready.len() <= slot {
                self.ready.resize_with(slot.saturating_add(1), || None);
            }
            if let Some(entry) = self.ready.get_mut(slot) {
                *entry = Some(done);
            }
        }
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        // Closing the job channel ends every worker's loop.
        self.jobs = None;
        for worker in self.workers.drain(..) {
            if worker.join().is_err() {
                tracing::warn!("a CRAM compression worker panicked");
            }
        }
    }
}
