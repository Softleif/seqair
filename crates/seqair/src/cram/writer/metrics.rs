//! Which codecs a block is tried with: every one on a trial slice, and on the
//! others the two that did best on earlier trial slices, as htslib's
//! `cram_metrics` learns them.

use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, VecDeque};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};

// r[impl cram_writer.compression.learning]

/// A trial slice every this many slices, the first included.
const TRIAL_SPAN: u64 = 70;

/// The runner-up is tried too while its score is within this factor of the
/// best's: a close second may win the next block, a distant one rarely does
/// (on a 3 M-read chr20 BAM, 5% kept the size to 0.02% and cut the
/// compression time by a quarter).
const RUNNER_UP_MARGIN: f64 = 1.05;

/// The most containers the pool may hold between submitting and writing.
/// Below [`TRIAL_SPAN`], so a slice still to plan is never more than one
/// trial behind the newest, and [`State::snapshots`] stays short.
pub(super) const MAX_IN_FLIGHT: u64 = 64;

/// A way to store a block. `Nx16` names an entry of
/// [`super::rans_nx16::normal_profile`], whose state count follows the
/// block's length.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Method {
    Raw,
    Gzip,
    Rans4x8Order0,
    Rans4x8Order1,
    Nx16(u8),
    Tok3,
    /// An entry of [`super::arith::PROFILE`].
    Arith(u8),
    /// An entry of [`super::fqzcomp::Strategy::ALL`].
    Fqz(u8),
}

/// A total order on methods, for ties: raw, then the rANS codecs, then
/// gzip and the name tokeniser, then the arithmetic coder and fqzcomp, each
/// slower to decode than the one before.
pub(super) fn rank(method: Method) -> u16 {
    match method {
        Method::Raw => 0,
        Method::Rans4x8Order0 => 1,
        Method::Rans4x8Order1 => 2,
        Method::Nx16(i) => 3_u16.saturating_add(u16::from(i)),
        Method::Gzip => 300,
        Method::Tok3 => 301,
        Method::Arith(i) => 400_u16.saturating_add(u16::from(i)),
        Method::Fqz(i) => 500_u16.saturating_add(u16::from(i)),
    }
}

/// One block's sizes on a trial slice: every method tried, and its length.
#[derive(Debug, Clone)]
pub(super) struct Trial {
    pub content_id: i32,
    pub len: usize,
    pub sizes: Vec<(Method, usize)>,
}

/// What a slice tries.
#[derive(Debug, Clone)]
pub(super) enum Plan {
    /// Every method, recording the sizes.
    Trial,
    /// The learned choices, by content id.
    Learned(Arc<Choices>),
}

pub(super) type Choices = FxHashMap<i32, LearnedChoice>;

/// The methods a non-trial slice tries for one content id.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct LearnedChoice {
    best: Method,
    runner_up: Option<Method>,
    /// The trials' (weighted) average length of the block.
    avg_len: f64,
}

impl LearnedChoice {
    /// The methods to try for a block of `len` bytes, or `None` when the
    /// length is far from what the trials saw: then every method.
    pub(super) fn methods(self, len: usize) -> Option<[Option<Method>; 2]> {
        let Self { best, runner_up, avg_len } = self;
        #[allow(clippy::cast_precision_loss, reason = "a block length")]
        let len = len as f64;
        let drifted =
            (len + 1000.0) / 4.0 > avg_len + 1000.0 || len + 1000.0 < (avg_len + 1000.0) / 4.0;
        (!drifted).then_some([Some(best), runner_up])
    }
}

/// The learned scores, shared by every thread that assembles containers.
#[derive(Debug, Default)]
pub(super) struct Metrics {
    exhaustive: bool,
    state: Mutex<State>,
    recorded: Condvar,
}

#[derive(Debug, Default)]
struct State {
    /// Trial slices recorded but not folded into `scores` yet, by index.
    pending: BTreeMap<u64, Vec<Trial>>,
    scores: FxHashMap<i32, Score>,
    /// The index of the last trial folded into `scores`.
    folded: Option<u64>,
    /// The choices after each of the last few folded trials, oldest first.
    snapshots: VecDeque<(u64, Arc<Choices>)>,
}

/// One content id's running scores: each method's size, and the block's
/// length, summed over the trials with the older ones halved per trial.
#[derive(Debug, Default)]
struct Score {
    sizes: FxHashMap<Method, f64>,
    len: f64,
    weight: f64,
}

impl Metrics {
    pub(super) fn new(exhaustive: bool) -> Self {
        Self { exhaustive, ..Self::default() }
    }

    pub(super) fn is_trial(&self, slice: u64) -> bool {
        self.exhaustive || slice.is_multiple_of(TRIAL_SPAN)
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        // Folding panics on nothing it could leave half done.
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// What slice `slice` tries. Its choices come from the trials up to the
    /// one a span before its own span's (the first trial, for the first two
    /// spans), so a worker waits only at the start of the file, and gets the
    /// choices the single-threaded writer makes.
    pub(super) fn plan(&self, slice: u64) -> Plan {
        if self.is_trial(slice) {
            return Plan::Trial;
        }
        let source = (slice / TRIAL_SPAN).saturating_sub(1).saturating_mul(TRIAL_SPAN);
        let mut state = self.lock();
        loop {
            if let Some((_, choices)) = state.snapshots.iter().find(|(at, _)| *at == source) {
                return Plan::Learned(Arc::clone(choices));
            }
            if !state.fold_next() {
                state = self.recorded.wait(state).unwrap_or_else(PoisonError::into_inner);
            }
        }
    }

    /// Record a trial slice's sizes, waking the slices that wait for it.
    /// Every trial slice must be recorded, however its assembly ends.
    pub(super) fn record(&self, slice: u64, trials: Vec<Trial>) {
        if self.exhaustive {
            return;
        }
        self.lock().pending.entry(slice).or_insert(trials);
        self.recorded.notify_all();
    }
}

impl State {
    /// Fold the next trial into the scores if it has been recorded.
    fn fold_next(&mut self) -> bool {
        let next = self.folded.map_or(0, |f| f.saturating_add(TRIAL_SPAN));
        let Some(trials) = self.pending.remove(&next) else { return false };
        if self.folded.is_some() {
            for score in self.scores.values_mut() {
                score.sizes.values_mut().for_each(|s| *s /= 2.0);
                score.len /= 2.0;
                score.weight /= 2.0;
            }
        }
        for trial in &trials {
            self.scores.entry(trial.content_id).or_default().add(trial);
        }
        let choices: Choices =
            self.scores.iter().filter_map(|(&id, score)| Some((id, score.choice()?))).collect();
        self.snapshots.push_back((next, Arc::new(choices)));
        if self.snapshots.len() > 3 {
            self.snapshots.pop_front();
        }
        self.folded = Some(next);
        true
    }
}

impl Score {
    fn add(&mut self, trial: &Trial) {
        // A method this trial lacks scores the raw size, as does one for the
        // trials before it that lacked it.
        let history = self.len;
        for &(method, _) in &trial.sizes {
            self.sizes.entry(method).or_insert(history);
        }
        for (method, total) in &mut self.sizes {
            let size = trial.sizes.iter().find(|(m, _)| m == method).map_or(trial.len, |&(_, s)| s);
            #[allow(clippy::cast_precision_loss, reason = "block sizes")]
            let size = size as f64;
            *total += size;
        }
        #[allow(clippy::cast_precision_loss, reason = "block sizes")]
        let len = trial.len as f64;
        self.len += len;
        self.weight += 1.0;
    }

    fn choice(&self) -> Option<LearnedChoice> {
        let mut ranked: Vec<(Method, f64)> = self.sizes.iter().map(|(&m, &s)| (m, s)).collect();
        ranked.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| rank(a.0).cmp(&rank(b.0))));
        let mut ranked = ranked.into_iter();
        let (best, best_score) = ranked.next()?;
        let runner_up = ranked
            .next()
            .filter(|&(_, score)| score <= best_score * RUNNER_UP_MARGIN)
            .map(|(m, _)| m);
        let avg_len = if self.weight > 0.0 { self.len / self.weight } else { 0.0 };
        Some(LearnedChoice { best, runner_up, avg_len })
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::cast_precision_loss,
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects,
        reason = "test code"
    )]
    use super::*;
    use hegel::prelude::*;

    const METHODS: [Method; 4] = [Method::Raw, Method::Gzip, Method::Nx16(0), Method::Nx16(1)];

    fn trial(content_id: i32, len: usize, sizes: &[(Method, usize)]) -> Trial {
        Trial { content_id, len, sizes: sizes.to_vec() }
    }

    fn learned(metrics: &Metrics, slice: u64) -> Arc<Choices> {
        match metrics.plan(slice) {
            Plan::Learned(choices) => choices,
            Plan::Trial => panic!("slice {slice} is a trial"),
        }
    }

    // r[verify cram_writer.compression.learning]
    /// The scores as a closed form: trial `i` of `n` weighs `2^-(n-1-i)`,
    /// and a method missing from it counts the block's length.
    #[hegel::test]
    fn the_scores_are_halving_weighted_sums(tc: TestCase) {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
        let mut trials = Vec::new();
        for _ in 0..n {
            let len = tc.draw(gs::integers::<usize>().min_value(32).max_value(5000));
            let mut sizes = Vec::new();
            for &m in &METHODS {
                if m == Method::Raw {
                    sizes.push((m, len));
                } else if tc.draw(gs::booleans()) {
                    sizes.push((m, tc.draw(gs::integers::<usize>().max_value(6000))));
                }
            }
            trials.push(trial(1, len, &sizes));
        }
        let metrics = Metrics::new(false);
        for (i, t) in trials.iter().enumerate().rev() {
            metrics.record(i as u64 * TRIAL_SPAN, vec![t.clone()]);
        }
        // The last trial, `n - 1`, serves the span after its own.
        let choice = learned(&metrics, n as u64 * TRIAL_SPAN + 1)[&1];

        let weight = |i: usize| 0.5_f64.powi(i32::try_from(n - 1 - i).unwrap());
        let mut oracle: Vec<(Method, f64)> = METHODS
            .iter()
            .filter(|&&m| trials.iter().any(|t| t.sizes.iter().any(|s| s.0 == m)))
            .map(|&m| {
                let total = trials.iter().enumerate().map(|(i, t)| {
                    let size = t.sizes.iter().find(|s| s.0 == m).map_or(t.len, |s| s.1);
                    weight(i) * size as f64
                });
                (m, total.sum())
            })
            .collect();
        oracle.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| rank(a.0).cmp(&rank(b.0))));
        // Sums of halvings of integers this small are exact in f64.
        assert_eq!(choice.best, oracle[0].0, "{oracle:?}");
        let runner_up = oracle.get(1).filter(|o| o.1 <= oracle[0].1 * 1.05).map(|o| o.0);
        assert_eq!(choice.runner_up, runner_up, "{oracle:?}");
        let avg: f64 =
            trials.iter().enumerate().map(|(i, t)| weight(i) * t.len as f64).sum::<f64>()
                / (0..n).map(weight).sum::<f64>();
        assert!((choice.avg_len - avg).abs() < 1e-6, "{} vs {avg}", choice.avg_len);
    }

    #[test]
    fn a_far_longer_or_shorter_block_and_an_unseen_one_try_everything() {
        let metrics = Metrics::new(false);
        let sizes = [(Method::Raw, 1000), (Method::Gzip, 310), (Method::Nx16(0), 300)];
        metrics.record(0, vec![trial(1, 1000, &sizes)]);
        let choices = learned(&metrics, 1);
        let c = choices[&1];
        assert_eq!(c.methods(1000), Some([Some(Method::Nx16(0)), Some(Method::Gzip)]));
        assert_eq!(c.methods(10), Some([Some(Method::Nx16(0)), Some(Method::Gzip)]));
        assert_eq!(c.methods(100_000), None);
        assert!(!choices.contains_key(&2));
    }

    /// The first two spans learn from the first trial, the third from the
    /// first two.
    #[test]
    fn a_slice_learns_from_the_trial_a_span_back() {
        let metrics = Metrics::new(false);
        let at = |size| vec![trial(1, 1000, &[(Method::Raw, 1000), (Method::Gzip, size)])];
        // Recorded out of order, as threads may.
        metrics.record(TRIAL_SPAN, at(100));
        metrics.record(0, at(900));
        let first = learned(&metrics, 1);
        assert!(Arc::ptr_eq(&first, &learned(&metrics, TRIAL_SPAN + 5)));
        let both = learned(&metrics, 2 * TRIAL_SPAN + 1);
        assert!(!Arc::ptr_eq(&first, &both));
        // A plan for an older slice still gets its own snapshot.
        assert!(Arc::ptr_eq(&first, &learned(&metrics, 3)));
        assert_eq!(first[&1].best, Method::Gzip);
        // Raw (1500) is far behind gzip (550): no runner-up.
        assert_eq!(both[&1].runner_up, None);
        // Raw (1000) is within 5% of gzip (900)? No: 1000 > 945.
        assert_eq!(first[&1].runner_up, None);
    }

    #[test]
    fn a_slice_waits_for_its_trial() {
        let metrics = Arc::new(Metrics::new(false));
        let waiter = {
            let metrics = Arc::clone(&metrics);
            std::thread::spawn(move || learned(&metrics, 5).len())
        };
        // Give the waiter time to block; either order must give trial 0.
        std::thread::sleep(std::time::Duration::from_millis(20));
        metrics.record(0, vec![trial(1, 100, &[(Method::Raw, 100)])]);
        assert_eq!(waiter.join().unwrap(), 1);
    }
}
