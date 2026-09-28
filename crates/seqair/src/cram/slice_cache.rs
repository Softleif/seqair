//! Decoded CRAM slices shared by a reader and its forks.
//!
//! Neighbouring region queries overlap the same slices at their edges, and a
//! slice is only useful once decoded whole, so without sharing every query
//! decodes each slice it touches again. The cache keeps decoded slices under
//! a byte budget, and makes threads that want the same missing slice at the
//! same time wait for one decode instead of each running their own.
// r[impl cram.slice_cache.shared]
// r[impl cram.slice_cache.budget+3]

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, OnceLock};

use super::reader::CramError;
use super::slice::DecodedSlice;

/// Where a decoded slice comes from: the slice's byte position in the file
/// and the reference it was decoded for (a multi-reference slice decodes
/// differently per reference).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct SliceKey {
    pub(crate) container_offset: u64,
    pub(crate) slice_offset: u64,
    pub(crate) tid: u32,
}

/// A slice being decoded or decoded: the first caller initialises it, and
/// every caller that finds it meanwhile blocks until it is set. `None`
/// means the decode failed.
type Cell = Arc<OnceLock<Option<Arc<DecodedSlice>>>>;

struct Entry {
    cell: Cell,
    /// When a fetch last asked for it, in [`Inner::tick`]s.
    last_use: u64,
    /// Its heap size once decoded and counted; 0 while being decoded.
    bytes: usize,
}

#[derive(Default)]
struct Inner {
    entries: HashMap<SliceKey, Entry>,
    tick: u64,
    /// Sum of the counted entries' `bytes`.
    used: usize,
    budget: SliceCacheBudget,
}

/// How much a reader's decoded-slice cache may keep
/// ([`IndexedCramReader::set_slice_cache_budget`](super::reader::IndexedCramReader::set_slice_cache_budget)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SliceCacheBudget {
    /// Room for the slices the open handles on the file (the reader and its
    /// live forks) are working through, at the largest slice's size: one
    /// for a lone handle, the one its last fetch ended in and a sweep's next
    /// fetch starts in; for several, the slices their fetches in flight
    /// span beyond the first — `handles` times the mean a fetch wants less
    /// one, rounded up — plus one.
    #[default]
    Auto,
    /// At most this many heap bytes; 0 turns caching off.
    Bytes(usize),
}

/// What a reader's decoded-slice cache has done so far
/// ([`IndexedCramReader::slice_cache_stats`](super::reader::IndexedCramReader::slice_cache_stats)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub struct SliceCacheStats {
    /// The budget in bytes as it stands now; [`SliceCacheBudget::Auto`]
    /// grows it with the handles and the largest slice.
    pub budget_bytes: usize,
    /// Slices decoded: one per cache miss, and every fetched slice when
    /// the budget is 0.
    pub decoded: u64,
    /// Slices the cache holds now.
    pub cached_slices: usize,
    /// Their heap bytes, which the budget bounds.
    pub cached_bytes: usize,
    /// The heap bytes of the largest slice decoded so far — what one slice
    /// of this file costs to keep.
    pub largest_slice_bytes: usize,
}

/// Decoded slices under a byte budget, least recently used evicted first.
pub(crate) struct SliceCache {
    inner: Mutex<Inner>,
    decodes: AtomicU64,
    largest: AtomicUsize,
    /// Readers sharing the cache: the one that opened the file and its live
    /// forks.
    pub(super) handles: AtomicUsize,
    /// Fetches that wanted a slice, and the slices they wanted in all.
    fetches: AtomicU64,
    fetched_slices: AtomicU64,
}

impl std::fmt::Debug for SliceCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SliceCache").field("stats", &self.stats()).finish_non_exhaustive()
    }
}

impl SliceCache {
    /// A cache for one open handle.
    pub(crate) fn new(budget: SliceCacheBudget) -> Self {
        Self {
            inner: Mutex::new(Inner { budget, ..Inner::default() }),
            decodes: AtomicU64::new(0),
            largest: AtomicUsize::new(0),
            handles: AtomicUsize::new(1),
            fetches: AtomicU64::new(0),
            fetched_slices: AtomicU64::new(0),
        }
    }

    /// A fork shares the cache.
    pub(crate) fn handle_opened(&self) {
        self.handles.fetch_add(1, Ordering::Relaxed);
    }

    /// A handle sharing the cache was dropped.
    pub(crate) fn handle_closed(&self) {
        self.handles.fetch_sub(1, Ordering::Relaxed);
    }

    /// A fetch wanted `slices` slices, which the automatic budget sizes by.
    pub(crate) fn record_fetch(&self, slices: usize) {
        if slices > 0 {
            self.fetches.fetch_add(1, Ordering::Relaxed);
            self.fetched_slices
                .fetch_add(u64::try_from(slices).unwrap_or(u64::MAX), Ordering::Relaxed);
        }
    }

    /// Slices [`SliceCacheBudget::Auto`] keeps room for.
    fn auto_slices(&self) -> usize {
        let handles = u64::try_from(self.handles.load(Ordering::Relaxed)).unwrap_or(u64::MAX);
        if handles <= 1 {
            return 1;
        }
        let fetches = self.fetches.load(Ordering::Relaxed);
        let beyond_first = self.fetched_slices.load(Ordering::Relaxed).saturating_sub(fetches);
        // Before any fetch, one slice each.
        let span = if fetches == 0 {
            handles
        } else {
            handles.saturating_mul(beyond_first).div_ceil(fetches)
        };
        usize::try_from(span.saturating_add(1)).unwrap_or(usize::MAX)
    }

    /// Change the budget, evicting down to it at once.
    pub(crate) fn set_budget(&self, budget: SliceCacheBudget) {
        let mut inner = self.lock();
        inner.budget = budget;
        let bytes = self.budget_bytes(&inner);
        evict(&mut inner, bytes);
    }

    /// What the budget allows now.
    fn budget_bytes(&self, inner: &Inner) -> usize {
        match inner.budget {
            SliceCacheBudget::Bytes(bytes) => bytes,
            // Never 0, which would turn the cache off before the first
            // decode tells it the slice size.
            SliceCacheBudget::Auto => {
                self.auto_slices().saturating_mul(self.largest.load(Ordering::Relaxed).max(1))
            }
        }
    }

    pub(crate) fn stats(&self) -> SliceCacheStats {
        let inner = self.lock();
        SliceCacheStats {
            budget_bytes: self.budget_bytes(&inner),
            decoded: self.decodes.load(Ordering::Relaxed),
            cached_slices: inner.entries.values().filter(|entry| entry.bytes > 0).count(),
            cached_bytes: inner.used,
            largest_slice_bytes: self.largest.load(Ordering::Relaxed),
        }
    }

    /// Decode one slice, counting it.
    fn decode(
        &self,
        decode: impl FnOnce() -> Result<DecodedSlice, CramError>,
    ) -> Result<Arc<DecodedSlice>, CramError> {
        self.decodes.fetch_add(1, Ordering::Relaxed);
        let slice = decode()?;
        self.largest.fetch_max(slice.heap_bytes(), Ordering::Relaxed);
        Ok(Arc::new(slice))
    }

    /// The slice for `key`, decoded with `decode` unless the cache holds it
    /// or another thread is decoding it (then this waits for that result).
    pub(crate) fn get_or_decode(
        &self,
        key: SliceKey,
        decode: impl FnOnce() -> Result<DecodedSlice, CramError>,
    ) -> Result<Arc<DecodedSlice>, CramError> {
        let cell = {
            let mut inner = self.lock();
            if self.budget_bytes(&inner) == 0 {
                drop(inner);
                return self.decode(decode);
            }
            inner.tick = inner.tick.wrapping_add(1);
            let tick = inner.tick;
            let entry = inner.entries.entry(key).or_insert_with(|| Entry {
                cell: Arc::default(),
                last_use: tick,
                bytes: 0,
            });
            entry.last_use = tick;
            Arc::clone(&entry.cell)
        };

        let mut decode = Some(decode);
        let mut error = None;
        let mut decoded_here = false;
        let slice = cell.get_or_init(|| {
            decoded_here = true;
            match decode.take().map(|decode| self.decode(decode)) {
                Some(Ok(slice)) => Some(slice),
                Some(Err(e)) => {
                    error = Some(e);
                    None
                }
                None => None,
            }
        });
        match slice {
            Some(slice) => {
                if decoded_here {
                    self.admit(key, &cell, slice.heap_bytes());
                }
                Ok(Arc::clone(slice))
            }
            None => {
                // Not cached, so a later fetch tries again.
                self.forget(key, &cell);
                match (error, decode.take()) {
                    (Some(e), _) => Err(e),
                    // Another thread's decode failed: decode here for this
                    // caller's own error.
                    (None, Some(decode)) => self.decode(decode),
                    // This thread ran the decode, which always leaves a slice
                    // or an error behind.
                    (None, None) => Err(CramError::Truncated { context: "slice cache decode" }),
                }
            }
        }
    }

    fn lock(&self) -> MutexGuard<'_, Inner> {
        // A panic while holding the lock leaves only bookkeeping behind:
        // entries and byte counts, never a half-written slice.
        self.inner.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Count a freshly decoded slice against the budget and evict down to
    /// it — possibly the slice itself, if it alone is over.
    fn admit(&self, key: SliceKey, cell: &Cell, bytes: usize) {
        let mut inner = self.lock();
        let budget = self.budget_bytes(&inner);
        let Some(entry) = inner.entries.get_mut(&key) else { return };
        if !Arc::ptr_eq(&entry.cell, cell) || entry.bytes != 0 {
            return;
        }
        // Counted as at least one byte, so a counted entry is never 0.
        entry.bytes = bytes.max(1);
        inner.used = inner.used.saturating_add(bytes.max(1));
        evict(&mut inner, budget);
    }

    fn forget(&self, key: SliceKey, cell: &Cell) {
        let mut inner = self.lock();
        if inner.entries.get(&key).is_some_and(|entry| Arc::ptr_eq(&entry.cell, cell))
            && let Some(entry) = inner.entries.remove(&key)
        {
            inner.used = inner.used.saturating_sub(entry.bytes);
        }
    }
}

/// Drop least recently used, counted entries until `used <= budget`. Slices
/// still being decoded are not counted yet and stay.
fn evict(inner: &mut Inner, budget: usize) {
    while inner.used > budget {
        let oldest = inner
            .entries
            .iter()
            .filter(|(_, entry)| entry.bytes > 0)
            .min_by_key(|(_, entry)| entry.last_use)
            .map(|(&key, _)| key);
        let Some(key) = oldest else { break };
        if let Some(entry) = inner.entries.remove(&key) {
            inner.used = inner.used.saturating_sub(entry.bytes);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(slice_offset: u64) -> SliceKey {
        SliceKey { container_offset: 0, slice_offset, tid: 0 }
    }

    fn slice() -> Result<DecodedSlice, CramError> {
        Ok(DecodedSlice::default())
    }

    #[test]
    fn a_cached_slice_is_decoded_once() {
        let cache = SliceCache::new(SliceCacheBudget::Bytes(1 << 20));
        let a = cache.get_or_decode(key(1), slice).unwrap();
        let b = cache.get_or_decode(key(1), || panic!("decoded twice")).unwrap();
        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(cache.stats().decoded, 1);
        cache.get_or_decode(key(2), slice).unwrap();
        assert_eq!(cache.stats().decoded, 2);
    }

    #[test]
    fn a_failed_decode_is_not_cached() {
        let cache = SliceCache::new(SliceCacheBudget::Bytes(1 << 20));
        let err = cache
            .get_or_decode(key(1), || Err(CramError::Truncated { context: "test" }))
            .unwrap_err();
        assert!(matches!(err, CramError::Truncated { context: "test" }));
        cache.get_or_decode(key(1), slice).unwrap();
        assert_eq!(cache.stats().decoded, 2);
    }

    #[test]
    fn budget_zero_decodes_every_time() {
        let cache = SliceCache::new(SliceCacheBudget::Bytes(0));
        cache.get_or_decode(key(1), slice).unwrap();
        cache.get_or_decode(key(1), slice).unwrap();
        assert_eq!(cache.stats().decoded, 2);
    }

    #[test]
    fn least_recently_used_goes_first() {
        // Every empty slice counts one byte: room for two.
        let cache = SliceCache::new(SliceCacheBudget::Bytes(2));
        cache.get_or_decode(key(1), slice).unwrap();
        cache.get_or_decode(key(2), slice).unwrap();
        cache.get_or_decode(key(1), || panic!("1 was evicted")).unwrap();
        cache.get_or_decode(key(3), slice).unwrap(); // evicts 2
        cache.get_or_decode(key(1), || panic!("1 was evicted")).unwrap();
        cache.get_or_decode(key(2), slice).unwrap();
        assert_eq!(cache.stats().decoded, 4);
        cache.set_budget(SliceCacheBudget::Bytes(0));
        assert_eq!(cache.lock().used, 0);
    }

    #[test]
    fn auto_budget_follows_handles_and_fetch_width() {
        let cache = SliceCache::new(SliceCacheBudget::Auto);
        let budget = || cache.stats().budget_bytes;
        assert_eq!(budget(), 1, "nothing decoded yet, but not off");

        let slice = 100 << 20;
        cache.largest.store(slice, Ordering::Relaxed);
        assert_eq!(budget(), slice, "a lone handle keeps its last slice");
        cache.record_fetch(4);
        assert_eq!(budget(), slice, "however wide its fetches");

        cache.handle_opened();
        cache.handle_opened();
        // Three handles, fetches four slices wide sharing one at the edge.
        assert_eq!(budget(), (3 * 3 + 1) * slice);
        // A mean of 1.25 beyond the first: 3.75 for three handles, rounded up.
        cache.record_fetch(1);
        cache.record_fetch(2);
        cache.record_fetch(0);
        cache.record_fetch(2);
        assert_eq!(budget(), (4 + 1) * slice);
        cache.handle_closed();
        assert_eq!(budget(), (3 + 1) * slice);
        cache.handle_closed();
        assert_eq!(budget(), slice);

        cache.set_budget(SliceCacheBudget::Bytes(7));
        assert_eq!(budget(), 7);
    }

    #[test]
    fn threads_wanting_one_slice_share_its_decode() {
        let cache = Arc::new(SliceCache::new(SliceCacheBudget::Bytes(1 << 20)));
        let barrier = Arc::new(std::sync::Barrier::new(8));
        let slices: Vec<_> = (0..8)
            .map(|_| {
                let (cache, barrier) = (Arc::clone(&cache), Arc::clone(&barrier));
                std::thread::spawn(move || {
                    barrier.wait();
                    cache
                        .get_or_decode(key(7), || {
                            std::thread::sleep(std::time::Duration::from_millis(50));
                            slice()
                        })
                        .unwrap()
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|t| t.join().unwrap())
            .collect();
        assert_eq!(cache.stats().decoded, 1);
        assert!(slices.iter().all(|s| Arc::ptr_eq(s, &slices[0])));
    }
}
