//! Decoded CRAM slices shared by a reader and its forks.
//!
//! Neighbouring region queries overlap the same slices at their edges, and a
//! slice is only useful once decoded whole, so without sharing every query
//! decodes each slice it touches again. The cache keeps decoded slices under
//! a byte budget, and makes threads that want the same missing slice at the
//! same time wait for one decode instead of each running their own.
// r[impl cram.slice_cache.shared]
// r[impl cram.slice_cache.budget]

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
}

/// Decoded slices under a byte budget, least recently used evicted first.
pub(crate) struct SliceCache {
    budget: AtomicUsize,
    inner: Mutex<Inner>,
    decodes: AtomicU64,
}

impl std::fmt::Debug for SliceCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SliceCache")
            .field("budget", &self.budget())
            .field("decodes", &self.decodes())
            .finish_non_exhaustive()
    }
}

impl SliceCache {
    pub(crate) fn new(budget: usize) -> Self {
        Self {
            budget: AtomicUsize::new(budget),
            inner: Mutex::new(Inner::default()),
            decodes: AtomicU64::new(0),
        }
    }

    pub(crate) fn budget(&self) -> usize {
        self.budget.load(Ordering::Relaxed)
    }

    /// Change the budget, evicting down to it at once.
    pub(crate) fn set_budget(&self, budget: usize) {
        self.budget.store(budget, Ordering::Relaxed);
        let mut inner = self.lock();
        evict(&mut inner, budget);
    }

    /// How many slices were decoded through the cache: one per miss.
    pub(crate) fn decodes(&self) -> u64 {
        self.decodes.load(Ordering::Relaxed)
    }

    /// The slice for `key`, decoded with `decode` unless the cache holds it
    /// or another thread is decoding it (then this waits for that result).
    pub(crate) fn get_or_decode(
        &self,
        key: SliceKey,
        decode: impl FnOnce() -> Result<DecodedSlice, CramError>,
    ) -> Result<Arc<DecodedSlice>, CramError> {
        if self.budget() == 0 {
            self.decodes.fetch_add(1, Ordering::Relaxed);
            return decode().map(Arc::new);
        }
        let cell = {
            let mut inner = self.lock();
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
            self.decodes.fetch_add(1, Ordering::Relaxed);
            match decode.take().map(|decode| decode()) {
                Some(Ok(slice)) => Some(Arc::new(slice)),
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
                    (None, Some(decode)) => {
                        self.decodes.fetch_add(1, Ordering::Relaxed);
                        decode().map(Arc::new)
                    }
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
        let budget = self.budget();
        let mut inner = self.lock();
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
        let cache = SliceCache::new(1 << 20);
        let a = cache.get_or_decode(key(1), slice).unwrap();
        let b = cache.get_or_decode(key(1), || panic!("decoded twice")).unwrap();
        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(cache.decodes(), 1);
        cache.get_or_decode(key(2), slice).unwrap();
        assert_eq!(cache.decodes(), 2);
    }

    #[test]
    fn a_failed_decode_is_not_cached() {
        let cache = SliceCache::new(1 << 20);
        let err = cache
            .get_or_decode(key(1), || Err(CramError::Truncated { context: "test" }))
            .unwrap_err();
        assert!(matches!(err, CramError::Truncated { context: "test" }));
        cache.get_or_decode(key(1), slice).unwrap();
        assert_eq!(cache.decodes(), 2);
    }

    #[test]
    fn budget_zero_decodes_every_time() {
        let cache = SliceCache::new(0);
        cache.get_or_decode(key(1), slice).unwrap();
        cache.get_or_decode(key(1), slice).unwrap();
        assert_eq!(cache.decodes(), 2);
    }

    #[test]
    fn least_recently_used_goes_first() {
        // Every empty slice counts one byte: room for two.
        let cache = SliceCache::new(2);
        cache.get_or_decode(key(1), slice).unwrap();
        cache.get_or_decode(key(2), slice).unwrap();
        cache.get_or_decode(key(1), || panic!("1 was evicted")).unwrap();
        cache.get_or_decode(key(3), slice).unwrap(); // evicts 2
        cache.get_or_decode(key(1), || panic!("1 was evicted")).unwrap();
        cache.get_or_decode(key(2), slice).unwrap();
        assert_eq!(cache.decodes(), 4);
        cache.set_budget(0);
        assert_eq!(cache.lock().used, 0);
    }

    #[test]
    fn threads_wanting_one_slice_share_its_decode() {
        let cache = Arc::new(SliceCache::new(1 << 20));
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
        assert_eq!(cache.decodes(), 1);
        assert!(slices.iter().all(|s| Arc::ptr_eq(s, &slices[0])));
    }
}
