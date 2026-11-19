//! Walking a BCF INFO or FORMAT section entry by entry, and the lazily built
//! index of one (`r[vcf_fields.lazy]`, `r[vcf_values.index+2]`).

use super::bcf_value::{BcfValueError, FormatColumn, TypedValue, read_typed_int};
use super::error::{KeyKind, RecordError, RecordPart};
use super::file::FileShared;
use super::index::{DictId, Entry, FieldIndex};
use crate::io::text_scan::{for_each_info_entry, for_each_of_two};
use std::cell::UnsafeCell;
use std::ops::{ControlFlow, Range};
use std::sync::atomic::{AtomicU8, Ordering};

/// The entries of a BCF INFO section (`samples: None`) or FORMAT section
/// (`Some(n)` samples): each key and the byte span of its value.
#[derive(Debug, Clone)]
pub(crate) struct BcfEntries<'a> {
    bytes: &'a [u8],
    at: usize,
    left: usize,
    samples: Option<usize>,
}

impl<'a> BcfEntries<'a> {
    /// The `n` entries of `bytes`, from `cursor` on.
    pub(crate) fn new(bytes: &'a [u8], n: usize, samples: Option<usize>, cursor: Cursor) -> Self {
        Self { bytes, at: cursor.at, left: n.saturating_sub(cursor.walked), samples }
    }

    /// Where the entries walked so far end.
    pub(crate) fn end(&self) -> usize {
        self.at
    }

    #[inline]
    fn entry(&mut self) -> Result<(i32, Range<usize>), BcfValueError> {
        match self.quick() {
            Some((key, span)) => {
                self.at = span.end;
                Ok((key, span))
            }
            None => self.checked_entry(),
        }
    }

    /// The common entry, read without building an error: a non-negative
    /// int8 or int16 key, then a value of a known type whose count is inline
    /// or a non-negative int8/int16, ending within the bytes. `None` sends
    /// everything else, malformed entries included, to [`Self::checked_entry`],
    /// which finds the same span or the error.
    #[inline]
    fn quick(&self) -> Option<(i32, Range<usize>)> {
        let rest = self.bytes.get(self.at..)?;
        let (key, key_len) = quick_int(rest)?;
        let value = rest.get(key_len..)?;
        let &descriptor = value.first()?;
        let width = WIDTH.get(usize::from(descriptor & 0x0F)).copied().flatten()?;
        let (count, header) = match descriptor >> 4 {
            15 => {
                let (count, used) = quick_int(value.get(1..)?)?;
                (usize::try_from(count).ok()?, used.checked_add(1)?)
            }
            n => (usize::from(n), 1),
        };
        let per_sample = count.checked_mul(usize::from(width))?;
        let len = match self.samples {
            None => per_sample,
            Some(n) => per_sample.checked_mul(n)?,
        };
        let start = self.at.checked_add(key_len)?;
        let end = start.checked_add(header)?.checked_add(len)?;
        (end <= self.bytes.len()).then_some((key, start..end))
    }

    #[cold]
    #[inline(never)]
    fn checked_entry(&mut self) -> Result<(i32, Range<usize>), BcfValueError> {
        let rest = self.bytes.get(self.at..).unwrap_or_default();
        let (key, key_len) = read_typed_int(rest)?;
        let value = rest.get(key_len..).unwrap_or_default();
        let len = match self.samples {
            None => TypedValue::parse(value)?.1,
            Some(n) => FormatColumn::parse(value, n)?.1,
        };
        let start = self.at.saturating_add(key_len);
        let end = start.saturating_add(len);
        self.at = end;
        Ok((key, start..end))
    }
}

/// Bytes per element by BCF type code; `None` for a code with no type.
const WIDTH: [Option<u8>; 16] = [
    Some(0), // missing
    Some(1), // int8
    Some(2), // int16
    Some(4), // int32
    None,
    Some(4), // float
    None,
    Some(1), // char
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
];

/// A typed scalar int8 or int16 that is not negative (so neither missing
/// nor end-of-vector), and its encoded length.
#[inline]
fn quick_int(bytes: &[u8]) -> Option<(i32, usize)> {
    match *bytes {
        [0x11, b, ..] if b < 0x80 => Some((i32::from(b), 2)),
        [0x12, lo, hi, ..] if hi < 0x80 => Some((i32::from(i16::from_le_bytes([lo, hi])), 3)),
        _ => None,
    }
}

impl Iterator for BcfEntries<'_> {
    type Item = Result<(i32, Range<usize>), BcfValueError>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.left = self.left.checked_sub(1)?;
        let entry = self.entry();
        if entry.is_err() {
            self.left = 0;
        }
        Some(entry)
    }
}

/// How far an index got: it stopped at the key asked for, or indexed the
/// whole section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scan {
    Found,
    End,
}

/// Where an incremental index stopped: the byte after the last entry
/// indexed, and how many entries that was.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Cursor {
    pub(crate) at: usize,
    pub(crate) walked: usize,
}

// r[impl vcf_values.index+2]
/// Index a BCF INFO or FORMAT section from `cursor` on, stopping after the
/// first entry of `stop`: every key must be a declared field of its kind
/// (`r[vcf_reader.field_table]`), and the entries must fill the section.
pub(crate) fn index_bcf(
    index: &mut FieldIndex,
    mut entries: BcfEntries<'_>,
    cursor: &mut Cursor,
    stop: Option<DictId>,
    kind: KeyKind,
    file: &FileShared,
) -> Result<Scan, RecordError> {
    let part = if kind == KeyKind::Format { RecordPart::Format } else { RecordPart::Info };
    let left = entries.left;
    for entry in &mut entries {
        let (key, span) = entry.map_err(|source| RecordError::Value { part, source })?;
        let id = u32::try_from(key)
            .ok()
            .and_then(DictId::new)
            .filter(|&id| file.fields.declared(kind, id).is_some())
            .ok_or(RecordError::UnknownKey { index: key, kind })?;
        let (Ok(start), Ok(end)) = (u32::try_from(span.start), u32::try_from(span.end)) else {
            return Err(RecordError::TooLarge { size: span.end as u64, max: u64::from(u32::MAX) });
        };
        if !index.push(Entry { id, start, end }) {
            warn_duplicate_id(file, id);
        } else if stop == Some(id) {
            cursor.at = entries.at;
            cursor.walked = cursor.walked.saturating_add(left.saturating_sub(entries.left));
            return Ok(Scan::Found);
        }
    }
    let len = entries.bytes.len().saturating_sub(entries.end());
    if len > 0 {
        return Err(RecordError::TrailingBytes { part, len });
    }
    Ok(Scan::End)
}

#[cold]
fn warn_duplicate_id(file: &FileShared, id: DictId) {
    file.warn_duplicate(file.fields.name(&file.header, id).unwrap_or_default());
}

/// An index a [`LazyIndex`] builds and keeps the buffers of.
pub(crate) trait Rebuild: Default {
    /// Start indexing a new section.
    fn restart(&mut self);
}

impl Rebuild for FieldIndex {
    fn restart(&mut self) {
        self.begin();
    }
}

/// The FORMAT fields of a record: which are present, and for VCF text
/// where each sample's value of each lies.
#[derive(Debug, Clone, Default)]
pub(crate) struct FormatIndex {
    /// For BCF, each field's column span; for VCF text, each field's
    /// position among the FORMAT keys (`k..k + 1`).
    pub(crate) keys: FieldIndex,
    /// VCF text only: sample `s`'s value of the `k`th key at
    /// `s * n_keys + k`, empty where the sample dropped it.
    pub(crate) values: Vec<Range<u32>>,
    pub(crate) n_keys: usize,
}

impl Rebuild for FormatIndex {
    fn restart(&mut self) {
        self.keys.begin();
        self.values.clear();
        self.n_keys = 0;
    }
}

/// [`LazyIndex`] states.
const EMPTY: u8 = 0;
const BUILDING: u8 = 1;
const PARTIAL: u8 = 2;
const DONE: u8 = 3;

/// Indexes a section from a cursor on, stopping after the first entry of
/// the key given, if any.
pub(crate) trait Step<T>:
    FnOnce(&mut T, &mut Cursor, Option<DictId>) -> Result<Scan, RecordError>
{
}
impl<T, F: FnOnce(&mut T, &mut Cursor, Option<DictId>) -> Result<Scan, RecordError>> Step<T> for F {}

// r[impl vcf_fields.lazy+2]
// r[impl vcf_values.index+2]
/// A section's index, built on first access, in place: its buffers stay
/// for the next record's build. A lookup ([`find`](LazyIndex::find)) indexes
/// no further than the key it asks for; [`get`](Self::get) indexes the rest.
/// A failed build is kept and returned on every access that reaches it.
///
/// One atomic byte guards the cells: `EMPTY`, `BUILDING` while the one
/// thread that won the swap from `EMPTY` or `PARTIAL` indexes, `PARTIAL`
/// between lookups that stopped early, and `DONE` once the section is
/// indexed, until [`reset`](Self::reset). A `std::sync::OnceLock` plus a
/// mutex for the buffers cost four locked atomics per record; this costs
/// one swap per access until `DONE`, and none after.
pub(crate) struct LazyIndex<T = FieldIndex> {
    state: AtomicU8,
    index: UnsafeCell<T>,
    cursor: UnsafeCell<Cursor>,
    /// The build's outcome, written by the builder before `DONE`.
    result: UnsafeCell<Result<(), RecordError>>,
}

// SAFETY: the cells are written only by the thread that moved `state` to
// `BUILDING` (or through `&mut self`), and read either by that thread or
// after an `Acquire` load saw the `Release` store of `DONE`, after which
// nothing writes them until `reset` takes `&mut self`. So `&LazyIndex` hands
// out `&T` and `&RecordError` across threads, which needs `T: Sync`, and the
// builder may be another thread than the owner, which needs `T: Send`.
unsafe impl<T: Send + Sync> Sync for LazyIndex<T> {}

impl<T: Default> Default for LazyIndex<T> {
    fn default() -> Self {
        Self {
            state: AtomicU8::new(EMPTY),
            index: UnsafeCell::new(T::default()),
            cursor: UnsafeCell::new(Cursor::default()),
            result: UnsafeCell::new(Ok(())),
        }
    }
}

impl<T: Clone + Default> Clone for LazyIndex<T> {
    /// A built index is cloned; one not built (or partly, or being built)
    /// is not, and the clone builds its own on first access.
    fn clone(&self) -> Self {
        match self.done() {
            Some((index, result)) => Self {
                state: AtomicU8::new(DONE),
                index: UnsafeCell::new(index.clone()),
                cursor: UnsafeCell::new(Cursor::default()),
                result: UnsafeCell::new(result.clone()),
            },
            None => Self::default(),
        }
    }
}

impl<T: std::fmt::Debug> std::fmt::Debug for LazyIndex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.done() {
            Some((index, Ok(()))) => f.debug_tuple("LazyIndex").field(index).finish(),
            Some((_, Err(e))) => f.debug_tuple("LazyIndex").field(e).finish(),
            None => f.write_str("LazyIndex(<not built>)"),
        }
    }
}

impl<T> LazyIndex<T> {
    /// The index and the build's outcome, once the whole section is indexed.
    #[inline]
    fn done(&self) -> Option<(&T, &Result<(), RecordError>)> {
        if self.state.load(Ordering::Acquire) != DONE {
            return None;
        }
        // SAFETY: `DONE` was stored (Release) after the last write to the
        // cells, and nothing writes them again before `reset(&mut self)`.
        Some(unsafe { (&*self.index.get(), &*self.result.get()) })
    }
}

/// Holds `BUILDING`; stores `next` when dropped: `PARTIAL` or `DONE` as the
/// build got, `EMPTY` (the default) if it panicked, so another access
/// builds again instead of waiting for ever.
struct Building<'a> {
    state: &'a AtomicU8,
    next: u8,
}

impl Drop for Building<'_> {
    fn drop(&mut self) {
        self.state.store(self.next, Ordering::Release);
    }
}

impl<T: Rebuild> LazyIndex<T> {
    /// Forget the index of the last record, keeping its buffers.
    #[inline]
    pub(crate) fn reset(&mut self) {
        *self.state.get_mut() = EMPTY;
    }

    /// The whole index, built by `step` on first access.
    #[inline]
    pub(crate) fn get(&self, step: impl Step<T>) -> Result<&T, &RecordError> {
        let (index, result) = match self.done() {
            Some(done) => done,
            None => self.build(step),
        };
        result.as_ref().map(|()| index)
    }

    #[cold]
    #[inline(never)]
    fn build(&self, step: impl Step<T>) -> (&T, &Result<(), RecordError>) {
        let Some(mut building) = self.lock() else { return self.wait_done() };
        // SAFETY: `lock` gave this thread `BUILDING`, so it is the only one
        // touching the cells until `building` is dropped.
        unsafe {
            let result = match step(&mut *self.index.get(), &mut *self.cursor.get(), None) {
                Ok(_) => Ok(()),
                Err(e) => Err(e),
            };
            *self.result.get() = result;
        }
        building.next = DONE;
        drop(building);
        self.wait_done()
    }

    /// `BUILDING`, taken from `EMPTY` (restarting the index) or `PARTIAL`;
    /// `None` once the index is `DONE`.
    fn lock(&self) -> Option<Building<'_>> {
        let mut from = EMPTY;
        loop {
            match self.state.compare_exchange_weak(
                from,
                BUILDING,
                Ordering::Acquire,
                Ordering::Acquire,
            ) {
                Ok(_) => break,
                Err(DONE) => return None,
                Err(now @ (EMPTY | PARTIAL)) => from = now,
                // Another thread indexes.
                Err(_) => std::thread::yield_now(),
            }
        }
        if from == EMPTY {
            // SAFETY: this thread holds `BUILDING` (see `build`).
            unsafe {
                (*self.index.get()).restart();
                *self.cursor.get() = Cursor::default();
            }
        }
        Some(Building { state: &self.state, next: EMPTY })
    }

    /// The index once `DONE`, which the caller saw or stored.
    fn wait_done(&self) -> (&T, &Result<(), RecordError>) {
        loop {
            if let Some(done) = self.done() {
                return done;
            }
            std::thread::yield_now();
        }
    }

    /// Whether anything has been indexed (or failed), for tests.
    #[cfg(test)]
    pub(crate) fn is_built(&self) -> bool {
        self.state.load(Ordering::Acquire) != EMPTY
    }
}

impl LazyIndex<FieldIndex> {
    /// The first entry of `id`, indexing the section no further than that
    /// (`r[vcf_values.index+2]`). Malformed data fails the lookup only if the
    /// index meets it before `id`.
    #[inline]
    pub(crate) fn find(
        &self,
        id: DictId,
        step: impl Step<FieldIndex>,
    ) -> Result<Option<Entry>, &RecordError> {
        match self.done() {
            Some((index, result)) => found(index, result, id),
            None => self.find_on(id, step),
        }
    }

    #[cold]
    #[inline(never)]
    fn find_on(
        &self,
        id: DictId,
        step: impl Step<FieldIndex>,
    ) -> Result<Option<Entry>, &RecordError> {
        let Some(mut building) = self.lock() else {
            let (index, result) = self.wait_done();
            return found(index, result, id);
        };
        // SAFETY: `lock` gave this thread `BUILDING`.
        let (index, cursor) = unsafe { (&mut *self.index.get(), &mut *self.cursor.get()) };
        if let Some(entry) = index.get(id) {
            building.next = PARTIAL;
            return Ok(Some(entry));
        }
        let (entry, result) = match step(index, cursor, Some(id)) {
            Ok(Scan::Found) => {
                building.next = PARTIAL;
                return Ok(index.get(id));
            }
            Ok(Scan::End) => (index.get(id), Ok(())),
            Err(e) => (None, Err(e)),
        };
        // SAFETY: as above.
        unsafe { *self.result.get() = result };
        building.next = DONE;
        drop(building);
        match entry {
            Some(entry) => Ok(Some(entry)),
            None => {
                let (index, result) = self.wait_done();
                found(index, result, id)
            }
        }
    }
}

/// What a lookup of `id` in a whole index reads: its entry if it has one
/// (indexed before any malformed data), else the build's error or absent.
fn found<'a>(
    index: &FieldIndex,
    result: &'a Result<(), RecordError>,
    id: DictId,
) -> Result<Option<Entry>, &'a RecordError> {
    match (index.get(id), result) {
        (Some(entry), _) => Ok(Some(entry)),
        (None, Ok(())) => Ok(None),
        (None, Err(e)) => Err(e),
    }
}

/// The span of a text entry or field, as `u32`s.
fn span32(r: Range<usize>) -> Result<(u32, u32), RecordError> {
    match (u32::try_from(r.start), u32::try_from(r.end)) {
        (Ok(s), Ok(e)) => Ok((s, e)),
        _ => Err(RecordError::TooLarge { size: r.end as u64, max: u64::from(u32::MAX) }),
    }
}

// r[impl vcf_reader.field_table]
// r[impl vcf_values.duplicates+2]
/// Index a VCF text INFO column: each key's span is what follows it (empty,
/// or `=` and the value). A key the field table does not hold yet is
/// added to it; `.` and empty entries hold nothing.
pub(crate) fn index_text_info(
    index: &mut FieldIndex,
    info: &[u8],
    cursor: &mut Cursor,
    stop: Option<DictId>,
    file: &FileShared,
) -> Result<Scan, RecordError> {
    let base = cursor.at;
    let Some(rest) = info.get(base..).filter(|_| info != b".") else { return Ok(Scan::End) };
    let flow = for_each_info_entry(rest, |entry| {
        let (key_start, key_end) =
            (base.saturating_add(entry.key.start), base.saturating_add(entry.key.end));
        let key = info.get(key_start..key_end).unwrap_or_default();
        if key.is_empty() {
            return ControlFlow::Continue(());
        }
        let Some(id) = file.fields.text_id(KeyKind::Info, key) else {
            return ControlFlow::Continue(());
        };
        let end = entry.value.map_or(key_end, |v| base.saturating_add(v.end));
        match span32(key_end..end) {
            Ok((start, span_end)) => {
                if !index.push(Entry { id, start, end: span_end }) {
                    file.warn_duplicate(std::str::from_utf8(key).unwrap_or_default());
                } else if stop == Some(id) {
                    return ControlFlow::Break(Ok(end.saturating_add(1)));
                }
                ControlFlow::Continue(())
            }
            Err(e) => ControlFlow::Break(Err(e)),
        }
    });
    match flow {
        ControlFlow::Break(Ok(next)) => {
            cursor.at = next;
            Ok(Scan::Found)
        }
        ControlFlow::Break(Err(e)) => Err(e),
        ControlFlow::Continue(()) => Ok(Scan::End),
    }
}

// r[impl vcf_fields.samples]
// r[impl vcf_reader.field_table]
/// Index VCF text FORMAT data: `keys` (the FORMAT column) and `samples`
/// (`n_samples` TAB-separated columns), in one pass over the samples. A
/// sample may drop trailing fields; one with more fields than FORMAT has
/// keys, or a column count other than `n_samples`, is malformed.
pub(crate) fn index_text_format(
    index: &mut FormatIndex,
    keys: &[u8],
    samples: &[u8],
    n_samples: usize,
    file: &FileShared,
) -> Result<(), RecordError> {
    if keys.is_empty() || keys == b"." {
        return Ok(());
    }
    let mut n_keys = 0usize;
    for (k, key) in keys.split(|&b| b == b':').enumerate() {
        n_keys = k.saturating_add(1);
        let Some(id) = file.fields.text_id(KeyKind::Format, key) else {
            continue;
        };
        let (start, end) = span32(k..n_keys)?;
        if !index.keys.push(Entry { id, start, end }) {
            file.warn_duplicate(std::str::from_utf8(key).unwrap_or_default());
        }
    }
    index.n_keys = n_keys;
    let total = n_samples.saturating_mul(n_keys);
    index.values.resize(total, 0..0);
    let (mut sample, mut field, mut start) = (0usize, 0usize, 0usize);
    let flow = for_each_of_two(samples, b'\t', b':', |pos, next_sample| {
        let Some(slot) = index.values.get_mut(sample.saturating_mul(n_keys).saturating_add(field))
        else {
            return ControlFlow::Break(too_many(sample, field, n_samples, n_keys));
        };
        match span32(start..pos) {
            Ok((s, e)) => *slot = s..e,
            Err(e) => return ControlFlow::Break(e),
        }
        start = pos.saturating_add(1);
        if next_sample {
            sample = sample.saturating_add(1);
            field = 0;
        } else {
            field = field.saturating_add(1);
            if field >= n_keys {
                return ControlFlow::Break(too_many(sample, field, n_samples, n_keys));
            }
        }
        ControlFlow::Continue(())
    });
    if let ControlFlow::Break(e) = flow {
        return Err(e);
    }
    let found = sample.saturating_add(1);
    if found != n_samples {
        return Err(RecordError::SampleCount {
            found: u32::try_from(found).unwrap_or(u32::MAX),
            expected: n_samples,
        });
    }
    let (s, e) = span32(start..samples.len())?;
    if let Some(slot) = index.values.get_mut(sample.saturating_mul(n_keys).saturating_add(field)) {
        *slot = s..e;
    }
    Ok(())
}

/// Sample `sample` holds a field past the last key, or there are more
/// samples than `n_samples`.
fn too_many(sample: usize, field: usize, n_samples: usize, n_keys: usize) -> RecordError {
    if sample >= n_samples {
        RecordError::SampleCount {
            found: u32::try_from(sample.saturating_add(1)).unwrap_or(u32::MAX),
            expected: n_samples,
        }
    } else {
        RecordError::SampleFields { sample, found: field.saturating_add(1), keys: n_keys }
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::arithmetic_side_effects, reason = "test sizes below 1000")]

    use super::*;
    use crate::vcf::bcf_encoding::{encode_type_byte, encode_typed_int_key};
    use hegel::prelude::*;

    /// Walk `bytes` the quick way where it applies, and every entry the
    /// checked way too: both must agree on each span, and the walk on the
    /// first error.
    fn walk(bytes: &[u8], n: usize, samples: Option<usize>) -> (usize, usize) {
        let mut entries = BcfEntries::new(bytes, n, samples, Cursor::default());
        let mut quick = 0;
        for _ in 0..n {
            let mut checked = entries.clone();
            let expected = checked.checked_entry();
            if let Some(got) = entries.quick() {
                assert_eq!(Ok(got), expected, "at {}", entries.at);
                quick += 1;
            }
            match entries.entry() {
                Ok(_) => assert_eq!(entries.at, checked.at),
                Err(e) => {
                    assert_eq!(Err(e), expected);
                    break;
                }
            }
        }
        (quick, entries.at)
    }

    // r[verify vcf_values.bcf_sentinels]
    /// Generated entries (keys and counts at every width, every type code,
    /// cut short or with a byte flipped) read the same through the quick
    /// path as through the checked one.
    #[hegel::test]
    fn the_quick_path_agrees_with_the_checked_one(tc: TestCase) {
        let samples = tc.draw(gs::optional(gs::integers::<usize>().min_value(1).max_value(4)));
        let n = tc.draw(gs::integers::<usize>().max_value(6));
        let mut bytes = Vec::new();
        for _ in 0..n {
            let key = tc.draw(gs::sampled_from(vec![0u32, 1, 127, 128, 32767, 32768, 70000]));
            encode_typed_int_key(&mut bytes, key);
            let code = tc.draw(gs::integers::<u8>().max_value(15));
            let count = tc.draw(gs::sampled_from(vec![0usize, 1, 3, 14, 15, 200, 40000]));
            encode_type_byte(&mut bytes, count, code);
            let width = WIDTH.get(usize::from(code)).copied().flatten().unwrap_or(1);
            let len = count * usize::from(width) * samples.unwrap_or(1);
            bytes.extend(std::iter::repeat_n(7u8, len.min(1000)));
        }
        if !bytes.is_empty() && tc.draw(gs::booleans()) {
            let at = tc.draw(gs::integers::<usize>().max_value(bytes.len() - 1));
            let flip = tc.draw(gs::integers::<u8>());
            if let Some(b) = bytes.get_mut(at) {
                *b ^= flip;
            }
        }
        let cut = tc.draw(gs::integers::<usize>().max_value(bytes.len()));
        walk(bytes.get(..cut).unwrap_or_default(), n, samples);
    }

    // r[verify vcf_fields.lazy+2]
    /// Threads racing for a fresh record's index see one build, its index or
    /// its error, record after record on the same buffers.
    #[test]
    fn racing_first_accesses_build_once() {
        use std::sync::atomic::AtomicUsize;
        let rounds: u32 = if cfg!(miri) { 4 } else { 300 };
        let threads = if cfg!(miri) { 3 } else { 8 };
        let mut lazy = LazyIndex::<FieldIndex>::default();
        for round in 0..rounds {
            lazy.reset();
            let builds = AtomicUsize::new(0);
            let id = DictId::new(round % 50).unwrap();
            let entry = Entry { id, start: round, end: round + 1 };
            let fails = round % 3 == 0;
            std::thread::scope(|s| {
                for _ in 0..threads {
                    s.spawn(|| {
                        let got = lazy.get(|index: &mut FieldIndex, _: &mut Cursor, _| {
                            builds.fetch_add(1, Ordering::Relaxed);
                            assert!(index.entries().is_empty(), "restarted");
                            index.push(entry);
                            if fails { Err(RecordError::Truncated) } else { Ok(Scan::End) }
                        });
                        match got {
                            Ok(index) => {
                                assert!(!fails);
                                assert_eq!(index.get(id), Some(entry));
                            }
                            Err(e) => assert!(fails && matches!(e, RecordError::Truncated)),
                        }
                    });
                }
            });
            assert_eq!(builds.load(Ordering::Relaxed), 1, "round {round}");
            assert!(lazy.is_built());
            if !fails {
                assert_eq!(
                    lazy.clone()
                        .get(|_: &mut FieldIndex, _: &mut Cursor, _| unreachable!())
                        .unwrap()
                        .get(id),
                    Some(entry)
                );
            }
        }
    }

    /// A build that panics leaves the index unbuilt; the next access builds it.
    #[test]
    fn a_panicking_build_is_retried() {
        let lazy = LazyIndex::<FieldIndex>::default();
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = lazy.get(|_: &mut FieldIndex, _: &mut Cursor, _| panic!("build"));
        }));
        assert!(panicked.is_err());
        assert!(!lazy.is_built());
        assert!(lazy.get(|_: &mut FieldIndex, _: &mut Cursor, _| Ok(Scan::End)).is_ok());
    }

    /// The common entries take the quick path, a flag ending the section
    /// (three bytes, nothing after them) included.
    #[test]
    fn common_entries_take_the_quick_path() {
        // DP=140 as int16, an int8 vector of 20 (overflow count), a flag.
        let info = [0x11, 3, 0x12, 140, 0, 0x11, 9, 0xF1, 0x11, 20];
        let mut bytes = info.to_vec();
        bytes.extend(1..=20u8);
        bytes.extend([0x11, 4, 0x00]);
        assert_eq!(walk(&bytes, 3, None), (3, bytes.len()));
        // A FORMAT column of two samples, int16 key.
        let format = [0x12, 0x00, 0x01, 0x21, 1, 2, 3, 4];
        assert_eq!(walk(&format, 1, Some(2)), (1, format.len()));
    }
}

#[cfg(test)]
mod incremental {
    #![allow(
        clippy::arithmetic_side_effects,
        clippy::cast_possible_truncation,
        reason = "test sizes below 100"
    )]

    use super::*;
    use crate::vcf::bcf_encoding::{BCF_BT_INT8, encode_type_byte, encode_typed_int_key};
    use crate::vcf::reader::VcfFile;
    use hegel::prelude::*;

    const KEYS: usize = 8;

    /// A file whose header declares INFO keys `K0`..`K7`.
    fn file() -> VcfFile {
        use std::fmt::Write as _;
        let mut text = String::from("##fileformat=VCFv4.3\n##contig=<ID=chr1>\n");
        for k in 0..KEYS {
            writeln!(text, "##INFO=<ID=K{k},Number=.,Type=Integer,Description=\"\">").unwrap();
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        VcfFile::from_bytes(text.into_bytes(), None).unwrap()
    }

    /// The id of `Kk`, or of a key the header does not declare.
    fn id(file: &FileShared, k: usize) -> DictId {
        file.fields.get(format!("K{k}").as_bytes()).unwrap_or(DictId::new(60_000).unwrap())
    }

    /// One access: a lookup of key `k`, or (`None`) the whole index.
    fn accesses(tc: &TestCase) -> Vec<Option<usize>> {
        tc.draw(gs::vecs(gs::optional(gs::integers::<usize>().max_value(KEYS + 1))).max_size(6))
    }

    /// Lookups (and whole-index accesses) on one lazy index, against the
    /// index built whole at once: a lookup reads what the whole index has
    /// for its key, or its error if it has none.
    fn check(
        file: &FileShared,
        accesses: &[Option<usize>],
        step: impl Fn(&mut FieldIndex, &mut Cursor, Option<DictId>) -> Result<Scan, RecordError>,
    ) {
        let mut whole = FieldIndex::default();
        whole.begin();
        let result = step(&mut whole, &mut Cursor::default(), None).map(drop);
        let mut lazy = LazyIndex::<FieldIndex>::default();
        for _ in 0..2 {
            lazy.reset();
            for &access in accesses {
                match access {
                    Some(k) => {
                        let id = id(file, k);
                        let expected = found(&whole, &result, id).map_err(Clone::clone);
                        assert_eq!(lazy.find(id, &step).map_err(Clone::clone), expected, "K{k}");
                    }
                    None => {
                        let got = lazy.get(&step).map(|i| i.entries().to_vec());
                        let expected = result.as_ref().map(|()| whole.entries().to_vec());
                        assert_eq!(got, expected);
                    }
                }
            }
        }
    }

    // r[verify vcf_values.index+2]
    // r[verify vcf_values.duplicates+2]
    /// BCF INFO sections of declared keys (repeats included), maybe with an
    /// undeclared key, a value cut short, a wrong entry count or trailing
    /// bytes.
    #[hegel::test]
    fn bcf_lookups_read_what_the_whole_index_reads(tc: TestCase) {
        let file = file();
        let file = &*file.0;
        let keys = tc.draw(gs::vecs(gs::integers::<usize>().max_value(KEYS)).max_size(12));
        let mut bytes = Vec::new();
        for &k in &keys {
            encode_typed_int_key(&mut bytes, id(file, k).get().min(100));
            let len = tc.draw(gs::integers::<usize>().max_value(3));
            encode_type_byte(&mut bytes, len, BCF_BT_INT8);
            bytes.extend(std::iter::repeat_n(1u8, len));
        }
        let cut = tc.draw(
            gs::integers::<usize>().min_value(bytes.len().saturating_sub(2)).max_value(bytes.len()),
        );
        bytes.truncate(cut);
        bytes.extend(tc.draw(gs::binary().max_size(2)));
        let n = (keys.len() + tc.draw(gs::integers::<usize>().max_value(2))).saturating_sub(1);
        check(file, &accesses(&tc), |index, cursor, stop| {
            index_bcf(
                index,
                BcfEntries::new(&bytes, n, None, *cursor),
                cursor,
                stop,
                KeyKind::Info,
                file,
            )
        });
    }

    // r[verify vcf_values.index+2]
    // r[verify vcf_values.duplicates+2]
    /// VCF text INFO columns of declared keys, flags and values, repeats
    /// and empty entries included.
    #[hegel::test]
    fn text_lookups_read_what_the_whole_index_reads(tc: TestCase) {
        let file = file();
        let file = &*file.0;
        let entries = tc.draw(gs::vecs(gs::integers::<usize>().max_value(KEYS + 1)).max_size(12));
        let info: Vec<String> = entries
            .iter()
            .map(|&k| match k {
                KEYS => String::new(),
                k if k % 3 == 0 => format!("K{k}"),
                k => format!("K{k}={k},{k}"),
            })
            .collect();
        let info = info.join(";");
        check(file, &accesses(&tc), |index, cursor, stop| {
            index_text_info(index, info.as_bytes(), cursor, stop, file)
        });
    }
}
