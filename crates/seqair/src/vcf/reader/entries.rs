//! Walking a BCF INFO or FORMAT section entry by entry, and the lazily built
//! index of one (`r[vcf_fields.lazy]`, `r[vcf_values.index]`).

use super::bcf_value::{BcfValueError, FormatColumn, TypedValue, read_typed_int};
use super::error::{KeyKind, RecordError, RecordPart};
use super::file::FileShared;
use super::index::{DictId, Entry, FieldIndex};
use crate::io::text_scan::{for_each_info_entry, for_each_of_two};
use std::ops::{ControlFlow, Range};
use std::sync::{Mutex, OnceLock, PoisonError};

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
    pub(crate) fn new(bytes: &'a [u8], n: usize, samples: Option<usize>) -> Self {
        Self { bytes, at: 0, left: n, samples }
    }

    /// Where the entries walked so far end.
    pub(crate) fn end(&self) -> usize {
        self.at
    }

    fn entry(&mut self) -> Result<(i32, Range<usize>), BcfValueError> {
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

impl Iterator for BcfEntries<'_> {
    type Item = Result<(i32, Range<usize>), BcfValueError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.left = self.left.checked_sub(1)?;
        let entry = self.entry();
        if entry.is_err() {
            self.left = 0;
        }
        Some(entry)
    }
}

/// Index a BCF INFO or FORMAT section: every key must be a declared field of
/// its kind (`r[vcf_reader.undeclared]`), and the entries must fill it.
pub(crate) fn index_bcf(
    index: &mut FieldIndex,
    mut entries: BcfEntries<'_>,
    kind: KeyKind,
    file: &FileShared,
) -> Result<(), RecordError> {
    let (part, declared) = match kind {
        KeyKind::Format => (RecordPart::Format, &file.format_ids),
        KeyKind::Info | KeyKind::Filter => (RecordPart::Info, &file.info_ids),
    };
    for entry in &mut entries {
        let (key, span) = entry.map_err(|source| RecordError::Value { part, source })?;
        let id = u32::try_from(key)
            .ok()
            .filter(|&k| declared.contains(k))
            .and_then(DictId::new)
            .ok_or(RecordError::UnknownKey { index: key, kind })?;
        let (Ok(start), Ok(end)) = (u32::try_from(span.start), u32::try_from(span.end)) else {
            return Err(RecordError::TooLarge { size: span.end as u64, max: u64::from(u32::MAX) });
        };
        if !index.push(Entry { id, start, end }) {
            file.warn_duplicate(
                file.header.string_map().name(id.get() as usize).unwrap_or_default(),
            );
        }
    }
    let len = entries.bytes.len().saturating_sub(entries.end());
    if len > 0 {
        return Err(RecordError::TrailingBytes { part, len });
    }
    Ok(())
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

// r[impl vcf_fields.lazy]
/// A section's index, built on first access. Filling a record takes the
/// last one out and keeps its buffers for the next build; a failed build is
/// kept and returned on every access.
#[derive(Debug)]
pub(crate) struct LazyIndex<T = FieldIndex> {
    built: OnceLock<Result<T, RecordError>>,
    spare: Mutex<Option<T>>,
}

impl<T> Default for LazyIndex<T> {
    fn default() -> Self {
        Self { built: OnceLock::new(), spare: Mutex::new(None) }
    }
}

impl<T: Clone> Clone for LazyIndex<T> {
    fn clone(&self) -> Self {
        Self { built: self.built.clone(), spare: Mutex::default() }
    }
}

impl<T: Rebuild> LazyIndex<T> {
    /// Forget the index of the last record, keeping its buffers.
    pub(crate) fn reset(&mut self) {
        if let Some(Ok(index)) = self.built.take() {
            *self.spare.get_mut().unwrap_or_else(PoisonError::into_inner) = Some(index);
        }
    }

    /// The index, built by `build` on first access.
    pub(crate) fn get(
        &self,
        build: impl FnOnce(&mut T) -> Result<(), RecordError>,
    ) -> Result<&T, &RecordError> {
        self.built
            .get_or_init(|| {
                let spare = self.spare.lock().unwrap_or_else(PoisonError::into_inner).take();
                let mut index = spare.unwrap_or_default();
                index.restart();
                match build(&mut index) {
                    Ok(()) => Ok(index),
                    Err(e) => {
                        *self.spare.lock().unwrap_or_else(PoisonError::into_inner) = Some(index);
                        Err(e)
                    }
                }
            })
            .as_ref()
    }

    /// Whether the index has been built (or failed), for tests.
    #[cfg(test)]
    pub(crate) fn is_built(&self) -> bool {
        self.built.get().is_some()
    }
}

/// The span of a text entry or field, as `u32`s.
fn span32(r: Range<usize>) -> Result<(u32, u32), RecordError> {
    match (u32::try_from(r.start), u32::try_from(r.end)) {
        (Ok(s), Ok(e)) => Ok((s, e)),
        _ => Err(RecordError::TooLarge { size: r.end as u64, max: u64::from(u32::MAX) }),
    }
}

// r[impl vcf_reader.undeclared]
// r[impl vcf_values.duplicates]
/// Index a VCF text INFO column: each declared key's span is what follows
/// it (empty, or `=` and the value). Undeclared keys are ignored with one
/// warning per file; `.` and empty entries hold nothing.
pub(crate) fn index_text_info(
    index: &mut FieldIndex,
    info: &[u8],
    file: &FileShared,
) -> Result<(), RecordError> {
    if info == b"." {
        return Ok(());
    }
    let flow = for_each_info_entry(info, |entry| {
        let key = info.get(entry.key.clone()).unwrap_or_default();
        if key.is_empty() {
            return ControlFlow::Continue(());
        }
        let Some(id) = file.keys.get(key).filter(|id| file.info_ids.contains(id.get())) else {
            file.warn_undeclared(KeyKind::Info, key);
            return ControlFlow::Continue(());
        };
        let end = entry.value.map_or(entry.key.end, |v| v.end);
        match span32(entry.key.end..end) {
            Ok((start, end)) => {
                if !index.push(Entry { id, start, end }) {
                    file.warn_duplicate(std::str::from_utf8(key).unwrap_or_default());
                }
                ControlFlow::Continue(())
            }
            Err(e) => ControlFlow::Break(e),
        }
    });
    match flow {
        ControlFlow::Break(e) => Err(e),
        ControlFlow::Continue(()) => Ok(()),
    }
}

// r[impl vcf_fields.samples]
// r[impl vcf_reader.undeclared]
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
        let Some(id) = file.keys.get(key).filter(|id| file.format_ids.contains(id.get())) else {
            file.warn_undeclared(KeyKind::Format, key);
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
