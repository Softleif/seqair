//! Walking a BCF INFO or FORMAT section entry by entry, and the lazily built
//! index of one (`r[vcf_fields.lazy]`, `r[vcf_values.index]`).

use super::bcf_value::{BcfValueError, FormatColumn, TypedValue, read_typed_int};
use super::error::{KeyKind, RecordError, RecordPart};
use super::file::FileShared;
use super::index::{DictId, Entry, FieldIndex};
use std::ops::Range;
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

/// The span of `id`'s first entry, walking the entries.
pub(crate) fn scan(
    mut entries: BcfEntries<'_>,
    id: DictId,
) -> Result<Option<Range<usize>>, BcfValueError> {
    let wanted = i64::from(id.get());
    entries
        .find_map(|e| match e {
            Ok((key, span)) => (i64::from(key) == wanted).then_some(Ok(span)),
            Err(e) => Some(Err(e)),
        })
        .transpose()
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

// r[impl vcf_fields.lazy]
/// A section's index, built on first access. Filling a record takes the
/// last one out and keeps its buffers for the next build; a failed build is
/// kept and returned on every access.
#[derive(Debug, Default)]
pub(crate) struct LazyIndex {
    built: OnceLock<Result<FieldIndex, RecordError>>,
    spare: Mutex<Option<FieldIndex>>,
}

impl Clone for LazyIndex {
    fn clone(&self) -> Self {
        Self { built: self.built.clone(), spare: Mutex::default() }
    }
}

impl LazyIndex {
    /// Forget the index of the last record, keeping its buffers.
    pub(crate) fn reset(&mut self) {
        if let Some(Ok(index)) = self.built.take() {
            *self.spare.get_mut().unwrap_or_else(PoisonError::into_inner) = Some(index);
        }
    }

    /// The index, built by `build` on first access.
    pub(crate) fn get(
        &self,
        build: impl FnOnce(&mut FieldIndex) -> Result<(), RecordError>,
    ) -> Result<&FieldIndex, &RecordError> {
        self.built
            .get_or_init(|| {
                let spare = self.spare.lock().unwrap_or_else(PoisonError::into_inner).take();
                let mut index = spare.unwrap_or_default();
                index.begin();
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
