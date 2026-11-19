//! A record section's index of present fields (`r[vcf_values.index+2]`):
//! [`SectionIndex`] for INFO (and the keys of FORMAT), [`FormatIndex`] for
//! FORMAT. A component holds one in a `RefCell` and resets it in `fill`;
//! a view calls [`SectionIndex::find`] / [`FormatIndex::finish`], which
//! scan on (`entries.rs`) only as far as the access needs
//! (`r[vcf_fields.lazy+4]`).

use super::entries::BcfEntries;
use super::error::{KeyKind, RecordError};
use super::file::FileShared;
use super::info::InfoSource;
use super::samples::FormatSource;
pub(crate) use crate::vcf::header::DictId;
use std::ops::Range;

/// One field present in a section: its id and where its value is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Entry {
    pub(crate) id: DictId,
    /// `start..end` in the section's bytes ([`Self::span`]): a BCF INFO
    /// typed value; a BCF FORMAT column (descriptor and every sample's
    /// elements); the VCF text after an INFO key, its `=` included
    /// (`InfoSource::value` strips it). A VCF text FORMAT entry holds the
    /// key's ordinal instead ([`Self::ordinal`]).
    pub(crate) start: u32,
    pub(crate) end: u32,
}

impl Entry {
    /// A VCF text FORMAT key's entry: its position among the FORMAT keys,
    /// where other entries hold a byte span.
    pub(crate) fn ordinal_of(id: DictId, ordinal: u32) -> Self {
        Self { id, start: ordinal, end: ordinal.saturating_add(1) }
    }

    /// The byte span of the value, in the section's bytes.
    pub(crate) fn span(self) -> Range<usize> {
        self.start as usize..self.end as usize
    }

    /// A VCF text FORMAT key's position among the FORMAT keys
    /// ([`Self::ordinal_of`]).
    pub(crate) fn ordinal(self) -> usize {
        self.start as usize
    }
}

/// How far a section's index got.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Progress {
    /// Indexed up to byte `at`, `walked` entries in (BCF counts its entries;
    /// VCF text leaves `walked` at 0). `u32`, as an [`Entry`]'s offsets are:
    /// the state stays one 8-byte word.
    Scanning { at: u32, walked: u32 },
    /// The whole section is indexed.
    Done,
    /// The scan met malformed data; what it indexed before stays readable.
    /// Boxed: the error is rare, and a small `Progress` keeps the per-record
    /// reset to a few stores.
    Failed(Box<RecordError>),
}

// r[impl vcf_values.index+2]
// r[impl vcf_fields.lazy+4]
/// The fields of one section in record order, a table with one `u32` per
/// dictionary id, `epoch << 16 | slot`, and how far the scan that fills
/// them got. A lookup is one table read, an epoch compare and one list
/// read, and an id whose stamp is from another epoch is absent. Starting a
/// record advances the epoch instead of clearing the table; on wrap the
/// table is zeroed. The table grows to the largest id seen, never to the
/// dictionary's size.
#[derive(Debug, Clone)]
pub(crate) struct SectionIndex {
    /// The fields found so far, in record order, first occurrences only.
    entries: Vec<Entry>,
    /// By `DictId`: `epoch << 16 | slot in entries`.
    table: Vec<u32>,
    /// The current record's stamp; never `0`, so a zeroed stamp is never current.
    epoch: u16,
    /// Where the scan stopped; the scans in `entries.rs` advance it.
    pub(crate) progress: Progress,
}

impl Default for SectionIndex {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
            table: Vec::new(),
            epoch: 1,
            progress: Progress::Scanning { at: 0, walked: 0 },
        }
    }
}

impl SectionIndex {
    /// Forget the last section's fields, keeping the buffers, to index a new one.
    #[inline]
    pub(crate) fn reset(&mut self) {
        self.entries.clear();
        self.progress = Progress::Scanning { at: 0, walked: 0 };
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.wrap();
        }
    }

    #[cold]
    fn wrap(&mut self) {
        self.table.fill(0);
        self.epoch = 1;
    }

    #[cold]
    fn grow(&mut self, id: usize) {
        self.table.resize(id.saturating_add(1), 0);
    }

    /// Add a field; `false`, leaving the index unchanged, when its id is
    /// already present (the first occurrence wins, `r[vcf_values.duplicates+2]`).
    #[inline]
    pub(crate) fn push(&mut self, entry: Entry) -> bool {
        let id = entry.id.get() as usize;
        if id >= self.table.len() {
            self.grow(id);
        }
        let Some(stamp) = self.table.get_mut(id) else { return false };
        if *stamp >> 16 == u32::from(self.epoch) {
            return false;
        }
        // At most one entry per id, and ids are 16 bits: the slot fits.
        let Ok(slot) = u16::try_from(self.entries.len()) else { return false };
        *stamp = (u32::from(self.epoch) << 16) | u32::from(slot);
        self.entries.push(entry);
        true
    }

    /// The entry of `id`, if the index has reached it.
    #[inline]
    pub(crate) fn get(&self, id: DictId) -> Option<Entry> {
        let stamp = *self.table.get(id.get() as usize)?;
        if stamp >> 16 != u32::from(self.epoch) {
            return None;
        }
        self.entries.get((stamp & 0xFFFF) as usize).copied()
    }

    /// The fields indexed, in record order, duplicates left out.
    pub(crate) fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// The first entry of `id` in the INFO section `info`, indexing no
    /// further than that (`r[vcf_values.index+2]`). Malformed data fails the
    /// lookup only if the scan meets it before `id`. The error is the
    /// index's, borrowed: the hot path moves nothing larger than an entry.
    #[inline]
    pub(crate) fn find(
        &mut self,
        id: DictId,
        info: InfoSource<'_>,
        file: &FileShared,
    ) -> Result<Option<Entry>, &RecordError> {
        if let Some(entry) = self.get(id) {
            return Ok(Some(entry));
        }
        if let Progress::Scanning { at, walked } = self.progress
            && let Some(entry) = self.scan(info, at, walked, Some(id), file)
        {
            return Ok(Some(entry));
        }
        self.outcome().map(|()| None)
    }

    /// Index the rest of the INFO section `info`.
    pub(crate) fn finish(
        &mut self,
        info: InfoSource<'_>,
        file: &FileShared,
    ) -> Result<(), &RecordError> {
        if let Progress::Scanning { at, walked } = self.progress {
            self.scan(info, at, walked, None, file);
        }
        self.outcome()
    }

    /// A settled scan's outcome: its error, if it failed.
    #[inline]
    fn outcome(&self) -> Result<(), &RecordError> {
        match &self.progress {
            Progress::Failed(e) => Err(e),
            Progress::Scanning { .. } | Progress::Done => Ok(()),
        }
    }

    /// Scan on from `at`, stopping after the first entry of `until`.
    #[cold]
    #[inline(never)]
    fn scan(
        &mut self,
        info: InfoSource<'_>,
        at: u32,
        walked: u32,
        until: Option<DictId>,
        file: &FileShared,
    ) -> Option<Entry> {
        let found = match info {
            InfoSource::Bcf { bytes, n } => {
                let entries = BcfEntries::new(bytes, n, None, at, walked);
                self.scan_bcf(entries, until, KeyKind::Info, file)
            }
            InfoSource::Text(bytes) => self.scan_text_info(bytes, at, until, file),
        };
        self.settle(found)
    }

    /// Record a scan's end: the entry it stopped at (still scanning), or the
    /// section's end (done) or its error (failed), which [`Self::outcome`]
    /// then returns.
    pub(crate) fn settle(&mut self, found: Result<Option<Entry>, RecordError>) -> Option<Entry> {
        match found {
            Ok(Some(entry)) => return Some(entry),
            Ok(None) => self.progress = Progress::Done,
            Err(e) => self.progress = Progress::Failed(Box::new(e)),
        }
        None
    }
}

/// The FORMAT fields of a record: which are present, and for VCF text
/// where each sample's value of each lies. Scanned whole on first access:
/// only INFO is scanned incrementally (`r[vcf_values.index+2]`).
#[derive(Debug, Clone, Default)]
pub(crate) struct FormatIndex {
    /// The fields present and the scan's progress. A BCF entry spans its
    /// column (every sample's values); a VCF text entry holds the key's
    /// position among the FORMAT keys ([`Entry::ordinal`]).
    pub(crate) keys: SectionIndex,
    /// VCF text only: sample `s`'s value of the `k`th key at
    /// `s * n_keys + k`, empty where the sample dropped it.
    pub(crate) cells: Vec<Range<u32>>,
    /// VCF text only: how many keys the FORMAT column names.
    pub(crate) n_keys: usize,
}

impl FormatIndex {
    /// Forget the last record's FORMAT fields, keeping the buffers.
    #[inline]
    pub(crate) fn reset(&mut self) {
        self.keys.reset();
        self.cells.clear();
        self.n_keys = 0;
    }

    /// Fail every access with `e`: the FORMAT data could not be copied.
    pub(crate) fn fail(&mut self, e: RecordError) {
        self.keys.progress = Progress::Failed(Box::new(e));
    }

    /// Where sample `sample`'s value of the `k`th key lies in the sample
    /// columns (VCF text).
    pub(crate) fn cell(&self, sample: usize, k: usize) -> Option<Range<usize>> {
        let r = self.cells.get(sample.checked_mul(self.n_keys)?.checked_add(k)?)?;
        Some(r.start as usize..r.end as usize)
    }

    /// Index the whole FORMAT data `source` (once; later calls return the
    /// first one's outcome).
    #[inline]
    pub(crate) fn finish(
        &mut self,
        source: FormatSource<'_>,
        file: &FileShared,
    ) -> Result<(), &RecordError> {
        if let Progress::Scanning { .. } = self.keys.progress {
            self.scan(source, file);
        }
        self.keys.outcome()
    }

    #[cold]
    #[inline(never)]
    fn scan(&mut self, source: FormatSource<'_>, file: &FileShared) {
        let found = match source {
            FormatSource::Bcf { indiv, n_fmt, n_samples } => {
                let entries = BcfEntries::new(indiv, n_fmt, Some(n_samples), 0, 0);
                self.keys.scan_bcf(entries, None, KeyKind::Format, file)
            }
            FormatSource::Text { keys, samples, n_samples } => {
                self.scan_text(keys, samples, n_samples, file).map(|()| None)
            }
        };
        self.keys.settle(found);
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        clippy::cast_possible_truncation,
        reason = "test indices below 30"
    )]

    use super::*;
    use hegel::prelude::*;
    use std::collections::HashMap;

    // r[verify vcf_values.index+2]
    // r[verify vcf_values.duplicates+2]
    /// Against a map that keeps each id's first entry, over many records on
    /// one index, starting next to the epoch's wrap.
    #[hegel::test]
    fn the_index_finds_each_ids_first_entry(tc: TestCase) {
        let mut index = SectionIndex {
            epoch: tc.draw(gs::sampled_from(vec![1u16, u16::MAX - 2])),
            ..SectionIndex::default()
        };
        let records = tc.draw(
            gs::vecs(gs::vecs(gs::integers::<u16>().max_value(300)).max_size(30)).max_size(8),
        );
        for ids in records {
            index.reset();
            let mut first = HashMap::new();
            for (i, &id) in ids.iter().enumerate() {
                let id = DictId::new(u32::from(id)).unwrap();
                let entry = Entry { id, start: i as u32, end: i as u32 + 1 };
                assert_eq!(index.push(entry), !first.contains_key(&id));
                first.entry(id).or_insert(entry);
            }
            for id in 0..=300u32 {
                let id = DictId::new(id).unwrap();
                assert_eq!(index.get(id), first.get(&id).copied(), "id {id:?}");
            }
            let mut expected: Vec<_> = first.into_values().collect();
            expected.sort_by_key(|e| e.start);
            assert_eq!(index.entries(), expected);
        }
    }
}
