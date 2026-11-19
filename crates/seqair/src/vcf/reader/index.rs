//! A record section's index of present fields (`r[vcf_values.index+2]`).

use super::entries::BcfEntries;
use super::error::{KeyKind, RecordError};
use super::file::FileShared;
use super::info::InfoSource;
pub(crate) use crate::vcf::header::DictId;

/// One field present in a section: its id and the byte span of its value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Entry {
    pub(crate) id: DictId,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

impl Entry {
    pub(crate) fn span(self) -> std::ops::Range<usize> {
        self.start as usize..self.end as usize
    }
}

/// How far a section's index got.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Progress {
    /// Indexed up to byte `at`, `walked` entries in (BCF counts its entries;
    /// VCF text leaves `walked` at 0).
    Scanning { at: usize, walked: usize },
    /// The whole section is indexed.
    Done,
    /// The scan met malformed data; what it indexed before stays readable.
    Failed(RecordError),
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
    entries: Vec<Entry>,
    table: Vec<u32>,
    /// Never `0`: a zeroed stamp is never current.
    epoch: u16,
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
    /// lookup only if the scan meets it before `id`.
    #[inline]
    pub(crate) fn find(
        &mut self,
        id: DictId,
        info: InfoSource<'_>,
        file: &FileShared,
    ) -> Result<Option<Entry>, RecordError> {
        match self.get(id) {
            Some(entry) => Ok(Some(entry)),
            None => self.resume(info, Some(id), file),
        }
    }

    /// Index the rest of the INFO section `info`.
    pub(crate) fn finish(
        &mut self,
        info: InfoSource<'_>,
        file: &FileShared,
    ) -> Result<(), RecordError> {
        self.resume(info, None, file).map(drop)
    }

    #[inline]
    fn resume(
        &mut self,
        info: InfoSource<'_>,
        until: Option<DictId>,
        file: &FileShared,
    ) -> Result<Option<Entry>, RecordError> {
        match self.progress {
            Progress::Done => Ok(None),
            Progress::Failed(ref e) => Err(e.clone()),
            Progress::Scanning { at, walked } => self.scan(info, at, walked, until, file),
        }
    }

    /// Scan on from `at`, stopping after the first entry of `until`.
    #[cold]
    #[inline(never)]
    fn scan(
        &mut self,
        info: InfoSource<'_>,
        at: usize,
        walked: usize,
        until: Option<DictId>,
        file: &FileShared,
    ) -> Result<Option<Entry>, RecordError> {
        let found = match info {
            InfoSource::Bcf { bytes, n } => {
                let entries = BcfEntries::new(bytes, n, None, at, walked);
                self.scan_bcf(entries, until, KeyKind::Info, file)
            }
            InfoSource::Text(bytes) => self.scan_text_info(bytes, at, until, file),
        };
        self.settle(found)
    }

    /// Record a scan's end: found (still scanning), done or failed.
    pub(crate) fn settle(
        &mut self,
        found: Result<Option<Entry>, RecordError>,
    ) -> Result<Option<Entry>, RecordError> {
        match found {
            Ok(Some(entry)) => Ok(Some(entry)),
            Ok(None) => {
                self.progress = Progress::Done;
                Ok(None)
            }
            Err(e) => {
                self.progress = Progress::Failed(e.clone());
                Err(e)
            }
        }
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
