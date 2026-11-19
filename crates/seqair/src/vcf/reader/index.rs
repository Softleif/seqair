//! A record section's index of present fields (`r[vcf_values.index+2]`).

/// A string-dictionary index: a header declares at most
/// [`MAX_DICT_ENTRIES`](crate::vcf::header::MAX_DICT_ENTRIES) = 2^16 ids.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct DictId(u16);

impl DictId {
    /// `None` for an index no header can declare.
    pub(crate) fn new(index: u32) -> Option<Self> {
        u16::try_from(index).ok().map(Self)
    }

    pub(crate) fn get(self) -> u32 {
        u32::from(self.0)
    }
}

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

// r[impl vcf_values.index+2]
/// The fields of one section in record order, and a table with one `u32` per
/// dictionary id, `epoch << 16 | slot`: a lookup is one table read, an epoch
/// compare and one list read, and an id whose stamp is from another epoch is
/// absent. Starting a record advances the epoch instead of clearing the
/// table; on wrap the table is zeroed. The table grows to the largest id
/// seen, never to the dictionary's size.
#[derive(Debug, Clone, Default)]
pub(crate) struct FieldIndex {
    entries: Vec<Entry>,
    table: Vec<u32>,
    /// `0` only before the first record: every stamp in the table is from
    /// an epoch `1..=u16::MAX`.
    epoch: u16,
}

impl FieldIndex {
    /// Start indexing a new section.
    #[inline]
    pub(crate) fn begin(&mut self) {
        self.entries.clear();
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
        let id = entry.id.0 as usize;
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

    #[inline]
    pub(crate) fn get(&self, id: DictId) -> Option<Entry> {
        let stamp = *self.table.get(id.0 as usize)?;
        if stamp >> 16 != u32::from(self.epoch) {
            return None;
        }
        self.entries.get((stamp & 0xFFFF) as usize).copied()
    }

    /// The fields present, in record order, duplicates left out.
    pub(crate) fn entries(&self) -> &[Entry] {
        &self.entries
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
        let mut index = FieldIndex {
            epoch: tc.draw(gs::sampled_from(vec![0u16, u16::MAX - 2])),
            ..FieldIndex::default()
        };
        let records = tc.draw(
            gs::vecs(gs::vecs(gs::integers::<u16>().max_value(300)).max_size(30)).max_size(8),
        );
        for ids in records {
            index.begin();
            let mut first = HashMap::new();
            for (i, &id) in ids.iter().enumerate() {
                let entry = Entry { id: DictId(id), start: i as u32, end: i as u32 + 1 };
                assert_eq!(index.push(entry), !first.contains_key(&id));
                first.entry(id).or_insert(entry);
            }
            for id in 0..=300u16 {
                assert_eq!(index.get(DictId(id)), first.get(&id).copied(), "id {id}");
            }
            let mut expected: Vec<_> = first.into_values().collect();
            expected.sort_by_key(|e| e.start);
            assert_eq!(index.entries(), expected);
        }
    }
}
