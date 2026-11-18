//! VCF text keys to dictionary ids, through a table built once per file
//! (`r[vcf_values.text_resolution]`).

use super::index::DictId;
use std::num::NonZeroU32;

/// One key: its first 8 bytes (zero-padded), its length, where its name
/// starts in the flat name buffer, and its id.
#[derive(Debug, Clone, Copy)]
struct Slot {
    prefix: u64,
    len: NonZeroU32,
    name_at: u32,
    id: DictId,
}

// r[impl vcf_values.text_resolution]
/// Every FILTER/INFO/FORMAT id of a header, by name: open addressing at
/// most half full, linear probing. A key of up to 8 bytes is told apart by
/// its prefix and length alone; a longer one is checked against the flat
/// name buffer. The hash is one 8-byte load, plus the last 8 bytes of a
/// longer key.
#[derive(Debug, Clone)]
pub(crate) struct TextKeys {
    slots: Vec<Option<Slot>>,
    names: Vec<u8>,
    /// `64 - log2(slots.len())`: the hash's top bits pick the slot.
    shift: u32,
}

/// `key`'s first 8 bytes and, for a longer key, its last 8 (else 0).
#[inline]
fn words(key: &[u8]) -> (u64, u64) {
    match (key.first_chunk::<8>(), key.last_chunk::<8>()) {
        (Some(first), Some(last)) if key.len() > 8 => {
            (u64::from_le_bytes(*first), u64::from_le_bytes(*last))
        }
        _ => {
            let mut padded = [0u8; 8];
            let n = key.len().min(8);
            if let (Some(dst), Some(src)) = (padded.get_mut(..n), key.get(..n)) {
                dst.copy_from_slice(src);
            }
            (u64::from_le_bytes(padded), 0)
        }
    }
}

#[inline]
fn hash(first: u64, last: u64, len: usize) -> u64 {
    (first ^ last.rotate_left(29) ^ len as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

impl TextKeys {
    /// The table of `names` (dictionary index, name). Empty names and
    /// indices past 16 bits, which no header holds, are left out.
    pub(crate) fn new<'a>(names: impl IntoIterator<Item = (usize, &'a str)>) -> Self {
        let names: Vec<(DictId, &str)> = names
            .into_iter()
            .filter_map(|(i, n)| Some((DictId::new(u32::try_from(i).ok()?)?, n)))
            .filter(|(_, n)| !n.is_empty())
            .collect();
        let capacity = names.len().saturating_mul(2).max(8).next_power_of_two();
        let mut table = Self {
            slots: vec![None; capacity],
            names: Vec::new(),
            shift: 64u32.saturating_sub(capacity.trailing_zeros()),
        };
        for (id, name) in names {
            table.insert(id, name.as_bytes());
        }
        table
    }

    fn insert(&mut self, id: DictId, name: &[u8]) {
        let (Ok(name_at), Some(len)) = (
            u32::try_from(self.names.len()),
            u32::try_from(name.len()).ok().and_then(NonZeroU32::new),
        ) else {
            return;
        };
        let (first, last) = words(name);
        let mask = self.slots.len().saturating_sub(1);
        let mut i = usize::try_from(hash(first, last, name.len()) >> self.shift).unwrap_or(0);
        // At most half full: a free slot is always found.
        for _ in 0..self.slots.len() {
            match self.slots.get(i) {
                Some(None) => break,
                // A header's names are distinct; keep the first if not.
                Some(Some(s)) if self.same(s, first, name) => return,
                _ => i = i.wrapping_add(1) & mask,
            }
        }
        if let Some(slot @ None) = self.slots.get_mut(i) {
            *slot = Some(Slot { prefix: first, len, name_at, id });
            self.names.extend_from_slice(name);
        }
    }

    #[inline]
    fn same(&self, slot: &Slot, first: u64, key: &[u8]) -> bool {
        if slot.prefix != first || slot.len.get() as usize != key.len() {
            return false;
        }
        key.len() <= 8 || {
            let start = slot.name_at as usize;
            self.names.get(start..start.saturating_add(key.len())) == Some(key)
        }
    }

    /// The id of `key`, if the header declares it.
    #[inline]
    pub(crate) fn get(&self, key: &[u8]) -> Option<DictId> {
        let (first, last) = words(key);
        let mask = self.slots.len().saturating_sub(1);
        let mut i = usize::try_from(hash(first, last, key.len()) >> self.shift).unwrap_or(0);
        loop {
            let slot = self.slots.get(i)?.as_ref()?;
            if self.same(slot, first, key) {
                return Some(slot.id);
            }
            i = i.wrapping_add(1) & mask;
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::cast_possible_truncation,
        clippy::arithmetic_side_effects,
        reason = "test indices below 60"
    )]

    use super::*;
    use hegel::prelude::*;
    use std::collections::HashMap;

    /// Names that share long prefixes, differ only in length or in a late
    /// byte, and include NUL, so the prefix compare alone cannot tell them
    /// apart.
    fn name(tc: &TestCase) -> String {
        let stem = tc.draw(gs::sampled_from(vec!["", "A", "AAAAAAA", "AAAAAAAA", "AAAAAAAAB"]));
        let n = tc.draw(gs::integers::<usize>().max_value(12));
        let tail: String =
            (0..n).map(|_| tc.draw(gs::sampled_from(vec!['A', 'B', '\0', '_']))).collect();
        format!("{stem}{tail}")
    }

    // r[verify vcf_values.text_resolution]
    /// Against a `HashMap`: every name finds its id, every other key none.
    #[hegel::test]
    fn the_table_finds_exactly_the_declared_names(tc: TestCase) {
        let mut oracle: HashMap<String, usize> = HashMap::new();
        let n = tc.draw(gs::integers::<usize>().max_value(60));
        for i in 0..n {
            oracle.entry(name(&tc)).or_insert(i);
        }
        let table = TextKeys::new(oracle.iter().map(|(n, &i)| (i, n.as_str())));
        assert!(table.slots.iter().filter(|s| s.is_some()).count() * 2 <= table.slots.len());
        for (n, &i) in &oracle {
            let expected = (!n.is_empty()).then(|| DictId::new(i as u32)).flatten();
            assert_eq!(table.get(n.as_bytes()), expected, "{n:?}");
        }
        for _ in 0..20 {
            let key = name(&tc);
            if !oracle.contains_key(&key) {
                assert_eq!(table.get(key.as_bytes()), None, "{key:?}");
            }
        }
    }

    #[test]
    fn an_empty_table_finds_nothing() {
        let table = TextKeys::new([]);
        assert_eq!(table.get(b"DP"), None);
        assert_eq!(table.get(b""), None);
    }
}
