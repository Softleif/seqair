//! A growable bitset of string-dictionary ids, for per-record "already seen"
//! tracking.

// r[impl vcf_header.field_set]
/// Set of BCF string-dictionary indices.
///
/// Starts empty without allocating and grows to fit the largest id inserted,
/// so it needs no header to be constructed. [`clear`](Self::clear) keeps the
/// words, so a set reused across records allocates only while it grows.
#[derive(Debug, Default, Clone)]
pub(crate) struct FieldSet {
    words: Box<[u64]>,
}

/// The word index and bit mask of `id`.
fn locate(id: u32) -> (usize, u64) {
    ((id >> 6) as usize, 1u64 << (id & 63))
}

impl FieldSet {
    /// Add `id`; `true` when it was not yet present.
    pub(crate) fn insert(&mut self, id: u32) -> bool {
        let (word, mask) = locate(id);
        if word >= self.words.len() {
            let mut words = std::mem::take(&mut self.words).into_vec();
            words.resize(word.saturating_add(1), 0);
            self.words = words.into_boxed_slice();
        }
        #[allow(clippy::indexing_slicing, reason = "the words were just grown past `word`")]
        let w = &mut self.words[word];
        let new = *w & mask == 0;
        *w |= mask;
        new
    }

    #[allow(dead_code, reason = "for the VCF/BCF reader's warn-once sets")]
    pub(crate) fn contains(&self, id: u32) -> bool {
        let (word, mask) = locate(id);
        self.words.get(word).is_some_and(|w| w & mask != 0)
    }

    /// Remove every id, keeping the allocation.
    pub(crate) fn clear(&mut self) {
        self.words.fill(0);
    }

    /// The ids in the set, ascending.
    #[allow(dead_code, reason = "for the VCF/BCF reader's warn-once sets")]
    pub(crate) fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        (0u32..).zip(self.words.iter()).flat_map(|(word, &bits)| {
            std::iter::successors((bits != 0).then_some(bits), |b| {
                let rest = b & b.wrapping_sub(1);
                (rest != 0).then_some(rest)
            })
            .map(move |b| {
                #[allow(
                    clippy::arithmetic_side_effects,
                    reason = "word < 2^26 because ids are u32, so word * 64 + 63 fits"
                )]
                let id = word * 64 + b.trailing_zeros();
                id
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::prelude::*;
    use std::collections::BTreeSet;

    // r[verify vcf_header.field_set]
    /// Every operation agrees with a `BTreeSet`, ids drawn across several
    /// word boundaries up to the dictionary limit so the set grows
    /// mid-sequence, and the iterator yields exactly the oracle's ids in order.
    #[hegel::test]
    fn field_set_matches_a_btreeset(tc: TestCase) {
        let mut set = FieldSet::default();
        let mut oracle = BTreeSet::new();
        let n_ops = tc.draw(gs::integers::<usize>().max_value(64));
        for _ in 0..n_ops {
            // Mostly small ids (dense words, repeat hits), sometimes far ones (growth).
            let id = tc.draw(hegel::one_of!(
                gs::integers::<u32>().max_value(200),
                gs::integers::<u32>().max_value(65_535),
                gs::sampled_from(vec![0u32, 63, 64, 127, 128, 65_535]),
            ));
            match tc.draw(gs::integers::<u8>().max_value(4)) {
                0 | 1 => assert_eq!(set.insert(id), oracle.insert(id), "insert {id}"),
                2 | 3 => assert_eq!(set.contains(id), oracle.contains(&id), "contains {id}"),
                _ => {
                    let words = set.words.len();
                    set.clear();
                    oracle.clear();
                    assert_eq!(set.words.len(), words, "clear keeps the allocation");
                }
            }
            assert!(set.iter().eq(oracle.iter().copied()), "iter disagrees with {oracle:?}");
        }
    }
}
