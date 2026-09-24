//! Vectorized byte scans for text formats: the next newline in a block, and
//! the TAB-separated fields of a line.
//!
//! Both are written once over `fearless_simd` (`r[io.simd_portable]`) and
//! compare one native-width vector per iteration. `find_byte` gates on
//! `any_true` and builds a bitmask only for the vector that holds the hit;
//! `split_tabs` walks the set bits of each vector's TAB mask.

use fearless_simd::{Level, dispatch, prelude::*};
use fearless_simd_macros::simd;

/// Index of the first `needle` in `hay`, if any.
pub(crate) fn find_byte(hay: &[u8], needle: u8) -> Option<usize> {
    find_byte_at(Level::new(), hay, needle)
}

// r[impl io.simd_portable]
fn find_byte_at(level: Level, hay: &[u8], needle: u8) -> Option<usize> {
    dispatch!(level, simd => find_byte_simd(simd, hay, needle))
}

#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn find_byte_simd<S: Simd>(simd: S, hay: &[u8], needle: u8) -> Option<usize> {
    let mut chunks = hay.chunks_exact(S::u8s::LEN);
    let mut offset = 0usize;
    for chunk in &mut chunks {
        let hit = S::u8s::from_slice(simd, chunk).simd_eq(needle);
        if hit.any_true() {
            // A set lane makes the bitmask nonzero, so this is < LEN.
            let lane = hit.to_bitmask().trailing_zeros() as usize;
            return Some(offset.saturating_add(lane));
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    let rel = chunks.remainder().iter().position(|&b| b == needle)?;
    Some(offset.saturating_add(rel))
}

/// At most this many fields: exactly what `line.splitn(MAX_FIELDS, b'\t')`
/// yields, the last one holding the unsplit rest of the line.
pub(crate) const MAX_FIELDS: usize = 12;

/// Split `line` on TAB into at most [`MAX_FIELDS`] fields, the last of which
/// keeps any further TABs — `line.splitn(MAX_FIELDS, |&b| b == b'\t')`
/// without allocating. Returns the fields and how many of `fields` are set.
pub(crate) fn split_tabs(line: &[u8]) -> ([&[u8]; MAX_FIELDS], usize) {
    split_tabs_at(Level::new(), line)
}

// r[impl io.simd_portable]
fn split_tabs_at(level: Level, line: &[u8]) -> ([&[u8]; MAX_FIELDS], usize) {
    dispatch!(level, simd => split_tabs_simd(simd, line))
}

/// Records one field boundary: the field `start..tab` goes to `fields[n]`.
struct Splitter<'a> {
    line: &'a [u8],
    fields: [&'a [u8]; MAX_FIELDS],
    n: usize,
    start: usize,
}

impl<'a> Splitter<'a> {
    /// All fields but the last are closed; the last takes the rest verbatim.
    fn full(&self) -> bool {
        self.n >= MAX_FIELDS.saturating_sub(1)
    }

    fn tab_at(&mut self, tab: usize) {
        if let Some(slot) = self.fields.get_mut(self.n) {
            *slot = self.line.get(self.start..tab).unwrap_or_default();
        }
        self.n = self.n.saturating_add(1);
        self.start = tab.saturating_add(1);
    }

    fn finish(mut self) -> ([&'a [u8]; MAX_FIELDS], usize) {
        if let Some(slot) = self.fields.get_mut(self.n) {
            *slot = self.line.get(self.start..).unwrap_or_default();
        }
        (self.fields, self.n.saturating_add(1))
    }
}

#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(clippy::arithmetic_side_effects, reason = "`bits - 1` runs only on a nonzero `bits`")]
fn split_tabs_simd<S: Simd>(simd: S, line: &[u8]) -> ([&[u8]; MAX_FIELDS], usize) {
    let mut s = Splitter { line, fields: [&[]; MAX_FIELDS], n: 0, start: 0 };
    let mut chunks = line.chunks_exact(S::u8s::LEN);
    let mut offset = 0usize;
    for chunk in &mut chunks {
        if s.full() {
            return s.finish();
        }
        let tabs = S::u8s::from_slice(simd, chunk).simd_eq(b'\t');
        if tabs.any_true() {
            let mut bits = tabs.to_bitmask();
            while bits != 0 && !s.full() {
                s.tab_at(offset.saturating_add(bits.trailing_zeros() as usize));
                bits &= bits - 1;
            }
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    for (i, &b) in chunks.remainder().iter().enumerate() {
        if s.full() {
            break;
        }
        if b == b'\t' {
            s.tab_at(offset.saturating_add(i));
        }
    }
    s.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    /// Draws text dense in TABs and newlines, so both scans have many hits,
    /// runs of adjacent hits, and hits at vector boundaries.
    fn text(tc: &TestCase, max: usize) -> Vec<u8> {
        let alphabet = gs::sampled_from(vec![b'\t', b'\n', b'A', b'3', 0x00, 0xFF, b'\r']);
        tc.draw(gs::vecs(alphabet).max_size(max))
    }

    // r[verify sam.perf.text_parsing]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn find_byte_every_level_matches_position(tc: TestCase) {
        let hay = text(&tc, 300);
        for needle in [b'\n', b'\t', 0xFF] {
            let expected = hay.iter().position(|&b| b == needle);
            for level in crate::simd_levels::levels() {
                assert_eq!(find_byte_at(level, &hay, needle), expected, "{level:?}");
            }
        }
    }

    // r[verify sam.perf.text_parsing]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn split_tabs_every_level_matches_splitn(tc: TestCase) {
        let line = text(&tc, 400);
        let expected: Vec<&[u8]> = line.splitn(MAX_FIELDS, |&b| b == b'\t').collect();
        for level in crate::simd_levels::levels() {
            let (fields, n) = split_tabs_at(level, &line);
            assert_eq!(fields.get(..n), Some(expected.as_slice()), "{level:?}");
        }
    }

    // r[verify sam.perf.text_parsing]
    #[test]
    fn split_tabs_boundary_lengths() {
        // Every TAB position across two vectors, and more TABs than fields.
        for len in [0usize, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 129] {
            for tab in 0..len {
                let mut line = vec![b'x'; len];
                line[tab] = b'\t';
                let expected: Vec<&[u8]> = line.splitn(MAX_FIELDS, |&b| b == b'\t').collect();
                let (fields, n) = split_tabs(&line);
                assert_eq!(&fields[..n], expected.as_slice(), "len={len} tab={tab}");
            }
            let line = vec![b'\t'; len];
            let expected: Vec<&[u8]> = line.splitn(MAX_FIELDS, |&b| b == b'\t').collect();
            let (fields, n) = split_tabs(&line);
            assert_eq!(&fields[..n], expected.as_slice(), "all tabs, len={len}");
        }
    }
}
