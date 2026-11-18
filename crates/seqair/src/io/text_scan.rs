//! Vectorized byte scans: the next occurrence of a byte (a newline in a SAM
//! block, a CRAM `BYTE_ARRAY_STOP` terminator), the TAB-separated fields of a
//! line, every occurrence of a separator, and the entries of a VCF INFO column.
//!
//! All are written once over `fearless_simd` (`r[io.simd_portable]`) and
//! compare one native-width vector per iteration. `find_byte` gates on
//! `any_true` and builds a bitmask only for the vector that holds the hit;
//! the others walk the set bits of each vector's separator mask.

use core::ops::{ControlFlow, Range};
use fearless_simd::{Level, dispatch, prelude::*, u8x64};
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
fn find_byte_simd<S: Simd>(simd: S, hay: &[u8], needle: u8) -> Option<usize> {
    find_byte_in(simd, hay, needle)
}

/// [`find_byte`] inside a kernel already running at level `S`: a scan per
/// field (each aux `Z` string, say) then costs no dispatch call.
#[inline(always)]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
pub(crate) fn find_byte_in<S: Simd>(simd: S, hay: &[u8], needle: u8) -> Option<usize> {
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

/// Every TAB offset of `line`, into `tabs` (cleared first), and the offset of
/// the first NUL, if any — what SAM encoding needs from a whole line, in one
/// pass instead of a NUL scan, a field split and a scan per optional field.
pub(crate) fn tabs_and_nul(line: &[u8], tabs: &mut Vec<usize>) -> Option<usize> {
    tabs_and_nul_at(Level::new(), line, tabs)
}

// r[impl io.simd_portable]
fn tabs_and_nul_at(level: Level, line: &[u8], tabs: &mut Vec<usize>) -> Option<usize> {
    dispatch!(level, simd => tabs_and_nul_simd(simd, line, tabs))
}

/// Every line goes through here, so it takes 64 bytes at a time rather than
/// `S::u8s`: one 64-bit TAB mask per block (cheap to build even on NEON), and
/// NUL folded into a running byte minimum checked once at the end.
#[simd]
#[allow(clippy::arithmetic_side_effects, reason = "`64 - tail.len()` on a nonempty tail")]
fn tabs_and_nul_simd<S: Simd>(simd: S, line: &[u8], tabs: &mut Vec<usize>) -> Option<usize> {
    tabs.clear();
    let mut min = u8x64::splat(simd, u8::MAX);
    let (blocks, tail) = line.as_chunks::<64>();
    let mut offset = 0usize;
    for chunk in blocks {
        scan_block(u8x64::from_slice(simd, chunk), offset, 0, &mut min, tabs);
        offset = offset.saturating_add(64);
    }
    if !tail.is_empty() {
        if let Some(last) = line.last_chunk::<64>() {
            // The last 64 bytes, overlapping the blocks before.
            scan_block(u8x64::from_slice(simd, last), offset, 64 - tail.len(), &mut min, tabs);
        } else {
            // Shorter than a block: padding that is neither TAB nor NUL.
            let mut last = [u8::MAX; 64];
            last.get_mut(..tail.len()).unwrap_or_default().copy_from_slice(tail);
            scan_block(u8x64::from_slice(simd, &last), offset, 0, &mut min, tabs);
        }
    }
    if min.simd_eq(0).any_true() { find_byte_in(simd, line, 0) } else { None }
}

/// [`tabs_and_nul_simd`]'s block `v` at `offset`, whose first `skip` (< 64)
/// bytes were scanned already.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "`bits - 1` runs only on a nonzero `bits`")]
fn scan_block<S: Simd>(
    v: u8x64<S>,
    offset: usize,
    skip: usize,
    min: &mut u8x64<S>,
    tabs: &mut Vec<usize>,
) {
    *min = min.min(v);
    let mut bits = v.simd_eq(b'\t').to_bitmask() >> skip;
    while bits != 0 {
        tabs.push(offset.saturating_add(bits.trailing_zeros() as usize));
        bits &= bits - 1;
    }
}

/// Calls `f` with the offset of every `needle` in `hay`, in order, until `f`
/// breaks — the sample columns of a VCF line (`\t`), the fields of a sample
/// (`:`), the elements of a value (`,`).
// r[impl vcf_reader.text_scan]
pub(crate) fn for_each_byte<B, F>(hay: &[u8], needle: u8, f: F) -> ControlFlow<B>
where
    F: FnMut(usize) -> ControlFlow<B>,
{
    for_each_byte_at(Level::new(), hay, needle, f)
}

// r[impl io.simd_portable]
fn for_each_byte_at<B, F>(level: Level, hay: &[u8], needle: u8, f: F) -> ControlFlow<B>
where
    F: FnMut(usize) -> ControlFlow<B>,
{
    dispatch!(level, simd => for_each_byte_simd(simd, hay, needle, f))
}

#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(clippy::arithmetic_side_effects, reason = "`bits - 1` runs only on a nonzero `bits`")]
fn for_each_byte_simd<S: Simd, B, F>(simd: S, hay: &[u8], needle: u8, mut f: F) -> ControlFlow<B>
where
    F: FnMut(usize) -> ControlFlow<B>,
{
    let mut chunks = hay.chunks_exact(S::u8s::LEN);
    let mut offset = 0usize;
    for chunk in &mut chunks {
        let hits = S::u8s::from_slice(simd, chunk).simd_eq(needle);
        if hits.any_true() {
            let mut bits = hits.to_bitmask();
            while bits != 0 {
                f(offset.saturating_add(bits.trailing_zeros() as usize))?;
                bits &= bits - 1;
            }
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    for (i, &b) in chunks.remainder().iter().enumerate() {
        if b == needle {
            f(offset.saturating_add(i))?;
        }
    }
    ControlFlow::Continue(())
}

/// Calls `f` with the offset of every `outer` or `inner` byte in `hay`, in
/// order, and whether it is `outer`, until `f` breaks: the sample columns
/// (`\t`) and their fields (`:`) of a VCF line in one pass.
// r[impl vcf_reader.text_scan]
pub(crate) fn for_each_of_two<B, F>(hay: &[u8], outer: u8, inner: u8, f: F) -> ControlFlow<B>
where
    F: FnMut(usize, bool) -> ControlFlow<B>,
{
    for_each_of_two_at(Level::new(), hay, outer, inner, f)
}

// r[impl io.simd_portable]
fn for_each_of_two_at<B, F>(level: Level, hay: &[u8], outer: u8, inner: u8, f: F) -> ControlFlow<B>
where
    F: FnMut(usize, bool) -> ControlFlow<B>,
{
    dispatch!(level, simd => for_each_of_two_simd(simd, hay, outer, inner, f))
}

/// One mask of both bytes per vector, walked in order; the `outer` mask
/// alone says which each set bit is.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(clippy::arithmetic_side_effects, reason = "`bits - 1` runs only on a nonzero `bits`")]
fn for_each_of_two_simd<S: Simd, B, F>(
    simd: S,
    hay: &[u8],
    outer: u8,
    inner: u8,
    mut f: F,
) -> ControlFlow<B>
where
    F: FnMut(usize, bool) -> ControlFlow<B>,
{
    let mut chunks = hay.chunks_exact(S::u8s::LEN);
    let mut offset = 0usize;
    for chunk in &mut chunks {
        let v = S::u8s::from_slice(simd, chunk);
        let outers = v.simd_eq(outer);
        let both = outers | v.simd_eq(inner);
        if both.any_true() {
            let outer_bits = outers.to_bitmask();
            let mut bits = both.to_bitmask();
            while bits != 0 {
                let lane = bits.trailing_zeros();
                f(offset.saturating_add(lane as usize), (outer_bits >> lane) & 1 == 1)?;
                bits &= bits - 1;
            }
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    for (i, &b) in chunks.remainder().iter().enumerate() {
        if b == outer || b == inner {
            f(offset.saturating_add(i), b == outer)?;
        }
    }
    ControlFlow::Continue(())
}

/// One `KEY` or `KEY=VALUE` entry of an INFO column, as byte ranges into it.
/// The first `=` splits key from value; later ones belong to the value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct InfoEntry {
    pub(crate) key: Range<usize>,
    /// `None` for a flag (no `=`); `Some` and possibly empty otherwise.
    pub(crate) value: Option<Range<usize>>,
}

/// Calls `f` with every `;`-separated entry of `info`, in order, until `f`
/// breaks: exactly the entries of `info.split(|&b| b == b';')`, empty ones
/// included (`""` is one empty entry), so the caller decides what an empty
/// key or a `.` column means.
// r[impl vcf_reader.text_scan]
pub(crate) fn for_each_info_entry<B, F>(info: &[u8], f: F) -> ControlFlow<B>
where
    F: FnMut(InfoEntry) -> ControlFlow<B>,
{
    for_each_info_entry_at(Level::new(), info, f)
}

// r[impl io.simd_portable]
fn for_each_info_entry_at<B, F>(level: Level, info: &[u8], f: F) -> ControlFlow<B>
where
    F: FnMut(InfoEntry) -> ControlFlow<B>,
{
    dispatch!(level, simd => for_each_info_entry_simd(simd, info, f))
}

/// The entry being read: where it starts and its first `=`, if seen.
struct InfoSplitter<F> {
    start: usize,
    eq: Option<usize>,
    f: F,
}

impl<B, F: FnMut(InfoEntry) -> ControlFlow<B>> InfoSplitter<F> {
    #[inline(always)]
    fn eq_at(&mut self, pos: usize) {
        self.eq.get_or_insert(pos);
    }

    #[inline(always)]
    fn semi_at(&mut self, pos: usize) -> ControlFlow<B> {
        let entry = self.entry(pos);
        self.start = pos.saturating_add(1);
        self.eq = None;
        (self.f)(entry)
    }

    #[inline(always)]
    fn entry(&self, end: usize) -> InfoEntry {
        match self.eq {
            Some(eq) => InfoEntry { key: self.start..eq, value: Some(eq.saturating_add(1)..end) },
            None => InfoEntry { key: self.start..end, value: None },
        }
    }
}

/// One mask of `;` and `=` per vector; its set bits are walked in order, and
/// the `;` mask alone says which of the two each one is.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(clippy::arithmetic_side_effects, reason = "`bits - 1` runs only on a nonzero `bits`")]
fn for_each_info_entry_simd<S: Simd, B, F>(simd: S, info: &[u8], f: F) -> ControlFlow<B>
where
    F: FnMut(InfoEntry) -> ControlFlow<B>,
{
    let mut s = InfoSplitter { start: 0, eq: None, f };
    let mut chunks = info.chunks_exact(S::u8s::LEN);
    let mut offset = 0usize;
    for chunk in &mut chunks {
        let v = S::u8s::from_slice(simd, chunk);
        let semis = v.simd_eq(b';');
        let seps = semis | v.simd_eq(b'=');
        if seps.any_true() {
            let semi_bits = semis.to_bitmask();
            let mut bits = seps.to_bitmask();
            while bits != 0 {
                let lane = bits.trailing_zeros();
                let pos = offset.saturating_add(lane as usize);
                if (semi_bits >> lane) & 1 == 1 {
                    s.semi_at(pos)?;
                } else {
                    s.eq_at(pos);
                }
                bits &= bits - 1;
            }
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    for (i, &b) in chunks.remainder().iter().enumerate() {
        let pos = offset.saturating_add(i);
        match b {
            b';' => s.semi_at(pos)?,
            b'=' => s.eq_at(pos),
            _ => {}
        }
    }
    let last = s.entry(info.len());
    (s.f)(last)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::arithmetic_side_effects, reason = "test arithmetic on small offsets")]

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

    // r[verify io.simd_portable]
    #[hegel::test]
    fn tabs_and_nul_every_level_matches_scalar(tc: TestCase) {
        let line = text(&tc, 400);
        let tabs: Vec<usize> = (0..line.len()).filter(|&i| line[i] == b'\t').collect();
        let nul = line.iter().position(|&b| b == 0);
        let mut got = vec![7];
        for level in crate::simd_levels::levels() {
            assert_eq!(tabs_and_nul_at(level, &line, &mut got), nul, "{level:?}");
            assert_eq!(got, tabs, "{level:?}");
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

    /// Draws text dense in INFO and sample separators.
    fn vcf_text(tc: &TestCase, max: usize) -> Vec<u8> {
        let alphabet = gs::sampled_from(vec![b';', b'=', b',', b':', b'\t', b'A', b'.']);
        tc.draw(gs::vecs(alphabet).max_size(max))
    }

    // r[verify vcf_reader.text_scan]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn for_each_byte_every_level_matches_scalar(tc: TestCase) {
        let hay = vcf_text(&tc, 400);
        let stop_after = tc.draw(gs::integers::<usize>().max_value(400));
        for needle in *b"\t:," {
            let all: Vec<usize> = (0..hay.len()).filter(|&i| hay[i] == needle).collect();
            for level in crate::simd_levels::levels() {
                let mut got = Vec::new();
                let flow = for_each_byte_at(level, &hay, needle, |pos| {
                    got.push(pos);
                    if got.len() > stop_after {
                        ControlFlow::Break(pos)
                    } else {
                        ControlFlow::Continue(())
                    }
                });
                let expected = &all[..all.len().min(stop_after + 1)];
                assert_eq!(got, expected, "{level:?}");
                assert_eq!(flow.break_value(), all.get(stop_after).copied(), "{level:?}");
            }
        }
    }

    // r[verify vcf_reader.text_scan]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn for_each_of_two_every_level_matches_scalar(tc: TestCase) {
        let hay = vcf_text(&tc, 400);
        let stop_after = tc.draw(gs::integers::<usize>().max_value(400));
        let all: Vec<(usize, bool)> = (0..hay.len())
            .filter(|&i| hay[i] == b'\t' || hay[i] == b':')
            .map(|i| (i, hay[i] == b'\t'))
            .collect();
        for level in crate::simd_levels::levels() {
            let mut got = Vec::new();
            let flow = for_each_of_two_at(level, &hay, b'\t', b':', |pos, outer| {
                got.push((pos, outer));
                if got.len() > stop_after {
                    ControlFlow::Break(pos)
                } else {
                    ControlFlow::Continue(())
                }
            });
            assert_eq!(got, &all[..all.len().min(stop_after + 1)], "{level:?}");
            assert_eq!(flow.break_value(), all.get(stop_after).map(|&(p, _)| p), "{level:?}");
        }
    }

    /// The splitter's oracle: `split` and `position`, nothing shared with it.
    fn scalar_info_entries(info: &[u8]) -> Vec<InfoEntry> {
        let mut start = 0;
        info.split(|&b| b == b';')
            .map(|entry| {
                let e = match entry.iter().position(|&b| b == b'=') {
                    Some(eq) => InfoEntry {
                        key: start..start + eq,
                        value: Some(start + eq + 1..start + entry.len()),
                    },
                    None => InfoEntry { key: start..start + entry.len(), value: None },
                };
                start += entry.len() + 1;
                e
            })
            .collect()
    }

    // r[verify vcf_reader.text_scan]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn info_entries_every_level_match_scalar(tc: TestCase) {
        let info = vcf_text(&tc, 400);
        let expected = scalar_info_entries(&info);
        let stop_after = tc.draw(gs::integers::<usize>().max_value(expected.len()));
        for level in crate::simd_levels::levels() {
            let mut got = Vec::new();
            let flow = for_each_info_entry_at(level, &info, |e| {
                got.push(e);
                if got.len() > stop_after {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            });
            assert_eq!(got, expected[..expected.len().min(stop_after + 1)], "{level:?}");
            assert_eq!(flow.is_break(), stop_after < expected.len(), "{level:?}");
        }
    }

    // r[verify vcf_reader.text_scan]
    #[test]
    fn info_entries_flags_values_and_empty_entries() {
        let info = b"DP=10;DB;AF=0.5,0.25;;X=a=b;E=";
        let mut got = Vec::new();
        let _ = for_each_info_entry(info, |e| {
            got.push((&info[e.key], e.value.map(|v| &info[v])));
            ControlFlow::<()>::Continue(())
        });
        let expected: Vec<(&[u8], Option<&[u8]>)> = vec![
            (b"DP", Some(b"10")),
            (b"DB", None),
            (b"AF", Some(b"0.5,0.25")),
            (b"", None),
            (b"X", Some(b"a=b")),
            (b"E", Some(b"")),
        ];
        assert_eq!(got, expected);
    }
}
