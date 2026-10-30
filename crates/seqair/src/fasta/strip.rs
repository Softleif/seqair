//! Vectorized newline-strip + ASCII-uppercase for fetched FASTA/FASTQ spans.
//!
//! Line terminators are sparse (one or two per line), so instead of a general
//! byte-compress the pass is split into three vectorized steps: a memchr-style
//! scan for the next `\n`/`\r`, a `memcpy` of the run between delimiters, and
//! an in-place uppercase of the copied run. Uppercasing is range-masked to
//! `a..=z` (0x61..=0x7A) — every other byte, including `>= 0x80`, passes
//! through unchanged, matching `u8::to_ascii_uppercase` exactly.
//!
//! The `.fai` line layout says where the terminators are, so the scan is the
//! fallback: [`strip_uppercase`] copies whole lines at the stride the index
//! claims, checks the claim on the way, and scans only when the file disagrees.

use super::index::{LineLayout, LineTerminator};
use core::num::NonZeroUsize;
use fearless_simd::{Level, dispatch, prelude::*};
use fearless_simd_macros::simd;

// r[impl fasta.fetch.line_stride]
/// Strip `\n`/`\r` from `raw`, the bytes of a span starting at base `first` of
/// a record laid out as `layout`, appending the rest ASCII-uppercased to
/// `out`. The output is the scan's whatever the layout claims.
pub(super) fn strip_uppercase(raw: &[u8], layout: LineLayout, first: u64, out: &mut Vec<u8>) {
    stride_at(Level::new(), raw, Stride::new(layout, first), out);
}

/// Where a span's line terminators sit, as its record's layout claims.
#[derive(Debug, Clone, Copy)]
struct Stride {
    /// Bases left on the span's first line.
    first_line: usize,
    /// Nonzero, so every copied line makes progress.
    line_bases: NonZeroUsize,
    terminator: LineTerminator,
}

impl Stride {
    fn new(layout: LineLayout, first: u64) -> Self {
        let line_bases = layout.linebases();
        // The remainder is below `line_bases`, so this never saturates.
        let first_line = line_bases.get().saturating_sub(first % line_bases);
        // A line longer than the address space is the rest of any span.
        Self {
            first_line: usize::try_from(first_line).unwrap_or(usize::MAX),
            line_bases: NonZeroUsize::try_from(line_bases).unwrap_or(NonZeroUsize::MAX),
            terminator: layout.terminator(),
        }
    }
}

// r[impl io.simd_portable]
fn stride_at(level: Level, raw: &[u8], stride: Stride, out: &mut Vec<u8>) {
    dispatch!(level, simd => stride_simd(simd, raw, stride, out));
}

/// A skipped run that is not the terminator catches a line too long; a
/// terminator among the copied bytes catches a line too short or a stray
/// `\r`. Either way the index does not describe the file and the scan decides.
#[simd]
fn stride_simd<S: Simd>(simd: S, raw: &[u8], stride: Stride, out: &mut Vec<u8>) {
    let start = out.len();
    let clean = copy_lines(raw, stride, out)
        && out.get_mut(start..).is_some_and(|copied| uppercase_in_place(simd, copied));
    if !clean {
        out.truncate(start);
        strip_simd(simd, raw, out);
    }
}

/// Append `raw`'s lines to `out`, skipping the terminator between them;
/// `false` at the first skipped run that is not the terminator. A correct
/// span ends on a base, so one ending inside a terminator is a mismatch too.
fn copy_lines(raw: &[u8], stride: Stride, out: &mut Vec<u8>) -> bool {
    let terminator = stride.terminator.bytes();
    let mut rest = raw;
    let mut line = stride.first_line;
    loop {
        let (bases, after) = rest.split_at_checked(line).unwrap_or((rest, &[]));
        out.extend_from_slice(bases);
        if after.is_empty() {
            return true;
        }
        let Some(next) = after.strip_prefix(terminator) else {
            return false;
        };
        rest = next;
        line = stride.line_bases.get();
    }
}

/// Strip `\n`/`\r` from `raw` by scanning for them, appending the rest
/// ASCII-uppercased to `out`: a byte equal to `b'\n'` or `b'\r'` is dropped,
/// anything else goes through `to_ascii_uppercase`. Consecutive delimiters
/// (`\r\n`) yield empty runs, which the copy turns into a no-op.
///
/// The helpers below are `#[inline(always)]` rather than `#[simd]`: they run
/// inside this function's SIMD context, and a `#[simd]` helper that stays out
/// of line is a call (through its own dispatcher) per line.
#[simd]
fn strip_simd<S: Simd>(simd: S, raw: &[u8], out: &mut Vec<u8>) {
    let mut search = 0usize;
    while let Some(hit) = find_newline(simd, raw, search) {
        if let Some(run) = raw.get(search..hit) {
            append_uppercased(simd, run, out);
        }
        search = hit.saturating_add(1);
    }
    if let Some(tail) = raw.get(search..) {
        append_uppercased(simd, tail, out);
    }
}

/// Append `run` to `out`, uppercased. The bytes are copied first (a
/// `memcpy`) so the vector stores below land on initialised memory, but every
/// vector is *loaded from `run`*: reading back what the copy just wrote would
/// stall on store-to-load forwarding, and the ragged tail — one vector
/// overlapping the last whole one — would stall on the kernel's own store.
#[inline(always)]
fn append_uppercased<S: Simd>(simd: S, run: &[u8], out: &mut Vec<u8>) {
    let start = out.len();
    out.extend_from_slice(run);
    if let Some(appended) = out.get_mut(start..) {
        uppercase_into(simd, run, appended);
    }
}

/// Index of the first `\n` or `\r` at or after `from`, if any: memchr2, one
/// native-width vector per iteration. `any_true` gates the loop because it is
/// one instruction everywhere, while a bitmask costs four on NEON — so the
/// bitmask is built only for the vector that has the hit.
#[inline(always)]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn find_newline<S: Simd>(simd: S, raw: &[u8], from: usize) -> Option<usize> {
    let mut chunks = raw.get(from..)?.chunks_exact(S::u8s::LEN);
    let mut offset = from;
    for chunk in &mut chunks {
        let v = S::u8s::from_slice(simd, chunk);
        let hit = v.simd_eq(b'\n') | v.simd_eq(b'\r');
        if hit.any_true() {
            // A set lane makes the bitmask nonzero, so this is < LEN.
            let lane = hit.to_bitmask().trailing_zeros() as usize;
            return Some(offset.saturating_add(lane));
        }
        offset = offset.saturating_add(S::u8s::LEN);
    }
    find_newline_scalar(raw, offset)
}

/// `dst = src` with `a..=z` uppercased, every other byte unchanged; the two
/// have the same length. Only that range earns the 0x20 flip — NOT the
/// `& 0xDF` blanket uppercase, which would also mangle `[`, `{`, `~` and
/// bytes `>= 0x80`. The flip is `xor`ed in from a `select` against zero rather than
/// `select`ing between the byte and its flip: the latter is a blend on x86,
/// the former an `and`.
#[inline(always)]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn uppercase_into<S: Simd>(simd: S, src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    let n = S::u8s::LEN;
    let flip = S::u8s::splat(simd, 0x20);
    let keep = S::u8s::splat(simd, 0);
    let mut to = dst.chunks_exact_mut(n);
    for (from, to) in src.chunks_exact(n).zip(&mut to) {
        uppercase_vector(simd, from, to, flip, keep);
    }
    if to.into_remainder().is_empty() {
        return;
    }
    // Uppercasing is idempotent, so a ragged tail is one more vector that
    // overlaps the last whole one rather than a byte loop: FASTA lines are
    // 60–80 bytes, so the tail is a large share of every run.
    let tail = src.len().checked_sub(n);
    match tail.and_then(|t| Some((src.get(t..)?, dst.get_mut(t..)?))) {
        Some((from, to)) => uppercase_vector(simd, from, to, flip, keep),
        None => {
            dst.copy_from_slice(src);
            uppercase_scalar(dst);
        }
    }
}

/// One vector of [`uppercase_into`]. A function, not a closure: a closure
/// does not get the enclosing `#[simd]` level's target features, so on x86
/// every vector op in it becomes an out-of-line call.
#[inline(always)]
fn uppercase_vector<S: Simd>(simd: S, from: &[u8], to: &mut [u8], flip: S::u8s, keep: S::u8s) {
    uppercased::<S>(S::u8s::from_slice(simd, from), flip, keep).store_slice(to);
}

#[inline(always)]
fn uppercased<S: Simd>(v: S::u8s, flip: S::u8s, keep: S::u8s) -> S::u8s {
    let lower = v.simd_ge(b'a') & v.simd_le(b'z');
    v ^ lower.select(flip, keep)
}

/// [`uppercase_into`] in place, and `true` when no byte is `\n` or `\r`: the
/// check rides on the loads the uppercase makes anyway. This second pass over
/// the copied span measured no slower than uppercasing each line as it is
/// copied (which pays a ragged tail per line) and 3–4% faster on 1–10 kb spans.
#[inline(always)]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn uppercase_in_place<S: Simd>(simd: S, buf: &mut [u8]) -> bool {
    let n = S::u8s::LEN;
    let flip = S::u8s::splat(simd, 0x20);
    let keep = S::u8s::splat(simd, 0);
    let mut terminators = S::mask8s::splat(simd, false);
    let mut chunks = buf.chunks_exact_mut(n);
    for chunk in &mut chunks {
        terminators |= uppercase_vector_in_place(simd, chunk, flip, keep);
    }
    if chunks.into_remainder().is_empty() {
        return !terminators.any_true();
    }
    let tail = buf.len().checked_sub(n);
    match tail.and_then(|t| buf.get_mut(t..)) {
        Some(tail) => terminators |= uppercase_vector_in_place(simd, tail, flip, keep),
        None => {
            uppercase_scalar(buf);
            return !buf.iter().any(|&b| b == b'\n' || b == b'\r');
        }
    }
    !terminators.any_true()
}

#[inline(always)]
fn uppercase_vector_in_place<S: Simd>(
    simd: S,
    bytes: &mut [u8],
    flip: S::u8s,
    keep: S::u8s,
) -> S::mask8s {
    let v = S::u8s::from_slice(simd, bytes);
    uppercased::<S>(v, flip, keep).store_slice(bytes);
    v.simd_eq(b'\n') | v.simd_eq(b'\r')
}

/// Scalar fallback: the exact predicate of the loop this replaced.
fn find_newline_scalar(raw: &[u8], from: usize) -> Option<usize> {
    let rel = raw.get(from..)?.iter().position(|&b| b == b'\n' || b == b'\r')?;
    Some(from.saturating_add(rel))
}

/// Scalar fallback: exact `u8::to_ascii_uppercase` semantics.
fn uppercase_scalar(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        b.make_ascii_uppercase();
    }
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::unwrap_used,
    clippy::indexing_slicing,
    reason = "test code with small, known-valid values"
)]
mod tests {
    use super::*;
    use crate::fasta::FaiEntry;
    use hegel::prelude::*;

    /// Independent scalar oracle: the exact semantics of the loop this
    /// replaced, written as plainly as possible.
    fn strip_upper_reference(raw: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(raw.len());
        for &b in raw {
            if b != b'\n' && b != b'\r' {
                out.push(b.to_ascii_uppercase());
            }
        }
        out
    }

    fn scan(level: Level, raw: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        dispatch!(level, simd => strip_simd(simd, raw, &mut out));
        out
    }

    fn run(raw: &[u8]) -> Vec<u8> {
        scan(Level::new(), raw)
    }

    // r[verify fasta.fetch.newline_stripping]
    // r[verify fasta.fetch.uppercase]
    #[test]
    fn targeted_vectors_match_reference() {
        let mut cases: Vec<Vec<u8>> = vec![
            Vec::new(),
            b"\n".to_vec(),
            b"\r".to_vec(),
            b"\r\n".to_vec(),
            b"\n\r".to_vec(),
            b"\nacgt".to_vec(), // delimiter at offset 0
            b"acgt\n".to_vec(), // delimiter as the final byte
            b"acgt".to_vec(),   // no delimiter at all
            b"acgtnryswkmbdhv".to_vec(),
            b"ACGTNRYSWKMBDHV".to_vec(),
            (0u8..=255).collect(),
            vec![b'\n'; 100],
            b"\r\n".repeat(50),
            b"ACGT\n".repeat(40), // realistic FASTA line structure
            vec![0x80, 0xFF, 0xFE, b'[', b'{', b'~', b'@', b'+', b'`'],
        ];
        // Straddle the 16/32-byte vector bodies and their scalar tails, with
        // the delimiter leading, trailing, and absent.
        for len in [1usize, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 1024] {
            cases.push(vec![b'g'; len]);
            let mut lead = vec![b'g'; len];
            if let Some(first) = lead.first_mut() {
                *first = b'\n';
            }
            cases.push(lead);
            let mut trail = vec![b'g'; len];
            if let Some(last) = trail.last_mut() {
                *last = b'\r';
            }
            cases.push(trail);
        }
        for case in &cases {
            assert_eq!(run(case), strip_upper_reference(case), "input: {case:?}");
        }
    }

    // r[verify fasta.fetch.newline_stripping]
    // r[verify fasta.fetch.uppercase]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn every_level_matches_scalar_reference(tc: TestCase) {
        let raw = tc.draw(gs::binary().max_size(1024));
        let expected = strip_upper_reference(&raw);
        for level in crate::simd_levels::levels() {
            assert_eq!(scan(level, &raw), expected, "{level:?}");
        }
    }

    // r[verify fasta.fetch.newline_stripping]
    // r[verify fasta.fetch.uppercase]
    /// Every byte value at every length, covering the vector bodies, their
    /// scalar tails, and delimiters at arbitrary positions: the dispatched
    /// SIMD path must equal the oracle bit-for-bit.
    #[hegel::test]
    fn simd_matches_scalar_reference(tc: TestCase) {
        let raw = tc.draw(gs::binary().max_size(4096));
        assert_eq!(run(&raw), strip_upper_reference(&raw));
    }

    /// Soft-masked bases and IUPAC codes, both cases.
    const BASES: &[u8] = b"ACGTNRYSWKMBDHVacgtnryswkmbdhv";

    /// One record's sequence block as `samtools faidx` would index it, each
    /// line terminated, with its entry (offset 0). An empty terminator is the
    /// single-line record with no trailing newline.
    #[hegel::composite]
    fn record(tc: &TestCase) -> (Vec<u8>, FaiEntry) {
        use LineTerminator::{CrLf, Empty, Lf};
        let terminator = tc.draw(gs::sampled_from(vec![Lf, CrLf, Empty]));
        let len = tc.draw(gs::integers::<usize>().min_value(1).max_value(1500));
        let bases = tc.draw(gs::vecs(gs::sampled_from(BASES.to_vec())).min_size(len).max_size(len));
        let linebases = if terminator == Empty {
            len
        } else {
            tc.draw(gs::integers::<usize>().min_value(1).max_value(200))
        };
        let mut file = Vec::new();
        for line in bases.chunks(linebases) {
            file.extend_from_slice(line);
            file.extend_from_slice(terminator.bytes());
        }
        let linewidth = linebases as u64 + terminator.width();
        let entry = FaiEntry {
            name: "r".into(),
            length: len as u64,
            offset: 0,
            layout: LineLayout::new(linebases as u64, linewidth).unwrap(),
            qual_offset: None,
        };
        (file, entry)
    }

    /// A span of `entry` anywhere, mid-line ends and the short last line
    /// included, as `(first base, file bytes start..end)`.
    fn draw_span(tc: &TestCase, entry: &FaiEntry) -> (u64, usize, usize) {
        let first = tc.draw(gs::integers::<u64>().max_value(entry.length - 1));
        let last = tc.draw(gs::integers::<u64>().min_value(first).max_value(entry.length - 1));
        let start = entry.byte_offset(first).unwrap() as usize;
        let end = entry.byte_offset(last).unwrap() as usize + 1;
        (first, start, end)
    }

    fn stride(level: Level, raw: &[u8], entry: &FaiEntry, first: u64) -> Vec<u8> {
        let mut out = Vec::new();
        stride_at(level, raw, Stride::new(entry.layout, first), &mut out);
        out
    }

    // r[verify fasta.fetch.line_stride]
    // r[verify io.simd_portable]
    /// On a file its index describes, the stride gives the scan's output and
    /// never needs the scan to get it.
    #[hegel::test]
    fn stride_matches_the_scan_on_an_honest_index(tc: TestCase) {
        let (file, entry) = tc.draw(record());
        let (first, start, end) = draw_span(&tc, &entry);
        let raw = &file[start..end];
        let mut copied = Vec::new();
        assert!(
            copy_lines(raw, Stride::new(entry.layout, first), &mut copied),
            "skipped runs must be terminators"
        );
        assert!(!copied.iter().any(|&b| b == b'\n' || b == b'\r'), "copied a terminator");
        for level in crate::simd_levels::levels() {
            assert_eq!(stride(level, raw, &entry, first), strip_upper_reference(raw), "{level:?}");
        }
    }

    // r[verify fasta.fetch.line_stride]
    /// An index that does not describe its file — a line one byte longer or
    /// shorter, two lines joined, a stray `\r` or `\n` — still gives exactly
    /// the scan's output. The span crosses the corrupted line's end, where
    /// the stride would go wrong.
    #[hegel::test]
    fn stride_matches_the_scan_when_the_index_lies(tc: TestCase) {
        let (mut file, entry) = tc.draw(record());
        let (lb, lw, len) =
            (entry.layout.linebases().get(), entry.layout.linewidth(), entry.length);
        let lines = len.div_ceil(lb);
        // A line with another after it, so its terminator lies inside a span.
        let line = tc.draw(gs::integers::<u64>().max_value(lines.saturating_sub(2)));
        let line_end = (line * lw + lb.min(len - line * lb)) as usize;
        let anywhere = tc.draw(gs::integers::<usize>().max_value(file.len() - 1));
        match tc.draw(gs::integers::<u8>().max_value(4)) {
            0 => file.insert(line_end, tc.draw(gs::sampled_from(BASES.to_vec()))),
            1 => drop(file.remove(line_end - 1)),
            2 if lw > lb => file[line_end] = tc.draw(gs::sampled_from(BASES.to_vec())),
            3 => file.insert(anywhere, b'\r'),
            _ => file[anywhere] = b'\n',
        }
        let first = tc.draw(gs::integers::<u64>().max_value(line * lb));
        let last = tc.draw(
            gs::integers::<u64>().min_value(((line + 1) * lb).min(len - 1)).max_value(len - 1),
        );
        let start = entry.byte_offset(first).unwrap() as usize;
        let end = entry.byte_offset(last).unwrap() as usize + 1;
        let raw = &file[start.min(file.len())..end.min(file.len())];
        let expected = strip_upper_reference(raw);
        let mut copied = Vec::new();
        if !copy_lines(raw, Stride::new(entry.layout, first), &mut copied) {
            tc.event("a skipped run is not the terminator");
        } else if copied.contains(&b'\n') || copied.contains(&b'\r') {
            tc.event("a terminator among the copied bases");
        }
        for level in crate::simd_levels::levels() {
            assert_eq!(stride(level, raw, &entry, first), expected, "{level:?}");
        }
    }
}
