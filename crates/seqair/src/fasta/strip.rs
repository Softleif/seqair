//! Vectorized newline-strip + ASCII-uppercase for fetched FASTA/FASTQ spans.
//!
//! Line terminators are sparse (one or two per line), so instead of a general
//! byte-compress the pass is split into three vectorized steps: a memchr-style
//! scan for the next `\n`/`\r`, a `memcpy` of the run between delimiters, and
//! an in-place uppercase of the copied run. Uppercasing is range-masked to
//! `a..=z` (0x61..=0x7A) — every other byte, including `>= 0x80`, passes
//! through unchanged, matching `u8::to_ascii_uppercase` exactly.

use fearless_simd::{Level, dispatch, prelude::*};
use fearless_simd_macros::simd;

/// Strip `\n`/`\r` from `raw`, appending the rest ASCII-uppercased to `out`.
///
/// Byte-for-byte equivalent to the scalar loop this replaced: a byte equal to
/// `b'\n'` or `b'\r'` is dropped, anything else goes through
/// `to_ascii_uppercase`. Consecutive delimiters (`\r\n`) yield empty runs,
/// which the copy turns into the no-op the old loop's skip was.
pub(super) fn strip_newlines_uppercase(raw: &[u8], out: &mut Vec<u8>) {
    strip_at(Level::new(), raw, out);
}

/// One dispatch per fetch, not per line: the whole pass runs inside the
/// kernel the level selects.
// r[impl io.simd_portable]
fn strip_at(level: Level, raw: &[u8], out: &mut Vec<u8>) {
    dispatch!(level, simd => strip_simd(simd, raw, out));
}

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
    let v = S::u8s::from_slice(simd, from);
    let lower = v.simd_ge(b'a') & v.simd_le(b'z');
    (v ^ lower.select(flip, keep)).store_slice(to);
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
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
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

    fn run(raw: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(raw.len());
        strip_newlines_uppercase(raw, &mut out);
        out
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
            let mut out = Vec::new();
            strip_at(level, &raw, &mut out);
            assert_eq!(out, expected, "{level:?}");
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
}
