use crate::Base;
use fearless_simd::{Level, dispatch, prelude::*, u8x16};
use fearless_simd_macros::simd;

/// One BAM sequence base: the sixteen IUPAC codes a BAM `seq` field can
/// store, each variant's discriminant being its 4-bit code. Lossless where
/// [`Base`] is not; use `Base` for calling and `IupacBase` for records that
/// pass through unchanged.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, bytemuck::NoUninit)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u8)]
#[must_use]
pub enum IupacBase {
    /// `=`: identical to the reference.
    Eq = 0,
    /// Adenine
    A = 1,
    /// Cytosine
    C = 2,
    /// A or C
    M = 3,
    /// Guanine
    G = 4,
    /// A or G
    R = 5,
    /// C or G
    S = 6,
    /// A, C or G
    V = 7,
    /// Thymine
    T = 8,
    /// A or T
    W = 9,
    /// C or T
    Y = 10,
    /// A, C or T
    H = 11,
    /// G or T
    K = 12,
    /// A, G or T
    D = 13,
    /// C, G or T
    B = 14,
    /// Any base
    #[default]
    N = 15,
}

const _: () = assert!(size_of::<IupacBase>() == 1);

impl IupacBase {
    /// Every variant, indexed by code.
    pub const ALL: [IupacBase; 16] = [
        Self::Eq,
        Self::A,
        Self::C,
        Self::M,
        Self::G,
        Self::R,
        Self::S,
        Self::V,
        Self::T,
        Self::W,
        Self::Y,
        Self::H,
        Self::K,
        Self::D,
        Self::B,
        Self::N,
    ];

    /// The upper-case letter of each code, indexed by code (htslib's
    /// `seq_nt16_str`).
    pub const ASCII: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

    /// ASCII → 4-bit code, htslib's `seq_nt16_table`: either case, `U` as
    /// `T`, the colour-space digits `0`..`3` as A/C/G/T, anything else `N`.
    /// The scalar reference for vectorized encoders.
    #[allow(
        clippy::indexing_slicing,
        clippy::cast_possible_truncation,
        reason = "const indices: i < 16, bytes < 256"
    )]
    pub const CODE_OF_ASCII: [u8; 256] = {
        let mut table = [15u8; 256];
        let mut i = 0;
        while i < 16 {
            let letter = Self::ASCII[i];
            table[letter as usize] = i as u8;
            table[letter.to_ascii_lowercase() as usize] = i as u8;
            i += 1;
        }
        table[b'U' as usize] = 8;
        table[b'u' as usize] = 8;
        table[b'0' as usize] = 1;
        table[b'1' as usize] = 2;
        table[b'2' as usize] = 4;
        table[b'3' as usize] = 8;
        table
    };

    /// The base for a BAM 4-bit code, or `None` above 15.
    #[inline]
    pub fn from_code(code: u8) -> Option<Self> {
        Self::ALL.get(usize::from(code)).copied()
    }

    /// The base for the low four bits of `nibble`; the high bits are ignored.
    #[inline]
    pub fn from_nibble(nibble: u8) -> Self {
        Self::ALL.get(usize::from(nibble & 0x0F)).copied().unwrap_or_default()
    }

    /// The BAM 4-bit code.
    #[inline]
    pub fn code(self) -> u8 {
        self as u8
    }

    /// The base for an ASCII byte, as htslib reads SAM text (see
    /// [`Self::CODE_OF_ASCII`]).
    #[inline]
    pub fn from_ascii(byte: u8) -> Self {
        Self::from_nibble(Self::CODE_OF_ASCII.get(usize::from(byte)).copied().unwrap_or(15))
    }

    /// The upper-case letter, `=` for [`IupacBase::Eq`].
    #[inline]
    pub fn to_ascii(self) -> u8 {
        Self::ASCII.get(usize::from(self.code())).copied().unwrap_or(b'N')
    }

    /// The IUPAC complement.
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::M => Self::K,
            Self::K => Self::M,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::V => Self::B,
            Self::B => Self::V,
            Self::H => Self::D,
            Self::D => Self::H,
            Self::W | Self::S | Self::N | Self::Eq => self,
        }
    }

    /// A/C/G/T as [`Base`]; every other code (ambiguity codes, `=`) is
    /// [`Base::Unknown`].
    pub fn to_base(self) -> Base {
        match self {
            Self::A => Base::A,
            Self::C => Base::C,
            Self::G => Base::G,
            Self::T => Base::T,
            _ => Base::Unknown,
        }
    }

    // r[impl types.iupac_base.slices]
    /// [`Self::from_ascii`] over a slice.
    pub fn from_ascii_slice(bytes: &[u8]) -> Vec<Self> {
        let mut out = Vec::with_capacity(bytes.len());
        Self::from_ascii_into(bytes, &mut out);
        out
    }

    // r[impl types.iupac_base.slices]
    /// [`Self::to_ascii`] over a slice.
    pub fn to_ascii_slice(bases: &[Self]) -> Vec<u8> {
        let mut out = Vec::with_capacity(bases.len());
        Self::to_ascii_into(bases, &mut out);
        out
    }

    // r[impl types.iupac_base.bulk]
    // r[impl types.iupac_base.simd]
    /// Append [`Self::from_ascii`] of every byte to `out`.
    pub fn from_ascii_into(bytes: &[u8], out: &mut Vec<Self>) {
        from_ascii_at(Level::new(), bytes, out);
    }

    // r[impl types.iupac_base.bulk]
    // r[impl types.iupac_base.simd]
    /// Append [`Self::to_ascii`] of every base to `out`.
    pub fn to_ascii_into(bases: &[Self], out: &mut Vec<u8>) {
        to_ascii_at(Level::new(), bases, out);
    }

    // r[impl types.iupac_base.bulk]
    // r[impl types.iupac_base.simd]
    /// [`Self::complement`] every base.
    pub fn complement_in_place(bases: &mut [Self]) {
        complement_at(Level::new(), bases);
    }

    // r[impl types.iupac_base.bulk]
    // r[impl types.iupac_base.simd]
    /// Reverse `bases` and [`Self::complement`] every one: the other strand,
    /// read 5' to 3'.
    pub fn reverse_complement_in_place(bases: &mut [Self]) {
        reverse_complement_at(Level::new(), bases);
    }

    // r[impl types.iupac_base.as_codes_mut]
    /// The codes of `bases`, to write in place: how a bulk kernel (here and
    /// in `seqair`) fills an `IupacBase` buffer without an intermediate one.
    /// Hidden and `unsafe` rather than a safe fill that masks afterwards:
    /// that mask pass measured 2x on a 150-base complement and 1.5x at 100 kb.
    ///
    /// # Safety
    /// Every byte MUST be below 16 when the returned borrow ends, unwinding
    /// included: the slice is `bases`, and a larger byte is not an
    /// `IupacBase`.
    #[doc(hidden)]
    pub unsafe fn as_codes_mut(bases: &mut [Self]) -> &mut [u8] {
        // SAFETY: `IupacBase` is `repr(u8)`, so the layouts match and every
        // byte is a code now; the caller leaves only codes behind.
        unsafe { std::slice::from_raw_parts_mut(bases.as_mut_ptr().cast::<u8>(), bases.len()) }
    }
}

// r[impl io.simd_portable]
fn from_ascii_at(level: Level, bytes: &[u8], out: &mut Vec<IupacBase>) {
    let start = out.len();
    out.resize(start.saturating_add(bytes.len()), IupacBase::N);
    if let Some(dst) = out.get_mut(start..) {
        // SAFETY: `FromAscii` writes only `CODE_OF_ASCII` values, all < 16.
        let dst = unsafe { IupacBase::as_codes_mut(dst) };
        dispatch!(level, simd => map_into(simd, FromAscii, bytes, dst));
    }
}

// r[impl io.simd_portable]
fn to_ascii_at(level: Level, bases: &[IupacBase], out: &mut Vec<u8>) {
    let start = out.len();
    out.resize(start.saturating_add(bases.len()), 0);
    if let Some(dst) = out.get_mut(start..) {
        let codes = bytemuck::cast_slice(bases);
        dispatch!(level, simd => map_into(simd, Nibble(&IupacBase::ASCII), codes, dst));
    }
}

// r[impl io.simd_portable]
fn complement_at(level: Level, bases: &mut [IupacBase]) {
    // SAFETY: `COMPLEMENT` maps a code to a code.
    let codes = unsafe { IupacBase::as_codes_mut(bases) };
    dispatch!(level, simd => map_in_place(simd, Nibble(&COMPLEMENT), codes));
}

// r[impl io.simd_portable]
fn reverse_complement_at(level: Level, bases: &mut [IupacBase]) {
    // SAFETY: `COMPLEMENT` maps a code to a code; reversing only moves codes.
    let codes = unsafe { IupacBase::as_codes_mut(bases) };
    dispatch!(level, simd => reverse_complement(simd, codes));
}

/// `complement()` by code.
#[allow(clippy::indexing_slicing, reason = "i < 16 = ALL.len()")]
const COMPLEMENT: [u8; 16] = {
    let mut table = [0u8; 16];
    let mut i = 0;
    while i < 16 {
        table[i] = IupacBase::ALL[i].complement() as u8;
        i += 1;
    }
    table
};

/// `CODE_OF_ASCII ^ 15` for the bytes `0x30..0x60`, one row per high nibble:
/// with `0x6_`/`0x7_`, the lower case of `0x4_`/`0x5_`, the only rows where a
/// byte is not `N`. Stored `^ 15` so that a byte no row picks comes out
/// `0 ^ 15`, `N`.
#[allow(clippy::indexing_slicing, reason = "row < 3, i < 16: indices < 0x60")]
const ASCII_ROWS: [[u8; 16]; 3] = {
    let mut rows = [[0u8; 16]; 3];
    let mut row = 0;
    while row < 3 {
        let mut i = 0;
        while i < 16 {
            rows[row][i] = IupacBase::CODE_OF_ASCII[((row + 3) << 4) | i] ^ 15;
            i += 1;
        }
        row += 1;
    }
    rows
};

// The rows hold all of `CODE_OF_ASCII`: every byte decodes as
// `FromAscii::vector` picks it.
#[allow(clippy::indexing_slicing, reason = "b < 256, row < 3, b & 15 < 16")]
const _: () = {
    let mut b = 0;
    while b < 256 {
        let picked = match b >> 4 {
            3 => ASCII_ROWS[0][b & 15],
            4 | 6 => ASCII_ROWS[1][b & 15],
            5 | 7 => ASCII_ROWS[2][b & 15],
            _ => 0,
        };
        assert!(IupacBase::CODE_OF_ASCII[b] == picked ^ 15);
        b += 1;
    }
};

/// A byte map with a vector form and a scalar form for tails; the two agree.
trait ByteMap: Copy {
    fn vector<S: Simd>(self, simd: S, v: S::u8s) -> S::u8s;
    fn scalar(self, byte: u8) -> u8;
}

/// A code through a 16-entry table. Defined on codes (bytes below 16) only:
/// the shuffle is not masked, so above 15 the two forms may disagree.
#[derive(Clone, Copy)]
struct Nibble(&'static [u8; 16]);

impl ByteMap for Nibble {
    #[inline(always)]
    fn vector<S: Simd>(self, simd: S, v: S::u8s) -> S::u8s {
        S::u8s::block_splat(u8x16::from_slice(simd, self.0)).swizzle_dyn_within_blocks(v)
    }

    #[inline(always)]
    #[allow(clippy::indexing_slicing, reason = "a nibble is < 16 = table.len()")]
    fn scalar(self, byte: u8) -> u8 {
        self.0[usize::from(byte & 15)]
    }
}

/// ASCII to code: the low nibble through each of [`ASCII_ROWS`], keeping the
/// row the high nibble names.
#[derive(Clone, Copy)]
struct FromAscii;

impl ByteMap for FromAscii {
    #[inline(always)]
    #[allow(
        clippy::arithmetic_side_effects,
        reason = "a lanewise shift by 4 of a `u8` cannot overflow"
    )]
    fn vector<S: Simd>(self, simd: S, v: S::u8s) -> S::u8s {
        let zero = S::u8s::splat(simd, 0);
        let lo = v & 0x0F;
        let hi = v >> 4;
        // Folds `0x6_`/`0x7_` onto `0x4_`/`0x5_`; no other high nibble
        // becomes 4 or 5.
        let letters = hi & 0x0D;
        let [digits, upper_a, upper_p] = &ASCII_ROWS;
        let digits = S::u8s::block_splat(u8x16::from_slice(simd, digits));
        let upper_a = S::u8s::block_splat(u8x16::from_slice(simd, upper_a));
        let upper_p = S::u8s::block_splat(u8x16::from_slice(simd, upper_p));
        let picked = hi.simd_eq(3).select(digits.swizzle_dyn_within_blocks(lo), zero)
            | letters.simd_eq(4).select(upper_a.swizzle_dyn_within_blocks(lo), zero)
            | letters.simd_eq(5).select(upper_p.swizzle_dyn_within_blocks(lo), zero);
        picked ^ 0x0F
    }

    #[inline(always)]
    #[allow(clippy::indexing_slicing, reason = "a byte is < 256 = table.len()")]
    fn scalar(self, byte: u8) -> u8 {
        IupacBase::CODE_OF_ASCII[usize::from(byte)]
    }
}

/// `map` from `src` into `dst`, which is as long. A ragged tail is one more
/// vector overlapping the last whole one; shorter than a vector, all scalar.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn map_into<S: Simd, M: ByteMap>(simd: S, map: M, src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    let n = S::u8s::LEN;
    let mut chunks = src.chunks_exact(n);
    for (s, d) in (&mut chunks).zip(dst.chunks_exact_mut(n)) {
        map.vector(simd, S::u8s::from_slice(simd, s)).store_slice(d);
    }
    if chunks.remainder().is_empty() {
        return;
    }
    let last = src.len().checked_sub(n);
    match last.and_then(|t| Some((src.get(t..)?, dst.get_mut(t..)?))) {
        Some((s, d)) => map.vector(simd, S::u8s::from_slice(simd, s)).store_slice(d),
        None => {
            for (d, &s) in dst.iter_mut().zip(src) {
                *d = map.scalar(s);
            }
        }
    }
}

/// `map` over `bytes` in place. A ragged tail is the last vector, loaded
/// before anything is stored, so the overlap is mapped from the input once
/// (and loading it first spares a store-forwarding stall).
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn map_in_place<S: Simd, M: ByteMap>(simd: S, map: M, bytes: &mut [u8]) {
    let n = S::u8s::LEN;
    let tail_start = (!bytes.len().is_multiple_of(n)).then(|| bytes.len().checked_sub(n)).flatten();
    let tail = tail_start.and_then(|t| bytes.get(t..)).map(|t| S::u8s::from_slice(simd, t));
    let mut chunks = bytes.chunks_exact_mut(n);
    for chunk in &mut chunks {
        map.vector(simd, S::u8s::from_slice(simd, chunk)).store_slice(chunk);
    }
    let remainder = chunks.into_remainder();
    match (tail, tail_start) {
        (Some(v), Some(t)) => {
            if let Some(to) = bytes.get_mut(t..) {
                map.vector(simd, v).store_slice(to);
            }
        }
        _ => {
            for b in remainder {
                *b = map.scalar(*b);
            }
        }
    }
}

/// Swap a vector from the front with one from the back, each reversed and
/// complemented, until fewer than two vectors are left. One to two vectors
/// left are two overlapping ones, both loaded before either is stored.
#[simd]
#[allow(clippy::arithmetic_side_effects, reason = "len >= 2 * LEN or >= LEN where subtracted")]
fn reverse_complement<S: Simd>(simd: S, mut bytes: &mut [u8]) {
    let n = S::u8s::LEN;
    let table = Nibble(&COMPLEMENT);
    while bytes.len() >= 2 * n {
        let (front, rest) = std::mem::take(&mut bytes).split_at_mut(n);
        let (middle, back) = rest.split_at_mut(rest.len() - n);
        let a = S::u8s::from_slice(simd, front);
        let b = S::u8s::from_slice(simd, back);
        table.vector(simd, b.reverse()).store_slice(front);
        table.vector(simd, a.reverse()).store_slice(back);
        bytes = middle;
    }
    let len = bytes.len();
    if len >= n {
        // `n <= len`, so neither split can fail.
        let a = S::u8s::from_slice(simd, bytes.split_at(n).0);
        let b = S::u8s::from_slice(simd, bytes.split_at(len - n).1);
        table.vector(simd, b.reverse()).store_slice(bytes.split_at_mut(n).0);
        table.vector(simd, a.reverse()).store_slice(bytes.split_at_mut(len - n).1);
    } else {
        bytes.reverse();
        for b in bytes {
            *b = table.scalar(*b);
        }
    }
}

impl From<Base> for IupacBase {
    fn from(base: Base) -> Self {
        match base {
            Base::A => Self::A,
            Base::C => Self::C,
            Base::G => Self::G,
            Base::T => Self::T,
            Base::Unknown => Self::N,
        }
    }
}

impl std::fmt::Display for IupacBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", char::from(self.to_ascii()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::TestCase;
    use hegel::generators as gs;

    // Hand-written, deliberately not derived from `complement`'s match.
    const COMPLEMENT_PAIRS: [(u8, u8); 10] = [
        (b'A', b'T'),
        (b'C', b'G'),
        (b'M', b'K'),
        (b'R', b'Y'),
        (b'W', b'W'),
        (b'S', b'S'),
        (b'V', b'B'),
        (b'H', b'D'),
        (b'N', b'N'),
        (b'=', b'='),
    ];

    // r[verify types.iupac_base.type]
    #[test]
    fn discriminants_are_the_sam_spec_codes() {
        for (code, letter) in b"=ACMGRSVTWYHKDBN".iter().enumerate() {
            let base = IupacBase::ALL[code];
            assert_eq!(usize::from(base.code()), code);
            assert_eq!(base.to_ascii(), *letter);
        }
    }

    // r[verify types.iupac_base.code]
    #[test]
    fn code_round_trips_and_rejects_above_15() {
        for code in 0..=255u8 {
            match IupacBase::from_code(code) {
                Some(b) => {
                    assert!(code < 16);
                    assert_eq!(b.code(), code);
                }
                None => assert!(code >= 16),
            }
            assert_eq!(IupacBase::from_nibble(code).code(), code & 15);
        }
    }

    // r[verify types.iupac_base.from_ascii]
    // r[verify types.iupac_base.to_ascii]
    #[test]
    fn every_byte_round_trips_through_ascii() {
        for byte in 0..=255u8 {
            let base = IupacBase::from_ascii(byte);
            assert_eq!(IupacBase::from_ascii(base.to_ascii()), base);
            let expected = match byte {
                b'U' | b'u' | b'3' => IupacBase::T,
                b'0' => IupacBase::A,
                b'1' => IupacBase::C,
                b'2' => IupacBase::G,
                b if b"=ACMGRSVTWYHKDBN".contains(&b.to_ascii_uppercase()) => {
                    IupacBase::from_ascii(b.to_ascii_uppercase())
                }
                _ => IupacBase::N,
            };
            assert_eq!(base, expected, "byte {byte}");
        }
    }

    // r[verify types.iupac_base.slices]
    #[hegel::test]
    fn slices_equal_elementwise(tc: TestCase) {
        let bytes: Vec<u8> = tc.draw(gs::vecs(gs::integers::<u8>()).max_size(100));
        let bases = IupacBase::from_ascii_slice(&bytes);
        assert_eq!(bases, bytes.iter().map(|&b| IupacBase::from_ascii(b)).collect::<Vec<_>>());
        assert_eq!(
            IupacBase::to_ascii_slice(&bases),
            bases.iter().map(|b| b.to_ascii()).collect::<Vec<_>>()
        );
    }

    /// Bytes of every value, in a random order with random repeats, so the
    /// 256 byte values meet every vector lane and tail position.
    fn any_bytes(tc: &TestCase) -> Vec<u8> {
        let len = tc.draw(gs::integers::<usize>().max_value(300));
        (0..len).map(|_| tc.draw(gs::integers::<u8>())).collect()
    }

    fn any_bases(tc: &TestCase) -> Vec<IupacBase> {
        let len = tc.draw(gs::integers::<usize>().max_value(300));
        (0..len).map(|_| IupacBase::from_nibble(tc.draw(gs::integers::<u8>()))).collect()
    }

    // r[verify types.iupac_base.bulk]
    #[hegel::test]
    fn bulk_appends_equal_elementwise(tc: TestCase) {
        let bytes = any_bytes(&tc);
        let prefix = tc.draw(gs::integers::<usize>().max_value(40));
        let mut out = vec![IupacBase::R; prefix];
        IupacBase::from_ascii_into(&bytes, &mut out);
        let (kept, new) = out.split_at(prefix);
        assert!(kept.iter().all(|&b| b == IupacBase::R));
        assert_eq!(new, bytes.iter().map(|&b| IupacBase::from_ascii(b)).collect::<Vec<_>>());

        let bases = any_bases(&tc);
        let mut out = vec![b'x'; prefix];
        IupacBase::to_ascii_into(&bases, &mut out);
        let (kept, new) = out.split_at(prefix);
        assert!(kept.iter().all(|&b| b == b'x'));
        assert_eq!(new, bases.iter().map(|b| b.to_ascii()).collect::<Vec<_>>());
    }

    // r[verify types.iupac_base.bulk]
    #[hegel::test]
    fn complement_in_place_and_reverse_complement(tc: TestCase) {
        let bases = any_bases(&tc);
        let mut c = bases.clone();
        IupacBase::complement_in_place(&mut c);
        assert_eq!(c, bases.iter().map(|b| b.complement()).collect::<Vec<_>>());
        let mut rc = bases.clone();
        IupacBase::reverse_complement_in_place(&mut rc);
        for (b, original) in rc.iter().zip(bases.iter().rev()) {
            assert_eq!(*b, original.complement());
        }
    }

    /// `bytes` placed `offset` bytes into a buffer, so a kernel sees it at
    /// every alignment.
    fn at_offset<T: Copy>(fill: T, offset: usize, items: &[T]) -> Vec<T> {
        std::iter::repeat_n(fill, offset).chain(items.iter().copied()).collect()
    }

    /// Every level's `from_ascii_at` and `to_ascii_at` on `bytes`, appended
    /// after a prefix, against the element-wise scalar functions.
    fn check_ascii_at_every_level(bytes: &[u8], prefix: usize) {
        let expected: Vec<IupacBase> = bytes.iter().map(|&b| IupacBase::from_ascii(b)).collect();
        let ascii: Vec<u8> = expected.iter().map(|b| b.to_ascii()).collect();
        for level in crate::simd_levels::levels() {
            let mut out = vec![IupacBase::R; prefix];
            from_ascii_at(level, bytes, &mut out);
            assert!(out[..prefix].iter().all(|&b| b == IupacBase::R), "{level:?}");
            assert_eq!(&out[prefix..], expected, "{level:?}, len={}", bytes.len());
            let mut out = vec![b'x'; prefix];
            to_ascii_at(level, &expected, &mut out);
            assert!(out[..prefix].iter().all(|&b| b == b'x'), "{level:?}");
            assert_eq!(&out[prefix..], ascii, "{level:?}, len={}", bytes.len());
        }
    }

    /// Every level's complement and reverse complement of `bases`.
    fn check_complement_at_every_level(bases: &mut [IupacBase]) {
        let complement: Vec<IupacBase> = bases.iter().map(|b| b.complement()).collect();
        let reverse_complement: Vec<IupacBase> = complement.iter().rev().copied().collect();
        let original = bases.to_vec();
        for level in crate::simd_levels::levels() {
            bases.copy_from_slice(&original);
            complement_at(level, bases);
            assert_eq!(bases, complement, "{level:?}, len={}", bases.len());
            bases.copy_from_slice(&original);
            reverse_complement_at(level, bases);
            assert_eq!(bases, reverse_complement, "{level:?}, len={}", bases.len());
        }
    }

    // r[verify types.iupac_base.simd]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn ascii_kernels_match_scalar_at_every_level(tc: TestCase) {
        let bytes = any_bytes(&tc);
        let offset = tc.draw(gs::integers::<usize>().max_value(63));
        let prefix = tc.draw(gs::integers::<usize>().max_value(63));
        let buf = at_offset(b'q', offset, &bytes);
        check_ascii_at_every_level(&buf[offset..], prefix);
    }

    // r[verify types.iupac_base.simd]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn complement_kernels_match_scalar_at_every_level(tc: TestCase) {
        let bases = any_bases(&tc);
        let offset = tc.draw(gs::integers::<usize>().max_value(63));
        let mut buf = at_offset(IupacBase::Eq, offset, &bases);
        check_complement_at_every_level(&mut buf[offset..]);
    }

    // r[verify types.iupac_base.simd]
    /// All 256 bytes (and all 16 codes, repeated) at every lane and tail
    /// position of the widest vector, at every level.
    #[test]
    fn kernels_cover_every_byte_at_every_position() {
        let all: Vec<u8> = (0..=255).collect();
        let codes: Vec<IupacBase> = IupacBase::ALL.iter().copied().cycle().take(256).collect();
        for pad in 0..64 {
            let bytes = at_offset(b'A', pad, &all);
            check_ascii_at_every_level(&bytes, 0);
            check_ascii_at_every_level(&bytes[pad..], pad);
            let mut bases = at_offset(IupacBase::T, pad, &codes);
            check_complement_at_every_level(&mut bases);
            check_complement_at_every_level(&mut bases[pad..]);
        }
    }

    // r[verify types.iupac_base.as_codes_mut]
    /// The contract `as_codes_mut` relies on, checked on plain bytes, where a
    /// stray byte is a failure rather than an invalid `IupacBase`: from any
    /// ASCII, and from any codes, every kernel leaves only codes.
    #[hegel::test]
    fn kernels_write_only_codes(tc: TestCase) {
        let bytes = any_bytes(&tc);
        let codes: Vec<u8> = bytes.iter().map(|b| b & 0x0F).collect();
        for level in crate::simd_levels::levels() {
            let mut out = vec![0u8; bytes.len()];
            dispatch!(level, simd => map_into(simd, FromAscii, &bytes, &mut out));
            assert!(out.iter().all(|&b| b < 16), "{level:?}");
            let mut out = codes.clone();
            dispatch!(level, simd => map_in_place(simd, Nibble(&COMPLEMENT), &mut out));
            assert!(out.iter().all(|&b| b < 16), "{level:?}");
            let mut out = codes.clone();
            dispatch!(level, simd => reverse_complement(simd, &mut out));
            assert!(out.iter().all(|&b| b < 16), "{level:?}");
        }
    }

    // r[verify types.iupac_base.complement]
    #[test]
    fn complement_is_an_involution_matching_the_table() {
        for b in IupacBase::ALL {
            assert_eq!(b.complement().complement(), b);
        }
        for (x, y) in COMPLEMENT_PAIRS {
            let (x, y) = (IupacBase::from_ascii(x), IupacBase::from_ascii(y));
            assert_eq!(x.complement(), y);
            assert_eq!(y.complement(), x);
        }
    }

    // r[verify types.iupac_base.base_conversions]
    #[test]
    fn base_conversions() {
        for base in [Base::A, Base::C, Base::G, Base::T, Base::Unknown] {
            let iupac = IupacBase::from(base);
            assert_eq!(iupac.to_base(), base);
            assert_eq!(iupac.to_ascii(), char::from(base.as_str().as_bytes()[0]) as u8);
        }
        for b in IupacBase::ALL {
            let known = matches!(b, IupacBase::A | IupacBase::C | IupacBase::G | IupacBase::T);
            assert_eq!(b.to_base() == Base::Unknown, !known);
        }
    }
}
