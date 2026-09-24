use crate::SmolStr;
use fearless_simd::{Level, dispatch, prelude::*};
use fearless_simd_macros::simd;
use thiserror::Error;

/// Represents a DNA base (A, C, G, T, or Unknown)
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, bytemuck::NoUninit)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u8)]
#[must_use]
pub enum Base {
    /// Adenine
    A = b'A',
    /// Cytosine
    C = b'C',
    /// Guanine
    G = b'G',
    /// Thymine
    T = b'T',
    /// Unknown base
    #[default]
    Unknown = b'N',
}

const _: () = {
    assert!(size_of::<Base>() == 1);
    assert!(size_of::<Option<Base>>() == 1);
};

#[cfg_attr(coverage_nightly, coverage(off))]
impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Table entry for a byte that is not one of the four canonical bases.
const NOT_KNOWN: u8 = u8::MAX;

/// `Base`'s discriminants are ASCII, so a 256-byte table indexed by the
/// discriminant covers every representable value. Entries other than A/C/G/T
/// are unreachable through the enum and map to [`NOT_KNOWN`] regardless.
const KNOWN_INDEX: [u8; 256] = {
    let mut table = [NOT_KNOWN; 256];
    table[Base::A as usize] = 0;
    table[Base::C as usize] = 1;
    table[Base::G as usize] = 2;
    table[Base::T as usize] = 3;
    table
};

impl Base {
    /// The four canonical DNA bases: A, C, G, T.
    pub const KNOWN: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

    /// Returns the index of this base in [`Base::KNOWN`], or `None` for `Unknown`.
    ///
    /// A table lookup rather than a `match`, because callers index per-base
    /// arrays with this in their innermost loop and the base is data: five
    /// match arms compile to a branch the predictor gets wrong every time the
    /// base changes, where the table is one always-resident load.
    #[inline]
    pub fn known_index(&self) -> Option<usize> {
        let slot = KNOWN_INDEX.get(*self as usize).copied().unwrap_or(NOT_KNOWN);
        (slot != NOT_KNOWN).then_some(slot as usize)
    }

    /// Get the inverse base (complementary base)
    pub fn inverse(&self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::Unknown => Base::Unknown,
        }
    }

    /// Get the string representation of the base
    pub fn as_str(&self) -> &'static str {
        match self {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
            Base::Unknown => "N",
        }
    }

    /// Get the character representation of the base
    pub fn as_char(&self) -> char {
        (*self) as u8 as char
    }

    // r[impl base_decode.ascii_batch]
    // r[impl base_decode.ascii_simd+2]
    /// Convert ASCII bytes to `Base` values in-place using SIMD acceleration.
    ///
    /// Reuses the input allocation — no new vector is allocated.
    pub fn from_ascii_vec(mut bytes: Vec<u8>) -> Vec<Base> {
        from_ascii_dispatch(&mut bytes);
        // Safety: from_ascii_dispatch has written only valid Base discriminants
        // (A=65, C=67, G=71, T=84, Unknown=78) into every byte.
        unsafe { Self::vec_u8_into_vec_base(bytes) }
    }

    /// Convert ASCII bytes to `Base` discriminants in-place (SIMD-accelerated).
    ///
    /// After this call, every byte in `bytes` is a valid `Base` discriminant.
    /// Use this when you need buffer reuse (the caller retains the Vec's capacity).
    pub fn convert_ascii_in_place(bytes: &mut [u8]) {
        from_ascii_dispatch(bytes);
    }

    /// Convert ASCII bytes to `Base` discriminants in-place and return a
    /// borrowed `&[Base]` view of the same memory.
    ///
    /// Unlike [`from_ascii_vec`], this does not consume the buffer: the
    /// caller keeps ownership of the underlying `Vec<u8>` (and its capacity).
    /// Use it when you want one allocation for repeated fetches and need
    /// to hand a `&[Base]` to downstream code (e.g. for `Rc::<[Base]>::from`).
    ///
    /// [`from_ascii_vec`]: Self::from_ascii_vec
    pub fn convert_ascii_in_place_as_slice(bytes: &mut [u8]) -> &[Base] {
        from_ascii_dispatch(bytes);
        // Safety: from_ascii_dispatch has written only valid Base discriminants
        // (A=65, C=67, G=71, T=84, Unknown=78). Base is repr(u8), so the layout
        // matches `[u8]` exactly. The returned reference borrows `bytes` for
        // the same lifetime, so the caller can't use the &mut [u8] view to
        // mutate the slice contents while the &[Base] is alive.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr().cast::<Base>(), bytes.len()) }
    }

    /// Reinterpret a `Vec<u8>` whose bytes are all valid `Base` discriminants.
    ///
    /// # Safety
    /// Every byte in `bytes` must be a valid `Base` discriminant:
    /// 65 (A), 67 (C), 71 (G), 84 (T), or 78 (Unknown).
    pub unsafe fn vec_u8_into_vec_base(bytes: Vec<u8>) -> Vec<Base> {
        let mut bytes = std::mem::ManuallyDrop::new(bytes);
        let ptr = bytes.as_mut_ptr() as *mut Base;
        let len = bytes.len();
        let cap = bytes.capacity();
        // Safety: Base is repr(u8) and all bytes are valid discriminants (caller
        // invariant). ptr came from a Vec<u8> so it satisfies the allocator,
        // alignment (1 == 1), and capacity requirements for Vec::from_raw_parts.
        unsafe { Vec::from_raw_parts(ptr, len, cap) }
    }
}

/// Dispatch to the best available ASCII→Base converter.
// r[impl base_decode.ascii_scalar_equivalence]
fn from_ascii_dispatch(bytes: &mut [u8]) {
    from_ascii_at(Level::new(), bytes);
}

// r[impl io.simd_portable]
fn from_ascii_at(level: Level, bytes: &mut [u8]) {
    dispatch!(level, simd => from_ascii_simd(simd, bytes));
}

/// One native-width vector per iteration (16 bytes on SSE/NEON, 32 on AVX2,
/// 64 on AVX-512):
/// 1. Uppercase: `byte & 0xDF` (clears bit 5, mapping a–z → A–Z)
/// 2. Compare uppercased byte against A(65), C(67), G(71), T(84)
/// 3. OR masks → `valid`
/// 4. Select: valid bytes keep their uppercased value, others become N(78)
///
/// The scalar tail goes through `BASE_LUT`, the same table `Base::from(u8)`
/// uses, so the two halves cannot disagree on a byte.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn from_ascii_simd<S: Simd>(simd: S, bytes: &mut [u8]) {
    let n = S::u8s::splat(simd, b'N');
    let mut chunks = bytes.chunks_exact_mut(S::u8s::LEN);
    for chunk in &mut chunks {
        let upper = S::u8s::from_slice(simd, chunk) & 0xDF;
        let valid =
            upper.simd_eq(b'A') | upper.simd_eq(b'C') | upper.simd_eq(b'G') | upper.simd_eq(b'T');
        valid.select(upper, n).store_slice(chunk);
    }
    from_ascii_scalar(chunks.into_remainder());
}

/// Scalar fallback: uses the existing `BASE_LUT` per byte.
#[allow(clippy::indexing_slicing, reason = "byte < 256 = BASE_LUT.len()")]
fn from_ascii_scalar(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        *b = BASE_LUT[*b as usize] as u8;
    }
}

impl std::ops::Deref for Base {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        match self {
            Base::A => &b'A',
            Base::C => &b'C',
            Base::G => &b'G',
            Base::T => &b'T',
            Base::Unknown => &b'N',
        }
    }
}

impl AsRef<str> for Base {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl PartialEq<str> for Base {
    fn eq(&self, other: &str) -> bool {
        *self.as_str() == *other
    }
}

impl PartialEq<char> for Base {
    fn eq(&self, other: &char) -> bool {
        self.as_char() == *other
    }
}

impl PartialEq<u8> for Base {
    fn eq(&self, other: &u8) -> bool {
        match self {
            Base::A => *other == b'A' || *other == b'a',
            Base::C => *other == b'C' || *other == b'c',
            Base::G => *other == b'G' || *other == b'g',
            Base::T => *other == b'T' || *other == b't',
            Base::Unknown => false,
        }
    }
}

impl PartialEq<Option<SmolStr>> for Base {
    fn eq(&self, other: &Option<SmolStr>) -> bool {
        if let Some(other) = other { other == self } else { false }
    }
}

impl PartialEq<Option<Base>> for Base {
    fn eq(&self, other: &Option<Base>) -> bool {
        if let Some(other) = other { other == self } else { false }
    }
}

impl PartialEq<Base> for Option<Base> {
    fn eq(&self, other: &Base) -> bool {
        if let Some(me) = self { other == me } else { false }
    }
}

impl PartialEq<Base> for SmolStr {
    fn eq(&self, other: &Base) -> bool {
        self.as_str() == other.as_str()
    }
}

impl From<Option<Base>> for Base {
    fn from(value: Option<Base>) -> Self {
        match value {
            Some(base) => base,
            None => Base::Unknown,
        }
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl std::fmt::Debug for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

// r[impl types.base.from_str_validation+2]
impl std::str::FromStr for Base {
    type Err = BaseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        let bytes = s.as_bytes();
        let Some(&first) = bytes.first() else {
            return Err(BaseError::Empty);
        };
        if bytes.len() > 1 {
            return Err(BaseError::MultipleChars);
        }
        match first | 0x20 {
            b'a' => Ok(Base::A),
            b'c' => Ok(Base::C),
            b'g' => Ok(Base::G),
            b't' => Ok(Base::T),
            b'n' => Ok(Base::Unknown),
            _ => Err(BaseError::InvalidBase(first)),
        }
    }
}

/// Lookup table for branchless u8 → Base conversion.
/// A single indexed load instead of a chain of comparisons.
#[allow(clippy::indexing_slicing, reason = "indices are byte literals, always < 256 = table len")]
static BASE_LUT: [Base; 256] = {
    let mut table = [Base::Unknown; 256];
    table[b'A' as usize] = Base::A;
    table[b'a' as usize] = Base::A;
    table[b'C' as usize] = Base::C;
    table[b'c' as usize] = Base::C;
    table[b'G' as usize] = Base::G;
    table[b'g' as usize] = Base::G;
    table[b'T' as usize] = Base::T;
    table[b't' as usize] = Base::T;
    table
};

impl From<u8> for Base {
    #[inline]
    fn from(value: u8) -> Self {
        #[allow(clippy::indexing_slicing, reason = "value is u8, always < 256 = BASE_LUT.len()")]
        BASE_LUT[value as usize]
    }
}

impl From<&str> for Base {
    fn from(val: &str) -> Self {
        match val {
            "A" | "a" => Base::A,
            "C" | "c" => Base::C,
            "G" | "g" => Base::G,
            "T" | "t" => Base::T,
            _ => Base::Unknown,
        }
    }
}

impl From<SmolStr> for Base {
    fn from(val: SmolStr) -> Self {
        val.as_str().into()
    }
}

impl From<&SmolStr> for Base {
    fn from(val: &SmolStr) -> Self {
        val.as_str().into()
    }
}

impl From<&u8> for Base {
    fn from(value: &u8) -> Self {
        Base::from(*value)
    }
}

impl From<Base> for SmolStr {
    fn from(val: Base) -> Self {
        SmolStr::new_inline(val.into())
    }
}

impl From<Base> for &'static str {
    fn from(val: Base) -> Self {
        val.as_str()
    }
}

// r[impl types.base.error_display]
/// Errors that can occur when parsing a [`Base`] from a string.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum BaseError {
    /// Input string was empty.
    #[error("Empty")]
    Empty,
    /// Input string contained more than one character after trimming.
    #[error("Multiple characters in base string")]
    MultipleChars,
    /// First byte is not a valid base (A/C/G/T/N, case-insensitive).
    #[error("Invalid base: 0x{0:02x}")]
    InvalidBase(u8),
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "test code with known small values")]
mod tests {
    use super::*;
    use hegel::prelude::*;
    use std::str::FromStr;

    /// The characters `Base::from_str` accepts.
    const BASE_CHARS: [char; 10] = ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'];

    // The index is a lookup table keyed on the ASCII discriminant, so the
    // property that matters is that *every* byte lands where the match arms
    // used to put it — including the ones that are not a `Base` at all and can
    // only be reached through a mis-transmute.
    #[test]
    fn known_index_agrees_with_the_variants_for_every_byte() {
        for byte in 0u8..=255 {
            let expected = match byte {
                b'A' | b'a' => Some(0),
                b'C' | b'c' => Some(1),
                b'G' | b'g' => Some(2),
                b'T' | b't' => Some(3),
                _ => None,
            };
            assert_eq!(Base::from(byte).known_index(), expected, "byte {byte}");
        }
        for (i, base) in Base::KNOWN.iter().enumerate() {
            assert_eq!(base.known_index(), Some(i), "{base} must index its own slot in KNOWN");
        }
        assert_eq!(Base::Unknown.known_index(), None);
    }

    #[hegel::test]
    fn from_u8_maps_correctly(tc: TestCase) {
        let input = tc.draw(gs::integers::<u8>());
        let base = Base::from(input);
        let expected = match input {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            _ => Base::Unknown,
        };
        assert_eq!(base, expected);
    }

    #[hegel::test]
    fn from_str_maps_to_correct_variant(tc: TestCase) {
        let input = tc.draw(gs::sampled_from(&BASE_CHARS));
        let base = Base::from_str(&input.to_string()).expect("valid single base char");
        let expected = match input {
            'A' | 'a' => Base::A,
            'C' | 'c' => Base::C,
            'G' | 'g' => Base::G,
            'T' | 't' => Base::T,
            'N' | 'n' => Base::Unknown,
            _ => unreachable!("BASE_CHARS holds only ACGTNacgtn"),
        };
        assert_eq!(base, expected);
    }

    #[hegel::test]
    fn from_str_rejects_multi_char(tc: TestCase) {
        let a = tc.draw(gs::sampled_from(&BASE_CHARS));
        let b = tc.draw(gs::sampled_from(&BASE_CHARS));
        let multi = format!("{a}{b}");
        assert!(Base::from_str(&multi).is_err());
    }

    // r[verify types.base.from_str_validation+2]
    /// The whole grammar the rule states, on arbitrary text: trim, then
    /// nothing is `Empty`, one byte is that base or `InvalidBase(byte)`, and
    /// anything longer is `MultipleChars`. The expected base comes from the
    /// rule's own list rather than the implementation's `| 0x20` case fold.
    #[hegel::test]
    fn from_str_matches_the_documented_grammar(tc: TestCase) {
        let input = tc.draw(gs::text().max_size(10));
        let parsed = Base::from_str(&input);

        match input.trim().as_bytes() {
            [] => assert!(matches!(parsed, Err(BaseError::Empty)), "{input:?} -> {parsed:?}"),
            [byte] => {
                let expected = match *byte {
                    b'A' | b'a' => Some(Base::A),
                    b'C' | b'c' => Some(Base::C),
                    b'G' | b'g' => Some(Base::G),
                    b'T' | b't' => Some(Base::T),
                    b'N' | b'n' => Some(Base::Unknown),
                    _ => None,
                };
                match expected {
                    Some(base) => {
                        assert!(matches!(parsed, Ok(b) if b == base), "{input:?} -> {parsed:?}");
                    }
                    None => assert!(
                        matches!(parsed, Err(BaseError::InvalidBase(b)) if b == *byte),
                        "{input:?} -> {parsed:?}, expected InvalidBase({byte:#04x})"
                    ),
                }
            }
            _ => {
                assert!(matches!(parsed, Err(BaseError::MultipleChars)), "{input:?} -> {parsed:?}")
            }
        }
    }

    // r[verify base_decode.ascii_simd+2]
    // r[verify base_decode.ascii_scalar_equivalence]
    #[hegel::test]
    fn from_ascii_vec_equivalence(tc: TestCase) {
        let input = tc.draw(gs::binary());
        let expected: Vec<Base> = input.iter().map(|&b| Base::from(b)).collect();
        let actual = Base::from_ascii_vec(input);
        assert_eq!(actual, expected);
    }

    /// Every SIMD level this CPU can run, the scalar fallback included, so
    /// each one is tested — not only the one `Level::new()` picks.
    fn levels() -> Vec<Level> {
        let best = Level::new();
        #[allow(unused_mut, reason = "only x86 has levels below the best one")]
        let mut levels = vec![Level::fallback(), best];
        #[cfg(target_arch = "x86_64")]
        {
            levels.extend(best.as_sse2().map(Level::Sse2));
            levels.extend(best.as_sse4_2().map(Level::Sse4_2));
            levels.extend(best.as_avx2().map(Level::Avx2));
        }
        levels
    }

    // r[verify base_decode.ascii_simd+2]
    // r[verify base_decode.ascii_scalar_equivalence]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn from_ascii_every_level_matches_scalar(tc: TestCase) {
        let input = tc.draw(gs::binary().max_size(300));
        let expected: Vec<u8> = input.iter().map(|&b| Base::from(b) as u8).collect();
        for level in levels() {
            let mut actual = input.clone();
            from_ascii_at(level, &mut actual);
            assert_eq!(actual, expected, "{level:?}");
        }
    }

    // r[verify base_decode.ascii_scalar_equivalence]
    #[hegel::test]
    fn from_ascii_vec_all_byte_values(tc: TestCase) {
        // All 256 byte values, but with varying-length padding to exercise
        // different SIMD alignment boundaries and tail lengths.
        let pad = || gs::integers::<usize>().min_value(0).max_value(31);
        let prefix_len = tc.draw(pad());
        let suffix_len = tc.draw(pad());
        let all_bytes: Vec<u8> = (0u16..=255).map(|b| b as u8).collect();
        let prefix: Vec<u8> = vec![b'N'; prefix_len];
        let suffix: Vec<u8> = vec![b'A'; suffix_len];
        let input: Vec<u8> = [&prefix[..], &all_bytes, &suffix[..]].concat();
        let expected: Vec<Base> = input.iter().map(|&b| Base::from(b)).collect();
        let actual = Base::from_ascii_vec(input);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_u8_to_base_valid() {
        let valid_bases = *b"ACGT";
        for &base in &valid_bases {
            let parsed = Base::from(base);
            assert_eq!(*parsed, base);
            // display
            assert_eq!(parsed.to_string(), (base as char).to_string());
            // debug -- same, actually
            assert_eq!(format!("{parsed:#?}"), (base as char).to_string());
        }
    }

    // r[verify base_decode.ascii_batch]
    #[test]
    fn test_from_ascii_vec_basic() {
        let input = b"ACGTacgtNnMR".to_vec();
        let result = Base::from_ascii_vec(input);
        assert_eq!(
            result,
            vec![
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::Unknown,
                Base::Unknown,
                Base::Unknown,
                Base::Unknown,
            ]
        );
    }

    #[test]
    fn test_from_ascii_vec_empty() {
        assert_eq!(Base::from_ascii_vec(vec![]), Vec::<Base>::new());
    }

    // r[verify base_decode.ascii_scalar_equivalence]
    // r[verify base_decode.ascii_simd+2]
    #[test]
    fn test_from_ascii_vec_simd_boundary_lengths() {
        for len in [0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 128] {
            let input: Vec<u8> = (0..len).map(|i| b"ACGTNacgtn"[i % 10]).collect();
            let expected: Vec<Base> = input.iter().map(|&b| Base::from(b)).collect();
            let actual = Base::from_ascii_vec(input);
            assert_eq!(actual, expected, "mismatch at len={len}");
        }
    }

    // r[verify base_decode.ascii_batch]
    #[test]
    fn test_convert_ascii_in_place_matches_from_ascii_vec() {
        let input = b"ACGTacgtNnMRWSYKVHDB\x00\xff".to_vec();
        let expected = Base::from_ascii_vec(input.clone());

        let mut buf = input;
        Base::convert_ascii_in_place(&mut buf);
        // After in-place conversion, raw bytes must match from_ascii_vec output
        assert_eq!(buf.len(), expected.len());
        for (i, (&byte, &base)) in buf.iter().zip(expected.iter()).enumerate() {
            assert_eq!(byte, base as u8, "mismatch at index {i}: byte={byte:#x}, base={base:?}");
        }
    }

    // r[verify base_decode.ascii_batch]
    /// `convert_ascii_in_place_as_slice` must yield the same bases as
    /// `from_ascii_vec`, while leaving the caller's `Vec<u8>` capacity intact
    /// (the whole point of the helper).
    #[test]
    fn test_convert_ascii_in_place_as_slice_matches_from_ascii_vec() {
        let input = b"ACGTacgtNnMRWSYKVHDB\x00\xff".to_vec();
        let expected = Base::from_ascii_vec(input.clone());

        let mut buf = input;
        let cap_before = buf.capacity();
        let view: &[Base] = Base::convert_ascii_in_place_as_slice(&mut buf);
        assert_eq!(view, expected.as_slice());
        // Caller still owns the buffer with its capacity preserved.
        let cap_after = buf.capacity();
        assert_eq!(cap_before, cap_after);
        // Length is unchanged — only the byte values were rewritten.
        assert_eq!(buf.len(), expected.len());
    }

    // r[verify types.base.from_str_validation+2]
    #[test]
    fn test_from_str_rejects_invalid_base() {
        assert!("Z".parse::<Base>().is_err());
        assert!("1".parse::<Base>().is_err());
        assert!("M".parse::<Base>().is_err());
        assert!("".parse::<Base>().is_err());
    }

    // r[verify types.base.from_str_validation+2]
    #[test]
    fn test_from_str_rejects_multi_char() {
        assert!(matches!("AT".parse::<Base>(), Err(BaseError::MultipleChars)));
        assert!(matches!("AC".parse::<Base>(), Err(BaseError::MultipleChars)));
        assert!(matches!("NN".parse::<Base>(), Err(BaseError::MultipleChars)));
        // Trimmed single char still works
        assert_eq!(" A ".parse::<Base>().unwrap(), Base::A);
    }

    // r[verify types.base.from_str_validation+2]
    #[test]
    fn test_from_str_accepts_valid_bases() {
        assert_eq!("A".parse::<Base>().unwrap(), Base::A);
        assert_eq!("a".parse::<Base>().unwrap(), Base::A);
        assert_eq!("C".parse::<Base>().unwrap(), Base::C);
        assert_eq!("c".parse::<Base>().unwrap(), Base::C);
        assert_eq!("G".parse::<Base>().unwrap(), Base::G);
        assert_eq!("g".parse::<Base>().unwrap(), Base::G);
        assert_eq!("T".parse::<Base>().unwrap(), Base::T);
        assert_eq!("t".parse::<Base>().unwrap(), Base::T);
        assert_eq!("N".parse::<Base>().unwrap(), Base::Unknown);
        assert_eq!("n".parse::<Base>().unwrap(), Base::Unknown);
    }

    // r[verify types.base.error_display]
    #[test]
    fn test_base_error_display_shows_actual_byte() {
        let err = "Z".parse::<Base>().unwrap_err();
        assert_eq!(err.to_string(), "Invalid base: 0x5a");

        let err = "".parse::<Base>().unwrap_err();
        assert_eq!(err.to_string(), "Empty");

        let err = "AT".parse::<Base>().unwrap_err();
        assert_eq!(err.to_string(), "Multiple characters in base string");
    }

    // r[verify base_decode.ascii_batch]
    #[test]
    fn test_vec_u8_into_vec_base_roundtrip() {
        let discriminants: Vec<u8> = vec![b'A', b'C', b'G', b'T', b'N', b'A', b'N', b'T'];
        // Safety: this is the test
        let bases = unsafe { Base::vec_u8_into_vec_base(discriminants) };
        assert_eq!(
            bases,
            vec![
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::Unknown,
                Base::A,
                Base::Unknown,
                Base::T
            ]
        );
    }
}
