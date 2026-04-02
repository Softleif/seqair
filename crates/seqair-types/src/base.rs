use crate::SmolStr;
use thiserror::Error;

/// Represents a DNA base (A, C, G, T, or Unknown)
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
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

impl Base {
    /// The four canonical DNA bases: A, C, G, T.
    pub const KNOWN: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

    /// Returns the index of this base in [`Base::KNOWN`], or `None` for `Unknown`.
    pub fn known_index(&self) -> Option<usize> {
        match self {
            Base::A => Some(0),
            Base::C => Some(1),
            Base::G => Some(2),
            Base::T => Some(3),
            Base::Unknown => None,
        }
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
    // r[impl base_decode.ascii_simd]
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
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // Safety: AVX2 verified by feature detection.
            unsafe { from_ascii_avx2(bytes) };
            return;
        }
        if is_x86_feature_detected!("ssse3") {
            // Safety: SSSE3 verified by feature detection.
            unsafe { from_ascii_ssse3(bytes) };
            return;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON is always available on aarch64.
        unsafe { from_ascii_neon(bytes) };
        return;
    }

    #[allow(unreachable_code)]
    from_ascii_scalar(bytes);
}

/// Scalar fallback: uses the existing BASE_LUT per byte.
#[allow(clippy::indexing_slicing, reason = "byte < 256 = BASE_LUT.len()")]
fn from_ascii_scalar(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        *b = BASE_LUT[*b as usize] as u8;
    }
}

/// AVX2 implementation: processes 32 bytes per iteration.
///
/// Algorithm per chunk (same as SSSE3, using 256-bit registers):
/// 1. Uppercase: `byte & 0xDF` (clears bit 5, mapping a–z → A–Z)
/// 2. Compare uppercased byte against A(65), C(67), G(71), T(84)
/// 3. OR masks → `valid` (0xFF where ACGT, 0x00 elsewhere)
/// 4. Select: valid bytes keep their uppercased value, others become N(78)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn from_ascii_avx2(bytes: &mut [u8]) {
    use std::arch::x86_64::*;

    // Safety: AVX2 intrinsics require AVX2 availability (ensured by #[target_feature]).
    let (upper_mask, v_a, v_c, v_g, v_t, v_n) = (
        _mm256_set1_epi8(0xDFu8 as i8),
        _mm256_set1_epi8(b'A' as i8),
        _mm256_set1_epi8(b'C' as i8),
        _mm256_set1_epi8(b'G' as i8),
        _mm256_set1_epi8(b'T' as i8),
        _mm256_set1_epi8(b'N' as i8),
    );

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i = 0;

    while i + 32 <= len {
        // Safety: pointer ops below are valid because `i + 32 <= len` guarantees
        // we never read/write past the end of the slice.
        unsafe {
            let chunk = _mm256_loadu_si256(ptr.add(i) as *const __m256i);
            let upper = _mm256_and_si256(chunk, upper_mask);

            let is_a = _mm256_cmpeq_epi8(upper, v_a);
            let is_c = _mm256_cmpeq_epi8(upper, v_c);
            let is_g = _mm256_cmpeq_epi8(upper, v_g);
            let is_t = _mm256_cmpeq_epi8(upper, v_t);

            let valid = _mm256_or_si256(_mm256_or_si256(is_a, is_c), _mm256_or_si256(is_g, is_t));
            let result =
                _mm256_or_si256(_mm256_and_si256(valid, upper), _mm256_andnot_si256(valid, v_n));

            _mm256_storeu_si256(ptr.add(i) as *mut __m256i, result);
        }
        i += 32;
    }

    debug_assert!(i <= bytes.len(), "AVX2 loop overshot: i={i}, len={}", bytes.len());
    #[allow(
        clippy::indexing_slicing,
        reason = "i..len is in bounds, byte value < 256 = BASE_LUT.len()"
    )]
    for b in &mut bytes[i..] {
        *b = BASE_LUT[*b as usize] as u8;
    }
}

/// SSSE3 implementation: processes 16 bytes per iteration.
///
/// Algorithm per chunk:
/// 1. Uppercase: `byte & 0xDF` (clears bit 5, mapping a–z → A–Z)
/// 2. Compare uppercased byte against A(65), C(67), G(71), T(84)
/// 3. OR masks → `valid` (0xFF where ACGT, 0x00 elsewhere)
/// 4. Select: valid bytes keep their uppercased value, others become N(78)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn from_ascii_ssse3(bytes: &mut [u8]) {
    use std::arch::x86_64::*;

    // Safety: SSE2 intrinsics require SSSE3 availability (ensured by #[target_feature]).
    let (upper_mask, v_a, v_c, v_g, v_t, v_n) = (
        _mm_set1_epi8(0xDFu8 as i8),
        _mm_set1_epi8(b'A' as i8),
        _mm_set1_epi8(b'C' as i8),
        _mm_set1_epi8(b'G' as i8),
        _mm_set1_epi8(b'T' as i8),
        _mm_set1_epi8(b'N' as i8),
    );

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i = 0;

    while i + 16 <= len {
        // Safety: pointer ops below are valid because `i + 16 <= len` guarantees
        // we never read/write past the end of the slice.
        unsafe {
            let chunk = _mm_loadu_si128(ptr.add(i) as *const __m128i);
            let upper = _mm_and_si128(chunk, upper_mask);

            let is_a = _mm_cmpeq_epi8(upper, v_a);
            let is_c = _mm_cmpeq_epi8(upper, v_c);
            let is_g = _mm_cmpeq_epi8(upper, v_g);
            let is_t = _mm_cmpeq_epi8(upper, v_t);

            let valid = _mm_or_si128(_mm_or_si128(is_a, is_c), _mm_or_si128(is_g, is_t));
            let result = _mm_or_si128(_mm_and_si128(valid, upper), _mm_andnot_si128(valid, v_n));

            _mm_storeu_si128(ptr.add(i) as *mut __m128i, result);
        }
        i += 16;
    }

    debug_assert!(i <= bytes.len(), "SSSE3 loop overshot: i={i}, len={}", bytes.len());
    #[allow(
        clippy::indexing_slicing,
        reason = "i..len is in bounds, byte value < 256 = BASE_LUT.len()"
    )]
    for b in &mut bytes[i..] {
        *b = BASE_LUT[*b as usize] as u8;
    }
}

/// NEON implementation: processes 16 bytes per iteration.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn from_ascii_neon(bytes: &mut [u8]) {
    use std::arch::aarch64::*;

    let upper_mask = vdupq_n_u8(0xDF);
    let v_a = vdupq_n_u8(b'A');
    let v_c = vdupq_n_u8(b'C');
    let v_g = vdupq_n_u8(b'G');
    let v_t = vdupq_n_u8(b'T');
    let v_n = vdupq_n_u8(b'N');

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i: usize = 0;

    while i.wrapping_add(16) <= len {
        // Safety: pointer ops below are valid because `i + 16 <= len` guarantees
        // we never read/write past the end of the slice.
        unsafe {
            let chunk = vld1q_u8(ptr.add(i));
            let upper = vandq_u8(chunk, upper_mask);

            let is_a = vceqq_u8(upper, v_a);
            let is_c = vceqq_u8(upper, v_c);
            let is_g = vceqq_u8(upper, v_g);
            let is_t = vceqq_u8(upper, v_t);

            let valid = vorrq_u8(vorrq_u8(is_a, is_c), vorrq_u8(is_g, is_t));
            let result = vbslq_u8(valid, upper, v_n);

            vst1q_u8(ptr.add(i), result);
        }
        i = i.wrapping_add(16);
    }

    debug_assert!(i <= bytes.len(), "NEON loop overshot: i={i}, len={}", bytes.len());
    #[allow(
        clippy::indexing_slicing,
        reason = "i..len is in bounds, byte value < 256 = BASE_LUT.len()"
    )]
    for b in &mut bytes[i..] {
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
mod tests {
    use super::*;
    use std::str::FromStr;

    proptest::proptest! {
        #[test]
        fn proptest_from_u8_maps_correctly(input: u8) {
            let base = Base::from(input);
            let expected = match input {
                b'A' | b'a' => Base::A,
                b'C' | b'c' => Base::C,
                b'G' | b'g' => Base::G,
                b'T' | b't' => Base::T,
                _ => Base::Unknown,
            };
            proptest::prop_assert_eq!(base, expected);
        }

        #[test]
        fn proptest_from_str_maps_to_correct_variant(input in "[ACGTNacgtn]") {
            let base = Base::from_str(&input).expect("valid single base char");
            let ch = input.chars().next().expect("non-empty string");
            let expected = match ch {
                'A' | 'a' => Base::A,
                'C' | 'c' => Base::C,
                'G' | 'g' => Base::G,
                'T' | 't' => Base::T,
                'N' | 'n' => Base::Unknown,
                _ => unreachable!("regex guarantees only ACGTNacgtn"),
            };
            proptest::prop_assert_eq!(base, expected);
        }

        #[test]
        fn proptest_from_str_rejects_multi_char(a in "[ACGTNacgtn]", b in "[ACGTNacgtn]") {
            let multi = format!("{a}{b}");
            proptest::prop_assert!(Base::from_str(&multi).is_err());
        }

        #[test]
        fn proptest_from_str_never_panics(input in r"\PC{0,10}") {
            // Just ensure no panic; errors are fine
            let _ = Base::from_str(&input);
        }

        // r[verify base_decode.ascii_simd]
        // r[verify base_decode.ascii_scalar_equivalence]
        #[test]
        fn proptest_from_ascii_vec_equivalence(input: Vec<u8>) {
            let expected: Vec<Base> = input.iter().map(|&b| Base::from(b)).collect();
            let actual = Base::from_ascii_vec(input);
            proptest::prop_assert_eq!(actual, expected);
        }

        // r[verify base_decode.ascii_scalar_equivalence]
        #[test]
        fn proptest_from_ascii_vec_all_byte_values(prefix_len in 0usize..32, suffix_len in 0usize..32) {
            // All 256 byte values, but with varying-length padding to exercise
            // different SIMD alignment boundaries and tail lengths.
            let all_bytes: Vec<u8> = (0u16..=255).map(|b| b as u8).collect();
            let prefix: Vec<u8> = vec![b'N'; prefix_len];
            let suffix: Vec<u8> = vec![b'A'; suffix_len];
            let input: Vec<u8> = [&prefix[..], &all_bytes, &suffix[..]].concat();
            let expected: Vec<Base> = input.iter().map(|&b| Base::from(b)).collect();
            let actual = Base::from_ascii_vec(input);
            proptest::prop_assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_u8_to_base_valid() {
        let valid_bases = [b'A', b'C', b'G', b'T'];
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
    // r[verify base_decode.ascii_simd]
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
