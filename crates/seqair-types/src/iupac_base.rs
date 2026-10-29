use crate::Base;

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
    #[allow(clippy::indexing_slicing, reason = "const indices: i < 16, bytes < 256")]
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
    pub fn complement(self) -> Self {
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

    /// [`Self::from_ascii`] over a slice.
    pub fn from_ascii_slice(bytes: &[u8]) -> Vec<Self> {
        bytes.iter().map(|&b| Self::from_ascii(b)).collect()
    }

    /// [`Self::to_ascii`] over a slice.
    pub fn to_ascii_slice(bases: &[Self]) -> Vec<u8> {
        bases.iter().map(|b| b.to_ascii()).collect()
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
