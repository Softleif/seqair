//! The fixture the bench and the profiling example share: one 150 bp read cut
//! from a 200 bp haplotype, with a realistic `CpG` density so the TAPS arm
//! exercises the conversion rows rather than falling through everywhere.

use compair::{Band, Base, BaseQuality, Haplotype, Probability, Read, Strand};

pub const HAPLOTYPE_LEN: usize = 200;
pub const READ_LEN: usize = 150;
pub const READ_OFFSET: usize = 25;

/// How many haplotypes a caller scores one read against.
pub const HAPLOTYPES_PER_READ: usize = 8;

/// A deterministic pseudo-random sequence.
pub fn sequence(len: usize, seed: u64) -> Vec<Base> {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            match state % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            }
        })
        .collect()
}

pub struct Fixture {
    /// The same locus on `HAPLOTYPES_PER_READ` candidate haplotypes: the
    /// first is the one the read was cut from, each next one carries one more
    /// substitution.
    pub haplotypes: Vec<Haplotype>,
    pub read: Read,
    pub betas: Vec<Probability>,
    pub band: Band,
}

pub fn fixture() -> Fixture {
    let bases = sequence(HAPLOTYPE_LEN, 0x9e37_79b9_7f4a_7c15);
    let read_bases: Vec<Base> =
        bases.get(READ_OFFSET..READ_OFFSET + READ_LEN).map(<[Base]>::to_vec).unwrap_or_default();
    let quals: Vec<BaseQuality> = (0..read_bases.len())
        .map(|index| BaseQuality::from_byte(25 + u8::try_from(index % 15).unwrap_or(0)))
        .collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("the bench fixture must build");
    let betas = (0..HAPLOTYPE_LEN)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    let band = Band::anchored(i32::try_from(READ_OFFSET).unwrap_or(0));
    let haplotypes = (0..HAPLOTYPES_PER_READ)
        .map(|variant| {
            let mut alternative = bases.clone();
            for edit in 0..variant {
                let at = READ_OFFSET + 10 + edit * 17;
                if let Some(base) = alternative.get_mut(at) {
                    *base = base.inverse();
                }
            }
            Haplotype::new(alternative)
        })
        .collect();
    Fixture { haplotypes, read, betas, band }
}
