//! The fixture the bench and the profiling example share: one 150 bp read cut
//! from a 200 bp haplotype, with a realistic `CpG` density so the TAPS arm
//! exercises the conversion rows rather than falling through everywhere.
//! And a shadow-scoring locus, the shape a variant caller scores many reads
//! against a few haplotypes in.

#![allow(dead_code, reason = "each consumer uses a different part of this")]

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

/// Reads per locus in the shadow fixture.
pub const SHADOW_READS: usize = 64;

/// The shadow fixture's band width: rastair's `max(46, 2 * longest allele +
/// 32)` rounded to eight, for a short indel.
pub const SHADOW_WIDTH: u32 = 48;

/// How many haplotypes the shadow fixture holds: the reference, an insertion
/// allele and a deletion allele.
pub const SHADOW_HAPLOTYPES: usize = 3;

/// One indel candidate's locus as a variant caller's shadow path scores it:
/// every read against every haplotype, each read through its own band.
pub struct Shadow {
    pub haplotypes: Vec<Haplotype>,
    /// 150 bp reads from both strands, each cut at the locus's anchor give or
    /// take ten bases and banded at the place it was cut from.
    pub reads: Vec<(Read, Band)>,
    pub betas: Vec<Probability>,
}

pub fn shadow() -> Shadow {
    let reference = sequence(HAPLOTYPE_LEN, 0x2545_f491_4f6c_dd1d);
    let mut insertion = reference.clone();
    insertion.splice(HAPLOTYPE_LEN / 2..HAPLOTYPE_LEN / 2, [Base::A, Base::C]);
    let mut deletion = reference.clone();
    deletion.drain(HAPLOTYPE_LEN / 2..HAPLOTYPE_LEN / 2 + 3);
    let haplotypes =
        [reference.clone(), insertion, deletion].into_iter().map(Haplotype::new).collect();

    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let reads = (0..SHADOW_READS)
        .map(|_| {
            let jitter = usize::try_from(next() % 21).unwrap_or(10);
            let start = READ_OFFSET - 10 + jitter;
            let mut bases: Vec<Base> =
                reference.get(start..start + READ_LEN).map(<[Base]>::to_vec).unwrap_or_default();
            let at = usize::try_from(next() % 150).unwrap_or(0);
            if let Some(base) = bases.get_mut(at) {
                *base = base.inverse();
            }
            let quals: Vec<BaseQuality> = (0..READ_LEN)
                .map(|_| BaseQuality::from_byte(20 + u8::try_from(next() % 21).unwrap_or(0)))
                .collect();
            let strand = if next() % 2 == 0 { Strand::OT } else { Strand::OB };
            let read = Read::uniform(
                bases,
                &quals,
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(10),
                strand,
            )
            .expect("the shadow fixture must build");
            let band = Band::new(SHADOW_WIDTH, i32::try_from(start).unwrap_or(0))
                .expect("the shadow width is valid");
            (read, band)
        })
        .collect();
    let betas = (0..HAPLOTYPE_LEN + 2)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    Shadow { haplotypes, reads, betas }
}
