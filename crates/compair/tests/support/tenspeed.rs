//! The "10s" PairHMM benchmark, parsed.
//!
//! This is the dataset the gpuPairHMM and Endeavor papers both report times on,
//! so it is the one place compair's numbers can be put next to a published
//! figure without re-deriving anyone's hardware. `tests/data/PROVENANCE.md` says
//! where the files came from.
//!
//! Shared by `tests/tenspeed.rs` and `benches/tenspeed.rs` — unlike the bench
//! fixtures, which are deliberately copied, this is a parser for a fixed file
//! and there is nothing here for one consumer to tune out from under the other.

#![allow(dead_code, reason = "each consumer uses a different part of this")]

use compair::{Band, Base, BaseQuality, Haplotype, Read, Strand};

/// GATK floors base qualities at 6 before the DP; the shipped file is already
/// clamped, so this only documents the convention.
const MIN_BASE_QUALITY: u8 = 6;

/// One read-haplotype pair.
pub struct Pair {
    pub haplotype: Haplotype,
    pub read: Read,
    /// The haplotype offset the read's first base is expected near, for the
    /// band. A variant caller gets this from the alignment; here it is seeded.
    pub offset: i32,
    /// The raw bytes, kept so the gkl arms can be handed exactly what gkl's own
    /// harness hands them.
    pub raw: Raw,
}

/// A pair as gkl's API wants it: ASCII bases and phred bytes.
pub struct Raw {
    pub haplotype: Vec<u8>,
    pub read: Vec<u8>,
    pub base_quals: Vec<u8>,
    pub insertion_quals: Vec<u8>,
    pub deletion_quals: Vec<u8>,
    pub gap_quals: Vec<u8>,
}

#[derive(Debug)]
pub enum Error {
    /// A group header that was not two counts.
    Header {
        line: usize,
    },
    /// A read line that was not bases plus four quality tracks.
    Read {
        line: usize,
    },
    /// The file ended while a group still wanted rows.
    Truncated {
        line: usize,
    },
    Compair(compair::Error),
}

impl From<compair::Error> for Error {
    fn from(error: compair::Error) -> Self {
        Self::Compair(error)
    }
}

const INPUT: &str = include_str!("../data/10s.in");

/// Every pair in the benchmark: group by group, and within a group read-major,
/// haplotype-minor.
pub fn pairs() -> Result<Vec<Pair>, Error> {
    let mut out = Vec::new();
    let mut lines = INPUT.lines().enumerate();
    while let Some((index, header)) = lines.next() {
        if header.trim().is_empty() {
            continue;
        }
        let mut counts = header.split_whitespace();
        let [Some(reads), Some(haplotypes)] = [counts.next(), counts.next()] else {
            return Err(Error::Header { line: index + 1 });
        };
        let [reads, haplotypes] = [reads, haplotypes]
            .map(|count| count.parse::<usize>().map_err(|_| Error::Header { line: index + 1 }));
        let (reads, haplotypes) = (reads?, haplotypes?);

        let mut group_reads = Vec::with_capacity(reads);
        for _ in 0..reads {
            let (number, line) = lines.next().ok_or(Error::Truncated { line: index + 1 })?;
            group_reads.push(parse_read(line, number + 1)?);
        }
        let mut group_haplotypes = Vec::with_capacity(haplotypes);
        for _ in 0..haplotypes {
            let (_, line) = lines.next().ok_or(Error::Truncated { line: index + 1 })?;
            let bases = line.split_whitespace().next().unwrap_or("").as_bytes().to_vec();
            group_haplotypes.push(bases);
        }

        for raw_read in &group_reads {
            for raw_haplotype in &group_haplotypes {
                let haplotype = Haplotype::from_ascii(raw_haplotype);
                let read = Read::new(
                    raw_read.read.iter().copied().map(Base::from).collect::<Vec<_>>(),
                    &quals(&raw_read.base_quals),
                    &quals(&raw_read.insertion_quals),
                    &quals(&raw_read.deletion_quals),
                    &quals(&raw_read.gap_quals),
                    Strand::OT,
                )?;
                let offset = seed_offset(&haplotype, &read);
                out.push(Pair {
                    haplotype,
                    read,
                    offset,
                    raw: Raw {
                        haplotype: raw_haplotype.clone(),
                        read: raw_read.read.clone(),
                        base_quals: raw_read.base_quals.clone(),
                        insertion_quals: raw_read.insertion_quals.clone(),
                        deletion_quals: raw_read.deletion_quals.clone(),
                        gap_quals: raw_read.gap_quals.clone(),
                    },
                });
            }
        }
    }

    Ok(out)
}

fn parse_read(line: &str, number: usize) -> Result<Raw, Error> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    let [read, base, insertion, deletion, gap] = fields[..] else {
        return Err(Error::Read { line: number });
    };
    Ok(Raw {
        haplotype: Vec::new(),
        read: read.as_bytes().to_vec(),
        base_quals: phred(base, MIN_BASE_QUALITY),
        insertion_quals: phred(insertion, 0),
        deletion_quals: phred(deletion, 0),
        gap_quals: phred(gap, 0),
    })
}

/// ASCII phred+33 to raw phred, floored at `min` the way GATK floors base
/// qualities.
fn phred(track: &str, min: u8) -> Vec<u8> {
    track.bytes().map(|byte| byte.saturating_sub(33).max(min)).collect()
}

fn quals(track: &[u8]) -> Vec<BaseQuality> {
    track.iter().copied().map(BaseQuality::from_byte).collect()
}

/// The haplotype offset at which the read has the most matching bases.
///
/// The band needs an anchor and a variant caller has one from the alignment;
/// this stands in for it, and is never inside anything being timed.
#[allow(
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "scaffolding over sequences of a few hundred bases"
)]
fn seed_offset(haplotype: &Haplotype, read: &Read) -> i32 {
    let (h, r) = (haplotype.len() as i64, read.len() as i64);
    let mut best = (0i64, -1i64);
    for offset in -(r - 1)..h {
        let mut hits = 0i64;
        for (index, base) in read.bases().iter().enumerate() {
            let column = offset + index as i64;
            if column >= 0 && haplotype.bases().get(column as usize) == Some(base) {
                hits += 1;
            }
        }
        if hits > best.1 {
            best = (offset, hits);
        }
    }
    i32::try_from(best.0).unwrap_or(0)
}

/// How many matrix cells the full DP visits: the sum of `read x haplotype` over
/// every pair. This is the denominator both papers divide by to get CUPS, so it
/// is the one that makes our numbers comparable to theirs.
pub fn matrix_cells(pairs: &[Pair]) -> u64 {
    pairs.iter().map(|pair| pair.read.len() as u64 * pair.haplotype.len() as u64).sum()
}

/// How many cells a band of this width actually visits, counted the way the
/// band defines itself: `|j - i - offset| <= width / 2`, clipped to the matrix.
///
/// The gap between this and [`matrix_cells`] is the whole point of banding, and
/// reporting CUPS against the wrong one of the two is how a banded kernel gets
/// to claim a speedup it did not earn.
pub fn band_cells(pairs: &[Pair], band: impl Fn(&Pair) -> Band) -> u64 {
    pairs
        .iter()
        .map(|pair| {
            let band = band(pair);
            let half = i64::from(band.width() / 2);
            let offset = i64::from(band.offset());
            let (rows, columns) = (pair.read.len() as i64, pair.haplotype.len() as i64);
            (0..rows)
                .map(|row| {
                    let low = (row + offset - half).max(0);
                    let high = (row + offset + half).min(columns - 1);
                    (high - low + 1).max(0) as u64
                })
                .sum::<u64>()
        })
        .sum()
}
