//! What the same 150 bp read against a 200 bp haplotype costs through
//! rust-bio's `PairHMM`, so the README's compair numbers have something
//! external to sit next to.
//!
//! The two arms are the same problem, not the same shape of work: rust-bio
//! runs the whole `201 x 150` matrix in log space, with a `ln_1p` and a
//! polynomial `exp` per state per cell, while `align_full` runs it in `f64`
//! with a power-of-two renormalisation per row. `tests/rust_bio.rs` proves
//! they compute the same number; this measures what each pays for it.
//!
//! The fixture is `benches/align.rs`'s, copied rather than shared so that
//! neither bench can quietly change the other's numbers.

use bio::stats::{
    LogProb, Prob,
    pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission},
};
use compair::Band;
use compair::{
    Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, align_banded_simd, align_full,
};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const HAPLOTYPE_LEN: usize = 200;
const READ_LEN: usize = 150;
const READ_OFFSET: usize = 25;
const INSERTION_QUAL: u8 = 45;
const DELETION_QUAL: u8 = 45;
const GAP_QUAL: u8 = 10;

/// A deterministic sequence with a realistic `CpG` density.
fn sequence(len: usize) -> Vec<Base> {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
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

fn fixture() -> (Haplotype, Read) {
    let bases = sequence(HAPLOTYPE_LEN);
    let read_bases: Vec<Base> = bases
        .get(READ_OFFSET..READ_OFFSET.saturating_add(READ_LEN))
        .map(<[Base]>::to_vec)
        .unwrap_or_default();
    let quals: Vec<BaseQuality> = (0..read_bases.len())
        .map(|index| BaseQuality::from_byte(25 + u8::try_from(index % 15).unwrap_or(0)))
        .collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(INSERTION_QUAL),
        BaseQuality::from_byte(DELETION_QUAL),
        BaseQuality::from_byte(GAP_QUAL),
        Strand::OT,
    )
    .expect("the bench fixture must build");
    (Haplotype::new(bases), read)
}

fn error_probability(phred: u8) -> f64 {
    10f64.powf(-f64::from(phred) / 10.0)
}

/// The mapping `tests/rust_bio.rs` documents: `x` is the haplotype, with one
/// unemittable base prepended, and `y` is the read.
struct Emissions<'a> {
    haplotype: &'a [Base],
    read: &'a [Base],
    epsilon: &'a [f64],
}

impl EmissionParameters for Emissions<'_> {
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        let Some(&hap) = i.checked_sub(1).and_then(|index| self.haplotype.get(index)) else {
            return XYEmission::Mismatch(LogProb::ln_zero());
        };
        let (Some(&base), Some(&eps)) = (self.read.get(j), self.epsilon.get(j)) else {
            return XYEmission::Mismatch(LogProb::ln_zero());
        };
        if hap == Base::Unknown || base == Base::Unknown || hap == base {
            XYEmission::Match(LogProb::from(Prob(1.0 - eps)))
        } else {
            XYEmission::Mismatch(LogProb::from(Prob(eps / 3.0)))
        }
    }

    fn prob_emit_x(&self, _i: usize) -> LogProb {
        LogProb::ln_one()
    }

    fn prob_emit_y(&self, _j: usize) -> LogProb {
        LogProb::ln_one()
    }

    fn len_x(&self) -> usize {
        self.haplotype.len().saturating_add(1)
    }

    fn len_y(&self) -> usize {
        self.read.len()
    }
}

struct Gaps {
    insertion: f64,
    deletion: f64,
    continuation: f64,
}

impl GapParameters for Gaps {
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(Prob(self.insertion))
    }
    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(Prob(self.deletion))
    }
    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::from(Prob(self.continuation))
    }
    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::from(Prob(self.continuation))
    }
}

struct Mode {
    start: LogProb,
}

impl StartEndGapParameters for Mode {
    fn prob_start_gap_x(&self, _i: usize) -> LogProb {
        self.start
    }
    fn free_start_gap_x(&self) -> bool {
        true
    }
    fn free_end_gap_x(&self) -> bool {
        true
    }
}

fn align(c: &mut Criterion) {
    let (haplotype, read) = fixture();
    let epsilon: Vec<f64> =
        read.base_quals().iter().map(|q| error_probability(q.as_byte())).collect();
    let gaps = Gaps {
        insertion: error_probability(INSERTION_QUAL),
        deletion: error_probability(DELETION_QUAL),
        continuation: error_probability(GAP_QUAL),
    };
    let match_to_match = 1.0 - (gaps.insertion + gaps.deletion);
    #[allow(clippy::cast_precision_loss, reason = "a 200 base haplotype")]
    let h = haplotype.len() as f64;
    let mode = Mode { start: LogProb(((1.0 - gaps.continuation) / (h * match_to_match)).ln()) };
    let emissions =
        Emissions { haplotype: haplotype.bases(), read: read.bases(), epsilon: &epsilon };
    let mut hmm = PairHMM::new(&gaps);
    let band = Band::anchored(i32::try_from(READ_OFFSET).unwrap_or(0));

    let mut group = c.benchmark_group("pairhmm/150x200");
    group.bench_function("rust-bio/logprob", |b| {
        b.iter(|| hmm.prob_related(black_box(&emissions), black_box(&mode), None));
    });
    group.bench_function("compair/reference", |b| {
        b.iter(|| {
            align_full(black_box(&haplotype), black_box(&read), &StandardEmission::default())
        });
    });
    group.bench_function("compair/simd", |b| {
        b.iter(|| {
            align_banded_simd(
                black_box(&haplotype),
                black_box(&read),
                &StandardEmission::default(),
                band,
            )
        });
    });
    group.finish();
}

criterion_group!(benches, align);
criterion_main!(benches);
