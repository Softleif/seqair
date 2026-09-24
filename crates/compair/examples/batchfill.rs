//! Where the batch kernel starts beating the strip kernel, as a function of how
//! many haplotypes the caller has.
//!
//! `align_batch` computes all eight lanes whatever the caller asked for, so a
//! group of `n < 8` throws `8 - n` away and the per-alignment cost falls only
//! with the fill. This finds the crossover, which is the number the dispatch in
//! `Workspace::align_candidates` has to know -- `BATCH_BREAK_EVEN`.
//!
//! Both sides race through the *same* lane, `align_strips_simd` against
//! `align_batch`: racing two different lanes would move the crossover by
//! whatever the lanes differ by rather than by what the fill costs, which is
//! not the number the dispatch needs.
//!
//! `cargo run -p compair --release --example batchfill`

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Log10Likelihood, Probability, Read,
    StandardEmission, Strand, TapsEmission, Workspace,
};
use std::time::Instant;

const HAPLOTYPE_LEN: usize = 200;
const READ_OFFSET: usize = 25;

fn main() {
    for read_len in [55usize, 150] {
        let (haplotypes, read, betas) = fixture(read_len, 24);
        let band = Band::anchored(i32::try_from(READ_OFFSET).unwrap_or(0));
        let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
        let standard = StandardEmission::default();
        println!(
            "\n=== {read_len} bp read, {HAPLOTYPE_LEN} bp haplotypes, band {} ===",
            band.width()
        );
        println!(
            "{:>5} {:>6} {:>10} {:>10} {:>10} {:>9} {:>9}",
            "haps", "fill", "strips", "batch", "dispatch", "batch/s", "disp/s"
        );
        for n in [1usize, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 24] {
            let refs: Vec<&Haplotype> = haplotypes.iter().take(n).collect();
            let chunks = n.div_ceil(8);
            #[allow(clippy::cast_precision_loss, reason = "small counts")]
            let fill = n as f64 / (chunks * 8) as f64 * 100.0;
            let (strips, batch, dispatch) = race(&refs, &read, &taps, band);
            println!(
                "{n:>5} {fill:>5.0}% {strips:>10.2} {batch:>10.2} {dispatch:>10.2} {:>8.2}x {:>8.2}x",
                strips / batch,
                strips / dispatch
            );
        }
        let _ = &standard;
    }
}

/// Microseconds for one pass of each kernel over the whole group, best of many.
fn race<E: compair::Emission>(
    haplotypes: &[&Haplotype],
    read: &Read,
    emission: &E,
    band: Band,
) -> (f64, f64, f64) {
    let mut workspace = Workspace::new();
    let mut out: Vec<Log10Likelihood> = Vec::new();
    let strips = time(|| {
        let mut sum = 0.0;
        for haplotype in haplotypes {
            sum += workspace.align_strips_simd(haplotype, read, emission, band).get();
        }
        sum
    });
    let batch = time(|| {
        workspace.align_batch(haplotypes, read, emission, band, &mut out);
        out.iter().map(|score| score.get()).sum::<f64>()
    });
    let dispatch = time(|| {
        workspace.align_candidates(haplotypes, read, emission, band, &mut out);
        out.iter().map(|score| score.get()).sum::<f64>()
    });
    (strips, batch, dispatch)
}

fn time(mut body: impl FnMut() -> f64) -> f64 {
    for _ in 0..200 {
        std::hint::black_box(body());
    }
    let mut best = f64::INFINITY;
    for _ in 0..400 {
        let start = Instant::now();
        let value = body();
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        std::hint::black_box(value);
        best = best.min(elapsed);
    }
    best
}

fn fixture(read_len: usize, count: usize) -> (Vec<Haplotype>, Read, Vec<Probability>) {
    let bases = sequence(HAPLOTYPE_LEN, 0x9e37_79b9_7f4a_7c15);
    let read_bases: Vec<Base> =
        bases.get(READ_OFFSET..READ_OFFSET + read_len).map(<[Base]>::to_vec).unwrap_or_default();
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
    .expect("the fixture must build");
    let betas = (0..HAPLOTYPE_LEN)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    let haplotypes = (0..count)
        .map(|variant| {
            let mut alternative = bases.clone();
            for edit in 0..variant {
                if let Some(base) = alternative.get_mut(READ_OFFSET + 10 + edit * 5) {
                    *base = base.inverse();
                }
            }
            Haplotype::new(alternative)
        })
        .collect();
    (haplotypes, read, betas)
}

fn sequence(len: usize, seed: u64) -> Vec<Base> {
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
