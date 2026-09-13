//! How much of a batch call is re-deriving the haplotype columns.
//!
//! `BatchPlan::fill` rebuilds two things per call: the read's row tracks, which
//! genuinely change per read, and the haplotypes' column tracks, which do not --
//! they depend only on the haplotypes, the emission and the strand. A caller
//! scoring many reads against one candidate set pays for the second every time.
//!
//! This isolates that cost. The band is fixed and anchored near the read's
//! start, so the *kernel*'s work is set by the read length and the band width
//! and does not grow with the haplotype. Any time that does grow with haplotype
//! length is fill.
//!
//! `cargo run -p compair --release --example fillcost`

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Log10Likelihood, Probability, Read,
    Strand, TapsEmission, Workspace,
};
use std::time::Instant;

const READ_LEN: usize = 55;

fn main() {
    println!("55 bp read, band 46 anchored at 25 -- kernel work is constant across rows here");
    println!(
        "{:>10} {:>14} {:>10} {:>14} {:>10}",
        "hap len", "batch/8 (us)", "growth", "strips/1 (us)", "growth"
    );
    let mut baseline = 0.0;
    let mut baseline_strips = 0.0;
    for hap_len in [200usize, 400, 800, 1600] {
        let bases = sequence(hap_len);
        let read_bases: Vec<Base> =
            bases.get(25..25 + READ_LEN).map(<[Base]>::to_vec).unwrap_or_default();
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
        .expect("fixture");
        let betas: Vec<Probability> = (0..hap_len)
            .map(|index| {
                Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                    .unwrap_or(Probability::ZERO)
            })
            .collect();
        let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
        let haplotypes: Vec<Haplotype> = (0..8)
            .map(|variant| {
                let mut alternative = bases.clone();
                for edit in 0..variant {
                    if let Some(base) = alternative.get_mut(35 + edit * 5) {
                        *base = base.inverse();
                    }
                }
                Haplotype::new(alternative)
            })
            .collect();
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let band = Band::anchored(25);
        let mut workspace = Workspace::new();
        let mut out: Vec<Log10Likelihood> = Vec::new();
        let micros = time(|| {
            workspace.align_batch(&refs, &read, &taps, band, &mut out);
            out.iter().map(|score| score.get()).sum::<f64>()
        });
        let first = refs.first().copied().expect("eight haplotypes");
        let strips = time(|| workspace.align_strips_simd(first, &read, &taps, band).get());
        if hap_len == 200 {
            baseline = micros;
            baseline_strips = strips;
            println!("{hap_len:>10} {micros:>14.2} {:>10} {strips:>14.2} {:>10}", "-", "-");
        } else {
            println!(
                "{hap_len:>10} {micros:>14.2} {:>9.2}x {strips:>14.2} {:>9.2}x",
                micros / baseline,
                strips / baseline_strips
            );
        }
    }
    println!(
        "\nBoth plans loop `for index in 0..h` over the whole haplotype, calling\n\
         `site_weights` per site, however narrow the band is. At constant kernel\n\
         work the growth above is all of it."
    );
}

fn time(mut body: impl FnMut() -> f64) -> f64 {
    for _ in 0..200 {
        std::hint::black_box(body());
    }
    let mut best = f64::INFINITY;
    for _ in 0..500 {
        let start = Instant::now();
        let value = body();
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        std::hint::black_box(value);
        best = best.min(elapsed);
    }
    best
}

fn sequence(len: usize) -> Vec<Base> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
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
