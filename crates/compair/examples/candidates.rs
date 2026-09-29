//! Many candidate haplotypes, reads arriving one at a time: genotyping the
//! length of a `CA` repeat.
//!
//! An aligner cannot say how many units a short tandem repeat has -- every
//! count aligns about as well -- so a caller builds one haplotype per
//! plausible length and lets the pair-HMM decide. `Workspace::candidates`
//! prepares the haplotypes once, and `Candidates::align` then scores each read
//! against all of them.
//!
//! `cargo run -p compair --release --example candidates`

use compair::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};

const LEFT_FLANK: &str = "CCGTAATGCCTTTCCCTAAC";
const RIGHT_FLANK: &str = "AGAGTTTTTCGAACTCGTGT";
/// The repeat counts to consider.
const UNITS: std::ops::RangeInclusive<usize> = 3..=9;

fn main() -> Result<(), compair::Error> {
    // One haplotype per repeat count. They differ in length, which is fine:
    // what makes their scores comparable is that each spans the same locus,
    // flank to flank.
    let haplotypes: Vec<Haplotype> = UNITS
        .map(|units| {
            let sequence = format!("{LEFT_FLANK}{}{RIGHT_FLANK}", "CA".repeat(units));
            Haplotype::from_ascii(sequence.as_bytes())
        })
        .collect();

    // Reads from a sample heterozygous for 5 and 6 units, as `(SEQ, start)`.
    // Every read starts in the left flank, so its start is the same on every
    // haplotype and one band serves them all.
    let records: [(&[u8], i32); 5] = [
        (b"GTAATGCCTTTCCCTAACCACACACACAAGAGTTTTTCGAACTC", 2),
        (b"TGCCTTTCCCTAACCACACACACAAGAGTTTTTCGAACTCGTGT", 6),
        (b"TAATGCCTTTCCCTAACCACACACACACAAGAGTTTTTCGAACTCGTG", 3),
        (b"CCTTTCCCTAACCACACACACACAAGAGTTTTTCGAACTCGTGT", 8),
        (b"TTTCCCTAACCACACACACAAGAGTTTTTCGAACTCGTGT", 10),
    ];

    // `candidates` borrows the workspace and the haplotypes for as long as it
    // lives, and derives each haplotype's per-column terms once rather than
    // once per read. `align` then picks the fastest kernel for the number of
    // haplotypes; every kernel gives the same score to the last bit.
    let mut workspace = Workspace::new();
    let mut candidates = workspace.candidates(&haplotypes, StandardEmission::default());
    let mut scores = Vec::new();
    let mut per_read: Vec<Vec<f64>> = Vec::with_capacity(records.len());

    println!("read  best  runner-up  margin (log10)");
    for (index, (seq, start)) in records.into_iter().enumerate() {
        let read = Read::uniform(
            Base::from_ascii_vec(seq.to_vec()),
            &vec![BaseQuality::from_byte(30); seq.len()],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )?;
        // One score per haplotype, in the haplotypes' order.
        candidates.align(&read, Band::anchored(start), &mut scores);
        let log10s: Vec<f64> = scores.iter().map(|score| score.get()).collect();

        let mut ranked: Vec<(usize, f64)> = UNITS.zip(log10s.iter().copied()).collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
        if let [(best, top), (second, next), ..] = ranked.as_slice() {
            println!("{index:>4}  {best:>4}  {second:>9}  {:>14.2}", top - next);
        }
        per_read.push(log10s);
    }

    // Every unordered pair of haplotypes is a diploid genotype; each read
    // comes from either of its two with probability one half.
    let mut best: Option<((usize, usize), f64)> = None;
    for (first, a) in UNITS.enumerate() {
        for (second, b) in UNITS.enumerate().skip(first) {
            let likelihood: f64 = per_read
                .iter()
                .filter_map(|scores| Some(log10_mean(*scores.get(first)?, *scores.get(second)?)))
                .sum();
            if best.is_none_or(|(_, top)| likelihood > top) {
                best = Some(((a, b), likelihood));
            }
        }
    }
    if let Some(((a, b), likelihood)) = best {
        println!("\nmost likely genotype: {a}/{b} units (log10 P(reads|G) = {likelihood:.2})");
    }
    Ok(())
}

/// `log10((10^a + 10^b) / 2)` without leaving log space.
fn log10_mean(a: f64, b: f64) -> f64 {
    let (high, low) = if a >= b { (a, b) } else { (b, a) };
    if high == f64::NEG_INFINITY {
        return high;
    }
    high + (1.0 + 10f64.powf(low - high)).log10() - 2f64.log10()
}
