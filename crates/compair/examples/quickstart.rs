//! Score a few reads against the reference and an alternate haplotype, then
//! turn the scores into per-read evidence and diploid genotype likelihoods.
//!
//! `cargo run -p compair --release --example quickstart`

use compair::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};

/// One aligned read as a BAM record gives it: SEQ, QUAL as SAM prints it
/// (Phred+33), and where the aligner placed the first base, counted from the
/// start of the haplotype window.
struct Record {
    seq: &'static [u8],
    qual: &'static [u8],
    start: i32,
}

fn main() -> Result<(), compair::Error> {
    // The locus: a 60 bp reference window and the same window carrying a
    // candidate G>T SNV at position 30. Both haplotypes must cover the same
    // window so that their scores are comparable.
    let reference =
        Haplotype::from_ascii(b"GCTAAAGACAATTACATAACATACACATCAGCACAAAACTTGTTGGCCCAGTGTGAATCA");
    let alternate =
        Haplotype::from_ascii(b"GCTAAAGACAATTACATAACATACACATCATCACAAAACTTGTTGGCCCAGTGTGAATCA");
    let haplotypes = [&reference, &alternate];

    let records = [
        // Three reads carrying the reference G...
        Record {
            seq: b"AAGACAATTACATAACATACACATCAGCACAAAACTTGTT",
            qual: b"??????????????????????????????????????55",
            start: 4,
        },
        Record {
            seq: b"TACATAACATACACATCAGCACAAAACTTGTTGGCCCAGT",
            qual: b"????????????????????????????????????????",
            start: 12,
        },
        Record {
            seq: b"ACATACACATCAGCACAAAACTTGTTGGCCCAGTGTGAAT",
            qual: b"????????????????????????????????????????",
            start: 18,
        },
        // ...and three carrying the T, the last at Q10 on that base, so it
        // counts for less.
        Record {
            seq: b"GACAATTACATAACATACACATCATCACAAAACTTGTTGG",
            qual: b"????????????????????????????????????????",
            start: 6,
        },
        Record {
            seq: b"ATAACATACACATCATCACAAAACTTGTTGGCCCAGTGTG",
            qual: b"????????????????????????????????????????",
            start: 15,
        },
        Record {
            seq: b"ATTACATAACATACACATCATCACAAAACTTGTTGGCCCA",
            qual: b"????????????????????+???????????????????",
            start: 10,
        },
    ];

    let mut reads = Vec::with_capacity(records.len());
    for record in &records {
        let quals: Vec<BaseQuality> =
            record.qual.iter().map(|q| BaseQuality::from_byte(q.saturating_sub(33))).collect();
        // Without per-base indel qualities (BAM's BI/BD tags), GATK's
        // defaults: Q45 to open a gap on either side, Q10 to extend one. The
        // strand only matters to a conversion-aware emission; see `taps.rs`.
        let read = Read::uniform(
            Base::from_ascii_vec(record.seq.to_vec()),
            &quals,
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )?;
        // The band is the strip of the DP matrix around where the aligner put
        // the read; `anchored` centres the default width (46 columns) on it,
        // which leaves room for indels of up to ~20 bp either way.
        reads.push((read, Band::anchored(record.start)));
    }

    // Keep one workspace and reuse it: it owns the kernels' buffers, so
    // scoring the next locus allocates nothing.
    let mut workspace = Workspace::new();
    let mut scores = Vec::new();
    workspace.align_reads(&haplotypes, &reads, &StandardEmission::default(), &mut scores);

    // `scores` is read-major: read `r` against haplotype `h` is at
    // `r * haplotypes.len() + h`. A perfect match scores about -1.8, not 0:
    // GATK's model puts a uniform prior on where along the haplotype the read
    // starts, log10(1/60) here. Haplotypes of equal length share it, so it
    // cancels in every ratio below.
    println!("read  log10 P(r|ref)  log10 P(r|alt)  log10 LR(alt)  P(alt|r)");
    let mut genotypes = [0.0f64; 3];
    for (index, per_read) in scores.chunks_exact(haplotypes.len()).enumerate() {
        let &[on_ref, on_alt] = per_read else { continue };
        let (on_ref, on_alt) = (on_ref.get(), on_alt.get());
        // With equal priors on the two haplotypes, the posterior is the
        // likelihood ratio squashed into [0, 1].
        let posterior = 1.0 / (1.0 + 10f64.powf(on_ref - on_alt));
        println!(
            "{index:>4}  {on_ref:>14.3}  {on_alt:>14.3}  {:>13.3}  {posterior:>8.4}",
            on_alt - on_ref
        );

        // A diploid genotype explains each read by either of its two
        // haplotypes with probability one half, and the reads independently.
        genotypes[0] += on_ref;
        genotypes[1] += log10_mean(on_ref, on_alt);
        genotypes[2] += on_alt;
    }

    // As a VCF's PL: Phred-scaled and relative to the best genotype.
    let best = genotypes.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("\ngenotype  log10 P(reads|G)    PL");
    for (name, likelihood) in ["0/0", "0/1", "1/1"].iter().zip(genotypes) {
        println!("{name:>8}  {likelihood:>16.3}  {:>4.0}", 10.0 * (best - likelihood));
    }
    Ok(())
}

/// `log10((10^a + 10^b) / 2)` without leaving log space, where a poor read's
/// likelihood can be smaller than an `f64` holds. A pair the band cannot
/// relate at all scores `Log10Likelihood::IMPOSSIBLE`, `-inf`, which this
/// carries through.
fn log10_mean(a: f64, b: f64) -> f64 {
    let (high, low) = if a >= b { (a, b) } else { (b, a) };
    if high == f64::NEG_INFINITY {
        return high;
    }
    high + (1.0 + 10f64.powf(low - high)).log10() - 2f64.log10()
}
