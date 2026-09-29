//! Why the emission model matters for TAPS data: the same reads scored with
//! GATK's plain emission and with the conversion-aware one.
//!
//! TAPS converts *methylated* cytosines, so a read from a methylated `CpG`
//! shows a `T` where the reference has a `C`. To a plain pair-HMM that `T` is
//! a mismatch, and five such reads look like a confident `C>T` SNV.
//! `TapsEmission` scores a `T` over a `C` at the site's methylation level
//! instead, so the same reads say "methylated `CpG`", unless the caller has
//! reason to believe the site is unmethylated.
//!
//! `cargo run -p compair --release --example taps`

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Emission, Haplotype, Probability, Read,
    StandardEmission, Strand, TapsEmission, Workspace,
};

/// The `C` of the window's one `CpG`.
const CPG: usize = 25;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The reference carries a `CpG` at 25-26; the alternate haplotype is the
    // `C>T` SNV there, which also destroys the `CpG`. Each haplotype's `CpG`s
    // are read off its own sequence, so nothing marks the site by hand.
    let reference =
        Haplotype::from_ascii(b"ATAAGGGTTATATCTGTCTAAGTGGCGGATAGTTAGAAGGCACATAAGATCATATTAGTG");
    let alternate =
        Haplotype::from_ascii(b"ATAAGGGTTATATCTGTCTAAGTGGTGGATAGTTAGAAGGCACATAAGATCATATTAGTG");
    let haplotypes = [&reference, &alternate];

    // Five original-top-strand (OT) reads of the *reference*: four show the
    // `CpG`'s `C` converted to `T`, one shows it unconverted. `(SEQ, start)`.
    let records: [(&[u8], usize); 5] = [
        (b"AGGGTTATATCTGTCTAAGTGGTGGATAGTTAGAAG", 3),
        (b"TATATCTGTCTAAGTGGTGGATAGTTAGAAGGCACA", 8),
        (b"ATCTGTCTAAGTGGTGGATAGTTAGAAGGCACATAA", 11),
        (b"TGTCTAAGTGGCGGATAGTTAGAAGGCACATAAGAT", 14),
        (b"CTAAGTGGTGGATAGTTAGAAGGCACATAAGATCAT", 17),
    ];
    let mut reads = Vec::with_capacity(records.len());
    let mut shown = Vec::with_capacity(records.len());
    for (seq, start) in records {
        // The strand is what makes a conversion visible: on OT reads a
        // converted `C` shows as `T`, on original-bottom (OB) reads a
        // converted `G` shows as `A` (the bottom strand's `C`, in reference
        // orientation). A `T` over a `C` on an OB read is a plain mismatch.
        let read = Read::uniform(
            Base::from_ascii_vec(seq.to_vec()),
            &vec![BaseQuality::from_byte(30); seq.len()],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )?;
        reads.push((read, Band::anchored(i32::try_from(start)?)));
        shown.push(
            CPG.checked_sub(start).and_then(|at| seq.get(at)).map_or('-', |b| char::from(*b)),
        );
    }
    let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();

    // Methylation levels (beta) per haplotype position; a caller estimates
    // them, from the reads or from a prior. Positions without an estimate,
    // or past the end, count as unmethylated. `Betas::Uniform` is the same
    // level at every `CpG`.
    let methylated = betas_with(0.8)?;
    let unmethylated = betas_with(0.0)?;
    // TAPS converts ~96% of methylated and ~0.4% of unmethylated cytosines
    // in this project's spike-in controls; `ConversionModel::new` takes a
    // library's own rates.
    let chemistry = ConversionModel::taps_default();

    let mut workspace = Workspace::new();
    let columns = [
        ratios(&mut workspace, &haplotypes, &reads, &StandardEmission::default()),
        ratios(
            &mut workspace,
            &haplotypes,
            &reads,
            &TapsEmission::new(chemistry, Betas::PerSite(&methylated)),
        ),
        ratios(
            &mut workspace,
            &haplotypes,
            &reads,
            &TapsEmission::new(chemistry, Betas::PerSite(&unmethylated)),
        ),
    ];

    println!("log10 P(read | C>T haplotype) - log10 P(read | reference)\n");
    println!("read  base at CpG  standard  TAPS beta=0.8  TAPS beta=0");
    for (index, base) in shown.iter().enumerate() {
        let [standard, methylated, unmethylated] =
            columns.each_ref().map(|column| column.get(index).copied().unwrap_or(f64::NAN));
        println!(
            "{index:>4}  {base:>11}  {standard:>8.2}  {methylated:>13.2}  {unmethylated:>11.2}"
        );
    }
    let [standard, methylated, unmethylated] =
        columns.each_ref().map(|column| column.iter().sum::<f64>());
    println!(" sum  {:>11}  {standard:>8.2}  {methylated:>13.2}  {unmethylated:>11.2}", "");
    println!(
        "\nThe plain emission calls the SNV with odds of 10^{standard:.0} to one. At beta = 0.8\n\
         the converted T's are expected and the one unconverted C decides for the\n\
         reference. At beta = 0 a T over a C is a 0.4% event, so it is evidence for\n\
         the SNV again: the model explains T's away only where methylation can."
    );
    Ok(())
}

/// Every read's log10 likelihood ratio, the alternate haplotype over the
/// reference, under one emission.
fn ratios<E: Emission>(
    workspace: &mut Workspace,
    haplotypes: &[&Haplotype],
    reads: &[(&Read, Band)],
    emission: &E,
) -> Vec<f64> {
    let mut scores = Vec::new();
    workspace.align_reads(haplotypes, reads, emission, &mut scores);
    scores
        .chunks_exact(haplotypes.len())
        .map(|pair| match pair {
            [on_ref, on_alt] => on_alt.get() - on_ref.get(),
            _ => f64::NAN,
        })
        .collect()
}

/// Unmethylated everywhere but the window's `CpG`. The betas are indexed by
/// haplotype position, and one emission serves every haplotype, which is
/// simple for an SNV; around an indel the positions after it shift.
fn betas_with(beta: f64) -> Result<Vec<Probability>, Box<dyn std::error::Error>> {
    let mut betas = vec![Probability::ZERO; 60];
    for position in [CPG, CPG + 1] {
        if let Some(slot) = betas.get_mut(position) {
            *slot = Probability::new(beta)?;
        }
    }
    Ok(betas)
}
