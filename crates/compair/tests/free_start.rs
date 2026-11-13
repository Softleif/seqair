//! The strip, batch, pairs and candidates kernels at the two ends of a read:
//! it may begin with an insertion or with a mismatch, and end at the
//! haplotype's last column. The `f64` oracles (`align_full`, and
//! `align_banded_f64_rows` over a band) are the reference; the kernels agree
//! with one another bit for bit elsewhere.

#![allow(clippy::arithmetic_side_effects, clippy::cast_possible_truncation, reason = "test code")]

#[path = "../src/pinned.rs"]
mod pinned;
#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    BATCH, Band, Base, BaseQuality, Haplotype, Pair, Read, StandardEmission, Workspace,
    align_banded_f64_rows, align_full, align_strips, align_strips_simd,
};
use hegel::TestCase;
use hegel::generators as gs;
use support::{any_base, any_strand, length};

#[derive(Debug)]
struct Ends {
    haplotype: Haplotype,
    read: Read,
    start: usize,
}

/// A read cut from a random haplotype, flush with its last column half the
/// time, and then given a mismatching first base, an inserted one, or neither.
/// The indel qualities range widely, so that an insertion opens both more and
/// less easily than it continues.
#[hegel::composite]
fn any_ends(tc: &TestCase) -> Ends {
    let h = tc.draw(length(2, 40));
    let bases = tc.draw(gs::vecs(any_base()).min_size(h).max_size(h));
    let n = tc.draw(length(1, h));
    let start = if tc.draw(gs::booleans()) { h - n } else { tc.draw(length(0, h - n)) };
    let mut read: Vec<Base> = bases.get(start..start + n).unwrap_or_default().to_vec();
    match tc.draw(length(0, 2)) {
        1 => {
            let shift = tc.draw(length(1, 3));
            if let Some(first) = read.first_mut() {
                let order = [Base::A, Base::C, Base::G, Base::T];
                let at = order.iter().position(|base| base == first).unwrap_or(0);
                *first = order.get((at + shift) % 4).copied().unwrap_or(Base::A);
            }
        }
        2 => read.insert(0, tc.draw(any_base())),
        _ => {}
    }
    let quality =
        |lo: u8| BaseQuality::from_byte(tc.draw(gs::integers::<u8>().min_value(lo).max_value(45)));
    let base_quals: Vec<BaseQuality> = read.iter().map(|_| quality(10)).collect();
    let read =
        Read::uniform(read, &base_quals, quality(5), quality(5), quality(1), tc.draw(any_strand()))
            .expect("a valid read");
    Ends { haplotype: Haplotype::new(bases), read, start }
}

fn close(got: f64, want: f64) -> bool {
    (got - want).abs() < 1e-4 * (1.0 + want.abs())
}

/// Each kernel against `align_full` through a band that holds the whole
/// matrix.
#[hegel::test(test_cases = pinned::cases(512))]
fn every_strip_family_kernel_is_the_full_recurrence_at_the_read_ends(tc: TestCase) {
    let Ends { haplotype, read, .. } = tc.draw(any_ends());
    let emission = StandardEmission::default();
    let band = Band::new(Band::MAX_WIDTH, 0).expect("a legal width");
    let want = align_full(&haplotype, &read, &emission).get();

    let mut workspace = Workspace::new();
    let (mut batch, mut pairs, mut prepared) = (Vec::new(), Vec::new(), Vec::new());
    workspace.align_batch(&[&haplotype; BATCH], &read, &emission, band, &mut batch);
    workspace.align_reads(&[&haplotype], &[(&read, band); BATCH], &emission, &mut pairs);
    workspace.candidates(&[&haplotype; BATCH], emission).align(&read, band, &mut prepared);

    let mut all = vec![
        ("strips", align_strips(&haplotype, &read, &emission, band).get()),
        ("strips/simd", align_strips_simd(&haplotype, &read, &emission, band).get()),
    ];
    all.extend(batch.iter().map(|score| ("batch", score.get())));
    all.extend(pairs.iter().map(|score| ("pairs", score.get())));
    all.extend(prepared.iter().map(|score| ("candidates", score.get())));
    for (kernel, got) in all {
        assert!(close(got, want), "{kernel}: {got} against the full recurrence {want}");
    }
}

/// Pairs and a batch whose lanes differ in haplotype length, so that column
/// `h` is a different column in each, and in band offset, so that column 0 is
/// shifted into the band, against the `f64` banded recurrence.
#[hegel::test(test_cases = pinned::cases(512))]
fn lanes_with_their_own_ends_are_the_banded_recurrence(tc: TestCase) {
    let Ends { haplotype, read, start } = tc.draw(any_ends());
    let width = 2 * tc.draw(length(2, 20)) as u32;
    let emission = StandardEmission::default();
    let cut = tc.draw(length(0, 3));
    // Cut no shorter than the read's end on the haplotype, where it ends at the
    // last column: a read that cannot align is rescored by `align_banded_f64`.
    let floor = (start + read.len()).min(haplotype.len());
    let haplotypes: Vec<Haplotype> = (0..BATCH)
        .map(|k| {
            let bases = haplotype.bases();
            let len = bases.len().saturating_sub(k * cut).max(floor);
            Haplotype::new(bases.get(..len).unwrap_or(bases).to_vec())
        })
        .collect();
    let spread = tc.draw(length(0, 4));
    let bands: Vec<Band> = (0..BATCH)
        .map(|k| {
            let offset = i32::try_from(start + (k * spread) % 5).expect("a small offset");
            Band::new(width, offset).expect("a legal width")
        })
        .collect();
    let pairs: Vec<Pair<'_>> = haplotypes
        .iter()
        .zip(&bands)
        .map(|(haplotype, &band)| Pair::new(haplotype, &read, band))
        .collect();

    let mut workspace = Workspace::new();
    let mut scores = Vec::new();
    workspace.align_pairs(&pairs, &emission, &mut scores);
    for (lane, (pair, got)) in pairs.iter().zip(&scores).enumerate() {
        let want = align_banded_f64_rows(pair.haplotype, &read, &emission, pair.band).get();
        assert!(
            got.get() == want || close(got.get(), want),
            "pairs lane {lane}: {} against the banded recurrence {want}",
            got.get()
        );
    }

    let band = bands.first().copied().expect("a lane");
    let refs: Vec<&Haplotype> = haplotypes.iter().collect();
    workspace.align_batch(&refs, &read, &emission, band, &mut scores);
    for (lane, (haplotype, got)) in haplotypes.iter().zip(&scores).enumerate() {
        let want = align_banded_f64_rows(haplotype, &read, &emission, band).get();
        assert!(
            got.get() == want || close(got.get(), want),
            "batch lane {lane}: {} against the banded recurrence {want}",
            got.get()
        );
        let strips = align_strips(haplotype, &read, &emission, band).get();
        assert!(got.get().to_bits() == strips.to_bits(), "batch lane {lane}: not the strip kernel");
    }
}
