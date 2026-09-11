#![allow(dead_code, reason = "shared with the tests")]
#![allow(clippy::print_stdout, reason = "this example exists to print them")]
include!("../tests/support/mod.rs");

use compair::{StandardEmission, align_banded, align_full};

fn main() {
    let vectors = gatk_vectors().expect("the bundled test data parses");
    let mut worst_ref = 0.0f64;
    let mut worst_band = 0.0f64;
    let mut narrow_lower = 0usize;
    let mut offsets = Vec::new();
    for v in &vectors {
        let full = align_full(&v.haplotype, &v.read, &StandardEmission::default()).get();
        worst_ref = worst_ref.max((full - v.expected_log10).abs());
        let o = seed_offset(&v.haplotype, &v.read);
        offsets.push(o);
        let b = Band::anchored(o);
        let banded = align_banded(&v.haplotype, &v.read, &StandardEmission::default(), b).get();
        worst_band = worst_band.max((banded - full).abs());
        let narrow = Band::new(4, o).expect("w");
        let n = align_banded(&v.haplotype, &v.read, &StandardEmission::default(), narrow).get();
        if n < full - 1e-6 {
            narrow_lower += 1;
        }
    }
    offsets.sort_unstable();
    println!("vectors            : {}", vectors.len());
    println!("worst |full - GATK|: {worst_ref:.3e}");
    println!("worst |band - full|: {worst_band:.3e}");
    println!("offsets min/max    : {:?} / {:?}", offsets.first(), offsets.last());
    println!("width-4 band lower : {narrow_lower} of {}", vectors.len());
}
