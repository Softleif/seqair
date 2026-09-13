//! How wrong the band's anchor may be before the band starts losing paths.
//!
//! Every band measurement in this repo feeds `Band::anchored` a seed offset
//! computed by scanning for the best ungapped match -- test scaffolding standing
//! in for an aligner. A caller has a real anchor from its alignment, and how
//! good that anchor is decides whether the band is free or expensive. This
//! perturbs the anchor and reports what it costs, which is the requirement to
//! hand the caller.
//!
//! `cargo run -p compair --release --example anchorsweep`

use compair::{Band, StandardEmission, align_full, align_strips_simd};

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

/// A loss below this is rounding between an `f32` band and an `f64` matrix.
const KEPT: f64 = 1e-3;

fn main() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let standard = StandardEmission::default();
    let full: Vec<f64> =
        pairs.iter().map(|pair| align_full(&pair.haplotype, &pair.read, &standard).get()).collect();

    println!(
        "10s, {} pairs. Rows are band width, columns are anchor error in bases.\n\
         Each cell is the share of pairs whose score still matches the full matrix.\n",
        pairs.len()
    );
    let errors = [0i32, 2, 5, 10, 15, 20, 30, 45];
    print!("{:>6}", "width");
    for error in errors {
        print!("{:>9}", format!("+/-{error}"));
    }
    println!();

    for width in [14u32, 30, 46, 94] {
        print!("{width:>6}");
        for error in errors {
            // The worse of the two directions: an anchor is as likely to be
            // early as late, and a caller has to survive both.
            let mut worst_rate = 1.0f64;
            for signed in if error == 0 { vec![0] } else { vec![-error, error] } {
                let mut kept = 0usize;
                for (pair, reference) in pairs.iter().zip(&full) {
                    let band = Band::new(width, pair.offset.saturating_add(signed))
                        .unwrap_or(Band::anchored(pair.offset));
                    let got = align_strips_simd(&pair.haplotype, &pair.read, &standard, band).get();
                    if reference - got < KEPT {
                        kept += 1;
                    }
                }
                #[allow(clippy::cast_precision_loss, reason = "a few thousand pairs")]
                let rate = kept as f64 / pairs.len() as f64;
                worst_rate = worst_rate.min(rate);
            }
            print!("{:>8.1}%", worst_rate * 100.0);
        }
        println!();
    }
    println!(
        "\nA band of width w centred perfectly tolerates an anchor error of about\n\
         w/2 before the read's own diagonal leaves the strip, so these should fall\n\
         off a cliff near half the width -- and the useful question is how much\n\
         margin the default buys over the narrowest band that works when centred."
    );
}
