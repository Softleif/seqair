//! Profiling/A-B harness for single CRAM codec streams (not a benchmark).
//!
//! `codec_decode_loop <nx16|r4x8|arith|fqz|tok3> <compressed file> <uncompressed len> <iters>`
//! (`arith`, `fqz` and `tok3` read the length from the stream; pass 0.)
//! Prints the best and median MB/s (uncompressed bytes) over `iters` decodes.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::print_stdout, reason = "example")]
#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::panic,
    reason = "example"
)]
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let codec = &args[1];
    let src = std::fs::read(&args[2]).unwrap();
    let mut len: usize = args[3].parse().unwrap();
    let iters: usize = args[4].parse().unwrap();
    let expected = args.get(5).map(|p| std::fs::read(p).unwrap());
    let mut times = Vec::with_capacity(iters);
    let mut check = 0u64;
    // Tables kept between blocks, as the CRAM reader keeps them.
    let mut tok3 = seqair::cram::tok3::Decoder::new();
    for _ in 0..iters {
        let t = Instant::now();
        let out = match codec.as_str() {
            "nx16" => seqair::cram::rans_nx16::decode(&src, len).unwrap(),
            "r4x8" => seqair::cram::rans::decode(&src).unwrap(),
            "arith" => seqair::cram::arith::decode(&src, len).unwrap(),
            "fqz" => seqair::cram::fqzcomp::decode(&src).unwrap(),
            "tok3" => tok3.decode(&src).unwrap(),
            _ => panic!("codec"),
        };
        times.push(t.elapsed().as_secs_f64());
        if len == 0 {
            len = out.len();
        }
        assert_eq!(out.len(), len);
        if let Some(e) = &expected {
            assert!(out == *e, "decoded bytes differ from the expected file");
        }
        check = check.wrapping_add(out.iter().map(|&b| u64::from(b)).sum::<u64>());
    }
    times.sort_by(f64::total_cmp);
    let mbs = |t: f64| len as f64 / t / 1e6;
    println!(
        "{codec} {}: best {:.0} MB/s, median {:.0} MB/s (check {check})",
        args[2].rsplit('/').next().unwrap(),
        mbs(times[0]),
        mbs(times[iters / 2])
    );
}
