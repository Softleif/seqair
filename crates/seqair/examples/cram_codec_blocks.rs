//! seqair against htscodecs on a CRAM file's own blocks (not a benchmark
//! harness; `benches/cram.rs` has the criterion groups).
//!
//! `cram_codec_blocks <file.cram> [--codec arith|fqzcomp|tok3] [--dump <dir>]`
//!
//! Every block of the file coded with the arithmetic coder, fqzcomp or tok3
//! is decoded by both, once, and must come out the same; then each codec's
//! blocks are decoded as a set, alternating seqair and htscodecs rounds so
//! machine drift hits both alike, and the best round of each is reported in
//! MiB/s of decoded data. `--dump` also writes each block's payload to
//! `<dir>/<codec><n>.bin` for `codec_decode_loop` and profilers.
#![allow(
    clippy::print_stdout,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    reason = "example"
)]

#[path = "../tests/support/cram_raw_blocks.rs"]
#[allow(dead_code, reason = "only some block methods are timed")]
mod cram_raw_blocks;
#[path = "../benches/support/htscodecs.rs"]
mod htscodecs;

use std::hint::black_box;
use std::time::{Duration, Instant};

use cram_raw_blocks::{ARITH, FQZCOMP, RawBlock, TOK3, raw_blocks};

/// Rounds per decoder; each round decodes the codec's blocks `REPEATS` times.
const ROUNDS: usize = 7;
/// Enough repeats that a round of the smallest set takes a few milliseconds.
const REPEATS: usize = 20;

/// What seqair's decoders keep between blocks, as the CRAM reader keeps it.
#[derive(Default)]
struct Seqair {
    tok3: seqair::cram::tok3::Decoder,
}

struct Codec {
    name: &'static str,
    method: u8,
    seqair: fn(&mut Seqair, &RawBlock) -> Vec<u8>,
    htscodecs: fn(&RawBlock) -> htscodecs::MallocBuf,
}

const CODECS: [Codec; 3] = [
    Codec {
        name: "arith",
        method: ARITH,
        seqair: |_, b| seqair::cram::arith::decode(&b.data, b.uncompressed_len).unwrap(),
        htscodecs: |b| htscodecs::arith_uncompress(&b.data).unwrap(),
    },
    Codec {
        name: "fqzcomp",
        method: FQZCOMP,
        seqair: |_, b| seqair::cram::fqzcomp::decode(&b.data).unwrap(),
        htscodecs: |b| htscodecs::fqz_decompress(&b.data).unwrap(),
    },
    Codec {
        name: "tok3",
        method: TOK3,
        seqair: |s, b| s.tok3.decode(&b.data).unwrap(),
        htscodecs: |b| htscodecs::tok3_decode_names(&b.data).unwrap(),
    },
];

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: cram_codec_blocks <file.cram> [--codec C] [--dump DIR]");
    let (mut only, mut dump) = (None, None);
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--codec" => only = args.next(),
            "--dump" => dump = args.next(),
            _ => panic!("unknown argument {flag}"),
        }
    }
    let blocks = raw_blocks(&std::fs::read(&path).unwrap());
    let mut state = Seqair::default();

    for codec in CODECS.iter().filter(|c| only.as_deref().is_none_or(|o| o == c.name)) {
        let set: Vec<&RawBlock> = blocks.iter().filter(|b| b.method == codec.method).collect();
        if set.is_empty() {
            continue;
        }
        for (i, block) in set.iter().enumerate() {
            let ours = (codec.seqair)(&mut state, block);
            assert_eq!(ours.len(), block.uncompressed_len, "{} block {i}: length", codec.name);
            let theirs = (codec.htscodecs)(block);
            assert!(
                ours == theirs.as_slice(),
                "{} block {i}: seqair and htscodecs differ",
                codec.name
            );
            if let Some(dir) = &dump {
                std::fs::write(format!("{dir}/{}{i}.bin", codec.name), &block.data).unwrap();
            }
        }

        let (mut best_ours, mut best_theirs) = (Duration::MAX, Duration::MAX);
        for _ in 0..ROUNDS {
            best_ours = best_ours.min(time(|| {
                for block in &set {
                    black_box((codec.seqair)(&mut state, black_box(block)));
                }
            }));
            best_theirs = best_theirs.min(time(|| {
                for block in &set {
                    black_box((codec.htscodecs)(black_box(block)));
                }
            }));
        }
        let decoded: usize = set.iter().map(|b| b.uncompressed_len).sum();
        let mib_s = |t: Duration| (decoded * REPEATS) as f64 / t.as_secs_f64() / f64::from(1 << 20);
        println!(
            "{:8} {:3} blocks {:>10} B  seqair {:7.1} MiB/s  htscodecs {:7.1} MiB/s  ratio {:.2}",
            codec.name,
            set.len(),
            decoded,
            mib_s(best_ours),
            mib_s(best_theirs),
            best_theirs.as_secs_f64() / best_ours.as_secs_f64(),
        );
    }
}

/// How long `REPEATS` runs of `f` take.
fn time(mut f: impl FnMut()) -> Duration {
    let start = Instant::now();
    for _ in 0..REPEATS {
        f();
    }
    start.elapsed()
}
