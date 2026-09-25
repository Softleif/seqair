//! hts-specs' official vectors for the arithmetic coder (`range/`) and the
//! name tokeniser (`tok3/`, whose `.11`–`.19` files code their streams with
//! the arithmetic coder), checked against the gzip originals as the vectors'
//! README describes. Every file in both directories is decoded, so a vector
//! added upstream cannot be missed.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

use std::io::Read;
use std::path::{Path, PathBuf};

use seqair::cram::{arith, tok3};

const CODECS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs");

fn gunzip(name: &str) -> Vec<u8> {
    let path = Path::new(CODECS).join("gzip").join(format!("{name}.gz"));
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// The vector files of a codec directory as `(stem, suffix, path)`, where
/// `q40+dir.65` is `("q40+dir", 65)`.
fn vectors(dir: &str) -> Vec<(String, u32, PathBuf)> {
    let mut out: Vec<_> = std::fs::read_dir(Path::new(CODECS).join(dir))
        .unwrap()
        .map(|e| e.unwrap().path())
        .map(|path| {
            let name = path.file_name().unwrap().to_str().unwrap().to_owned();
            let (stem, suffix) = name.rsplit_once('.').unwrap();
            (stem.to_owned(), suffix.parse().unwrap(), path)
        })
        .collect();
    out.sort();
    out
}

/// The README's processing: the u32 file as-is; the quality files' first
/// whitespace-separated column of each line, newlines dropped.
fn expected_range_output(stem: &str) -> Vec<u8> {
    let original = gunzip(stem);
    if stem == "u32" {
        return original;
    }
    original
        .split(|&b| b == b'\n')
        .filter_map(|line| line.split(u8::is_ascii_whitespace).find(|f| !f.is_empty()))
        .flatten()
        .copied()
        .collect()
}

// r[verify cram.codec.arith+2]
// r[verify cram.codec.arith.wrapper]
// r[verify cram.codec.arith.order]
// r[verify cram.codec.arith.rle]
// r[verify cram.codec.arith.pack]
// r[verify cram.codec.arith.stripe]
// r[verify cram.codec.arith.ext]
#[test]
fn range_vectors_decode_to_the_originals() {
    let files = vectors("range");
    assert!(files.len() >= 32, "expected the 32 hts-specs range vectors, found {}", files.len());
    let mut flags_seen = 0u8;
    for (stem, suffix, path) in &files {
        let src = std::fs::read(path).unwrap();
        flags_seen |= src[0];
        let decoded = arith::decode(&src, 0).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let expected = expected_range_output(stem);
        assert_eq!(decoded.len(), expected.len(), "{stem}.{suffix}: length");
        assert!(decoded == expected, "{stem}.{suffix}: content differs");
    }
    // Their top-level flags reach order-1, EXT, STRIPE, RLE and PACK (NOSZ and
    // CAT occur in the stripe substreams).
    assert_eq!(flags_seen, 0x01 | 0x04 | 0x08 | 0x40 | 0x80, "{flags_seen:#x}");
}

// r[verify cram.codec.tok3_arith]
// r[verify cram.codec.tok3]
#[test]
fn tok3_vectors_decode_to_the_originals() {
    let files = vectors("tok3");
    assert!(files.len() >= 100, "expected the hts-specs tok3 vectors, found {}", files.len());
    let mut arith_files = 0;
    for (stem, suffix, path) in &files {
        let src = std::fs::read(path).unwrap();
        // `.11`–`.19` use the arithmetic coder; the header says so in byte 8.
        let use_arith = src[8] != 0;
        assert_eq!(use_arith, *suffix > 10, "{stem}.{suffix}: use_arith byte");
        arith_files += usize::from(use_arith);

        let decoded = tok3::decode(&src).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        // tok3 ends every name with NUL; the originals end every line with
        // a newline.
        let mut expected = gunzip(stem);
        for b in &mut expected {
            if *b == b'\n' {
                *b = 0;
            }
        }
        assert_eq!(decoded.len(), expected.len(), "{stem}.{suffix}: length");
        assert!(decoded == expected, "{stem}.{suffix}: content differs");
    }
    assert!(arith_files >= 50, "only {arith_files} arith-coded tok3 vectors");
}
