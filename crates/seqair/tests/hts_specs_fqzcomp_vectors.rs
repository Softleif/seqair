//! Decode the hts-specs fqzcomp test vectors (`tests/data/hts-specs/cram/codecs`)
//! and compare them to the gzipped originals, processed as the vectors' README
//! says: the first whitespace-separated column of each line, newlines removed.
//! The originals are ASCII phred+33; htscodecs' `fqzcomp_qual` test tool,
//! which made the vectors, subtracts 33 before compressing (as CRAM stores
//! qualities), so the decoded bytes are the originals minus 33.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

use std::{fs, io::Read, path::PathBuf};

use seqair::cram::fqzcomp;

fn codecs_dir() -> PathBuf {
    PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs"))
}

/// The original of a vector: first column of every line, concatenated, as
/// binary qualities.
fn expected(name: &str) -> Vec<u8> {
    let gz = fs::read(codecs_dir().join("gzip").join(format!("{name}.gz"))).unwrap();
    let mut text = Vec::new();
    flate2::read::MultiGzDecoder::new(gz.as_slice()).read_to_end(&mut text).unwrap();
    text.split(|&b| b == b'\n')
        .flat_map(|line| line.split(|b| b.is_ascii_whitespace()).next().unwrap_or_default())
        .map(|&ascii| ascii.checked_sub(33).expect("phred+33"))
        .collect()
}

// r[verify cram.codec.fqzcomp]
// r[verify cram.codec.fqzcomp.params]
// r[verify cram.codec.fqzcomp.param_block]
// r[verify cram.codec.fqzcomp.array]
// r[verify cram.codec.fqzcomp.selector]
// r[verify cram.codec.fqzcomp.record]
// r[verify cram.codec.fqzcomp.context]
// r[verify cram.codec.fqzcomp.reverse]
#[test]
fn every_hts_specs_fqzcomp_vector_decodes_to_its_original() {
    let mut vectors: Vec<PathBuf> = fs::read_dir(codecs_dir().join("fqzcomp"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect();
    vectors.sort();
    assert_eq!(vectors.len(), 16, "4 inputs x 4 parameter sets: {vectors:?}");

    for path in &vectors {
        let file_name = path.file_name().unwrap().to_str().unwrap();
        let (name, _variant) = file_name.rsplit_once('.').unwrap();
        let want = expected(name);
        let got = fqzcomp::decode(&fs::read(path).unwrap())
            .unwrap_or_else(|e| panic!("{file_name}: {e}"));
        assert_eq!(got.len(), want.len(), "{file_name}: length");
        if let Some(i) = got.iter().zip(&want).position(|(g, w)| g != w) {
            panic!("{file_name}: first difference at byte {i}");
        }
    }
}

#[test]
fn truncated_vectors_error_instead_of_panicking() {
    // qvar.0: variable lengths, dedup, a selector — every per-record step.
    let src = fs::read(codecs_dir().join("fqzcomp").join("qvar.0")).unwrap();
    for cut in (0..src.len()).step_by(src.len() / 40).chain([src.len() - 1]) {
        let prefix = src.get(..cut).unwrap();
        assert!(fqzcomp::decode(prefix).is_err(), "decoded a stream cut at {cut}");
    }
}
