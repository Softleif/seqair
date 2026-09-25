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

/// ITF8 as in CRAM3 §2.3, for values below 2^28.
fn itf8(v: u32) -> Vec<u8> {
    let [b0, b1, b2, b3] = v.to_be_bytes();
    match v {
        0..0x80 => vec![b3],
        0x80..0x4000 => vec![0x80 | b2, b3],
        0x4000..0x20_0000 => vec![0xC0 | b1, b2, b3],
        0x20_0000..0x1000_0000 => vec![0xE0 | b0, b1, b2, b3],
        _ => panic!("test values stay below 2^28"),
    }
}

/// An external-data CRAM block holding `data` compressed with method 7.
fn fqzcomp_block(data: &[u8], uncompressed_size: u32) -> Vec<u8> {
    let mut block = vec![7, 4];
    block.extend(itf8(11)); // content id
    block.extend(itf8(u32::try_from(data.len()).unwrap()));
    block.extend(itf8(uncompressed_size));
    block.extend_from_slice(data);
    let mut crc = libdeflater::Crc::new();
    crc.update(&block);
    block.extend(crc.sum().to_le_bytes());
    block
}

// r[verify cram.codec.fqzcomp]
#[test]
fn block_method_7_decodes_and_checks_the_header_size() {
    use seqair::cram::{CramError, block::parse_block};

    let data = fs::read(codecs_dir().join("fqzcomp").join("q8.2")).unwrap();
    let want = expected("q8");
    let size = u32::try_from(want.len()).unwrap();

    let bytes = fqzcomp_block(&data, size);
    let (block, used) = parse_block(&bytes).unwrap();
    assert_eq!(block.content_id, 11);
    assert_eq!(block.data, want);
    assert_eq!(used, bytes.len());

    let err = parse_block(&fqzcomp_block(&data, size - 1)).unwrap_err();
    assert!(
        matches!(err, CramError::FqzcompSizeMismatch { expected, found }
            if expected == want.len() - 1 && found == want.len()),
        "{err:?}"
    );
}
