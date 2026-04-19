//! Generate seeds for `fuzz_reader_indexed` from real test data.
//!
//! Uses `Input::encode` for exact round-trip fidelity — every byte of the seed
//! is used by the fuzz target, no waste.
//!
//! Run: cd crates/seqair/fuzz && cargo run --example gen_indexed_seeds

use seqair_fuzz::indexed_reader::Input;
use std::fs;
use std::io::Read;
use std::path::Path;

fn main() {
    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/data");
    let seed_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("seeds/fuzz_reader_indexed");
    fs::create_dir_all(&seed_dir).unwrap();

    let fai = fs::read(data_dir.join("test.fasta.gz.fai")).unwrap();
    let gzi = fs::read(data_dir.join("test.fasta.gz.gzi")).unwrap();

    // === BAM seeds ===
    let bam = fs::read(data_dir.join("test.bam")).unwrap();
    let bai = fs::read(data_dir.join("test.bam.bai")).unwrap();

    // Small BAM: 5 BGZF blocks — exercises header + first few records
    let small_bam = truncate_bgzf(&bam, 5);
    write_seed(&seed_dir, "bam_small", Input::encode(0, &small_bam, &bai, &[], &fai, &gzi));

    // Larger BAM: 20 blocks — deeper record coverage
    let medium_bam = truncate_bgzf(&bam, 20);
    write_seed(&seed_dir, "bam_medium", Input::encode(0, &medium_bam, &bai, &[], &fai, &gzi));

    // Header-only BAM: just enough blocks for the header, no records — exercises empty fetch
    let header_bam = truncate_bgzf(&bam, 2);
    write_seed(
        &seed_dir,
        "bam_header_only",
        Input::encode(0, &header_bam, &bai, &[], &fai, &gzi),
    );

    // BAM with empty BAI (magic + n_ref=0) — exercises index-mismatch path
    let empty_bai = b"BAI\x01\x00\x00\x00\x00";
    write_seed(
        &seed_dir,
        "bam_empty_bai",
        Input::encode(0, &small_bam, empty_bai, &[], &fai, &gzi),
    );

    // BAM without FASTA — exercises the no-reference path
    write_seed(&seed_dir, "bam_no_fasta", Input::encode(0, &small_bam, &bai, &[], b"", &[]));

    // === SAM-BGZF seed (format=1) ===
    let sam_text = b"@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t**********\n\
        read2\t0\tchr1\t200\t40\t5M1I4M\t*\t0\t0\tACGTACGTAC\t**********\n\
        read3\t0\tchr1\t300\t50\t3M2D7M\t*\t0\t0\tACGTACGTAC\t**********\n";
    let bgzf_sam = make_bgzf_block(sam_text);
    // TBI: magic + tabix header (n_ref=1, format=SAM=1, columns, meta=#, skip=0, l_nm=5, "chr1\0")
    // + one ref with zero bins (exercises empty-query path)
    let mut tbi = Vec::new();
    tbi.extend_from_slice(b"TBI\x01");
    tbi.extend_from_slice(&1i32.to_le_bytes()); // n_ref
    tbi.extend_from_slice(&1i32.to_le_bytes()); // format (SAM=1)
    tbi.extend_from_slice(&1i32.to_le_bytes()); // col_seq
    tbi.extend_from_slice(&2i32.to_le_bytes()); // col_beg
    tbi.extend_from_slice(&0i32.to_le_bytes()); // col_end (0 = compute from col_beg)
    tbi.extend_from_slice(&(b'#' as i32).to_le_bytes()); // meta
    tbi.extend_from_slice(&0i32.to_le_bytes()); // skip
    tbi.extend_from_slice(&5i32.to_le_bytes()); // l_nm
    tbi.extend_from_slice(b"chr1\0"); // names
    tbi.extend_from_slice(&0i32.to_le_bytes()); // n_bin for ref 0
    write_seed(&seed_dir, "sam_bgzf", Input::encode(1, &bgzf_sam, &tbi, &[], &fai, &gzi));

    // === Plain SAM seed (format=3) ===
    let plain_sam = b"@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n";
    write_seed(
        &seed_dir,
        "plain_sam",
        Input::encode(3, plain_sam, &[], &[], &fai, &gzi),
    );

    // Plain SAM with multiple contigs
    let plain_sam_multi = b"@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:2000\n\
        read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n\
        read2\t0\tchr2\t500\t40\t10M\t*\t0\t0\tACGTACGTAC\t*\n";
    write_seed(
        &seed_dir,
        "plain_sam_multi_contig",
        Input::encode(3, plain_sam_multi, &[], &[], &fai, &gzi),
    );

    // Plain SAM without FASTA
    write_seed(
        &seed_dir,
        "plain_sam_no_fasta",
        Input::encode(3, plain_sam, &[], &[], b"", &[]),
    );

    // === CRAM seeds ===
    let crai_v31_gz = fs::read(data_dir.join("test.cram.crai")).unwrap();
    let crai_v31 = decompress_gz(&crai_v31_gz);

    // CRAM v3.1
    let cram_v31 = fs::read(data_dir.join("test.cram")).unwrap();
    let small_cram_v31 = &cram_v31[..cram_v31.len().min(50_000)];
    write_seed(
        &seed_dir,
        "cram_v31",
        Input::encode(2, small_cram_v31, &crai_v31, &[], &fai, &gzi),
    );

    // CRAM v3.0
    let cram_v30 = fs::read(data_dir.join("test_v30.cram")).unwrap();
    let crai_v30_gz = fs::read(data_dir.join("test_v30.cram.crai")).unwrap();
    let crai_v30 = decompress_gz(&crai_v30_gz);
    let small_cram_v30 = &cram_v30[..cram_v30.len().min(50_000)];
    write_seed(
        &seed_dir,
        "cram_v30",
        Input::encode(2, small_cram_v30, &crai_v30, &[], &fai, &gzi),
    );

    // CRAM gzip codec
    let cram_gzip = fs::read(data_dir.join("test_gzip.cram")).unwrap();
    let crai_gzip_gz = fs::read(data_dir.join("test_gzip.cram.crai")).unwrap();
    let crai_gzip = decompress_gz(&crai_gzip_gz);
    let small_cram_gzip = &cram_gzip[..cram_gzip.len().min(50_000)];
    write_seed(
        &seed_dir,
        "cram_gzip",
        Input::encode(2, small_cram_gzip, &crai_gzip, &[], &fai, &gzi),
    );

    // CRAM without FASTA
    write_seed(
        &seed_dir,
        "cram_no_fasta",
        Input::encode(2, small_cram_v31, &crai_v31, &[], b"", &[]),
    );
}

fn write_seed(seed_dir: &Path, name: &str, data: Vec<u8>) {
    let parsed = Input::parse(&data).unwrap_or_else(|| panic!("{name} seed must round-trip parse"));
    let total = parsed.data1.len() + parsed.data2.len() + parsed.fasta_gz.len() + parsed.fai.len()
        + parsed.gzi.len();
    fs::write(seed_dir.join(name), &data).unwrap();
    println!("  {name}: {} bytes ({total} payload)", data.len());
}

fn truncate_bgzf(data: &[u8], max_blocks: usize) -> Vec<u8> {
    let mut pos = 0;
    let mut blocks = 0;
    while pos < data.len() && blocks < max_blocks {
        if data.get(pos..pos + 2) != Some(&[0x1f, 0x8b]) {
            break;
        }
        let xlen = u16::from_le_bytes([data[pos + 10], data[pos + 11]]) as usize;
        let extra = &data[pos + 12..pos + 12 + xlen];
        let mut ep = 0;
        let mut bsize = None;
        while ep + 4 <= extra.len() {
            let (si1, si2) = (extra[ep], extra[ep + 1]);
            let slen = u16::from_le_bytes([extra[ep + 2], extra[ep + 3]]) as usize;
            if si1 == b'B' && si2 == b'C' && slen == 2 {
                bsize = Some(u16::from_le_bytes([extra[ep + 4], extra[ep + 5]]));
                break;
            }
            ep += 4 + slen;
        }
        let Some(bsize) = bsize else { break };
        pos += bsize as usize + 1;
        blocks += 1;
    }
    data[..pos].to_vec()
}

fn decompress_gz(data: &[u8]) -> Vec<u8> {
    let mut decoder = flate2::read::GzDecoder::new(data);
    let mut out = Vec::new();
    decoder.read_to_end(&mut out).unwrap_or(0);
    out
}

/// Build a single BGZF block from uncompressed data, followed by an EOF block.
fn make_bgzf_block(data: &[u8]) -> Vec<u8> {
    use flate2::write::DeflateEncoder;
    use std::io::Write;

    let mut encoder = DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(data).unwrap();
    let compressed = encoder.finish().unwrap();

    let bsize = (18 + compressed.len() + 8 - 1) as u16;
    let mut block = Vec::with_capacity(18 + compressed.len() + 8 + 28);

    // BGZF header (18 bytes)
    block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]); // gzip magic + DEFLATE + FEXTRA
    block.extend_from_slice(&[0; 4]); // MTIME
    block.push(0); // XFL
    block.push(0xff); // OS
    block.extend_from_slice(&6u16.to_le_bytes()); // XLEN = 6
    block.extend_from_slice(&[b'B', b'C', 2, 0]); // BC subfield, SLEN=2
    block.extend_from_slice(&bsize.to_le_bytes()); // BSIZE

    // Compressed data
    block.extend_from_slice(&compressed);

    // Footer: CRC32 + ISIZE
    let crc = crc32(data);
    block.extend_from_slice(&crc.to_le_bytes());
    block.extend_from_slice(&(data.len() as u32).to_le_bytes());

    // EOF block
    block.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
        0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ]);

    block
}

fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in data {
        crc ^= byte as u32;
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ 0xEDB8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}
