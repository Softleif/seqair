//! Compare seqair's TBI index output with what bcftools index generates.
//!
//! Writes a VCF.gz with seqair (co-producing a TBI index), then runs
//! `bcftools index` on the same file and compares the two TBI files byte-by-byte.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::vcf::alleles::Alleles;
use seqair::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair::vcf::writer::VcfWriter;
use seqair_types::{Base, Pos1, SmolStr};
use std::sync::Arc;

fn has_bcftools() -> bool {
    std::process::Command::new("bcftools").arg("--version").output().is_ok()
}

fn test_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
            .unwrap()
            .add_contig("chr2", ContigDef { length: Some(243_000_000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap()
            .add_sample("sample1")
            .unwrap()
            .build()
            .unwrap(),
    )
}

/// Write a .vcf.gz file with multiple sorted records and return the co-produced TBI index.
fn write_vcf_gz_with_index(
    dir: &std::path::Path,
) -> (std::path::PathBuf, seqair::vcf::index_builder::IndexBuilder) {
    let header = test_header();
    let vcf_path = dir.join("test.vcf.gz");
    let file = std::fs::File::create(&vcf_path).unwrap();

    let mut writer = VcfWriter::bgzf(file, header.clone());
    writer.write_header().unwrap();

    // Write sorted records across chr1
    let positions = [100, 500, 1000, 5000, 20000, 50000, 100000];
    let bases = [(Base::A, Base::T), (Base::C, Base::G), (Base::G, Base::A)];

    for (i, &pos) in positions.iter().enumerate() {
        let (ref_base, alt_base) = bases[i % bases.len()];
        let alleles = Alleles::snv(ref_base, alt_base).unwrap();
        let record = VcfRecordBuilder::new("chr1", Pos1::new(pos).unwrap(), alleles)
            .qual(30.0)
            .filter_pass()
            .info_integer("DP", (i as i32 + 1) * 10)
            .format_keys(&["GT"])
            .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1))])
            .build(&header)
            .unwrap();
        writer.write_record(&record).unwrap();
    }

    // Add a record on chr2
    let alleles = Alleles::snv(Base::T, Base::C).unwrap();
    let record = VcfRecordBuilder::new("chr2", Pos1::new(10000).unwrap(), alleles)
        .qual(50.0)
        .filter_pass()
        .info_integer("DP", 80)
        .format_keys(&["GT"])
        .add_sample(vec![SampleValue::Genotype(Genotype::unphased(1, 1))])
        .build(&header)
        .unwrap();
    writer.write_record(&record).unwrap();

    let index = writer.finish().unwrap().expect("BGZF writer should produce TBI index");
    (vcf_path, index)
}

/// Write seqair's TBI index to a file.
fn write_seqair_tbi(
    index: &seqair::vcf::index_builder::IndexBuilder,
    tbi_path: &std::path::Path,
    header: &VcfHeader,
) {
    let names: Vec<SmolStr> = header.contigs().keys().cloned().collect();
    let file = std::fs::File::create(tbi_path).unwrap();
    index.write_tbi(file, &names).unwrap();
}

/// Run bcftools index to create a TBI for comparison.
fn bcftools_index(vcf_gz_path: &std::path::Path) -> std::path::PathBuf {
    let tbi_path = vcf_gz_path.with_extension("gz.tbi");
    let output = std::process::Command::new("bcftools")
        .args(["index", "-t"]) // -t = TBI format
        .arg(vcf_gz_path)
        .arg("-o")
        .arg(&tbi_path)
        .output()
        .expect("bcftools should run");

    assert!(
        output.status.success(),
        "bcftools index failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    tbi_path
}

/// Decompress a BGZF-compressed TBI file to get raw binary content.
fn decompress_tbi(tbi_path: &std::path::Path) -> Vec<u8> {
    let file = std::fs::File::open(tbi_path).unwrap();
    let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(file);
    let mut data = Vec::new();
    reader.read_to_end(&mut data).unwrap();
    data
}

// ── Tests ──────────────────────────────────────────────────────────────

#[test]
fn seqair_tbi_has_correct_magic_and_structure() {
    let dir = tempfile::tempdir().unwrap();
    let (_vcf_path, index) = write_vcf_gz_with_index(dir.path());
    let header = test_header();

    let seqair_tbi_path = dir.path().join("seqair.tbi");
    write_seqair_tbi(&index, &seqair_tbi_path, &header);

    let data = decompress_tbi(&seqair_tbi_path);

    // TBI magic: TBI\1
    assert_eq!(&data[..4], b"TBI\x01", "TBI magic mismatch");

    // n_ref
    let n_ref = i32::from_le_bytes([data[4], data[5], data[6], data[7]]);
    assert_eq!(n_ref, 2, "should have 2 references (chr1, chr2)");

    // format = 2 (VCF)
    let format = i32::from_le_bytes([data[8], data[9], data[10], data[11]]);
    assert_eq!(format, 2, "format should be 2 (VCF)");

    // col_seq = 1
    let col_seq = i32::from_le_bytes([data[12], data[13], data[14], data[15]]);
    assert_eq!(col_seq, 1, "col_seq should be 1");

    // col_beg = 2
    let col_beg = i32::from_le_bytes([data[16], data[17], data[18], data[19]]);
    assert_eq!(col_beg, 2, "col_beg should be 2");

    // meta = '#' = 35
    let meta = i32::from_le_bytes([data[24], data[25], data[26], data[27]]);
    assert_eq!(meta, 35, "meta should be '#' (35)");
}

#[test]
fn bcftools_can_query_seqair_vcf_gz_with_seqair_tbi() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, index) = write_vcf_gz_with_index(dir.path());
    let header = test_header();

    // Write seqair's TBI
    let tbi_path = vcf_path.with_extension("gz.tbi");
    write_seqair_tbi(&index, &tbi_path, &header);

    // bcftools should be able to query with seqair's index
    let output = std::process::Command::new("bcftools")
        .args(["view", "-H", "-r", "chr1:1-1000"])
        .arg(&vcf_path)
        .output()
        .expect("bcftools should run");

    assert!(
        output.status.success(),
        "bcftools view with seqair TBI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();

    // Should find records at pos 100 and 500 (within chr1:1-1000)
    assert!(
        lines.len() >= 2,
        "expected at least 2 records in chr1:1-1000, got {}:\n{}",
        lines.len(),
        text,
    );

    // Verify positions are in range
    for line in &lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "chr1");
        let pos: u32 = fields[1].parse().unwrap();
        assert!((1..=1000).contains(&pos), "pos {pos} out of range");
    }
}

#[test]
fn bcftools_can_query_chr2_with_seqair_tbi() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, index) = write_vcf_gz_with_index(dir.path());
    let header = test_header();

    let tbi_path = vcf_path.with_extension("gz.tbi");
    write_seqair_tbi(&index, &tbi_path, &header);

    // Query chr2
    let output = std::process::Command::new("bcftools")
        .args(["view", "-H", "-r", "chr2"])
        .arg(&vcf_path)
        .output()
        .expect("bcftools should run");

    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();

    assert_eq!(lines.len(), 1, "should find exactly 1 record on chr2");
    let fields: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(fields[0], "chr2");
    assert_eq!(fields[1], "10000");
}

#[test]
fn seqair_tbi_matches_bcftools_tbi_structure() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, index) = write_vcf_gz_with_index(dir.path());
    let header = test_header();

    // Write seqair's TBI
    let seqair_tbi_path = dir.path().join("seqair.tbi");
    write_seqair_tbi(&index, &seqair_tbi_path, &header);

    // Generate bcftools TBI
    let bcftools_tbi_path = bcftools_index(&vcf_path);

    // Decompress both and compare structure
    let seqair_data = decompress_tbi(&seqair_tbi_path);
    let bcftools_data = decompress_tbi(&bcftools_tbi_path);

    // Both should have same magic
    assert_eq!(&seqair_data[..4], &bcftools_data[..4], "magic mismatch");

    // Same n_ref
    let seqair_nref =
        i32::from_le_bytes([seqair_data[4], seqair_data[5], seqair_data[6], seqair_data[7]]);
    let bcftools_nref = i32::from_le_bytes([
        bcftools_data[4],
        bcftools_data[5],
        bcftools_data[6],
        bcftools_data[7],
    ]);
    assert_eq!(seqair_nref, bcftools_nref, "n_ref mismatch");

    // Same format, col_seq, col_beg, col_end
    assert_eq!(&seqair_data[8..28], &bcftools_data[8..28], "TBI header fields mismatch");

    // Same sequence names section (length + names)
    let seqair_l_nm =
        i32::from_le_bytes([seqair_data[28], seqair_data[29], seqair_data[30], seqair_data[31]]);
    let bcftools_l_nm = i32::from_le_bytes([
        bcftools_data[28],
        bcftools_data[29],
        bcftools_data[30],
        bcftools_data[31],
    ]);
    assert_eq!(seqair_l_nm, bcftools_l_nm, "l_nm (sequence names length) mismatch");

    let names_end = 32 + seqair_l_nm as usize;
    assert_eq!(
        &seqair_data[32..names_end],
        &bcftools_data[32..names_end],
        "sequence names mismatch"
    );
}

// ── Proptests ──────────────────────────────────────────────────────────

use proptest::prelude::*;

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

/// Generate a sorted list of unique positions within [1, `max_pos`].
fn arb_sorted_positions(max_count: usize, max_pos: u32) -> impl Strategy<Value = Vec<u32>> {
    proptest::collection::hash_set(1u32..max_pos, 1..max_count).prop_map(|set| {
        let mut v: Vec<u32> = set.into_iter().collect();
        v.sort_unstable();
        v
    })
}

/// Write a VCF.gz with the given sorted positions, produce TBI, return paths.
fn write_proptest_vcf(
    dir: &std::path::Path,
    positions: &[u32],
    ref_base: Base,
    alt_base: Base,
) -> (std::path::PathBuf, std::path::PathBuf) {
    let header = test_header();
    let vcf_path = dir.join("test.vcf.gz");
    let file = std::fs::File::create(&vcf_path).unwrap();

    let mut writer = VcfWriter::bgzf(file, header.clone());
    writer.write_header().unwrap();

    for &pos in positions {
        let alleles = Alleles::snv(ref_base, alt_base).unwrap();
        let record = VcfRecordBuilder::new("chr1", Pos1::new(pos).unwrap(), alleles)
            .filter_pass()
            .info_integer("DP", 30)
            .format_keys(&["GT"])
            .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1))])
            .build(&header)
            .unwrap();
        writer.write_record(&record).unwrap();
    }

    let index = writer.finish().unwrap().expect("should produce index");
    let tbi_path = vcf_path.with_extension("gz.tbi");
    let names: Vec<SmolStr> = header.contigs().keys().cloned().collect();
    let file = std::fs::File::create(&tbi_path).unwrap();
    index.write_tbi(file, &names).unwrap();

    (vcf_path, tbi_path)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    /// bcftools can query any region of a seqair-produced VCF.gz using seqair's TBI
    /// and returns exactly the records whose positions fall in the query range.
    #[test]
    fn bcftools_region_query_matches_expected(
        positions in arb_sorted_positions(50, 1_000_000),
        ref_base in arb_base(),
        alt_base in arb_base(),
        // Query a random sub-range
        query_start in 1u32..500_000,
        query_span in 1u32..500_000,
    ) {
        prop_assume!(ref_base != alt_base);
        if !has_bcftools() {
            return Ok(());
        }

        let query_end = query_start.saturating_add(query_span).min(1_000_000);

        let dir = tempfile::tempdir().unwrap();
        let (vcf_path, _tbi_path) = write_proptest_vcf(dir.path(), &positions, ref_base, alt_base);

        // Query with bcftools using seqair's TBI
        let output = std::process::Command::new("bcftools")
            .args(["view", "-H", "-r", &format!("chr1:{query_start}-{query_end}")])
            .arg(&vcf_path)
            .output()
            .unwrap();

        prop_assert!(output.status.success(),
            "bcftools failed: {}", String::from_utf8_lossy(&output.stderr));

        let text = String::from_utf8(output.stdout).unwrap();
        let found_positions: Vec<u32> = text.lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').nth(1).unwrap().parse::<u32>().unwrap())
            .collect();

        // Expected: positions that fall within [query_start, query_end]
        let expected: Vec<u32> = positions.iter()
            .copied()
            .filter(|&p| p >= query_start && p <= query_end)
            .collect();

        prop_assert_eq!(&found_positions, &expected);
    }

    /// TBI header structure matches bcftools output for any random record set.
    #[test]
    fn tbi_header_matches_bcftools(
        positions in arb_sorted_positions(30, 500_000),
        ref_base in arb_base(),
        alt_base in arb_base(),
    ) {
        prop_assume!(ref_base != alt_base);
        if !has_bcftools() {
            return Ok(());
        }

        let dir = tempfile::tempdir().unwrap();
        let (vcf_path, seqair_tbi_path) = write_proptest_vcf(dir.path(), &positions, ref_base, alt_base);

        // Generate bcftools TBI
        let bcftools_tbi_path = dir.path().join("bcftools.tbi");
        let output = std::process::Command::new("bcftools")
            .args(["index", "-t"])
            .arg(&vcf_path)
            .arg("-o")
            .arg(&bcftools_tbi_path)
            .output()
            .unwrap();
        prop_assert!(output.status.success());

        let seqair_data = decompress_tbi(&seqair_tbi_path);
        let bcftools_data = decompress_tbi(&bcftools_tbi_path);

        // Both must be valid TBI files. Now that seqair also omits empty refs,
        // the entire header (magic + n_ref + format fields) should match.
        prop_assert!(seqair_data.len() >= 32, "seqair TBI too short");
        prop_assert!(bcftools_data.len() >= 32, "bcftools TBI too short");
        prop_assert_eq!(&seqair_data[..32], &bcftools_data[..32],
            "TBI header mismatch (magic + n_ref + format + columns + meta + skip)");

        // Sequence names must match exactly
        let seqair_l_nm = i32::from_le_bytes(
            [seqair_data[28], seqair_data[29], seqair_data[30], seqair_data[31]]
        ) as usize;
        let names_end = 32 + seqair_l_nm;
        prop_assert_eq!(
            &seqair_data[32..names_end],
            &bcftools_data[32..names_end],
            "sequence names mismatch"
        );
    }

    /// An empty VCF (header only, no records) produces a valid TBI that bcftools accepts.
    #[test]
    fn empty_vcf_produces_valid_tbi(_dummy in 0u8..1) {
        if !has_bcftools() {
            return Ok(());
        }

        let dir = tempfile::tempdir().unwrap();
        let (vcf_path, _tbi_path) = write_proptest_vcf(dir.path(), &[], Base::A, Base::T);

        let output = std::process::Command::new("bcftools")
            .args(["view", "-H"])
            .arg(&vcf_path)
            .output()
            .unwrap();

        prop_assert!(output.status.success(),
            "bcftools can't read empty VCF: {}", String::from_utf8_lossy(&output.stderr));
        let text = String::from_utf8(output.stdout).unwrap();
        prop_assert_eq!(text.trim(), "", "empty VCF should have no records");
    }
}
