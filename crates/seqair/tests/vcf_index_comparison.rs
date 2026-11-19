//! Compare seqair's TBI index output with what bcftools index generates.
//!
//! Writes a VCF.gz with seqair unified Writer (co-producing a TBI index), then
//! runs `bcftools index` on the same file and compares the two TBI files.
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

#[path = "support/index_parse.rs"]
mod index_parse;
#[path = "../src/pinned.rs"]
mod pinned;
use index_parse::{assert_same_index, parse_index};
use seqair::vcf::record_encoder::{FormatFieldDef, FormatGt, Gt, InfoFieldDef, InfoInt, Scalar};
use seqair::vcf::{
    Alleles, ContigDef, CoordinateIndex, Genotype, IndexKind, OutputFormat, VcfError, VcfHeader,
    Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

fn has_bcftools() -> bool {
    std::process::Command::new("bcftools").arg("--version").output().is_ok()
}

struct IndexSetup {
    header: Arc<VcfHeader>,
    chr1: seqair::vcf::ContigId,
    chr2: seqair::vcf::ContigId,
    dp_info: InfoInt,
    gt_fmt: FormatGt,
}

fn make_index_setup() -> IndexSetup {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(250_000_000))).unwrap();
    let chr2 = builder.register_contig("chr2", ContigDef::new(Some(243_000_000))).unwrap();
    let mut builder = builder.infos();
    let dp_info: InfoInt =
        builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
    let mut builder = builder.formats();
    let gt_fmt: FormatGt =
        builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
    let mut builder = builder.samples();
    builder.add_sample("sample1").unwrap();
    let header = Arc::new(builder.build().unwrap());
    IndexSetup { header, chr1, chr2, dp_info, gt_fmt }
}

/// Write a .vcf.gz file with multiple sorted records and return the co-produced TBI index.
fn write_vcf_gz_with_index(
    dir: &std::path::Path,
) -> (std::path::PathBuf, seqair::io::IndexBuilder) {
    let setup = make_index_setup();
    let vcf_path = dir.join("test.vcf.gz");
    let file = std::fs::File::create(&vcf_path).unwrap();

    let writer = Writer::new(file, OutputFormat::VcfGz).index_kind(IndexKind::Tbi);
    let mut writer = writer.write_header(&setup.header).unwrap();

    let positions = [100u32, 500, 1000, 5000, 20000, 50000, 100000];
    let bases = [(Base::A, Base::T), (Base::C, Base::G), (Base::G, Base::A)];

    for (i, &pos) in positions.iter().enumerate() {
        let (ref_base, alt_base) = bases[i % bases.len()];
        let alleles = Alleles::snv(ref_base, alt_base).unwrap();
        let mut enc = writer
            .begin_record(&setup.chr1, Pos1::new(pos).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, (i as i32 + 1) * 10).unwrap();
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        enc.emit().unwrap();
    }

    // Add a record on chr2
    let alleles = Alleles::snv(Base::T, Base::C).unwrap();
    let mut enc = writer
        .begin_record(&setup.chr2, Pos1::new(10000).unwrap(), &alleles, Some(50.0))
        .unwrap()
        .filter_pass();
    setup.dp_info.encode(&mut enc, 80).unwrap();
    let mut enc = enc.begin_samples();
    setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(1, 1)]).unwrap();
    enc.emit().unwrap();

    let (_inner, index) = writer.finish().unwrap();
    let index = index.expect("VcfGz writer should produce TBI index");
    assert_eq!((index.kind(), index.suffix()), (IndexKind::Tbi, "tbi"));
    (vcf_path, index.into_builder())
}

/// Write seqair's TBI index to a file.
fn write_seqair_tbi(
    index: &seqair::io::IndexBuilder,
    tbi_path: &std::path::Path,
    header: &VcfHeader,
) {
    use seqair_types::SmolStr;
    let names: Vec<SmolStr> = header.contigs().map(|(c, _)| c.name().clone()).collect();
    let file = std::fs::File::create(tbi_path).unwrap();
    index.write_tbi(file, &names).unwrap();
}

/// Run bcftools index to create a TBI for comparison.
fn bcftools_index(vcf_gz_path: &std::path::Path) -> std::path::PathBuf {
    let tbi_path = vcf_gz_path.with_extension("gz.tbi");
    let output = std::process::Command::new("bcftools")
        .args(["index", "-t"])
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
    let setup = make_index_setup();

    let seqair_tbi_path = dir.path().join("seqair.tbi");
    write_seqair_tbi(&index, &seqair_tbi_path, &setup.header);

    let data = decompress_tbi(&seqair_tbi_path);

    assert_eq!(&data[..4], b"TBI\x01", "TBI magic mismatch");

    let n_ref = i32::from_le_bytes([data[4], data[5], data[6], data[7]]);
    assert_eq!(n_ref, 2, "should have 2 references (chr1, chr2)");

    let format = i32::from_le_bytes([data[8], data[9], data[10], data[11]]);
    assert_eq!(format, 2, "format should be 2 (VCF)");

    let col_seq = i32::from_le_bytes([data[12], data[13], data[14], data[15]]);
    assert_eq!(col_seq, 1, "col_seq should be 1");

    let col_beg = i32::from_le_bytes([data[16], data[17], data[18], data[19]]);
    assert_eq!(col_beg, 2, "col_beg should be 2");

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
    let setup = make_index_setup();

    let tbi_path = vcf_path.with_extension("gz.tbi");
    write_seqair_tbi(&index, &tbi_path, &setup.header);

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

    assert!(
        lines.len() >= 2,
        "expected at least 2 records in chr1:1-1000, got {}:\n{}",
        lines.len(),
        text,
    );

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
    let setup = make_index_setup();

    let tbi_path = vcf_path.with_extension("gz.tbi");
    write_seqair_tbi(&index, &tbi_path, &setup.header);

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
    let setup = make_index_setup();

    let seqair_tbi_path = dir.path().join("seqair.tbi");
    write_seqair_tbi(&index, &seqair_tbi_path, &setup.header);

    let bcftools_tbi_path = bcftools_index(&vcf_path);

    let seqair_data = decompress_tbi(&seqair_tbi_path);
    let bcftools_data = decompress_tbi(&bcftools_tbi_path);

    assert_eq!(&seqair_data[..4], &bcftools_data[..4], "magic mismatch");

    let seqair_nref =
        i32::from_le_bytes([seqair_data[4], seqair_data[5], seqair_data[6], seqair_data[7]]);
    let bcftools_nref = i32::from_le_bytes([
        bcftools_data[4],
        bcftools_data[5],
        bcftools_data[6],
        bcftools_data[7],
    ]);
    assert_eq!(seqair_nref, bcftools_nref, "n_ref mismatch");

    assert_eq!(&seqair_data[8..28], &bcftools_data[8..28], "TBI header fields mismatch");

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

use hegel::prelude::*;

/// The four concrete bases a REF or ALT allele can hold.
const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

fn arb_base() -> impl PrintableGenerator<Base> {
    gs::sampled_from(&ACGT).print_as_debug()
}

/// Generate a sorted list of unique positions within [1, `max_pos`].
fn arb_sorted_positions(max_count: usize, max_pos: u32) -> impl PrintableGenerator<Vec<u32>> {
    gs::vecs(gs::integers::<u32>().min_value(1).max_value(max_pos - 1))
        .min_size(1)
        .max_size(max_count - 1)
        .unique(true)
        .map(|mut v: Vec<u32>| {
            v.sort_unstable();
            v
        })
}

/// Write a VCF.gz with the given sorted positions, produce TBI, return paths.
fn write_generated_vcf(
    dir: &std::path::Path,
    positions: &[u32],
    ref_base: Base,
    alt_base: Base,
) -> (std::path::PathBuf, std::path::PathBuf) {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000_000))).unwrap();
    let mut builder = builder.infos();
    let dp_info: InfoInt =
        builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
    let mut builder = builder.formats();
    let gt_fmt: FormatGt =
        builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
    let mut builder = builder.samples();
    builder.add_sample("s1").unwrap();
    let header = Arc::new(builder.build().unwrap());

    let vcf_path = dir.join("test.vcf.gz");
    let file = std::fs::File::create(&vcf_path).unwrap();

    let writer = Writer::new(file, OutputFormat::VcfGz);
    let mut writer = writer.write_header(&header).unwrap();

    for &pos in positions {
        let alleles = Alleles::snv(ref_base, alt_base).unwrap();
        let mut enc = writer
            .begin_record(&chr1, Pos1::new(pos).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        dp_info.encode(&mut enc, 30).unwrap();
        let mut enc = enc.begin_samples();
        gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        enc.emit().unwrap();
    }

    let (_inner, index) = writer.finish().unwrap();
    let index = index.expect("should produce index").into_builder();
    let tbi_path = vcf_path.with_extension("gz.tbi");
    use seqair_types::SmolStr;
    let names: Vec<SmolStr> = header.contigs().map(|(c, _)| c.name().clone()).collect();
    let file = std::fs::File::create(&tbi_path).unwrap();
    index.write_tbi(file, &names).unwrap();

    (vcf_path, tbi_path)
}

/// bcftools can query any region of a seqair-produced VCF.gz using seqair's TBI
/// and returns exactly the records whose positions fall in the query range.
#[hegel::test(test_cases = crate::pinned::cases(20))]
fn bcftools_region_query_matches_expected(tc: TestCase) {
    let positions = tc.draw(arb_sorted_positions(50, 1_000_000));
    let ref_base = tc.draw(arb_base());
    let alt_base = tc.draw(arb_base());
    let query_start = tc.draw(gs::integers::<u32>().min_value(1).max_value(499999));
    let query_span = tc.draw(gs::integers::<u32>().min_value(1).max_value(499999));
    tc.assume(ref_base != alt_base);
    if !has_bcftools() {
        return;
    }

    let query_end = query_start.saturating_add(query_span).min(1_000_000);

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, _tbi_path) = write_generated_vcf(dir.path(), &positions, ref_base, alt_base);

    let output = std::process::Command::new("bcftools")
        .args(["view", "-H", "-r", &format!("chr1:{query_start}-{query_end}")])
        .arg(&vcf_path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "bcftools failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let text = String::from_utf8(output.stdout).unwrap();
    let found_positions: Vec<u32> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').nth(1).unwrap().parse::<u32>().unwrap())
        .collect();

    let expected: Vec<u32> =
        positions.iter().copied().filter(|&p| p >= query_start && p <= query_end).collect();

    assert_eq!(&found_positions, &expected);
}

/// TBI header structure matches bcftools output for any random record set.
#[hegel::test(test_cases = crate::pinned::cases(20))]
fn tbi_header_matches_bcftools(tc: TestCase) {
    let positions = tc.draw(arb_sorted_positions(30, 500_000));
    let ref_base = tc.draw(arb_base());
    let alt_base = tc.draw(arb_base());
    tc.assume(ref_base != alt_base);
    if !has_bcftools() {
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, seqair_tbi_path) =
        write_generated_vcf(dir.path(), &positions, ref_base, alt_base);

    let bcftools_tbi_path = dir.path().join("bcftools.tbi");
    let output = std::process::Command::new("bcftools")
        .args(["index", "-t"])
        .arg(&vcf_path)
        .arg("-o")
        .arg(&bcftools_tbi_path)
        .output()
        .unwrap();
    assert!(output.status.success());

    let seqair_data = decompress_tbi(&seqair_tbi_path);
    let bcftools_data = decompress_tbi(&bcftools_tbi_path);

    assert!(seqair_data.len() >= 32, "seqair TBI too short");
    assert!(bcftools_data.len() >= 32, "bcftools TBI too short");
    assert_eq!(
        &seqair_data[..32],
        &bcftools_data[..32],
        "TBI header mismatch (magic + n_ref + format + columns + meta + skip)"
    );

    let seqair_l_nm =
        i32::from_le_bytes([seqair_data[28], seqair_data[29], seqair_data[30], seqair_data[31]])
            as usize;
    let names_end = 32 + seqair_l_nm;
    assert_eq!(
        &seqair_data[32..names_end],
        &bcftools_data[32..names_end],
        "sequence names mismatch"
    );
}

/// An empty VCF (header only, no records) produces a valid TBI that bcftools accepts.
#[hegel::test(test_cases = crate::pinned::cases(20))]
fn empty_vcf_produces_valid_tbi(tc: TestCase) {
    let _dummy = tc.draw(gs::integers::<u8>().max_value(0));
    if !has_bcftools() {
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let (vcf_path, _tbi_path) = write_generated_vcf(dir.path(), &[], Base::A, Base::T);

    let output = std::process::Command::new("bcftools")
        .args(["view", "-H"])
        .arg(&vcf_path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "bcftools can't read empty VCF: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.trim(), "", "empty VCF should have no records");
}

// ── Whole-index comparison ─────────────────────────────────────────────

/// Write `n` records at sorted positions over a contig long enough that the
/// file spans many BGZF blocks, so `compress_binning` has bins to move; return
/// the VCF.gz path and the index builder the writer co-produced.
fn write_large_vcf(
    dir: &std::path::Path,
    positions: &[u32],
    format: OutputFormat,
    kind: IndexKind,
) -> (std::path::PathBuf, CoordinateIndex) {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(100_000_000))).unwrap();
    let mut builder = builder.infos();
    let dp: InfoInt =
        builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
    let header = Arc::new(builder.formats().samples().build().unwrap());
    let vcf_path = dir.join(if format == OutputFormat::Bcf { "large.bcf" } else { "large.vcf.gz" });
    let mut writer = Writer::new(std::fs::File::create(&vcf_path).unwrap(), format)
        .index_kind(kind)
        .write_header(&header)
        .unwrap();
    for (i, &pos) in positions.iter().enumerate() {
        let alleles = Alleles::snv(Base::A, Base::C).unwrap();
        let mut enc = writer
            .begin_record(&chr1, Pos1::new(pos).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        dp.encode(&mut enc, i as i32).unwrap();
        enc.emit().unwrap();
    }
    let (_inner, index) = writer.finish().unwrap();
    (vcf_path, index.expect("should produce index"))
}

// r[verify index_builder.tbi_format]
// r[verify index_builder.compress_bins]
// r[verify index_builder.no_coor_trailer]
// r[verify vcf_writer.index_amend_last]
// r[verify csi.write_format]
// r[verify record_encoder.index_kind]
// r[verify record_encoder.finish+2]
/// The indexes seqair co-produces hold exactly what `bcftools index` builds
/// for the same file: every bin with its chunks (after `compress_binning`),
/// the TBI linear index, the CSI `loffset`s, and the `n_no_coor` trailer —
/// TBI for VCF.gz (`-t`, the writer set to `IndexKind::Tbi`, the index
/// serialized as `Writer::finish` returns it) and CSI for BCF (`-c`). (A VCF.gz CSI is not
/// compared: htslib's tabix path gives it at least six levels at
/// `min_shift` 14, seqair the fewest that cover the contigs.)
#[hegel::test(test_cases = crate::pinned::cases(10))]
fn whole_tbi_and_csi_match_bcftools(tc: TestCase) {
    if !has_bcftools() {
        return;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(20_000));
    let span = tc.draw(gs::integers::<u32>().min_value(1).max_value(99_000_000));
    let mut positions: Vec<u32> =
        (0..n).map(|_| tc.draw(gs::integers::<u32>().min_value(1).max_value(span))).collect();
    positions.sort_unstable();

    for format in [OutputFormat::VcfGz, OutputFormat::Bcf] {
        let dir = tempfile::tempdir().unwrap();
        let csi = format == OutputFormat::Bcf;
        let kind = if csi { IndexKind::Csi } else { IndexKind::Tbi };
        let (path, index) = write_large_vcf(dir.path(), &positions, format, kind);
        assert_eq!(index.kind(), kind);
        let ours = dir.path().join("ours.idx");
        index.write_to_path(&ours).unwrap();
        let theirs = dir.path().join("theirs.idx");
        let out = std::process::Command::new("bcftools")
            .args(["index", if csi { "-c" } else { "-t" }, "-o"])
            .arg(&theirs)
            .arg(&path)
            .output()
            .unwrap();
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        let ours = parse_index(&decompress_tbi(&ours), csi);
        let theirs = parse_index(&decompress_tbi(&theirs), csi);
        assert_same_index(&ours, &theirs, &format!("{format:?}"));
    }
}

// ── Choosing the index ─────────────────────────────────────────────────

// r[verify record_encoder.index_kind]
// r[verify record_encoder.writer_new+2]
/// VCF.gz is indexed by CSI unless told otherwise, and BCF by CSI only: a
/// writer asked for a TBI of a BCF refuses at `write_header`, before writing.
#[test]
fn csi_is_the_default_and_bcf_has_no_tbi() {
    let setup = make_index_setup();
    let (_, index) = Writer::new(Vec::new(), OutputFormat::VcfGz)
        .write_header(&setup.header)
        .unwrap()
        .finish()
        .unwrap();
    let index = index.expect("VCF.gz is indexed");
    assert_eq!((index.kind(), index.suffix()), (IndexKind::Csi, "csi"));
    assert_eq!(IndexKind::Tbi.suffix(), "tbi");

    let mut out = Vec::new();
    let err = Writer::new(&mut out, OutputFormat::Bcf)
        .index_kind(IndexKind::Tbi)
        .write_header(&setup.header)
        .err();
    assert!(matches!(err, Some(VcfError::TbiForBcf)), "{err:?}");
    assert!(out.is_empty(), "nothing is written");

    let (_, index) = Writer::new(Vec::new(), OutputFormat::Vcf)
        .index_kind(IndexKind::Tbi)
        .write_header(&setup.header)
        .unwrap()
        .finish()
        .unwrap();
    assert!(index.is_none(), "plain VCF is not indexed, whatever the kind");
}

// r[verify record_encoder.index_kind]
/// The TBI the writer returns is read by `tabix` (which lists the contigs
/// with records and answers a region from it) and by seqair's own TBI
/// parser, which finds the same chunks in it as in the TBI `bcftools index
/// -t` builds for the same file.
#[test]
fn tabix_and_seqair_read_the_returned_tbi() {
    if !has_bcftools() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let positions: Vec<u32> = (1..4000).map(|i| i * 997).collect();
    let (path, index) =
        write_large_vcf(dir.path(), &positions, OutputFormat::VcfGz, IndexKind::Tbi);
    let ours = dir.path().join(format!("large.vcf.gz.{}", index.suffix()));
    index.write_to_path(&ours).unwrap();

    let tabix = |args: &[&str]| {
        let out = std::process::Command::new("tabix").args(args).arg(&path).output().unwrap();
        assert!(out.status.success(), "tabix: {}", String::from_utf8_lossy(&out.stderr));
        String::from_utf8(out.stdout).unwrap()
    };
    assert_eq!(tabix(&["-l"]), "chr1\n");
    let mut region = std::process::Command::new("tabix");
    let out = region.arg(&path).arg("chr1:100000-200000").output().unwrap();
    assert!(out.status.success(), "tabix: {}", String::from_utf8_lossy(&out.stderr));
    let found: Vec<u32> = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect();
    let expected: Vec<u32> =
        positions.iter().copied().filter(|p| (100_000..=200_000).contains(p)).collect();
    assert_eq!(found, expected);

    let theirs = dir.path().join("theirs.tbi");
    let out = std::process::Command::new("bcftools")
        .args(["index", "-t", "-o"])
        .arg(&theirs)
        .arg(&path)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let bam_header = seqair::bam::BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:100000000\n").unwrap();
    let parse =
        |p: &std::path::Path| seqair::bam::BamIndex::from_tabix_path(p, &bam_header).unwrap();
    let (ours, theirs) = (parse(&ours), parse(&theirs));
    for (start, end) in
        [(0, 100_000_000), (5_000, 6_000), (1_000_000, 3_000_000), (3_900_000, 4_000_000)]
    {
        let q = |i: &seqair::bam::BamIndex| {
            i.query(0, seqair::bam::Pos0::new(start).unwrap(), seqair::bam::Pos0::new(end).unwrap())
        };
        assert_eq!(format!("{:?}", q(&ours)), format!("{:?}", q(&theirs)), "{start}..{end}");
    }
}

// r[verify index_builder.range_check]
/// A TBI stops at 2^29: a record ending past it fails its `emit`, as htslib
/// refuses it, where a CSI (the default) indexes it.
#[test]
fn a_tbi_refuses_a_record_past_its_bins() {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1 << 30))).unwrap();
    let header = builder.build().unwrap();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();
    let far = Pos1::new((1 << 29) + 10).unwrap();
    for kind in [IndexKind::Csi, IndexKind::Tbi] {
        let mut writer = Writer::new(Vec::new(), OutputFormat::VcfGz)
            .index_kind(kind)
            .write_header(&header)
            .unwrap();
        let result = writer.begin_record(&chr1, far, &alleles, None).unwrap().filter_pass().emit();
        match kind {
            IndexKind::Csi => result.unwrap(),
            IndexKind::Tbi => assert!(
                matches!(
                    result,
                    Err(VcfError::Index(seqair::io::IndexError::BeyondIndexRange { .. }))
                ),
                "{result:?}"
            ),
        }
    }
}
