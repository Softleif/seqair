//! Pileup cross-validation: run seqair's pileup engine on htslib's mpileup
//! test SAMs and compare depth + base/deletion/refskip counts against
//! rust-htslib's `bam_plp_auto` (the reference per `r[pileup.htslib_compat]`).
//!
//! htslib's mpileup/ test files exercise specific CIGAR combinations:
//! - `mp_D`: deletions (D ops)
//! - `mp_I`: insertions (I ops)
//! - `mp_DI`: deletion then insertion
//! - `mp_ID`: insertion then deletion
//! - `mp_N`: ref-skips (N ops, RNA-seq introns)
//! - `mp_N2`: complex N+I+D+P combinations
//! - `mp_P`: padding (P ops) with insertions
//! - `mp_overlap1/2`: overlapping mate pairs
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

use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use seqair::bam::{IndexedBamReader, PileupEngine, Pos0, RecordStore};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn mpileup_dir() -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/mpileup/")).to_path_buf()
}

/// Convert SAM to sorted, indexed BAM.
fn sam_to_bam(dir: &Path, sam_path: &Path) -> PathBuf {
    let bam_path = dir.join("test.bam");
    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools sort failed");
    let status = Command::new("samtools")
        .arg("index")
        .arg(&bam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index failed");
    assert!(status.success());
    bam_path
}

struct HtsColumn {
    pos: u32,
    depth: u32,
    n_del: usize,
    n_refskip: usize,
    n_ins: usize,
}

/// Read pileup from rust-htslib (wraps `bam_plp_auto`).
fn htslib_pileup(bam_path: &Path, contig: &[u8], start: u64, end: u64) -> Vec<HtsColumn> {
    let mut reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let tid = reader.header().tid(contig).expect("tid");
    reader
        .fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64))
        .expect("htslib fetch");

    let mut columns = Vec::new();
    for p in reader.pileup() {
        let p = p.expect("htslib pileup");
        let pos = p.pos();
        if u64::from(pos) < start || u64::from(pos) >= end {
            continue;
        }
        let mut n_del = 0;
        let mut n_refskip = 0;
        let mut n_ins = 0;
        for a in p.alignments() {
            if a.is_del() {
                n_del += 1;
            }
            if a.is_refskip() {
                n_refskip += 1;
            }
            if matches!(a.indel(), Indel::Ins(_)) {
                n_ins += 1;
            }
        }
        columns.push(HtsColumn { pos, depth: p.depth(), n_del, n_refskip, n_ins });
    }
    columns
}

/// Read pileup from seqair.
fn seqair_pileup(bam_path: &Path, contig: &str, start: u32, end: u32) -> Vec<HtsColumn> {
    let mut reader = IndexedBamReader::open(bam_path).expect("seqair open");
    let tid = reader.header().tid(contig).expect("contig not found");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap(), &mut store)
        .expect("fetch");

    let engine = PileupEngine::new(store, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());

    engine
        .map(|col| {
            let mut n_del = 0;
            let mut n_refskip = 0;
            let mut n_ins = 0;
            for a in col.alignments() {
                match a.op {
                    seqair::bam::PileupOp::Deletion { .. } => n_del += 1,
                    // ComplexIndel counts as both a deletion and an insertion
                    // (htslib: is_del=true + indel>0). N→I also sets is_refskip.
                    seqair::bam::PileupOp::ComplexIndel { is_refskip, .. } => {
                        n_del += 1;
                        n_ins += 1;
                        if is_refskip {
                            n_refskip += 1;
                        }
                    }
                    // htslib sets is_del=true for BOTH D and N ops (meaning
                    // "no query base"). is_refskip is the additional flag for N.
                    seqair::bam::PileupOp::RefSkip => {
                        n_del += 1;
                        n_refskip += 1;
                    }
                    seqair::bam::PileupOp::Insertion { .. } => n_ins += 1,
                    seqair::bam::PileupOp::Match { .. } => {}
                }
            }
            HtsColumn { pos: *col.pos(), depth: col.depth() as u32, n_del, n_refskip, n_ins }
        })
        .collect()
}

/// Compare seqair and htslib pileup column-by-column.
// r[verify pileup.htslib_compat]
// r[verify pileup_indel.depth_includes_all]
// r[verify pileup_indel.deletions_included]
// r[verify pileup_indel.refskips_included]
fn assert_pileup_parity(sam_name: &str, contig: &str, contig_len: u32) {
    let sam_path = mpileup_dir().join(format!("{sam_name}.sam"));
    assert!(sam_path.exists(), "SAM not found: {}", sam_path.display());

    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_bam(dir.path(), &sam_path);

    let hts = htslib_pileup(&bam_path, contig.as_bytes(), 0, u64::from(contig_len));
    let ours = seqair_pileup(&bam_path, contig, 0, contig_len);

    assert_eq!(
        ours.len(),
        hts.len(),
        "{sam_name}: column count seqair={} htslib={}",
        ours.len(),
        hts.len()
    );

    for (i, (s, h)) in ours.iter().zip(&hts).enumerate() {
        assert_eq!(s.pos, h.pos, "{sam_name} col {i}: pos");
        assert_eq!(
            s.depth,
            h.depth,
            "{sam_name} pos {}: depth seqair={} htslib={}",
            s.pos + 1,
            s.depth,
            h.depth
        );
        assert_eq!(
            s.n_del,
            h.n_del,
            "{sam_name} pos {}: n_del seqair={} htslib={}",
            s.pos + 1,
            s.n_del,
            h.n_del
        );
        assert_eq!(
            s.n_refskip,
            h.n_refskip,
            "{sam_name} pos {}: n_refskip seqair={} htslib={}",
            s.pos + 1,
            s.n_refskip,
            h.n_refskip
        );
        assert_eq!(
            s.n_ins,
            h.n_ins,
            "{sam_name} pos {}: n_ins seqair={} htslib={}",
            s.pos + 1,
            s.n_ins,
            h.n_ins
        );
    }
}

// --- Tests for contig "z" (LN:13) SAMs ---

/// Deletions: 3 reads with M and D CIGAR ops.
#[test]
fn mpileup_deletions() {
    assert_pileup_parity("mp_D", "z", 13);
}

/// Insertions: reads with M and I CIGAR ops.
#[test]
fn mpileup_insertions() {
    assert_pileup_parity("mp_I", "z", 13);
}

/// Deletion then insertion: 7M2D2I3M and similar patterns.
#[test]
fn mpileup_del_then_ins() {
    assert_pileup_parity("mp_DI", "z", 13);
}

/// Insertion then deletion: 7M2I2D3M and similar patterns.
#[test]
fn mpileup_ins_then_del() {
    assert_pileup_parity("mp_ID", "z", 13);
}

/// Ref-skips (N ops): RNA-seq intron-like patterns.
#[test]
fn mpileup_refskips() {
    assert_pileup_parity("mp_N", "z", 13);
}

/// Complex N+I+D+P combinations.
#[test]
fn mpileup_complex_nidp() {
    assert_pileup_parity("mp_N2", "z", 13);
}

/// Padding (P ops) with insertions.
#[test]
fn mpileup_padding() {
    assert_pileup_parity("mp_P", "z", 13);
}

// --- Overlap tests (contig "1") ---
// Note: these tests exercise overlapping mate pairs. seqair does NOT do
// mate overlap dedup (per spec: 4-pileup-dedup.md says "not implemented").
// htslib's bam_plp_auto also does NOT do overlap dedup by default — that's
// done by samtools mpileup's -x flag. So seqair and htslib should match here.

/// Overlapping mate pairs (variant 1): read2 before read1 in file.
#[test]
fn mpileup_overlap1() {
    assert_pileup_parity("mp_overlap1", "1", 100_020);
}

/// Overlapping mate pairs (variant 2): read1 before read2 in file.
#[test]
fn mpileup_overlap2() {
    assert_pileup_parity("mp_overlap2", "1", 100_020);
}
