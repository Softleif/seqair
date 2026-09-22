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

/// One alignment's contribution to a column, in the shape both implementations
/// can describe: where it sits in the read (`None` inside a deletion or skip),
/// and which of htslib's three predicates hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Entry {
    qpos: Option<usize>,
    is_del: bool,
    is_refskip: bool,
    has_ins: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct HtsColumn {
    pos: u32,
    depth: u32,
    n_del: usize,
    n_refskip: usize,
    n_ins: usize,
    /// Sorted, so the two implementations' iteration orders do not matter.
    entries: Vec<Entry>,
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
        let mut entries = Vec::new();
        for a in p.alignments() {
            if a.is_del() {
                n_del += 1;
            }
            if a.is_refskip() {
                n_refskip += 1;
            }
            let has_ins = matches!(a.indel(), Indel::Ins(_));
            if has_ins {
                n_ins += 1;
            }
            entries.push(Entry {
                qpos: a.qpos(),
                is_del: a.is_del(),
                is_refskip: a.is_refskip(),
                has_ins,
            });
        }
        entries.sort_unstable();
        columns.push(HtsColumn { pos, depth: p.depth(), n_del, n_refskip, n_ins, entries });
    }
    columns
}

/// Read pileup from seqair.
fn seqair_pileup(bam_path: &Path, contig: &str, start: u32, end: u32) -> Vec<HtsColumn> {
    let mut reader = IndexedBamReader::open(bam_path).expect("seqair open");
    let tid = reader.header().tid(contig).expect("contig not found");
    let mut store = RecordStore::new();
    let span = (Pos0::new(start).unwrap()..=Pos0::new(end).unwrap()).into();
    reader.fetch_into(tid, span, &mut store).expect("fetch");

    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span);

    let mut columns = Vec::new();

    while let Some(col) = engine.pileups() {
        let mut n_del = 0;
        let mut n_refskip = 0;
        let mut n_ins = 0;
        let mut entries = Vec::new();
        for a in col.raw_alignments() {
            let entry = match a.op {
                seqair::bam::PileupOp::Deletion { .. } => {
                    n_del += 1;
                    Entry { qpos: None, is_del: true, is_refskip: false, has_ins: false }
                }
                seqair::bam::PileupOp::ComplexIndel { is_refskip, .. } => {
                    n_del += 1;
                    n_ins += 1;
                    if is_refskip {
                        n_refskip += 1;
                    }
                    Entry { qpos: None, is_del: true, is_refskip, has_ins: true }
                }
                seqair::bam::PileupOp::RefSkip => {
                    n_del += 1;
                    n_refskip += 1;
                    Entry { qpos: None, is_del: true, is_refskip: true, has_ins: false }
                }
                seqair::bam::PileupOp::Insertion { qpos, .. } => {
                    n_ins += 1;
                    Entry {
                        qpos: Some(qpos.get() as usize),
                        is_del: false,
                        is_refskip: false,
                        has_ins: true,
                    }
                }
                seqair::bam::PileupOp::Match { qpos, .. } => Entry {
                    qpos: Some(qpos.get() as usize),
                    is_del: false,
                    is_refskip: false,
                    has_ins: false,
                },
                seqair::bam::PileupOp::SoftClip { .. } => {
                    panic!("overhang 0 must not emit soft clips")
                }
            };
            entries.push(entry);
        }
        entries.sort_unstable();
        columns.push(HtsColumn {
            pos: u32::from(col.pos()),
            depth: col.depth() as u32,
            n_del,
            n_refskip,
            n_ins,
            entries,
        });
    }

    columns
}

/// Compare two column streams entry by entry. `context` names the input in
/// every failure message — a fixture name, or the SAM text itself.
fn assert_columns_match(ours: &[HtsColumn], hts: &[HtsColumn], context: &str) {
    assert_eq!(
        ours.len(),
        hts.len(),
        "{context}: column count seqair={} htslib={}",
        ours.len(),
        hts.len()
    );

    for (i, (s, h)) in ours.iter().zip(hts).enumerate() {
        assert_eq!(s.pos, h.pos, "{context} col {i}: pos");
        let at = s.pos + 1;
        assert_eq!(s.depth, h.depth, "{context} pos {at}: depth");
        assert_eq!(s.n_del, h.n_del, "{context} pos {at}: n_del");
        assert_eq!(s.n_refskip, h.n_refskip, "{context} pos {at}: n_refskip");
        assert_eq!(s.n_ins, h.n_ins, "{context} pos {at}: n_ins");
        // r[pileup.htslib_compat] also requires the same qpos per alignment,
        // which the counts above cannot see.
        assert_eq!(s.entries, h.entries, "{context} pos {at}: per-alignment entries");
    }
}

/// Compare seqair and htslib pileup column-by-column on a fixture SAM.
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

    assert_columns_match(&ours, &hts, sam_name);
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

// ── The same parity, on CIGARs nobody wrote by hand ─────────────────────
//
// The fixtures above are htslib's own mpileup test files: eight hand-built
// SAMs, each aimed at one CIGAR combination. They are the cases htslib's
// authors thought to write down. `r[pileup.htslib_compat]` claims parity for
// *any* BAM and region, so the generator below builds SAMs from the same op
// alphabet and compares the same way — the reference implementation is the
// oracle, and nothing about the input is fixed.

mod generated {
    use super::{assert_columns_match, htslib_pileup, sam_to_bam, seqair_pileup};
    use hegel::prelude::*;
    use std::fmt::Write as _;
    use std::io::Write as _;

    const CONTIG: &str = "ref";
    const CONTIG_LEN: u32 = 200;
    /// Reads start within this many bases of the contig start.
    const START_WINDOW: u32 = 40;

    /// A CIGAR as `(len, op_char)` pairs, valid by construction: matches at
    /// both ends of the aligned body, soft clips outside those, hard clips
    /// outside the soft clips. The body's interior draws from the ops the
    /// mpileup fixtures exercise — `I`, `D`, `N` and `P`.
    #[hegel::composite]
    fn arb_cigar(tc: &TestCase) -> Vec<(u32, char)> {
        let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);

        let mut body: Vec<(u32, char)> = vec![(tc.draw_silent(n(1, 10)), 'M')];
        let groups = tc.draw_silent(gs::integers::<usize>().max_value(3));
        for _ in 0..groups {
            // One or two consecutive non-match ops, then a match: `D` followed
            // by `I` is the `mp_DI` shape, `I` then `D` is `mp_ID`.
            let inner = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(2));
            for _ in 0..inner {
                let op = tc.draw_silent(gs::sampled_from(&['I', 'D', 'N', 'P']));
                let len = match op {
                    'N' => tc.draw_silent(n(1, 8)),
                    'P' => tc.draw_silent(n(1, 2)),
                    _ => tc.draw_silent(n(1, 4)),
                };
                body.push((len, op));
            }
            body.push((tc.draw_silent(n(1, 10)), 'M'));
        }

        let mut ops = Vec::new();
        if tc.draw_silent(gs::booleans()) {
            ops.push((tc.draw_silent(n(1, 5)), 'H'));
        }
        if tc.draw_silent(gs::booleans()) {
            ops.push((tc.draw_silent(n(1, 6)), 'S'));
        }
        ops.extend(body);
        if tc.draw_silent(gs::booleans()) {
            ops.push((tc.draw_silent(n(1, 6)), 'S'));
        }
        if tc.draw_silent(gs::booleans()) {
            ops.push((tc.draw_silent(n(1, 5)), 'H'));
        }
        ops
    }

    fn query_len(ops: &[(u32, char)]) -> u32 {
        ops.iter().filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X')).map(|(l, _)| l).sum()
    }

    fn ref_span(ops: &[(u32, char)]) -> u32 {
        ops.iter().filter(|(_, op)| matches!(op, 'M' | 'D' | 'N' | '=' | 'X')).map(|(l, _)| l).sum()
    }

    /// A whole single-contig SAM file: a header and one to ten reads whose
    /// alignments overlap heavily, so most columns carry several of them.
    #[hegel::composite]
    fn arb_sam(tc: &TestCase) -> String {
        let n_reads = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(10));
        let mut sam = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n");

        for i in 0..n_reads {
            let ops = tc.draw_silent(arb_cigar());
            let span = ref_span(&ops);
            // Keep the alignment inside the contig, and start it inside a
            // window much narrower than the contig so the reads pile up on
            // each other — a run whose columns are all depth 1 exercises none
            // of the active-set bookkeeping.
            let last_start = CONTIG_LEN.saturating_sub(span).clamp(1, START_WINDOW);
            let pos = tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(last_start));
            let qlen = query_len(&ops);
            let mapq = tc.draw_silent(gs::integers::<u32>().max_value(60));
            // Reverse strand and duplicate change nothing for either pileup,
            // but they are what a real file carries.
            let flag = tc.draw_silent(gs::sampled_from(&[0u32, 16, 1024]));

            let seq: String = (0..qlen)
                .map(|_| tc.draw_silent(gs::sampled_from(&['A', 'C', 'G', 'T'])))
                .collect();
            let qual: String = (0..qlen)
                .map(|_| {
                    let q = tc.draw_silent(gs::integers::<u8>().min_value(33).max_value(73));
                    q as char
                })
                .collect();
            let mut cigar = String::new();
            for (len, op) in &ops {
                write!(cigar, "{len}{op}").expect("writing to a String cannot fail");
            }

            writeln!(sam, "r{i}\t{flag}\t{CONTIG}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
                .expect("writing to a String cannot fail");
        }
        sam
    }

    // r[verify pileup.htslib_compat]
    // r[verify pileup_indel.depth_includes_all]
    // r[verify pileup_indel.deletions_included]
    // r[verify pileup_indel.refskips_included]
    /// For a generated SAM, every column seqair produces must match htslib's
    /// `bam_plp_auto` on position, depth, the del/refskip/ins counts, and the
    /// per-alignment `qpos`.
    #[hegel::test(test_cases = 64)]
    fn pileup_matches_htslib_on_generated_reads(tc: TestCase) {
        let sam = tc.draw(arb_sam().print_as_debug());

        let dir = tempfile::tempdir().expect("tempdir");
        let sam_path = dir.path().join("generated.sam");
        std::fs::File::create(&sam_path)
            .expect("create SAM")
            .write_all(sam.as_bytes())
            .expect("write SAM");
        let bam_path = sam_to_bam(dir.path(), &sam_path);

        let hts = htslib_pileup(&bam_path, CONTIG.as_bytes(), 0, u64::from(CONTIG_LEN));
        let ours = seqair_pileup(&bam_path, CONTIG, 0, CONTIG_LEN);

        // A generated SAM that produced no columns would make the comparison
        // vacuous; every read here has at least one `M`, so it cannot happen.
        assert!(!ours.is_empty(), "no columns for:\n{sam}");
        assert_columns_match(&ours, &hts, &sam);

        // What the run actually reached, for `HEGEL_STATISTICS=1`: a suite
        // that never generates a deletion is not testing deletions.
        tc.event_value("columns", ours.len() as f64);
        tc.event_value("max depth", ours.iter().map(|c| c.depth).max().unwrap_or(0).into());
        for (label, n) in [
            ("deletion columns", ours.iter().filter(|c| c.n_del > 0).count()),
            ("refskip columns", ours.iter().filter(|c| c.n_refskip > 0).count()),
            ("insertion columns", ours.iter().filter(|c| c.n_ins > 0).count()),
        ] {
            if n > 0 {
                tc.event(label);
            }
            tc.event_value(label, n as f64);
        }
    }
}
