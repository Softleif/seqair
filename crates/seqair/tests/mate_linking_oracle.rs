//! Mate fields and mate linking, refereed by samtools.
//!
//! `mate_linking.rs` checks seqair against seqair: it builds a store by hand,
//! links it, and asserts the links match what the same idea of "mate" predicts.
//! A systematic misreading of `next_ref_id`/`next_pos`/`template_len` would be
//! reproduced identically on both sides of that check and pass it.
//!
//! htslib can referee. `samtools fixmate` recomputes every mate field and every
//! mate flag bit from the alignments themselves, so the corpus here is built
//! the other way round: the generator decides only the *geometry* (where each
//! fragment lands, which strand, which flavour of alignment) and writes SAM
//! with the mate columns blank. `fixmate` fills them in, and that output is the
//! oracle — seqair must read back exactly what htslib computed, and a BAM
//! seqair writes from those records must give a second `fixmate` nothing to do.
//!
//! TLEN is the subtle field, and `tlen_is_the_signed_distance_between_five_prime_ends`
//! pins what htslib actually does rather than what the SAM prose suggests.
//! `cram_tlen.rs` is the neighbouring test; it replays htslib's own tlen/
//! fixtures through the CRAM and BAM readers and never invokes `fixmate`, so
//! the two do not overlap.
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
    clippy::cast_sign_loss,
    reason = "test code with known small values"
)]

use hegel::prelude::*;
use seqair::bam::header::BamHeader;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Two references, so a template can straddle them. Long enough that nothing
/// the generator draws runs off the end.
const CONTIGS: [(&str, u32); 2] = [("chr1", 20_000), ("chr2", 20_000)];

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

// ── running samtools ──────────────────────────────────────────────────────

fn samtools(args: &[&str]) -> String {
    let out = Command::new("samtools").args(args).output().expect("samtools on PATH");
    assert!(
        out.status.success(),
        "samtools {args:?} failed: {}",
        String::from_utf8(out.stderr).expect("samtools reports UTF-8")
    );
    String::from_utf8(out.stdout).expect("samtools prints UTF-8")
}

fn path_str(p: &Path) -> &str {
    p.to_str().expect("temp paths are UTF-8")
}

// ── addressing a record across the pipeline ───────────────────────────────

/// Which alignment of a template this is. `fixmate` rewrites POS (it places an
/// unmapped mate on its partner) and the mate flag bits, so neither can address
/// a record; the bits below are the ones it never touches, and the generator
/// emits at most one record per `(qname, role)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Role {
    FirstPrimary,
    SecondPrimary,
    Secondary,
    Supplementary,
}

fn role_of(flags: u16) -> Role {
    if flags & SECONDARY != 0 {
        Role::Secondary
    } else if flags & SUPPLEMENTARY != 0 {
        Role::Supplementary
    } else if flags & SECOND != 0 {
        Role::SecondPrimary
    } else {
        Role::FirstPrimary
    }
}

type Key = (String, Role);

/// The columns `fixmate` recomputes, and nothing else.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Mate {
    flags: u16,
    next_ref: i32,
    next_pos: i32,
    tlen: i32,
}

fn tid_of(name: &str) -> i32 {
    CONTIGS.iter().position(|(n, _)| *n == name).map_or(-1, |i| i as i32)
}

/// The mate columns of every alignment in a SAM text, keyed by `(qname, role)`.
fn parse_sam(text: &str) -> BTreeMap<Key, Mate> {
    let mut map = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 11, "SAM line has too few fields: {line}");
        let flags: u16 = f[1].parse().expect("FLAG is a number");
        let pnext: i64 = f[7].parse().expect("PNEXT is a number");
        let next_ref = match f[6] {
            "*" => -1,
            "=" => tid_of(f[2]),
            other => tid_of(other),
        };
        let mate = Mate {
            flags,
            next_ref,
            // SAM is 1-based and spells "unavailable" as 0; BAM is 0-based and
            // spells it -1.
            next_pos: (pnext - 1) as i32,
            tlen: f[8].parse().expect("TLEN is a number"),
        };
        let key = (f[0].to_owned(), role_of(flags));
        assert!(map.insert(key.clone(), mate).is_none(), "two records share the key {key:?}");
    }
    map
}

/// The same map, read out of a store seqair filled.
fn store_mates(store: &RecordStore, into: &mut BTreeMap<Key, Mate>) {
    for idx in store.indices() {
        let rec = store.record(idx).expect("index came from the store");
        let qname =
            std::str::from_utf8(rec.qname()).expect("generated qnames are ASCII").to_owned();
        let mate = Mate {
            flags: rec.flags.raw(),
            next_ref: rec.next_ref_id,
            next_pos: rec.next_pos,
            tlen: rec.template_len,
        };
        let key = (qname, role_of(rec.flags.raw()));
        assert!(into.insert(key.clone(), mate).is_none(), "two records share the key {key:?}");
    }
}

/// Report the first field that differs, naming the record — a bare
/// `assert_eq!` on two maps prints both in full and says nothing about which
/// of the four columns moved.
#[track_caller]
fn assert_same_mates(left: &BTreeMap<Key, Mate>, right: &BTreeMap<Key, Mate>, what: &str) {
    let left_keys: BTreeSet<&Key> = left.keys().collect();
    let right_keys: BTreeSet<&Key> = right.keys().collect();
    assert_eq!(left_keys, right_keys, "{what}: different records");
    for (key, l) in left {
        let r = &right[key];
        assert_eq!(l.flags, r.flags, "{what}: FLAG of {key:?} ({l:?} vs {r:?})");
        assert_eq!(l.next_ref, r.next_ref, "{what}: RNEXT of {key:?} ({l:?} vs {r:?})");
        assert_eq!(l.next_pos, r.next_pos, "{what}: PNEXT of {key:?} ({l:?} vs {r:?})");
        assert_eq!(l.tlen, r.tlen, "{what}: TLEN of {key:?} ({l:?} vs {r:?})");
    }
}

// ── the generated corpus ──────────────────────────────────────────────────

/// How a template's two fragments are placed. Everything else about a record
/// — its mate columns and its mate flag bits — is left for `fixmate` to decide.
#[derive(Clone, Copy, Debug)]
enum Shape {
    /// Both fragments on `chr1`. `gap` places the second relative to the
    /// first and may be zero (both start at the same position) or negative.
    SameRef { gap: i32, a_rev: bool, b_rev: bool },
    /// The second fragment on `chr2`: TLEN must be 0 and RNEXT must not be `=`.
    CrossRef { a_rev: bool, b_rev: bool },
    /// The second fragment unmapped, which `fixmate` places on the first.
    MateUnmapped { a_rev: bool },
    /// Neither fragment placed. Such records are outside any region query, so
    /// they test that seqair's reader reports exactly the placed ones.
    BothUnmapped,
}

#[derive(Clone, Copy, Debug)]
struct Template {
    shape: Shape,
    a_pos: u32,
    a_len: u32,
    b_len: u32,
    /// A secondary alignment of read 1, which is never a mate.
    secondary: bool,
    /// A supplementary alignment of read 1, likewise.
    supplementary: bool,
}

#[hegel::composite]
fn arb_template(tc: &TestCase) -> Template {
    let len = |tc: &TestCase| tc.draw_silent(gs::integers::<u32>().min_value(5).max_value(60));
    let a_len = len(tc);
    let shape = match tc.draw_silent(gs::integers::<u8>().max_value(7)) {
        0..=3 => Shape::SameRef {
            // Zero and the small offsets are drawn deliberately: they are
            // where the two fragments start together or contain one another,
            // which is where TLEN's sign rule and its span rule part company.
            gap: match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
                0 => 0,
                1 => tc.draw_silent(gs::integers::<i32>().min_value(-6).max_value(6)),
                _ => tc.draw_silent(gs::integers::<i32>().min_value(-400).max_value(400)),
            },
            a_rev: tc.draw_silent(gs::booleans()),
            b_rev: tc.draw_silent(gs::booleans()),
        },
        4 | 5 => Shape::CrossRef {
            a_rev: tc.draw_silent(gs::booleans()),
            b_rev: tc.draw_silent(gs::booleans()),
        },
        6 => Shape::MateUnmapped { a_rev: tc.draw_silent(gs::booleans()) },
        _ => Shape::BothUnmapped,
    };
    Template {
        shape,
        a_pos: tc.draw_silent(gs::integers::<u32>().min_value(500).max_value(5_000)),
        a_len,
        // Equal lengths as often as not: with `gap == 0` that is the exact tie
        // — same start, same end — where the sign has to come from somewhere
        // other than the coordinates.
        b_len: if tc.draw_silent(gs::booleans()) { a_len } else { len(tc) },
        secondary: tc.draw_silent(gs::booleans()),
        supplementary: tc.draw_silent(gs::booleans()),
    }
}

fn templates() -> impl PrintableGenerator<Vec<Template>> {
    gs::vecs(arb_template().print_as_debug()).min_size(1).max_size(6)
}

/// Record which mate cases a drawn corpus reached. Run with `HEGEL_STATISTICS=1`
/// to see the distribution: a generator that quietly stopped producing
/// same-start pairs or cross-reference templates would otherwise leave the
/// suite looking just as green.
fn note_shapes(tc: &TestCase, templates: &[Template]) {
    for tpl in templates {
        match tpl.shape {
            Shape::SameRef { gap, a_rev, b_rev } => {
                tc.event("same reference");
                if gap == 0 {
                    tc.event("both fragments start at the same position");
                }
                let b_pos = i64::from(tpl.a_pos) + i64::from(gap);
                let a_end = i64::from(tpl.a_pos) + i64::from(tpl.a_len);
                if b_pos >= i64::from(tpl.a_pos) && b_pos + i64::from(tpl.b_len) <= a_end {
                    tc.event("one fragment contains the other");
                }
                match (a_rev, b_rev) {
                    (false, true) | (true, false) => tc.event("FR or RF orientation"),
                    _ => tc.event("FF or RR orientation"),
                }
            }
            Shape::CrossRef { .. } => tc.event("mates on different references"),
            Shape::MateUnmapped { .. } => tc.event("one mate unmapped"),
            Shape::BothUnmapped => tc.event("both mates unmapped"),
        }
        if tpl.secondary {
            tc.event("a secondary alignment shares the qname");
        }
        if tpl.supplementary {
            tc.event("a supplementary alignment shares the qname");
        }
    }
}

fn seq_of(len: u32) -> String {
    "ACGT".repeat(len as usize).chars().take(len as usize).collect()
}

/// One alignment line with the mate columns blank (`*`, `0`, `0`).
fn sam_line(
    out: &mut String,
    qname: &str,
    flags: u16,
    rname: &str,
    pos1: u32,
    cigar: &str,
    len: u32,
) {
    let mapq = if flags & UNMAPPED == 0 { 60 } else { 0 };
    writeln!(
        out,
        "{qname}\t{flags}\t{rname}\t{pos1}\t{mapq}\t{cigar}\t*\t0\t0\t{}\t{}",
        seq_of(len),
        "I".repeat(len as usize)
    )
    .expect("writing to a String cannot fail");
}

fn corpus_sam(templates: &[Template]) -> String {
    let mut out = String::from("@HD\tVN:1.6\tSO:queryname\n");
    for (name, len) in CONTIGS {
        writeln!(out, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    for (t, tpl) in templates.iter().enumerate() {
        let q = format!("t{t}");
        let a_cigar = format!("{}M", tpl.a_len);
        let b_cigar = format!("{}M", tpl.b_len);
        let (a_mapped, a_rev) = match tpl.shape {
            Shape::SameRef { a_rev, .. } | Shape::CrossRef { a_rev, .. } => (true, a_rev),
            Shape::MateUnmapped { a_rev } => (true, a_rev),
            Shape::BothUnmapped => (false, false),
        };
        let a_flags = PAIRED | FIRST | if a_rev { REVERSE } else { 0 };

        match tpl.shape {
            Shape::SameRef { gap, b_rev, .. } => {
                let b_pos = (tpl.a_pos as i32 + gap) as u32;
                sam_line(&mut out, &q, a_flags, "chr1", tpl.a_pos + 1, &a_cigar, tpl.a_len);
                sam_line(
                    &mut out,
                    &q,
                    PAIRED | SECOND | if b_rev { REVERSE } else { 0 },
                    "chr1",
                    b_pos + 1,
                    &b_cigar,
                    tpl.b_len,
                );
            }
            Shape::CrossRef { b_rev, .. } => {
                sam_line(&mut out, &q, a_flags, "chr1", tpl.a_pos + 1, &a_cigar, tpl.a_len);
                sam_line(
                    &mut out,
                    &q,
                    PAIRED | SECOND | if b_rev { REVERSE } else { 0 },
                    "chr2",
                    tpl.a_pos + 1,
                    &b_cigar,
                    tpl.b_len,
                );
            }
            Shape::MateUnmapped { .. } => {
                sam_line(&mut out, &q, a_flags, "chr1", tpl.a_pos + 1, &a_cigar, tpl.a_len);
                sam_line(&mut out, &q, PAIRED | SECOND | UNMAPPED, "*", 0, "*", tpl.b_len);
            }
            Shape::BothUnmapped => {
                sam_line(&mut out, &q, PAIRED | FIRST | UNMAPPED, "*", 0, "*", tpl.a_len);
                sam_line(&mut out, &q, PAIRED | SECOND | UNMAPPED, "*", 0, "*", tpl.b_len);
            }
        }

        if a_mapped && tpl.secondary {
            sam_line(
                &mut out,
                &q,
                a_flags | SECONDARY,
                "chr1",
                tpl.a_pos + 18,
                &a_cigar,
                tpl.a_len,
            );
        }
        if a_mapped && tpl.supplementary {
            // Hard-clipped, which is what an aligner emits for a split read.
            sam_line(
                &mut out,
                &q,
                a_flags | SUPPLEMENTARY,
                "chr1",
                tpl.a_pos + 32,
                &format!("3H{}M", tpl.a_len),
                tpl.a_len,
            );
        }
    }
    out
}

/// The corpus after `samtools fixmate`: the mate fields htslib computed, plus
/// the same records sorted and indexed so seqair can read them back.
struct Oracle {
    _dir: tempfile::TempDir,
    dir: PathBuf,
    sorted_bam: PathBuf,
    /// Every alignment, including the ones no region query can reach.
    all: BTreeMap<Key, Mate>,
}

impl Oracle {
    fn build(templates: &[Template]) -> Self {
        let dir = tempfile::tempdir().expect("tempdir");
        let root = dir.path().to_owned();
        let raw = root.join("raw.sam");
        let fixed = root.join("fixed.sam");
        let sorted = root.join("sorted.bam");
        std::fs::write(&raw, corpus_sam(templates)).expect("write corpus");

        samtools(&["fixmate", "-O", "sam", path_str(&raw), path_str(&fixed)]);
        let all = parse_sam(&std::fs::read_to_string(&fixed).expect("read fixmate output"));
        samtools(&["sort", "-o", path_str(&sorted), path_str(&fixed)]);
        samtools(&["index", path_str(&sorted)]);

        Self { _dir: dir, dir: root, sorted_bam: sorted, all }
    }

    /// The alignments a region query can reach: everything `fixmate` left with
    /// a reference. An unplaced record is not in the index at all.
    fn placed(&self) -> BTreeMap<Key, Mate> {
        let text = samtools(&["view", path_str(&self.sorted_bam), "chr1", "chr2"]);
        parse_sam(&text)
    }
}

/// Fetch every placed record, one store per reference (`fetch_into` clears).
fn fetch_all(path: &Path) -> Vec<RecordStore> {
    let mut reader = IndexedBamReader::open(path).expect("open the BAM seqair must read");
    CONTIGS
        .iter()
        .map(|(name, len)| {
            let tid = reader.header().tid(name).expect("contig in header");
            let mut store = RecordStore::new();
            reader
                .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(*len).unwrap(), &mut store)
                .expect("fetch");
            store
        })
        .collect()
}

fn header() -> BamHeader {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    BamHeader::from_sam_text(&text).expect("two @SQ lines are a header")
}

/// Write every fetched record back out with seqair's own writer, in the
/// coordinate order the stores are already in.
fn write_back(dir: &Path, stores: &[RecordStore]) -> PathBuf {
    let path = dir.join("seqair.bam");
    let header = header();
    let mut writer =
        BamWriterBuilder::to_path(&path, &header).build().expect("build the BAM writer");
    for store in stores {
        for idx in store.indices() {
            writer.write_store_record(store, idx).expect("write a record seqair just read");
        }
    }
    writer.finish().expect("finish");
    path
}

// ── the TLEN convention ───────────────────────────────────────────────────

/// TLEN as htslib computes it: the signed distance from this fragment's 5' end
/// to its mate's. The 5' end of a forward fragment is its leftmost base; of a
/// reverse fragment, one past its rightmost. Note this is *not* the
/// leftmost-to-rightmost span the SAM prose describes — the two agree on a
/// canonical FR pair and disagree as soon as the orientation is FF/RR or one
/// fragment contains the other.
fn five_prime(pos: i64, span: i64, reverse: bool) -> i64 {
    if reverse { pos + span } else { pos }
}

// r[verify bam.record.mate_fields]
/// What `samtools fixmate` writes into TLEN, pinned case by case — including
/// the tie the generator cannot decide on its own: both fragments starting at
/// the same position. htslib breaks it by strand, not by order in the file and
/// not by the READ1 flag: the forward fragment takes the positive sign.
#[test]
fn tlen_is_the_signed_distance_between_five_prime_ends() {
    // (qname, a_pos, a_len, a_rev, b spec, a's expected TLEN)
    struct Case {
        name: &'static str,
        a: (u32, u32, bool),
        b: (&'static str, u32, u32, bool),
        a_tlen: i32,
    }
    let cases = [
        Case { name: "fr_apart", a: (100, 10, false), b: ("chr1", 190, 10, true), a_tlen: 100 },
        // Both start at 200. The forward fragment is the positive one.
        Case { name: "fr_same_start", a: (200, 10, false), b: ("chr1", 200, 10, true), a_tlen: 10 },
        Case {
            name: "rf_same_start",
            a: (250, 10, true),
            b: ("chr1", 250, 10, false),
            a_tlen: -10,
        },
        // Same 5' end on both sides, so no distance at all.
        Case { name: "ff_same_start", a: (300, 10, false), b: ("chr1", 300, 10, false), a_tlen: 0 },
        Case { name: "rr_same_start", a: (400, 10, true), b: ("chr1", 400, 10, true), a_tlen: 0 },
        // 90, not the 100 that leftmost-to-rightmost would give.
        Case { name: "ff_apart", a: (500, 10, false), b: ("chr1", 590, 10, false), a_tlen: 90 },
        // The first fragment contains the second: 9, not the span of 20.
        Case { name: "contained", a: (600, 20, false), b: ("chr1", 604, 5, true), a_tlen: 9 },
        // A template that spans two references has no length.
        Case { name: "cross_ref", a: (700, 10, false), b: ("chr2", 800, 10, true), a_tlen: 0 },
    ];

    let mut sam = String::from("@HD\tVN:1.6\tSO:queryname\n");
    for (name, len) in CONTIGS {
        writeln!(sam, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    for case in &cases {
        let (a_pos, a_len, a_rev) = case.a;
        let (b_ref, b_pos, b_len, b_rev) = case.b;
        sam_line(
            &mut sam,
            case.name,
            PAIRED | FIRST | if a_rev { REVERSE } else { 0 },
            "chr1",
            a_pos + 1,
            &format!("{a_len}M"),
            a_len,
        );
        sam_line(
            &mut sam,
            case.name,
            PAIRED | SECOND | if b_rev { REVERSE } else { 0 },
            b_ref,
            b_pos + 1,
            &format!("{b_len}M"),
            b_len,
        );
    }
    // And one pair whose mate is unmapped, which has no length either.
    sam_line(&mut sam, "mate_unmapped", PAIRED | FIRST, "chr1", 901, "10M", 10);
    sam_line(&mut sam, "mate_unmapped", PAIRED | SECOND | UNMAPPED, "*", 0, "*", 10);

    let dir = tempfile::tempdir().expect("tempdir");
    let raw = dir.path().join("tlen.sam");
    let fixed = dir.path().join("tlen.fixed.sam");
    std::fs::write(&raw, sam).expect("write");
    samtools(&["fixmate", "-O", "sam", path_str(&raw), path_str(&fixed)]);
    let fixed = parse_sam(&std::fs::read_to_string(&fixed).expect("read fixmate output"));

    for case in &cases {
        let first = fixed[&(case.name.to_owned(), Role::FirstPrimary)];
        let second = fixed[&(case.name.to_owned(), Role::SecondPrimary)];
        assert_eq!(first.tlen, case.a_tlen, "{}: read 1 TLEN", case.name);
        assert_eq!(second.tlen, -case.a_tlen, "{}: read 2 TLEN is the negation", case.name);

        // The same numbers from the formula, so the rule is pinned as a rule
        // and not only as eight literals.
        let (a_pos, a_len, a_rev) = case.a;
        let (b_ref, b_pos, b_len, b_rev) = case.b;
        let expected = if b_ref == "chr1" {
            five_prime(i64::from(b_pos), i64::from(b_len), b_rev)
                - five_prime(i64::from(a_pos), i64::from(a_len), a_rev)
        } else {
            0
        };
        assert_eq!(
            i64::from(first.tlen),
            expected,
            "{}: the five-prime formula disagrees with fixmate",
            case.name
        );
    }

    for role in [Role::FirstPrimary, Role::SecondPrimary] {
        assert_eq!(fixed[&("mate_unmapped".to_owned(), role)].tlen, 0, "unmapped mate, {role:?}");
    }
}

// ── the fixmate oracle ────────────────────────────────────────────────────

// r[verify bam.record.mate_fields]
// r[verify record_store.slim_record_fields]
/// seqair must read back exactly the mate fields htslib computed. Since the
/// generator never chose them, agreement here is agreement with htslib and not
/// with a second copy of seqair's own arithmetic.
#[hegel::test(test_cases = 24)]
fn seqair_reads_the_mate_fields_fixmate_computed(tc: TestCase) {
    let templates = tc.draw(templates());
    note_shapes(&tc, &templates);
    let oracle = Oracle::build(&templates);

    let mut seqair = BTreeMap::new();
    for store in fetch_all(&oracle.sorted_bam) {
        store_mates(&store, &mut seqair);
    }

    assert_same_mates(&oracle.placed(), &seqair, "fixmate vs seqair");

    // The records seqair did not report are exactly the unplaced ones: a
    // region query cannot reach them, and nothing else may go missing.
    let unreachable: BTreeSet<&Key> =
        oracle.all.keys().filter(|k| !seqair.contains_key(*k)).collect();
    for key in unreachable {
        assert_eq!(
            oracle.all[key].flags & UNMAPPED,
            UNMAPPED,
            "{key:?} is placed but seqair did not report it"
        );
        assert_eq!(oracle.all[key].next_ref, -1, "{key:?} was placed by fixmate");
    }
}

// r[verify bam_writer.mate_fields_verbatim]
/// The round trip in the other direction: seqair writes the records it read,
/// and `samtools fixmate` recomputes every mate field over that file. Anything
/// it changes is a disagreement between seqair and htslib about mates — either
/// seqair read a field wrongly or wrote it wrongly, and either way the file it
/// produced is not what htslib would have produced.
#[hegel::test(test_cases = 24)]
fn fixmate_finds_nothing_to_fix_in_a_bam_seqair_wrote(tc: TestCase) {
    let templates = tc.draw(templates());
    let oracle = Oracle::build(&templates);
    let stores = fetch_all(&oracle.sorted_bam);
    let written = write_back(&oracle.dir, &stores);

    // What seqair put in the file, read back by samtools rather than by
    // seqair, so the writer is not checked against the reader that fed it.
    let as_written = parse_sam(&samtools(&["view", path_str(&written)]));
    assert_same_mates(&oracle.placed(), &as_written, "fixmate vs the BAM seqair wrote");

    // `fixmate` needs the two halves of a template adjacent.
    let collated = oracle.dir.join("collated.bam");
    let refixed_path = oracle.dir.join("refixed.sam");
    samtools(&["sort", "-n", "-o", path_str(&collated), path_str(&written)]);
    samtools(&["fixmate", "-O", "sam", path_str(&collated), path_str(&refixed_path)]);
    let refixed = parse_sam(&std::fs::read_to_string(&refixed_path).expect("read fixmate output"));

    assert_same_mates(&as_written, &refixed, "seqair's BAM vs the same BAM re-fixmated");
}

// r[verify record_store.link_mates.htslib_agreement]
// r[verify record_store.link_mates+2]
/// The templates `link_mates` pairs must be the templates htslib calls pairs.
/// The expected set is derived from `samtools view`'s own output — flags and
/// RNEXT as samtools prints them — so a misreading of the mate fields inside
/// seqair cannot move both sides of the comparison at once.
#[hegel::test(test_cases = 24)]
fn link_mates_pairs_exactly_the_templates_samtools_calls_pairs(tc: TestCase) {
    let templates = tc.draw(templates());
    let oracle = Oracle::build(&templates);

    // Paired, mapped, mate mapped, primary — and the mate on this very
    // reference, which samtools spells `=`.
    let selected = samtools(&[
        "view",
        "-f",
        "0x1",
        "-F",
        "0x90c",
        path_str(&oracle.sorted_bam),
        "chr1",
        "chr2",
    ]);
    let mut same_ref_halves: BTreeMap<String, usize> = BTreeMap::new();
    for line in selected.lines().filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        if f[6] == "=" {
            *same_ref_halves.entry(f[0].to_owned()).or_default() += 1;
        }
    }
    let expected: BTreeSet<&String> =
        same_ref_halves.iter().filter(|(_, n)| **n == 2).map(|(q, _)| q).collect();

    let mut linked: BTreeSet<String> = BTreeSet::new();
    let mut pairs = 0u32;
    for mut store in fetch_all(&oracle.sorted_bam) {
        let stats = store.link_mates();
        assert!(stats.is_clean(), "generated qnames are unique per template");
        pairs += stats.pairs;
        for idx in store.indices() {
            let rec = store.record(idx).expect("index came from the store");
            if rec.mate_idx().is_some() {
                linked.insert(std::str::from_utf8(rec.qname()).expect("ASCII qnames").to_owned());
            }
        }
    }

    let linked_refs: BTreeSet<&String> = linked.iter().collect();
    assert_eq!(linked_refs, expected, "linked templates");
    assert_eq!(pairs as usize, expected.len(), "one pair per linked template");
}
