//! Mate linking in the record store and its surfacing in pileup columns.
//!
//! Covers `record_store.qname_hash`, `record_store.link_mates`,
//! `record_store.mate_overlap`, `record_store.link_mates.invalidated`,
//! `pileup.mate_link_cache`, `pileup.column_record_order` and
//! `pileup.column_find_record`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_lossless,
    reason = "test code with known small values"
)]

use proptest::prelude::*;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::{RecordStore, qname_hash};
use seqair_types::{BamFlags, Base, Pos0};

const PAIRED: u16 = 0x1;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// A read to push: the fields mate linking actually looks at.
#[derive(Clone, Copy, Debug)]
struct Read {
    pos: i32,
    len: u32,
    mate_pos: i32,
    flags: u16,
    tid: i32,
    mate_tid: i32,
}

impl Read {
    /// One half of a normal proper pair.
    fn mate(pos: i32, len: u32, mate_pos: i32, which: u16) -> Self {
        Self { pos, len, mate_pos, flags: PAIRED | which, tid: 0, mate_tid: 0 }
    }

    fn end(self) -> i32 {
        self.pos + self.len as i32
    }

    fn with_flags(mut self, flags: u16) -> Self {
        self.flags = flags;
        self
    }
}

fn push(store: &mut RecordStore, qname: &[u8], read: Read) -> u32 {
    let len = read.len as usize;
    let cigar = [CigarOp::new(CigarOpType::Match, read.len)];
    let bases = vec![Base::A; len];
    let qual = vec![30u8; len];
    store
        .push_fields(
            Pos0::new(read.pos as u32).unwrap(),
            Pos0::new(read.end() as u32).unwrap(),
            BamFlags::from(read.flags),
            60,
            read.len,
            0,
            qname,
            &cigar,
            &bases,
            &qual,
            &[],
            read.tid,
            read.mate_tid,
            read.mate_pos,
            0,
            &mut (),
        )
        .expect("push_fields accepts a well-formed record")
        .expect("no customizer rejects records")
}

/// Push both halves of a proper pair and link.
fn linked_pair(a: Read, b: Read) -> RecordStore {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", a);
    push(&mut store, b"frag", b);
    store.link_mates();
    store
}

fn overlap_of(store: &RecordStore, idx: u32) -> Option<(u64, u64)> {
    store.mate_overlap(idx).map(|r| (r.start.as_u64(), r.end.as_u64()))
}

// ── qname hash ────────────────────────────────────────────────────────────

// r[verify record_store.qname_hash]
#[test]
fn qname_hash_is_a_pure_function_of_the_bytes() {
    assert_eq!(
        qname_hash(b"A00123:45:HVWXY:1:1101:12345:6789"),
        qname_hash(b"A00123:45:HVWXY:1:1101:12345:6789")
    );
    assert_ne!(qname_hash(b"read1"), qname_hash(b"read2"));
    assert_ne!(qname_hash(b""), qname_hash(b"\0"));
}

// r[verify record_store.qname_hash]
/// Both mates share a qname, so they must share the hash — that is what makes
/// the hash a fragment id — and the store must expose it per record.
#[test]
fn mates_share_the_qname_hash() {
    let store = linked_pair(Read::mate(100, 50, 130, FIRST), Read::mate(130, 50, 100, SECOND));
    assert_eq!(store.record(0).qname_hash, store.record(1).qname_hash);
    assert_eq!(store.record(0).qname_hash, qname_hash(b"frag"));
}

// r[verify record_store.qname_hash]
/// Avalanche: flipping one bit of a realistic qname must change about half the
/// output bits. `FxHash` fails this badly on the low bits; the folded multiply
/// does not. The window is generous (a fair coin over 64 bits has sd ≈ 4) but
/// it is nowhere near what a weak-avalanche hash produces.
#[test]
fn qname_hash_avalanches() {
    let name = b"A00123:45:HVWXYDRXX:1:1101:12345:6789";
    let base = qname_hash(name);
    let mut total = 0u32;
    let mut trials = 0u32;
    let mut worst = u32::MAX;
    for byte in 0..name.len() {
        for bit in 0..8u32 {
            let mut flipped = name.to_vec();
            flipped[byte] ^= 1 << bit;
            let changed = (qname_hash(&flipped) ^ base).count_ones();
            total += changed;
            worst = worst.min(changed);
            trials += 1;
        }
    }
    let mean = f64::from(total) / f64::from(trials);
    assert!((24.0..40.0).contains(&mean), "mean flipped bits {mean} not near 32");
    assert!(worst >= 8, "a single input bit moved only {worst} output bits");
}

// ── linking ───────────────────────────────────────────────────────────────

// r[verify record_store.link_mates]
// r[verify record_store.mate_overlap]
#[test]
fn overlapping_mates_link_and_share_the_overlap() {
    let store = linked_pair(Read::mate(100, 50, 130, FIRST), Read::mate(130, 50, 100, SECOND));

    assert_eq!(store.record(0).mate_idx(), Some(1));
    assert_eq!(store.record(1).mate_idx(), Some(0));
    assert_eq!(overlap_of(&store, 0), Some((130, 150)));
    assert_eq!(overlap_of(&store, 1), Some((130, 150)));
}

// r[verify record_store.mate_overlap]
#[test]
fn disjoint_mates_link_with_an_empty_overlap() {
    let store = linked_pair(Read::mate(100, 50, 300, FIRST), Read::mate(300, 50, 100, SECOND));

    assert_eq!(store.record(0).mate_idx(), Some(1));
    let overlap = store.mate_overlap(0).expect("linked");
    assert!(overlap.start >= overlap.end, "disjoint mates must yield an empty range");
}

// r[verify record_store.mate_overlap]
#[test]
fn an_unlinked_record_has_no_overlap() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    store.link_mates();
    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.mate_overlap(0), None);
}

// r[verify record_store.link_mates]
/// A supplementary alignment shares the qname but is not a mate: it must not
/// link, and it must not steal the link from the real mate either.
#[test]
fn supplementary_alignments_do_not_link() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    push(
        &mut store,
        b"frag",
        Read::mate(100, 50, 130, FIRST).with_flags(PAIRED | FIRST | SUPPLEMENTARY),
    );
    store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(1));
    assert_eq!(store.record(1).mate_idx(), Some(0));
    assert_eq!(store.record(2).mate_idx(), None);
}

// r[verify record_store.link_mates]
#[test]
fn secondary_alignments_do_not_link() {
    let mut store = RecordStore::new();
    push(
        &mut store,
        b"frag",
        Read::mate(100, 50, 130, FIRST).with_flags(PAIRED | FIRST | SECONDARY),
    );
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    store.link_mates();

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates]
/// The mate lies outside the fetched region, so only one half is in the store.
#[test]
fn a_mate_outside_the_store_leaves_the_read_unlinked() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 9_000, FIRST));
    push(&mut store, b"other", Read::mate(120, 50, 8_000, FIRST));
    store.link_mates();

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates]
#[test]
fn mates_on_another_contig_do_not_link() {
    let a = Read { mate_tid: 1, ..Read::mate(100, 50, 130, FIRST) };
    let b = Read { tid: 1, ..Read::mate(130, 50, 100, SECOND) };
    let store = linked_pair(a, b);

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates]
#[test]
fn unpaired_reads_do_not_link() {
    let a = Read::mate(100, 50, 130, 0).with_flags(0);
    let b = Read::mate(130, 50, 100, 0).with_flags(0);
    let store = linked_pair(a, b);

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates]
/// Two templates whose reads sit at the same coordinates must not cross-link.
#[test]
fn distinct_qnames_at_identical_positions_do_not_cross_link() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag_a", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_b", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_a", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag_b", Read::mate(130, 50, 100, SECOND));
    store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(2));
    assert_eq!(store.record(2).mate_idx(), Some(0));
    assert_eq!(store.record(1).mate_idx(), Some(3));
    assert_eq!(store.record(3).mate_idx(), Some(1));
}

// r[verify record_store.link_mates.invalidated]
#[test]
fn sorting_and_dedup_drop_the_links() {
    let mut store = linked_pair(Read::mate(100, 50, 130, FIRST), Read::mate(130, 50, 100, SECOND));
    assert!(store.record(0).mate_idx().is_some());

    store.sort_by_pos();
    assert_eq!(store.record(0).mate_idx(), None, "sort_by_pos must invalidate mate indices");

    store.link_mates();
    assert!(store.record(0).mate_idx().is_some());
    store.dedup();
    assert_eq!(store.record(0).mate_idx(), None, "dedup must invalidate mate indices");
}

// ── pileup surfacing ──────────────────────────────────────────────────────

// r[verify pileup.mate_link_cache]
// r[verify pileup.column_find_record]
#[test]
fn columns_report_the_overlap_and_can_reach_the_mate() {
    let store = linked_pair(Read::mate(100, 50, 130, FIRST), Read::mate(130, 50, 100, SECOND));
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(200).unwrap());

    let mut checked = 0;
    while let Some(col) = engine.pileups() {
        let pos = col.pos().as_u64();
        for aln in col.raw_alignments() {
            assert_eq!(aln.mate_idx(), Some(1 - aln.record_idx()));
            assert_eq!(
                aln.in_mate_overlap(),
                (130..150).contains(&pos),
                "in_mate_overlap wrong at {pos}"
            );
            checked += 1;
        }
        if (130..150).contains(&pos) {
            let first = col.find_record(0).expect("record 0 covers the overlap");
            let second = col.find_record(1).expect("record 1 covers the overlap");
            assert_eq!(first.record_idx(), 0);
            assert_eq!(second.record_idx(), 1);
            assert_eq!(first.qname_hash(), second.qname_hash());
        } else {
            // Exactly one of the two records is in the column outside the overlap.
            assert_eq!(col.depth(), 1);
        }
    }
    assert!(checked > 0, "engine produced no columns");
}

// r[verify pileup.column_find_record]
#[test]
fn find_record_returns_none_for_a_record_outside_the_column() {
    let store = linked_pair(Read::mate(100, 50, 300, FIRST), Read::mate(300, 50, 100, SECOND));
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(120).unwrap());
    let col = engine.pileups().expect("a column at 100");
    assert!(col.find_record(0).is_some());
    assert!(col.find_record(1).is_none(), "the disjoint mate is not in this column");
    assert!(col.find_record(99).is_none());
}

// r[verify pileup.column_record_order]
#[test]
fn column_alignments_are_ordered_by_record_idx() {
    let mut store = RecordStore::new();
    for (i, pos) in [100, 101, 101, 105, 110].into_iter().enumerate() {
        let name = format!("frag{i}");
        push(&mut store, name.as_bytes(), Read::mate(pos, 50, pos + 20, FIRST));
    }
    store.link_mates();
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(160).unwrap());

    while let Some(col) = engine.pileups() {
        let idxs: Vec<u32> = col.raw_alignments().map(|a| a.record_idx()).collect();
        assert!(
            idxs.windows(2).all(|w| w[0] < w[1]),
            "column not ascending by record_idx: {idxs:?}"
        );
    }
}

// ── proptest ──────────────────────────────────────────────────────────────

/// One generated template: two mates, optionally with a supplementary copy of
/// the first mate, optionally with the second mate missing from the store.
#[derive(Clone, Copy, Debug)]
struct Template {
    a_pos: i32,
    a_len: u32,
    b_pos: i32,
    b_len: u32,
    supplementary: bool,
    mate_present: bool,
    same_contig: bool,
}

fn templates() -> impl Strategy<Value = Vec<Template>> {
    prop::collection::vec(
        (
            100i32..400,
            20u32..80,
            100i32..400,
            20u32..80,
            any::<bool>(),
            any::<bool>(),
            any::<bool>(),
        )
            .prop_map(
                |(a_pos, a_len, b_pos, b_len, supplementary, mate_present, same_contig)| Template {
                    a_pos,
                    a_len,
                    b_pos,
                    b_len,
                    supplementary,
                    mate_present,
                    same_contig,
                },
            ),
        0..8,
    )
}

proptest! {
    // r[verify record_store.link_mates]
    // r[verify record_store.mate_overlap]
    #[test]
    fn linking_is_symmetric_and_exact(templates in templates()) {
        let mut store = RecordStore::new();
        // (record idx, template idx, is the mate half, expected partner idx)
        let mut expected: Vec<(u32, Option<u32>)> = Vec::new();

        for (t, tpl) in templates.iter().enumerate() {
            let qname = format!("frag{t}");
            let mate_tid = if tpl.same_contig { 0 } else { 1 };
            let a = Read {
                mate_tid,
                ..Read::mate(tpl.a_pos, tpl.a_len, tpl.b_pos, FIRST)
            };
            let a_idx = push(&mut store, qname.as_bytes(), a);

            let b_idx = tpl.mate_present.then(|| {
                let b = Read {
                    tid: if tpl.same_contig { 0 } else { 1 },
                    ..Read::mate(tpl.b_pos, tpl.b_len, tpl.a_pos, SECOND)
                };
                push(&mut store, qname.as_bytes(), b)
            });

            if tpl.supplementary {
                // A supplementary copy at a position that is nobody's next_pos.
                let s = Read::mate(tpl.a_pos, tpl.a_len, tpl.b_pos, FIRST)
                    .with_flags(PAIRED | FIRST | SUPPLEMENTARY);
                let s_idx = push(&mut store, qname.as_bytes(), s);
                expected.push((s_idx, None));
            }

            // A pair links only when both halves are present on the same contig.
            match b_idx {
                Some(b_idx) if tpl.same_contig => {
                    expected.push((a_idx, Some(b_idx)));
                    expected.push((b_idx, Some(a_idx)));
                }
                Some(b_idx) => {
                    expected.push((a_idx, None));
                    expected.push((b_idx, None));
                }
                None => expected.push((a_idx, None)),
            }
        }

        store.link_mates();

        for (idx, partner) in expected {
            prop_assert_eq!(
                store.record(idx).mate_idx(),
                partner,
                "record {} linked to {:?}, expected {:?}",
                idx,
                store.record(idx).mate_idx(),
                partner
            );

            // Symmetry and the interval, computed independently from the two records.
            if let Some(partner) = partner {
                prop_assert_eq!(store.record(partner).mate_idx(), Some(idx));
                let a = store.record(idx);
                let b = store.record(partner);
                let start = a.pos.max(b.pos);
                let end = a.end_pos.min(b.end_pos).max(start);
                prop_assert_eq!(store.mate_overlap(idx), Some(start..end));
                prop_assert_eq!(store.mate_overlap(partner), Some(start..end));
            } else {
                prop_assert_eq!(store.mate_overlap(idx), None);
            }
        }
    }
}
