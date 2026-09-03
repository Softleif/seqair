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
use seqair::reader::{DepthLimit, Readers, SegmentOptions};
use seqair_types::{BamFlags, Base, Pos0};
use std::collections::HashSet;
use std::num::NonZeroU32;
use std::path::Path;

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

    /// Last covered position: `end_pos` is inclusive.
    fn end(self) -> i32 {
        self.pos + self.len as i32 - 1
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
    let _stats = store.link_mates();
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
    assert_eq!(store.record(0).qname_hash(), store.record(1).qname_hash());
    assert_eq!(store.record(0).qname_hash(), Some(qname_hash(b"frag")));
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

// r[verify record_store.qname_hash]
/// Avalanche over bit flips is not the whole story: what a hash table needs is
/// that *realistic* keys spread. Illumina qnames share a long prefix and differ
/// only in tile/x/y near the end, which is the shape that punishes a weak hash.
///
/// Both ends of the value are checked, because both are load-bearing in
/// `hashbrown`: the low bits pick the bucket, the top 7 bits are the SIMD
/// control byte compared before any key is touched.
#[test]
fn qname_hash_spreads_realistic_illumina_names() {
    const BUCKETS: usize = 4096;
    let mut hashes = Vec::new();
    for lane in 1..=4u32 {
        for tile in 1101..1126u32 {
            for x in 0..40u32 {
                for y in 0..5u32 {
                    let name = format!(
                        "A00123:45:HVWXYDRXX:{lane}:{tile}:{}:{}",
                        1000 + x * 137,
                        6000 + y * 11
                    );
                    hashes.push(qname_hash(name.as_bytes()));
                }
            }
        }
    }
    assert_eq!(hashes.len(), 20_000);

    let unique: HashSet<u64> = hashes.iter().copied().collect();
    assert_eq!(unique.len(), hashes.len(), "realistic qnames must not collide at 64 bits");

    // Bucket index: the low bits.
    let mut low = vec![0u32; BUCKETS];
    // Control byte: the top 7 bits.
    let mut top = vec![0u32; 128];
    for &h in &hashes {
        low[(h as usize) % BUCKETS] += 1;
        top[(h >> 57) as usize] += 1;
    }

    // 20k names over 4096 buckets averages 4.9 per bucket. For a hash that
    // behaves like a fair coin the Poisson tail puts the fullest bucket around
    // 15 and leaves about e^-4.9 * 4096 ~= 31 buckets empty. Both bounds are
    // generous multiples of that; what they rule out is clustering, where a
    // hash whose low bits barely move on this input piles hundreds into one
    // bucket and empties most of the table.
    let worst_low = low.iter().copied().max().unwrap_or(0);
    let empty_low = low.iter().filter(|&&n| n == 0).count();
    assert!(worst_low <= 20, "low bits cluster: worst bucket holds {worst_low} of 20000");
    assert!(empty_low <= 80, "low bits cluster: {empty_low} of {BUCKETS} buckets unused");

    // 20k over 128 control-byte values averages 156.
    let (worst_top, best_top) =
        (top.iter().copied().max().unwrap_or(0), top.iter().copied().min().unwrap_or(0));
    assert!(
        worst_top < 250 && best_top > 80,
        "control byte clusters: {best_top}..={worst_top} per value, expected ~156"
    );
}

// ── linking ───────────────────────────────────────────────────────────────

// r[verify record_store.link_mates+2]
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
    let _stats = store.link_mates();
    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.mate_overlap(0), None);
}

// r[verify record_store.link_mates+2]
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
    let _stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(1));
    assert_eq!(store.record(1).mate_idx(), Some(0));
    assert_eq!(store.record(2).mate_idx(), None);
}

// r[verify record_store.link_mates+2]
#[test]
fn secondary_alignments_do_not_link() {
    let mut store = RecordStore::new();
    push(
        &mut store,
        b"frag",
        Read::mate(100, 50, 130, FIRST).with_flags(PAIRED | FIRST | SECONDARY),
    );
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    let _stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates+2]
/// The mate lies outside the fetched region, so only one half is in the store.
#[test]
fn a_mate_outside_the_store_leaves_the_read_unlinked() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 9_000, FIRST));
    push(&mut store, b"other", Read::mate(120, 50, 8_000, FIRST));
    let _stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates+2]
#[test]
fn mates_on_another_contig_do_not_link() {
    let a = Read { mate_tid: 1, ..Read::mate(100, 50, 130, FIRST) };
    let b = Read { tid: 1, ..Read::mate(130, 50, 100, SECOND) };
    let store = linked_pair(a, b);

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates+2]
#[test]
fn unpaired_reads_do_not_link() {
    let a = Read::mate(100, 50, 130, 0).with_flags(0);
    let b = Read::mate(130, 50, 100, 0).with_flags(0);
    let store = linked_pair(a, b);

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
}

// r[verify record_store.link_mates+2]
/// Two templates whose reads sit at the same coordinates must not cross-link.
#[test]
fn distinct_qnames_at_identical_positions_do_not_cross_link() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag_a", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_b", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_a", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag_b", Read::mate(130, 50, 100, SECOND));
    let _stats = store.link_mates();

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

    let _stats = store.link_mates();
    assert!(store.record(0).mate_idx().is_some());
    store.dedup();
    assert_eq!(store.record(0).mate_idx(), None, "dedup must invalidate mate indices");
}

// r[verify record_store.qname_hash.no_name]
// r[verify record_store.link_mates+2]
/// A CRAM written with `RN=false` stores no read names. Such records have no
/// template identity and MUST NOT link — otherwise every nameless record in a
/// slice would compare equal to every other and unrelated reads at reciprocal
/// positions would be treated as one molecule.
#[test]
fn nameless_records_have_no_identity_and_never_link() {
    let mut store = RecordStore::new();
    push(&mut store, b"", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"", Read::mate(130, 50, 100, SECOND));
    let stats = store.link_mates();

    assert_eq!(store.record(0).qname_hash(), None, "an empty qname is no identity");
    assert_eq!(store.record(1).qname_hash(), None);
    assert_eq!(store.record(0).mate_idx(), None, "nameless records must not link");
    assert_eq!(store.record(1).mate_idx(), None);
    assert_eq!(stats.pairs, 0);
    assert!(stats.is_clean(), "nameless records are not a qname-uniqueness violation");
}

// r[verify record_store.qname_hash.no_name]
/// `0` is reserved for "no qname", so a real name must never produce it.
#[test]
fn the_empty_qname_hashes_to_the_reserved_zero() {
    assert_eq!(qname_hash(b""), 0);
    for name in [b"a".as_ref(), b"\0", b"frag", b"A00123:45:HVWXY:1:1101:12345:6789"] {
        assert_ne!(qname_hash(name), 0, "{name:?} must not take the reserved value");
    }
}

// r[verify record_store.link_mates+2]
// r[verify record_store.qname_hash.identity_is_probabilistic]
/// Two templates whose qnames collide in the hash must each still link to their
/// own mate: the hash only buckets candidates, the qname bytes decide. The
/// collision is forced, since a 64-bit hash with full avalanche cannot be
/// collided on purpose.
#[test]
fn a_hash_collision_costs_a_comparison_not_a_link() {
    let mut store = RecordStore::new();
    // Interleaved so the colliding record is scanned before the real mate.
    push(&mut store, b"frag_a", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_b", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag_a", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag_b", Read::mate(130, 50, 100, SECOND));
    for idx in 0..4 {
        store.set_qname_hash(idx, 0xC0111DEu64).expect("record exists");
    }

    let stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(2), "frag_a linked across the collision");
    assert_eq!(store.record(2).mate_idx(), Some(0));
    assert_eq!(store.record(1).mate_idx(), Some(3), "frag_b linked across the collision");
    assert_eq!(store.record(3).mate_idx(), Some(1));
    assert_eq!(stats.pairs, 2);
    assert!(stats.is_clean(), "a hash collision is not a qname-uniqueness violation");
}

// r[verify record_store.link_mates.qname_uniqueness]
// r[verify record_store.link_mates.stats]
/// Two templates sharing one qname — a merged or re-headered BAM. The first
/// free pair wins, the second pair links to itself, and the caller is told the
/// input broke the assumption.
#[test]
fn repeated_primary_qnames_pair_up_and_are_counted() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    let stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(2));
    assert_eq!(store.record(2).mate_idx(), Some(0));
    assert_eq!(store.record(1).mate_idx(), Some(3));
    assert_eq!(store.record(3).mate_idx(), Some(1));
    assert_eq!(stats.pairs, 2);
    assert_eq!(stats.ambiguous_qnames, 1, "the second copy of the left mate is reported");
    assert!(!stats.is_clean());
}

// r[verify record_store.link_mates.unique_records]
// r[verify record_store.link_mates.stats]
/// The same record twice — what a store would hold if a reader loaded it from
/// two overlapping index chunks. Indistinguishable from a repeated qname, and
/// reported the same way.
#[test]
fn a_duplicated_record_is_reported_as_ambiguous() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag", Read::mate(130, 50, 100, SECOND));
    let stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), Some(1));
    assert_eq!(store.record(1).mate_idx(), Some(0));
    assert_eq!(store.record(2).mate_idx(), None, "the third record has no free partner");
    assert_eq!(stats.pairs, 1);
    assert_eq!(stats.ambiguous_qnames, 1);
}

// r[verify record_store.link_mates+2]
/// Mate fields that do not point at each other are not a pair, however alike
/// the two records otherwise look.
#[test]
fn non_reciprocal_mate_positions_do_not_link() {
    let mut store = RecordStore::new();
    // `a` claims its mate is at 130, `b` sits at 130 but claims 99.
    push(&mut store, b"frag", Read::mate(100, 50, 130, FIRST));
    push(&mut store, b"frag", Read::mate(130, 50, 99, SECOND));
    let stats = store.link_mates();

    assert_eq!(store.record(0).mate_idx(), None);
    assert_eq!(store.record(1).mate_idx(), None);
    assert_eq!(stats.pairs, 0);
}

// r[verify record_store.link_mates+2]
// r[verify record_store.link_mates.invalidated]
/// Records pushed out of coordinate order link the same way once the store has
/// been sorted — which is the order the pileup engine requires anyway, and the
/// point at which links are dropped and must be rebuilt.
#[test]
fn linking_survives_a_sort_of_out_of_order_records() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag_b", Read::mate(300, 50, 260, SECOND));
    push(&mut store, b"frag_a", Read::mate(130, 50, 100, SECOND));
    push(&mut store, b"frag_b", Read::mate(260, 50, 300, FIRST));
    push(&mut store, b"frag_a", Read::mate(100, 50, 130, FIRST));

    store.sort_by_pos();
    let stats = store.link_mates();

    assert_eq!(stats.pairs, 2);
    assert!(stats.is_clean());
    for idx in 0..4 {
        let mate = store.record(idx).mate_idx().expect("every record is one half of a pair");
        assert_eq!(store.record(mate).mate_idx(), Some(idx), "links must be symmetric");
        assert_eq!(store.qname(idx), store.qname(mate), "linked records share a qname");
    }
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
    let _stats = store.link_mates();
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
    // r[verify record_store.link_mates+2]
    // r[verify record_store.link_mates.stats]
    /// Whatever the input, these hold: a link is symmetric, both ends carry the
    /// same qname and reciprocal positions, no nameless or non-primary record is
    /// ever linked, and every record is claimed by at most one partner. The
    /// generator deliberately produces repeated qnames, nameless records, and
    /// forced hash collisions.
    #[test]
    fn links_are_always_well_formed(
        reads in prop::collection::vec(
            (
                0usize..4,        // qname index, so qnames repeat
                100i32..300,      // pos
                20u32..60,        // len
                100i32..300,      // mate pos
                0usize..5,        // flag shape
                any::<bool>(),    // nameless
            ),
            0..12,
        ),
        collide in any::<bool>(),
    ) {
        const NAMES: [&[u8]; 4] = [b"a", b"bb", b"ccc", b"dddd"];
        let mut store = RecordStore::new();
        for &(name_idx, pos, len, mate_pos, shape, nameless) in &reads {
            let flags = match shape {
                0 => PAIRED | FIRST,
                1 => PAIRED | SECOND,
                2 => PAIRED | FIRST | SUPPLEMENTARY,
                3 => PAIRED | SECOND | SECONDARY,
                _ => 0,
            };
            let read = Read::mate(pos, len, mate_pos, 0).with_flags(flags);
            let name: &[u8] = if nameless { b"" } else { NAMES[name_idx] };
            push(&mut store, name, read);
        }
        if collide {
            // Every record hashes the same: linking must fall back entirely on
            // the qname bytes.
            for idx in 0..store.len() as u32 {
                store.set_qname_hash(idx, 1).expect("record exists");
            }
        }

        let stats = store.link_mates();

        let mut linked_count = 0u32;
        for idx in 0..store.len() as u32 {
            let rec = store.record(idx);
            let Some(partner) = rec.mate_idx() else { continue };
            linked_count += 1;
            prop_assert_ne!(partner, idx, "a record must not link to itself");
            prop_assert_eq!(store.record(partner).mate_idx(), Some(idx), "links are symmetric");
            prop_assert_eq!(store.qname(idx), store.qname(partner));
            prop_assert!(!store.qname(idx).is_empty(), "a nameless record must never link");
            prop_assert!(rec.qname_hash().is_some());
            prop_assert_eq!(rec.pos.as_i32(), store.record(partner).next_pos);
            prop_assert_eq!(store.record(partner).pos.as_i32(), rec.next_pos);
            prop_assert!(rec.flags.is_paired() && !rec.flags.is_unmapped());
            prop_assert!(!rec.flags.is_secondary() && !rec.flags.is_supplementary());
        }
        prop_assert_eq!(linked_count, stats.pairs * 2, "stats must count each pair once");
    }

    // r[verify record_store.link_mates+2]
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

        let _stats = store.link_mates();

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
                prop_assert_eq!(store.qname(idx), store.qname(partner), "links need equal qnames");
                prop_assert_eq!(
                    store.record(idx).pos.as_i32(),
                    store.record(partner).next_pos,
                    "links need reciprocal mate positions"
                );
                prop_assert_eq!(
                    store.record(partner).pos.as_i32(),
                    store.record(idx).next_pos
                );
                prop_assert_eq!(
                    store.record(idx).tid,
                    store.record(partner).tid,
                    "links never cross contigs"
                );
                let a = store.record(idx);
                let b = store.record(partner);
                let start = a.pos.max(b.pos);
                let last = a.end_pos.min(b.end_pos);
                let end = if last < start {
                    start
                } else {
                    Pos0::new(last.as_u64() as u32 + 1).unwrap()
                };
                prop_assert_eq!(store.mate_overlap(idx), Some(start..end));
                prop_assert_eq!(store.mate_overlap(partner), Some(start..end));
            } else {
                prop_assert_eq!(store.mate_overlap(idx), None);
            }
        }
    }
}

// ── through the readers, on real files ────────────────────────────────────

fn fixture(name: &str) -> &'static Path {
    // Leaked once per call site; these are test-lifetime paths.
    Box::leak(
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/"))
            .join(name)
            .into_boxed_path(),
    )
}

fn test_fasta() -> &'static Path {
    fixture("test.fasta.gz")
}

/// The region of `test.bam` / `test.cram` that actually has reads.
const REGION: (u32, u32) = (6_105_000, 6_140_000);

/// Pile up the fixture region through `Readers` and return, per column,
/// how many alignments carry a mate link and how many are inside the overlap.
fn link_counts_through_readers(path: &Path) -> (usize, usize, usize) {
    let mut readers = Readers::open(path, test_fasta()).expect("open readers");
    let tid = readers.header().tid("chr19").expect("chr19 in header");
    let contig = readers.header().target_name(tid).expect("contig name").to_string();
    let opts = SegmentOptions::new(NonZeroU32::new(2_000_000).unwrap());
    let start = Pos0::new(REGION.0).unwrap();
    let end = Pos0::new(REGION.1).unwrap();
    let segment = readers
        .segments((contig.as_str(), start, end), opts)
        .expect("segment plan")
        .next()
        .expect("one segment");

    let (mut alignments, mut linked, mut in_overlap) = (0, 0, 0);
    let mut engine = readers.pileup(&segment, DepthLimit::Unlimited).expect("pileup");
    while let Some(col) = engine.pileups() {
        for aln in col.raw_alignments() {
            alignments += 1;
            if let Some(mate_idx) = aln.mate_idx() {
                linked += 1;
                // The link must be usable: the mate is a real record, and it
                // agrees that this record is *its* mate.
                let mate = col.store().record(mate_idx);
                assert_eq!(
                    mate.mate_idx(),
                    Some(aln.record_idx()),
                    "mate link must be symmetric through the store"
                );
                assert_eq!(mate.qname_hash(), aln.qname_hash(), "linked records share a qname");
            }
            if aln.in_mate_overlap() {
                in_overlap += 1;
                assert!(aln.mate_idx().is_some(), "only a linked record can be in an overlap");
            }
        }
    }
    (alignments, linked, in_overlap)
}

/// The BAM re-encoded as bgzf SAM, so the reader tests cover all three formats.
/// Returns `None` when samtools is not on `PATH`.
fn sam_gz_of_the_bam(dir: &Path) -> Option<std::path::PathBuf> {
    let sam_gz = dir.join("test.sam.gz");
    let ok = std::process::Command::new("samtools")
        .args(["view", "-h", "--output-fmt", "SAM,level=6", "-o"])
        .arg(&sam_gz)
        .arg(fixture("test.bam"))
        .status()
        .ok()?
        .success();
    if !ok {
        return None;
    }
    let indexed = std::process::Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz)
        .status()
        .ok()?
        .success();
    indexed.then_some(sam_gz)
}

// r[verify record_store.link_mates.invalidated]
// r[verify pileup.mate_link_cache]
/// `Readers::pileup` must link mates itself — a consumer that never calls
/// `link_mates` still gets usable links, on every format.
#[test]
fn readers_pileup_links_mates_on_every_format() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mut paths: Vec<std::path::PathBuf> =
        ["test.bam", "test.cram", "test_v30.cram"].iter().map(|n| fixture(n).to_owned()).collect();
    // SAM.gz has to be built from the BAM; skip it where samtools is missing
    // rather than failing a test that is about mate links.
    paths.extend(sam_gz_of_the_bam(dir.path()));

    for path in paths {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("?").to_string();
        let (alignments, linked, in_overlap) = link_counts_through_readers(&path);
        assert!(alignments > 0, "{name}: fixture produced no alignments");
        assert!(
            linked > 0,
            "{name}: no alignment carried a mate link; Readers::pileup did not link"
        );
        assert!(
            in_overlap > 0,
            "{name}: no alignment fell inside a mate overlap, so dedup would never fire"
        );
    }
}

// r[verify record_store.link_mates.unique_records]
/// Linking assumes a record appears once. Every reader must deliver that, or
/// two copies of one alignment would pair up and be counted twice.
#[test]
fn a_fetched_store_holds_each_record_once() {
    for name in ["test.bam", "test.cram", "test_v30.cram"] {
        let path = fixture(name);
        let mut readers = Readers::open(path, test_fasta()).expect("open readers");
        let tid = readers.header().tid("chr19").expect("chr19 in header");
        let mut store = RecordStore::new();
        readers
            .fetch_into(tid, Pos0::new(REGION.0).unwrap(), Pos0::new(REGION.1).unwrap(), &mut store)
            .expect("fetch");
        assert!(!store.is_empty(), "{name}: no records fetched");
        let stats = store.link_mates();

        // Only the records that can link have to be unique: `test.bam` really
        // does contain one qname twice as a *placed-unmapped* mate (flag 133 at
        // chr19:6105793), which is exactly the kind of thing real files hold and
        // exactly why linking rejects unmapped records outright.
        let mut seen = HashSet::new();
        for idx in 0..store.len() as u32 {
            let rec = store.record(idx);
            let flags = rec.flags;
            let linkable = flags.is_paired()
                && !flags.is_unmapped()
                && !flags.is_secondary()
                && !flags.is_supplementary();
            if !linkable {
                assert_eq!(
                    rec.mate_idx(),
                    None,
                    "{name}: a non-primary or unmapped record must never be linked"
                );
                continue;
            }
            let identity = (rec.pos, flags.raw(), store.qname(idx).to_vec());
            assert!(
                seen.insert(identity),
                "{name}: linkable record {idx} (pos {}) appears twice in one fetched store",
                rec.pos.as_u64()
            );
        }

        assert!(
            stats.is_clean(),
            "{name}: {} records could not be linked despite a same-qname candidate",
            stats.ambiguous_qnames
        );
        assert!(stats.pairs > 0, "{name}: no pairs linked in a paired-end fixture");
    }
}
