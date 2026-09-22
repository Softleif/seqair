//! Per-record extras, carried from the push that computed them all the way
//! out to the alignment a column reports.
//!
//! `record_store_extras` covers the store side: an extra survives `sort_by_pos`
//! and `dedup` because `SlimRecord` addresses it through `extras_idx` rather
//! than by record position. The pileup side is where that indirection is
//! actually exercised — `prepare_for_pileup` sorts the store, so a read's store
//! index is no longer the index it was pushed at, and every column entry then
//! resolves its extra through two hops (`record_idx` → `SlimRecord` →
//! `extras_idx`). An off-by-one in either hop hands a caller another read's
//! payload, and nothing about the resulting pileup looks wrong.
//!
//! So the payload is chosen to make a mis-association visible: each read's
//! extra is its *push* ordinal plus a hash of its own (unique) qname. Push
//! order is deliberately not position order, so the ordinal and the store index
//! disagree for almost every read, and the oracle is the generated input list
//! keyed by qname — not the store, which is the thing under test.
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

mod helpers;

use hegel::prelude::*;
use helpers::SyntheticRead;
use seqair::bam::pileup::{PileupEngine, PileupOp};
use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};
use seqair::bam::{Pos0, RecordIdx};
use seqair_types::Base;
use std::collections::HashMap;

// ── the payload ─────────────────────────────────────────────────────────

/// What `compute` attaches to a record: where it was pushed, and which read it
/// was. Either half alone would be a weak witness — the ordinal repeats across
/// runs and the hash repeats across a read's columns — but together they name
/// exactly one generated read at exactly one point in the push sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Tag {
    push_ordinal: u32,
    qname_hash: u64,
}

/// FNV-1a. Deliberately not `RecordStore`'s own `qname_hash`: the oracle must
/// not be able to agree with the store by sharing its arithmetic.
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The qname of the `i`-th generated read. Unique, so no two reads can share a
/// `Tag` and no mate link can form to blur the two.
fn qname_for(i: usize) -> Vec<u8> {
    format!("q{i:03}").into_bytes()
}

/// Tags each pushed record with its push ordinal and its own name's hash.
#[derive(Clone, Default)]
struct TagByQname {
    next: u32,
}

impl CustomizeRecordStore for TagByQname {
    type Extra = Tag;

    fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<Tag>) -> Tag {
        let push_ordinal = self.next;
        self.next = self.next.wrapping_add(1);
        let qname = rec.qname(store).expect("a just-pushed record's name is in the slab");
        Tag { push_ordinal, qname_hash: fnv1a(qname) }
    }
}

// ── inputs ──────────────────────────────────────────────────────────────

/// Reads in *push* order, which is not position order: `arb_read_set` sorts,
/// and a sorted push would let `extras_idx` be a plain record index and still
/// pass. Every CIGAR class the repo generates is in here — clips, deletions,
/// insertions, ref-skips, `=`/`X` — so the columns include entries with no
/// query base at all.
#[hegel::composite]
fn arb_reads(tc: &TestCase) -> Vec<SyntheticRead> {
    tc.draw_silent(gs::vecs(helpers::arb_read()).min_size(1).max_size(12))
}

/// `helpers::make_record_with_cigar` names every record `read`; the oracle
/// needs the names to tell the reads apart. Same layout, different name.
fn raw_with_qname(qname: &[u8], read: &SyntheticRead) -> Vec<u8> {
    let mut name = qname.to_vec();
    name.push(0);
    let name_len = name.len();

    let packed: Vec<u32> =
        read.cigar_ops.iter().map(|&(len, op)| (len << 4) | u32::from(op)).collect();
    let cigar_bytes = packed.len() * 4;
    let seq_bytes = (read.seq_len as usize).div_ceil(2);
    let total = 32 + name_len + cigar_bytes + seq_bytes + read.seq_len as usize;

    let mut raw = vec![0u8; total];
    raw[0..4].copy_from_slice(&0i32.to_le_bytes());
    raw[4..8].copy_from_slice(&read.pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = read.mapq;
    raw[12..14].copy_from_slice(&(packed.len() as u16).to_le_bytes());
    raw[14..16].copy_from_slice(&read.flags.to_le_bytes());
    raw[16..20].copy_from_slice(&read.seq_len.to_le_bytes());
    raw[32..32 + name_len].copy_from_slice(&name);

    let cigar_start = 32 + name_len;
    for (i, op) in packed.iter().enumerate() {
        raw[cigar_start + i * 4..cigar_start + (i + 1) * 4].copy_from_slice(&op.to_le_bytes());
    }

    let seq_start = cigar_start + cigar_bytes;
    for i in 0..seq_bytes {
        raw[seq_start + i] = 0x11; // A, A
    }
    let qual_start = seq_start + seq_bytes;
    for i in 0..read.seq_len as usize {
        raw[qual_start + i] = 30;
    }
    raw
}

/// The region that holds every read, with a little slack past the last one so
/// the engine has to run out of records rather than out of region.
fn region(reads: &[SyntheticRead]) -> (Pos0, Pos0) {
    let start = reads.iter().map(|r| r.pos as u32).min().unwrap_or(0);
    let end = reads.iter().map(|r| r.pos as u32 + r.ref_span + 10).max().unwrap_or(start);
    (Pos0::new(start).unwrap(), Pos0::new(end).unwrap())
}

fn store_with<E: CustomizeRecordStore>(
    reads: &[SyntheticRead],
    customize: &mut E,
) -> RecordStore<E::Extra> {
    let mut store = RecordStore::new();
    for (i, read) in reads.iter().enumerate() {
        store.push_raw(&raw_with_qname(&qname_for(i), read), customize).unwrap();
    }
    store
}

/// One column entry, reduced to what does not depend on the extras type.
type Entry = (RecordIdx, PileupOp);

/// The whole column stream: position, reference base, and every entry in the
/// order the column reports them.
fn column_stream<U>(
    reads: &[SyntheticRead],
    store: RecordStore<U>,
) -> Vec<(Pos0, Base, Vec<Entry>)> {
    let (start, end) = region(reads);
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, (start..=end).into());
    let mut out = Vec::new();
    while let Some(col) = engine.pileups() {
        let entries =
            col.raw_alignments().map(|aln| (aln.record_idx(), aln.op)).collect::<Vec<Entry>>();
        out.push((col.pos(), col.reference_base(), entries));
    }
    out
}

/// A short name for the op, for the coverage events.
fn op_name(op: &PileupOp) -> &'static str {
    match op {
        PileupOp::Match { .. } => "match",
        PileupOp::Insertion { .. } => "insertion",
        PileupOp::Deletion { .. } => "deletion",
        PileupOp::ComplexIndel { .. } => "complex indel",
        PileupOp::RefSkip => "ref-skip",
        PileupOp::SoftClip { .. } => "soft clip",
    }
}

// ── properties ──────────────────────────────────────────────────────────

// r[verify pileup.extras.generic_param]
// r[verify pileup.alignment_view]
// r[verify record_store.extras.access]
// r[verify record_store.extras.sort_dedup_generic]
/// Every alignment in every column carries the extra of the record it came
/// from — including the alignments that carry no query base at all.
///
/// The oracle is the generated read list keyed by qname, so the assertion is
/// "the payload the pileup hands back is the one this read was tagged with at
/// push time", not "the payload agrees with wherever the store is currently
/// pointing".
#[hegel::test]
fn every_alignment_carries_its_own_records_extra(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());

    // What each read was tagged with, derived from the input alone.
    let expected: HashMap<Vec<u8>, Tag> = reads
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let name = qname_for(i);
            let hash = fnv1a(&name);
            (name, Tag { push_ordinal: i as u32, qname_hash: hash })
        })
        .collect();

    let store = store_with(&reads, &mut TagByQname::default());
    let (start, end) = region(&reads);
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, (start..=end).into());

    let mut checked = 0usize;
    let mut without_base = 0usize;
    let mut ops_seen: Vec<&'static str> = Vec::new();
    while let Some(col) = engine.pileups() {
        for aln in col.alignments() {
            let name = aln.qname();
            let want = expected.get(name).unwrap_or_else(|| {
                panic!("column at {} reported a qname no read was pushed with", col.pos())
            });
            assert_eq!(
                aln.extra(),
                want,
                "extra at {} for {}",
                col.pos(),
                String::from_utf8_lossy(name)
            );
            // The same answer through the store-facing spelling: `extra()` on
            // the view is documented to be the record's extra, and a column
            // entry that resolved to the wrong record would agree with itself
            // here while disagreeing with the oracle above.
            assert_eq!(aln.record().extra(), want, "record().extra() at {}", col.pos());

            if aln.qpos().is_none() {
                without_base += 1;
            }
            ops_seen.push(op_name(&aln.op));
            checked += 1;
        }
    }

    assert!(checked > 0, "a set of mapped reads must produce at least one alignment");

    for op in ops_seen {
        tc.event(op);
    }
    tc.event_value("alignments", checked as f64);
    if without_base > 0 {
        tc.event("an alignment with no query base carried its extra");
    }
    if reads.len() > 1 && reads.windows(2).any(|w| w[0].pos > w[1].pos) {
        tc.event("pushed out of position order");
    }
}

// r[verify pileup.extras.generic_param]
/// Attaching extras must not change the pileup itself. The same reads through a
/// `RecordStore<()>` and a `RecordStore<Tag>` must produce the same columns, in
/// the same order, with the same entries — the type parameter buys a payload,
/// not a different traversal.
#[hegel::test]
fn extras_do_not_change_the_column_stream(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());

    let unit = column_stream(&reads, store_with(&reads, &mut ()));
    let tagged = column_stream(&reads, store_with(&reads, &mut TagByQname::default()));

    assert_eq!(unit.len(), tagged.len(), "column count");
    for (a, b) in unit.iter().zip(&tagged) {
        assert_eq!(a, b, "column at {}", a.0);
    }

    tc.event_value("columns", unit.len() as f64);
    tc.event_value("entries", unit.iter().map(|c| c.2.len()).sum::<usize>() as f64);
}
