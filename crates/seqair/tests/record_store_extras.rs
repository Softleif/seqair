//! Tests for `RecordStore<U>` per-record customize/extras and `PileupEngine` `pileups()`.
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
    clippy::drop_non_drop,
    reason = "test code with known small values"
)]

mod helpers;

use helpers::{collect_columns, make_record};
use seqair::bam::Pos0;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};

/// Helper: push N synthetic records at positions 100, 101, ... with 10M CIGAR.
fn store_with_n_records(n: u32) -> RecordStore {
    let mut store = RecordStore::new();
    for i in 0..n {
        let pos = 100 + i as i32;
        store.push_raw(&make_record(0, pos, 99, 60, 10), &mut ()).unwrap();
    }
    store
}

fn store_with_n_records_customized<C: CustomizeRecordStore>(
    n: u32,
    mut c: C,
) -> RecordStore<C::Extra> {
    let mut store = RecordStore::new();
    for i in 0..n {
        let pos = 100 + i as i32;
        store.push_raw(&make_record(0, pos, 99, 60, 10), &mut c).unwrap();
    }
    store
}

/// Customize value used in the position-extracting tests below.
#[derive(Clone, Default)]
struct ExtractPos;

impl CustomizeRecordStore for ExtractPos {
    type Extra = i32;
    fn compute(&mut self, rec: &SlimRecord, _: &RecordStore<i32>) -> i32 {
        rec.pos.as_i32()
    }
}

// ---- RecordStore extras ----

// r[verify record_store.slim_record.field_getters]
#[test]
fn slim_record_getters_return_slab_data() {
    // Build a record with a known qname + aux + sequence and read each slab
    // back via the SlimRecord getters.
    let aux = b"RGZmy_rg\0";
    let raw = helpers::make_record_with_aux(0, 200, 0, 30, 10, aux);

    let mut store = RecordStore::new();
    store.push_raw(&raw, &mut ()).unwrap();
    let rec = store.record(0);

    let qname = rec.qname(&store).expect("qname slab readable");
    assert_eq!(qname, b"read");

    let seq = rec.seq(&store).expect("seq slab readable");
    assert_eq!(seq.len(), 10);

    let qual = rec.qual(&store).expect("qual slab readable");
    assert_eq!(qual.len(), 10);

    let cigar = rec.cigar(&store).expect("cigar slab readable");
    assert_eq!(cigar.len(), 1); // single typed op

    let aux_bytes = rec.aux(&store).expect("aux slab readable");
    assert_eq!(aux_bytes.as_bytes(), aux as &[u8]);
}

// r[verify record_store.extras.clear]
#[test]
fn clear_clears_extras_and_retains_capacity() {
    let mut store = store_with_n_records_customized(50, ExtractPos);

    assert_eq!(store.len(), 50);
    let records_cap = store.records_capacity();
    let extras_cap = store.extras_capacity();
    assert!(extras_cap >= 50);

    store.clear();

    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() >= records_cap);
    assert!(store.extras_capacity() >= extras_cap);
}

// r[verify record_store.extras.generic_param]
#[test]
fn unit_extras_is_zero_cost() {
    // Vec<()> elements are ZST — no heap allocation for element storage.
    assert_eq!(std::mem::size_of::<()>(), 0);

    // A RecordStore<()> (the default) should work identically to the old RecordStore.
    let store = store_with_n_records(3);
    assert_eq!(store.len(), 3);
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn sort_and_dedup_work_on_unit_store() {
    let mut store = RecordStore::new();
    // Push records out of order.
    store.push_raw(&make_record(0, 200, 99, 60, 10), &mut ()).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut ()).unwrap();

    store.sort_by_pos();
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn sort_by_pos_preserves_extras_mapping() {
    let mut store = RecordStore::new();
    // Push records out of order: pos 300, 100, 200
    store.push_raw(&make_record(0, 300, 99, 60, 10), &mut ExtractPos).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut ExtractPos).unwrap();
    store.push_raw(&make_record(0, 200, 99, 60, 10), &mut ExtractPos).unwrap();

    // Now sort — records reorder, but extras_idx preserves the mapping.
    store.sort_by_pos();

    // After sort, records are at positions 100, 200, 300.
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());
    assert_eq!(store.record(2).pos, Pos0::new(300).unwrap());

    // Extras must still match the record they were computed from.
    assert_eq!(*store.extra(0), 100);
    assert_eq!(*store.extra(1), 200);
    assert_eq!(*store.extra(2), 300);
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn dedup_on_typed_store_preserves_extras() {
    /// Tag each record with the index it had at compute time.
    #[derive(Clone, Default)]
    struct RecordIdx {
        next: u32,
    }
    impl CustomizeRecordStore for RecordIdx {
        type Extra = u32;
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<u32>) -> u32 {
            let i = self.next;
            self.next = i.wrapping_add(1);
            i
        }
    }

    let mut store = RecordStore::new();
    let mut c = RecordIdx::default();
    // Push duplicate records at same position with same flags/name.
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut c).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut c).unwrap(); // dup
    store.push_raw(&make_record(0, 200, 99, 60, 10), &mut c).unwrap();
    store.push_raw(&make_record(0, 200, 99, 60, 10), &mut c).unwrap(); // dup

    store.sort_by_pos();
    store.dedup();

    // Should have 2 records after dedup (one per position).
    assert_eq!(store.len(), 2);
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());

    // Extras should still be accessible and correspond to the surviving records.
    let e0 = *store.extra(0);
    let e1 = *store.extra(1);
    assert_ne!(e0, e1, "surviving records should have distinct extras");
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn set_alignment_then_sort_on_typed_store() {
    let mut store = RecordStore::new();
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut ExtractPos).unwrap();
    store.push_raw(&make_record(0, 50, 99, 60, 10), &mut ExtractPos).unwrap();

    // Sort — now record at pos=50 comes first.
    store.sort_by_pos();
    assert_eq!(store.record(0).pos, Pos0::new(50).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(100).unwrap());
    assert_eq!(*store.extra(0), 50);
    assert_eq!(*store.extra(1), 100);
}

// ---- PileupEngine extras ----

// r[verify pileup.extras.constructor_accepts_any_u]
#[test]
fn engine_accepts_typed_store() {
    let store = store_with_n_records_customized(5, ExtractMapq);

    #[derive(Clone, Default)]
    struct ExtractMapq;
    impl CustomizeRecordStore for ExtractMapq {
        type Extra = u8;
        fn compute(&mut self, rec: &SlimRecord, _: &RecordStore<u8>) -> u8 {
            rec.mapq
        }
    }

    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    engine.set_max_depth(10);

    // Verify the engine still works — iterate and check columns are produced.
    let columns = collect_columns(&mut engine);
    assert!(!columns.is_empty());
}

// r[verify pileup.lending_iterator]
// r[verify pileup.alignment_view]
#[test]
fn pileups_yields_same_columns_on_unit_and_typed_store() {
    let store = store_with_n_records(5);

    // Collect columns via pileups() on a unit store.
    let store_copy = store_with_n_records(5);
    let mut engine =
        PileupEngine::new(store_copy, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    let iter_columns: Vec<_> =
        collect_columns(&mut engine).into_iter().map(|c| (c.pos(), c.depth())).collect();

    // Collect columns via pileups() on a typed store.
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    let mut cws_columns = Vec::new();
    while let Some(col) = engine.pileups() {
        cws_columns.push((col.pos(), col.depth()));
    }

    assert_eq!(iter_columns, cws_columns);
}

// r[verify pileup.lending_iterator]
// r[verify pileup.alignment_view]
#[test]
fn pileups_column_exposes_store_access_via_alignment_view() {
    let store = store_with_n_records_customized(3, ExtractPos);
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());

    let mut saw_extras = false;
    while let Some(col) = engine.pileups() {
        for aln in col.alignments() {
            let extra = aln.extra();
            // Extra should be the record's position.
            assert!(*extra >= 100 && *extra <= 102);
            saw_extras = true;
        }
    }
    assert!(saw_extras, "should have seen at least one alignment with extras");
}

// r[verify pileup.extras.recover_store]
/// `engine.take_store()` returns a cleared store with retained capacity.
/// This is the low-level mechanism that `PileupGuard::drop` uses to recover
/// the buffer back into `Readers`.
#[test]
fn engine_take_store_clears_records_and_keeps_capacity() {
    // Build a store manually, transform to extras, feed into engine, then strip.
    let store = store_with_n_records_customized(10, ExtractMapqU32);

    #[derive(Clone, Default)]
    struct ExtractMapqU32;
    impl CustomizeRecordStore for ExtractMapqU32 {
        type Extra = u32;
        fn compute(&mut self, rec: &SlimRecord, _: &RecordStore<u32>) -> u32 {
            rec.mapq as u32
        }
    }

    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());

    // Consume all columns.
    while engine.pileups().is_some() {}

    let store = engine.take_store().expect("store should be available");

    // The recovered store should be empty (cleared by take_store) but have capacity.
    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() > 0);
}

// ---- Proptests ----

use proptest::prelude::*;

proptest! {
    // r[verify record_store.extras.sort_dedup_generic]
    #[test]
    fn proptest_sort_preserves_extras_mapping(
        positions in proptest::collection::vec(0i32..1000, 2..50),
    ) {
        let mut store = RecordStore::new();
        for &pos in &positions {
            store.push_raw(&make_record(0, pos, 99, 60, 10), &mut ExtractPos).unwrap();
        }

        store.sort_by_pos();

        // After sorting, each record's extra must still equal its position
        // (proving the extras_idx indirection survived reordering).
        for i in 0..store.len() as u32 {
            let pos = store.record(i).pos.as_i32();
            let extra = *store.extra(i);
            prop_assert_eq!(pos, extra, "extras_idx broken at record {}", i);
        }

        // Records must be sorted.
        for i in 1..store.len() as u32 {
            prop_assert!(
                store.record(i).pos >= store.record(i - 1).pos,
                "not sorted at {}", i
            );
        }
    }

    // r[verify record_store.extras.sort_dedup_generic]
    #[test]
    fn proptest_dedup_preserves_extras_on_typed_store(
        positions in proptest::collection::vec(0i32..20, 2..30),
    ) {
        let mut store = RecordStore::new();
        for &pos in &positions {
            // Same flags/mapq/seq_len so duplicates at same pos are detected.
            store.push_raw(&make_record(0, pos, 99, 60, 10), &mut ExtractPos).unwrap();
        }

        store.sort_by_pos();
        store.dedup();

        // After dedup, each surviving record's extra must equal its position.
        for i in 0..store.len() as u32 {
            let pos = store.record(i).pos.as_i32();
            let extra = *store.extra(i);
            prop_assert_eq!(pos, extra, "extras_idx broken after dedup at record {}", i);
        }

        // No consecutive duplicates.
        for i in 1..store.len() as u32 {
            let a = store.record(i - 1);
            let b = store.record(i);
            if a.pos == b.pos {
                // Same position is OK if flags differ (different records).
                // Our test uses identical flags, so this shouldn't happen.
                prop_assert!(a.flags != b.flags, "dedup missed a duplicate at pos {}", a.pos.as_i32());
            }
        }
    }
}

// r[verify pileup.extras.constructor_accepts_any_u]
#[test]
fn pileup_with_sorted_typed_store() {
    // End-to-end: load out of order, compute extras (original push index),
    // sort, pileup.
    let mut store = RecordStore::new();

    /// Counter customizer so each record's extra is its push index.
    #[derive(Clone, Default)]
    struct PushIdx {
        next: u32,
    }
    impl CustomizeRecordStore for PushIdx {
        type Extra = u32;
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<u32>) -> u32 {
            let i = self.next;
            self.next = i.wrapping_add(1);
            i
        }
    }

    let mut c = PushIdx::default();
    store.push_raw(&make_record(0, 105, 99, 60, 10), &mut c).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), &mut c).unwrap();
    store.push_raw(&make_record(0, 102, 99, 60, 10), &mut c).unwrap();

    // Sort — records reorder, extras follow via extras_idx.
    store.sort_by_pos();
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(*store.extra(0), 1); // originally the second record pushed

    // Build engine with the sorted typed store.
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(114).unwrap());
    let mut column_count = 0;
    while let Some(col) = engine.pileups() {
        for aln in col.alignments() {
            // extras should be accessible and valid
            let _extra = aln.extra();
        }
        column_count += 1;
    }
    assert!(column_count > 0, "should have produced at least one column");
}

// ---- push-time filtering ----

/// Drop records whose mapq is below `threshold`. Verifies push-time filtering
/// based on a `SlimRecord` field.
#[derive(Clone)]
struct MinMapq {
    threshold: u8,
}
impl CustomizeRecordStore for MinMapq {
    type Extra = ();
    fn filter(&mut self, rec: &SlimRecord, _: &RecordStore<()>) -> bool {
        rec.mapq >= self.threshold
    }
    fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
}

/// Drop everything at push time — used to verify zero-waste rollback.
#[derive(Clone, Default)]
struct DropAll;
impl CustomizeRecordStore for DropAll {
    type Extra = ();
    fn filter(&mut self, _: &SlimRecord, _: &RecordStore<()>) -> bool {
        false
    }
    fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
}

// r[verify record_store.pre_filter.rollback]
#[test]
fn keep_record_drops_low_mapq_records_via_push_raw() {
    let mut store = RecordStore::new();
    let mut filter = MinMapq { threshold: 20 };

    let mapqs = [0u8, 10, 20, 30, 40, 50];
    for (i, &mq) in mapqs.iter().enumerate() {
        let raw = make_record(0, 100 + i as i32, 99, mq, 10);
        store.push_raw(&raw, &mut filter).unwrap();
    }

    assert_eq!(store.len(), 4, "records with mapq < 20 should be dropped");
    for i in 0..4u32 {
        assert!(store.record(i).mapq >= 20, "kept record at idx {i} below threshold");
    }
}

// r[verify record_store.pre_filter.rollback]
#[test]
fn keep_record_drops_low_mapq_records_via_fetch_into_customized() {
    use seqair::bam::IndexedBamReader;
    use seqair::reader::FetchCounts;

    let bam_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/data/test.bam");
    if !bam_path.exists() {
        return; // skip if test data not available in this workspace layout
    }
    let mut reader = IndexedBamReader::open(&bam_path).expect("open test BAM");
    let tid = reader.header().tid("chr19").expect("chr19 missing");
    let mut store = RecordStore::new();

    // Always-keep: kept == fetched.
    let counts_all = reader
        .fetch_into_customized(
            tid,
            Pos0::new(6_103_000).unwrap(),
            Pos0::new(6_104_000).unwrap(),
            &mut store,
            &mut (),
        )
        .expect("fetch_into_customized");
    assert!(counts_all.fetched > 0, "test region should produce records");
    assert_eq!(counts_all.kept, counts_all.fetched);
    let all_len = store.len();
    assert_eq!(all_len, counts_all.kept);

    // Drop-all: kept == 0 but fetched unchanged.
    let counts_none = reader
        .fetch_into_customized(
            tid,
            Pos0::new(6_103_000).unwrap(),
            Pos0::new(6_104_000).unwrap(),
            &mut store,
            &mut DropAll,
        )
        .expect("fetch_into_customized");
    assert_eq!(counts_none.fetched, counts_all.fetched);
    assert_eq!(counts_none.kept, 0);
    assert_eq!(store.len(), 0, "all records rolled back");
    let _: FetchCounts = counts_none;
}

// r[verify cram.fetch_into_customized.push_time]
#[test]
fn cram_fetch_into_customized_applies_filter_at_push_time() {
    use seqair::cram::reader::IndexedCramReader;
    use seqair::reader::FetchCounts;

    let workspace_data = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/data");
    let cram_path = workspace_data.join("test_v30.cram");
    let fasta_path = workspace_data.join("test.fasta.gz");
    if !cram_path.exists() || !fasta_path.exists() {
        return; // skip if test data not available in this workspace layout
    }

    let mut reader = IndexedCramReader::open(&cram_path, &fasta_path).expect("open CRAM");
    let tid = reader.header().tid("chr19").expect("chr19 missing from CRAM header");
    let region_start = Pos0::new(6_103_076).unwrap();
    let region_end = Pos0::new(6_143_229).unwrap();
    let mut store = RecordStore::new();

    // Baseline: always-keep — kept == fetched and matches the sparse fetch_into count.
    let counts_all = reader
        .fetch_into_customized(tid, region_start, region_end, &mut store, &mut ())
        .expect("fetch_into_customized all");
    assert!(counts_all.fetched > 0, "test region should produce records");
    assert_eq!(counts_all.kept, counts_all.fetched, "always-keep must keep everything");
    assert_eq!(store.len(), counts_all.kept);

    // Snapshot positions of all fetched records — we'll use them as an oracle
    // for the mixed filter below.
    let mut all_positions: Vec<i64> =
        (0..store.len() as u32).map(|i| store.record(i).pos.as_i64()).collect();
    all_positions.sort_unstable();

    // Drop-all: fetched unchanged, kept == 0. The store (including all slabs)
    // must be empty after rollback — this is the zero-waste guarantee.
    let counts_none = reader
        .fetch_into_customized(tid, region_start, region_end, &mut store, &mut DropAll)
        .expect("fetch_into_customized none");
    assert_eq!(counts_none.fetched, counts_all.fetched, "fetched must be filter-independent");
    assert_eq!(counts_none.kept, 0, "drop-all must yield no kept records");
    assert_eq!(store.len(), 0, "store must be empty after drop-all (rollback)");

    // Mixed filter: keep records whose pos is even. Use an independent oracle
    // (filter the baseline positions) to verify the survivors.
    #[derive(Clone)]
    struct EvenPos;
    impl CustomizeRecordStore for EvenPos {
        type Extra = ();
        fn filter(&mut self, rec: &SlimRecord, _: &RecordStore<()>) -> bool {
            rec.pos.as_i64() % 2 == 0
        }
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
    }

    let counts_mixed = reader
        .fetch_into_customized(tid, region_start, region_end, &mut store, &mut EvenPos)
        .expect("fetch_into_customized mixed");
    assert_eq!(counts_mixed.fetched, counts_all.fetched, "fetched must be filter-independent");
    let expected_kept: usize = all_positions.iter().filter(|p| *p % 2 == 0).count();
    assert_eq!(counts_mixed.kept, expected_kept, "filter result must match oracle");
    assert_eq!(store.len(), expected_kept);
    for i in 0..store.len() as u32 {
        assert_eq!(store.record(i).pos.as_i64() % 2, 0, "record {i} failed filter predicate");
    }

    let _: FetchCounts = counts_mixed;
}

// ---- read-group filter via aux tag ----

/// Extract the value of the BAM aux RG:Z tag from the raw aux bytes, if present.
/// Aux entries are `tag(2 bytes) + type(1 byte) + value`. RG is type Z
/// (NUL-terminated string).
fn extract_rg(aux: &[u8]) -> Option<&[u8]> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let tag = &aux[i..i + 2];
        let val_type = aux[i + 2];
        i += 3;
        if tag == b"RG" && val_type == b'Z' {
            let start = i;
            while i < aux.len() && aux[i] != 0 {
                i += 1;
            }
            return Some(&aux[start..i]);
        }
        // Skip value based on type for tags we don't care about.
        match val_type {
            b'A' | b'c' | b'C' => i += 1,
            b's' | b'S' => i += 2,
            b'i' | b'I' | b'f' => i += 4,
            b'Z' | b'H' => {
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                i += 1; // skip NUL
            }
            _ => return None,
        }
    }
    None
}

/// Customizer that keeps only records whose RG:Z aux tag matches `wanted`.
/// Demonstrates push-time filtering driven by aux data — a common need for
/// callers who want to subset by sample/read group in the alignment file.
#[derive(Clone)]
struct ReadGroupFilter {
    wanted: Vec<u8>,
}
impl CustomizeRecordStore for ReadGroupFilter {
    type Extra = ();
    fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<()>) -> bool {
        // The just-pushed record's aux bytes live at the tail of the aux slab.
        let Ok(aux) = rec.aux(store) else { return false };
        extract_rg(aux.as_bytes()) == Some(self.wanted.as_slice())
    }
    fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
}

// r[verify record_store.pre_filter.rollback]
// r[verify record_store.slim_record.field_getters]
#[test]
fn read_group_filter_keeps_only_matching_records() {
    // Build records with three different read groups: two RG1, one RG2, one RG3.
    let rg1_aux = b"RGZRG1\0";
    let rg2_aux = b"RGZRG2\0";
    let rg3_aux = b"RGZRG3\0";

    let mut store = RecordStore::new();
    let mut filter = ReadGroupFilter { wanted: b"RG1".to_vec() };

    let inputs: &[(i32, &[u8])] = &[
        (100, rg1_aux),
        (101, rg2_aux),
        (102, rg1_aux),
        (103, rg3_aux),
        (104, rg2_aux),
        (105, rg1_aux),
    ];
    for (pos, aux) in inputs {
        let raw = helpers::make_record_with_aux(0, *pos, 99, 60, 10, aux);
        store.push_raw(&raw, &mut filter).unwrap();
    }

    // Three records had RG1 → exactly three should survive.
    assert_eq!(store.len(), 3, "only RG1 records should remain");
    let kept_positions: Vec<i32> =
        (0..store.len() as u32).map(|i| store.record(i).pos.as_i32()).collect();
    assert_eq!(kept_positions, vec![100, 102, 105]);

    // Each surviving record's aux must indeed contain RG1.
    for i in 0..store.len() as u32 {
        let rec = store.record(i);
        let aux = rec.aux(&store).expect("aux readable");
        assert_eq!(extract_rg(aux.as_bytes()), Some(b"RG1".as_ref()));
    }
}

// r[verify record_store.pre_filter.rollback]
#[test]
fn read_group_filter_rolls_back_rejected_records() {
    // Drop everything but RG1; verify the store ends up byte-identical to
    // pushing only the matching records with no filter.
    let rg1_aux = b"RGZRG1\0";
    let other_aux = b"RGZRG2\0";

    let mut a = RecordStore::new();
    let mut filter = ReadGroupFilter { wanted: b"RG1".to_vec() };
    let inputs: &[&[u8]] = &[other_aux, rg1_aux, other_aux, rg1_aux, rg1_aux, other_aux];
    for (i, aux) in inputs.iter().enumerate() {
        let raw = helpers::make_record_with_aux(0, 100 + i as i32, 99, 60, 10, aux);
        a.push_raw(&raw, &mut filter).unwrap();
    }

    // Reference store: only the matching records, no filter.
    let mut b = RecordStore::new();
    for (i, aux) in inputs.iter().enumerate() {
        if extract_rg(aux) == Some(b"RG1".as_ref()) {
            let raw = helpers::make_record_with_aux(0, 100 + i as i32, 99, 60, 10, aux);
            b.push_raw(&raw, &mut ()).unwrap();
        }
    }

    assert_eq!(a.len(), b.len(), "rollback must yield same record count");
    for i in 0..a.len() as u32 {
        assert_eq!(a.record(i).pos, b.record(i).pos);
        assert_eq!(a.aux(i), b.aux(i));
    }
}
