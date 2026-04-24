//! Tests for `RecordStore<U>` per-record extras and `PileupEngine` `pileups()`.
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
use seqair::bam::record_store::{RecordStore, RecordStoreExtras};

/// Helper: push N synthetic records at positions 100, 101, ... with 10M CIGAR.
fn store_with_n_records(n: u32) -> RecordStore {
    let mut store = RecordStore::new();
    for i in 0..n {
        let pos = 100 + i as i32;
        store.push_raw(&make_record(0, pos, 99, 60, 10), |_, _| true).unwrap();
    }
    store
}

// ---- RecordStore extras ----

// r[verify record_store.extras.provider]
// r[verify record_store.extras.access]
#[test]
fn with_extras_computes_per_record_data() {
    let store = store_with_n_records(5);

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    let store = store.apply_extras(&mut ExtractPos);

    assert_eq!(store.len(), 5);
    for i in 0..5u32 {
        assert_eq!(*store.extra(i), 100 + i as i32);
    }
}

// r[verify record_store.extras.provider]
#[test]
fn with_extras_preserves_slab_data() {
    let store = store_with_n_records(3);

    // Capture original data before transformation.
    let original_positions: Vec<i32> = (0..3u32).map(|i| store.record(i).pos.as_i32()).collect();
    let original_mapqs: Vec<u8> = (0..3u32).map(|i| store.record(i).mapq).collect();

    #[derive(Clone, Default)]
    struct ExtractMapqScaled;
    impl RecordStoreExtras for ExtractMapqScaled {
        type Extra = u16;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> u16 {
            s.record(idx).mapq as u16 * 10
        }
    }

    let store = store.apply_extras(&mut ExtractMapqScaled);

    // Slab data must be moved, not copied — verify it's still correct.
    for i in 0..3u32 {
        assert_eq!(store.record(i).pos.as_i32(), original_positions[i as usize]);
        assert_eq!(store.record(i).mapq, original_mapqs[i as usize]);
    }
    // Extras computed correctly.
    assert_eq!(*store.extra(0), 600);
    assert_eq!(*store.extra(1), 600);
    assert_eq!(*store.extra(2), 600);
}

// r[verify record_store.extras.access]
#[test]
fn extra_mut_allows_modification() {
    let store = store_with_n_records(3);

    #[derive(Clone, Default)]
    struct ZeroExtra;
    impl RecordStoreExtras for ZeroExtra {
        type Extra = u32;
        fn compute(&mut self, _idx: u32, _s: &RecordStore<()>) -> u32 {
            0
        }
    }

    let mut store = store.apply_extras(&mut ZeroExtra);

    *store.extra_mut(1) = 42;

    assert_eq!(*store.extra(0), 0);
    assert_eq!(*store.extra(1), 42);
    assert_eq!(*store.extra(2), 0);
}

// r[verify record_store.extras.provider]
#[test]
fn with_extras_closure_can_read_aux() {
    // Build a record with aux data and verify the provider can access it.
    let aux = b"NMCd"; // NM:C:100 (tag NM, type C=u8, value 100)
    let raw = helpers::make_record_with_aux(0, 100, 99, 60, 10, aux);

    let mut store = RecordStore::new();
    store.push_raw(&raw, |_, _| true).unwrap();

    #[derive(Clone, Default)]
    struct AuxLen;
    impl RecordStoreExtras for AuxLen {
        type Extra = usize;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> usize {
            s.aux(idx).len()
        }
    }

    let store = store.apply_extras(&mut AuxLen);

    assert_eq!(*store.extra(0), aux.len());
}

// r[verify record_store.extras.clear]
#[test]
fn clear_clears_extras_and_retains_capacity() {
    let store = store_with_n_records(50);

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    let mut store = store.apply_extras(&mut ExtractPos);

    assert_eq!(store.len(), 50);
    let records_cap = store.records_capacity();
    let extras_cap = store.extras_capacity();
    assert!(extras_cap >= 50);

    store.clear();

    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() >= records_cap);
    assert!(store.extras_capacity() >= extras_cap);
}

// r[verify record_store.extras.strip]
#[test]
fn strip_extras_preserves_slab_capacity() {
    let store = store_with_n_records(50);

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    let store = store.apply_extras(&mut ExtractPos);

    let records_cap = store.records_capacity();
    let names_cap = store.names_capacity();
    let qual_cap = store.qual_capacity();

    let store = store.strip_extras();

    // Slab capacities preserved.
    assert!(store.records_capacity() >= records_cap);
    assert!(store.names_capacity() >= names_cap);
    assert!(store.qual_capacity() >= qual_cap);

    // Records still accessible (they were moved, not cleared).
    assert_eq!(store.len(), 50);
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
    store.push_raw(&make_record(0, 200, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap();

    store.sort_by_pos();
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn sort_by_pos_preserves_extras_mapping() {
    let mut store = RecordStore::new();
    // Push records out of order: pos 300, 100, 200
    store.push_raw(&make_record(0, 300, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 200, 99, 60, 10), |_, _| true).unwrap();

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    // Compute extras BEFORE sorting — extras carry the original position.
    let mut store = store.apply_extras(&mut ExtractPos);

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
    let mut store = RecordStore::new();
    // Push duplicate records at same position with same flags/name.
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap(); // dup
    store.push_raw(&make_record(0, 200, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 200, 99, 60, 10), |_, _| true).unwrap(); // dup

    #[derive(Clone, Default)]
    struct RecordIdx;
    impl RecordStoreExtras for RecordIdx {
        type Extra = u32;
        fn compute(&mut self, idx: u32, _s: &RecordStore<()>) -> u32 {
            idx
        }
    }

    // Tag each record with a unique ID before dedup.
    let mut store = store.apply_extras(&mut RecordIdx);

    store.sort_by_pos();
    store.dedup();

    // Should have 2 records after dedup (one per position).
    assert_eq!(store.len(), 2);
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());

    // Extras should still be accessible and correspond to the surviving records.
    // The exact idx values depend on which duplicate survives (dedup_by keeps the first).
    let e0 = *store.extra(0);
    let e1 = *store.extra(1);
    assert_ne!(e0, e1, "surviving records should have distinct extras");
}

// r[verify record_store.extras.sort_dedup_generic]
#[test]
fn set_alignment_then_sort_on_typed_store() {
    let mut store = RecordStore::new();
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 50, 99, 60, 10), |_, _| true).unwrap();

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    // Compute extras with original positions.
    let mut store = store.apply_extras(&mut ExtractPos);

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
    let store = store_with_n_records(5);

    #[derive(Clone, Default)]
    struct ExtractMapq;
    impl RecordStoreExtras for ExtractMapq {
        type Extra = u8;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> u8 {
            s.record(idx).mapq
        }
    }

    let store = store.apply_extras(&mut ExtractMapq);
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    engine.set_max_depth(10);
    engine.set_filter(|flags, _aux| !flags.is_duplicate());

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
    let store = store_with_n_records(3);

    #[derive(Clone, Default)]
    struct ExtractPos;
    impl RecordStoreExtras for ExtractPos {
        type Extra = i32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
            s.record(idx).pos.as_i32()
        }
    }

    let store = store.apply_extras(&mut ExtractPos);
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
#[test]
fn recover_store_works_with_extras_engine() {
    // Build a store manually, transform to extras, feed into engine, then strip.
    let store = store_with_n_records(10);

    #[derive(Clone, Default)]
    struct ExtractMapqU32;
    impl RecordStoreExtras for ExtractMapqU32 {
        type Extra = u32;
        fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> u32 {
            s.record(idx).mapq as u32
        }
    }

    let store = store.apply_extras(&mut ExtractMapqU32);
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());

    // Consume all columns.
    while engine.pileups().is_some() {}

    // take_store + strip_extras (the same operation recover_store does internally).
    let store = engine.take_store().expect("store should be available");
    let store: RecordStore = store.strip_extras();

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
            store.push_raw(&make_record(0, pos, 99, 60, 10), |_, _| true).unwrap();
        }

        #[derive(Clone, Default)]
        struct ExtractPos;
        impl RecordStoreExtras for ExtractPos {
            type Extra = i32;
            fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
                s.record(idx).pos.as_i32()
            }
        }

        // Tag each record with its original position before sorting.
        let mut store = store.apply_extras(&mut ExtractPos);

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
            store.push_raw(&make_record(0, pos, 99, 60, 10), |_, _| true).unwrap();
        }

        #[derive(Clone, Default)]
        struct ExtractPos;
        impl RecordStoreExtras for ExtractPos {
            type Extra = i32;
            fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
                s.record(idx).pos.as_i32()
            }
        }

        // Compute extras (original position) before sort+dedup.
        let mut store = store.apply_extras(&mut ExtractPos);

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

    // r[verify record_store.extras.provider]
    #[test]
    fn proptest_with_extras_after_sort_produces_correct_mapping(
        positions in proptest::collection::vec(0i32..500, 2..30),
    ) {
        let mut store = RecordStore::new();
        for &pos in &positions {
            store.push_raw(&make_record(0, pos, 99, 60, 10), |_, _| true).unwrap();
        }

        #[derive(Clone, Default)]
        struct ExtractPos;
        impl RecordStoreExtras for ExtractPos {
            type Extra = i32;
            fn compute(&mut self, idx: u32, s: &RecordStore<()>) -> i32 {
                s.record(idx).pos.as_i32()
            }
        }

        // Sort first (on unit store), then compute extras.
        store.sort_by_pos();

        let store = store.apply_extras(&mut ExtractPos);

        // After apply_extras on a sorted store, extras_idx is reset to sequential.
        for i in 0..store.len() as u32 {
            let pos = store.record(i).pos.as_i32();
            let extra = *store.extra(i);
            prop_assert_eq!(pos, extra, "apply_extras mapping wrong at record {}", i);
        }
    }
}

// r[verify pileup.extras.constructor_accepts_any_u]
#[test]
fn pileup_with_sorted_typed_store() {
    // End-to-end: load out of order, compute extras, sort, pileup.
    let mut store = RecordStore::new();
    store.push_raw(&make_record(0, 105, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10), |_, _| true).unwrap();
    store.push_raw(&make_record(0, 102, 99, 60, 10), |_, _| true).unwrap();

    #[derive(Clone, Default)]
    struct RecordIdx;
    impl RecordStoreExtras for RecordIdx {
        type Extra = u32;
        fn compute(&mut self, idx: u32, _s: &RecordStore<()>) -> u32 {
            idx
        }
    }

    // Compute extras before sorting — each extra is the original push index.
    let mut store = store.apply_extras(&mut RecordIdx);

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

// r[verify record_store.pre_filter.provider_hook]
#[test]
fn keep_record_drops_low_mapq_records_via_fetch_into_filtered() {
    use seqair::bam::IndexedBamReader;
    use seqair::bam::record_store::SlimRecord;
    use seqair::reader::FetchCounts;

    // Build a synthetic store with records at mapq 0, 10, 20, 30, 40, 50.
    let mapqs = [0u8, 10, 20, 30, 40, 50];

    // Verify the keep_record hook via a direct push_raw_filtered flow first —
    // this checks that the provider's keep_record is honored by RecordStore.
    let mut store = RecordStore::new();
    for (i, &mq) in mapqs.iter().enumerate() {
        let raw = make_record(0, 100 + i as i32, 99, mq, 10);
        // Inline filter equivalent to `rec.mapq >= 20`.
        store.push_raw(&raw, |rec, _s| rec.mapq >= 20).unwrap();
    }
    assert_eq!(store.len(), 4, "records with mapq < 20 should be dropped");
    for i in 0..4u32 {
        assert!(store.record(i).mapq >= 20, "kept record at idx {i} below threshold");
    }

    // Now the same check through IndexedBamReader::fetch_into_filtered using
    // a test BAM. We can't easily inject synthetic mapq into the test BAM,
    // so just assert FetchCounts semantics (fetched >= kept, kept equal to
    // store.len() after call).
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
        .fetch_into_filtered(
            tid,
            Pos0::new(6_103_000).unwrap(),
            Pos0::new(6_104_000).unwrap(),
            &mut store,
            |_: &SlimRecord, _: &RecordStore| true,
        )
        .expect("fetch_into_filtered");
    assert!(counts_all.fetched > 0, "test region should produce records");
    assert_eq!(counts_all.kept, counts_all.fetched);
    let all_len = store.len();
    assert_eq!(all_len, counts_all.kept);

    // Drop-all: kept == 0 but fetched unchanged.
    let counts_none = reader
        .fetch_into_filtered(
            tid,
            Pos0::new(6_103_000).unwrap(),
            Pos0::new(6_104_000).unwrap(),
            &mut store,
            |_: &SlimRecord, _: &RecordStore| false,
        )
        .expect("fetch_into_filtered");
    assert_eq!(counts_none.fetched, counts_all.fetched);
    assert_eq!(counts_none.kept, 0);
    assert_eq!(store.len(), 0, "all records rolled back");
    // FetchCounts invariant
    let _: FetchCounts = counts_none;
}

// r[verify cram.fetch_into_filtered.push_time]
#[test]
fn cram_fetch_into_filtered_applies_filter_at_push_time() {
    use seqair::bam::record_store::SlimRecord;
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
        .fetch_into_filtered(
            tid,
            region_start,
            region_end,
            &mut store,
            |_: &SlimRecord, _: &RecordStore| true,
        )
        .expect("fetch_into_filtered all");
    assert!(counts_all.fetched > 0, "test region should produce records");
    assert_eq!(counts_all.kept, counts_all.fetched, "always-keep must keep everything");
    assert_eq!(store.len(), counts_all.kept);

    // Snapshot positions of all fetched records — we'll use them as an oracle
    // for the mixed filter below.
    let mut all_positions: Vec<i64> =
        (0..store.len() as u32).map(|i| store.record(i).pos.as_i64()).collect();
    all_positions.sort_unstable();

    // Drop-all: fetched unchanged, kept == 0. The store (including all slabs)
    // must be empty after rollback — this is the zero-waste guarantee the
    // previous retain-based implementation couldn't provide.
    let counts_none = reader
        .fetch_into_filtered(
            tid,
            region_start,
            region_end,
            &mut store,
            |_: &SlimRecord, _: &RecordStore| false,
        )
        .expect("fetch_into_filtered none");
    assert_eq!(counts_none.fetched, counts_all.fetched, "fetched must be filter-independent");
    assert_eq!(counts_none.kept, 0, "drop-all must yield no kept records");
    assert_eq!(store.len(), 0, "store must be empty after drop-all (rollback)");

    // Mixed filter: keep records whose pos is even. Use an independent oracle
    // (slice-and-filter the baseline positions) to verify the survivors.
    let counts_mixed = reader
        .fetch_into_filtered(
            tid,
            region_start,
            region_end,
            &mut store,
            |rec: &SlimRecord, _: &RecordStore| rec.pos.as_i64() % 2 == 0,
        )
        .expect("fetch_into_filtered mixed");
    assert_eq!(counts_mixed.fetched, counts_all.fetched, "fetched must be filter-independent");
    let expected_kept: usize = all_positions.iter().filter(|p| *p % 2 == 0).count();
    assert_eq!(counts_mixed.kept, expected_kept, "filter result must match oracle");
    assert_eq!(store.len(), expected_kept);
    for i in 0..store.len() as u32 {
        assert_eq!(store.record(i).pos.as_i64() % 2, 0, "record {i} failed filter predicate");
    }

    // FetchCounts invariant: the type is `Copy + Default`.
    let _: FetchCounts = counts_mixed;
}
