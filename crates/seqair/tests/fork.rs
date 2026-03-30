//! Tests for `IndexedBamReader::fork()` — the thread-safe shared-state pattern.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use seqair::bam::{IndexedBamReader, Pos, RecordStore, Zero};
use std::{path::Path, sync::Arc, thread};

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

fn fetch_record_positions(
    reader: &mut IndexedBamReader,
    contig: &str,
    start: u64,
    end: u64,
) -> Vec<(i64, i64)> {
    let mut store = RecordStore::new();
    let tid = reader.header().tid(contig).expect("tid lookup");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(start as u32).unwrap(),
            Pos::<Zero>::new(end as u32).unwrap(),
            &mut store,
        )
        .expect("fetch_into");
    (0..store.len() as u32)
        .map(|i| (store.record(i).pos.as_i64(), store.record(i).end_pos.as_i64()))
        .collect()
}

// r[verify bam.reader.fork]
// r[verify bam.reader.fork_equivalence]
#[test]
fn fork_produces_identical_results_to_independent_open() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open prototype");
    let mut forked = prototype.fork().expect("fork");
    let mut independent = IndexedBamReader::open(test_bam_path()).expect("open independent");

    for &(contig, start, end) in CONTIGS {
        let forked_records = fetch_record_positions(&mut forked, contig, start, end);
        let independent_records = fetch_record_positions(&mut independent, contig, start, end);

        assert_eq!(
            forked_records, independent_records,
            "{contig}: forked and independent readers produced different results"
        );
        assert!(!forked_records.is_empty(), "{contig}: expected records but got none");
    }
}

// r[verify bam.reader.fork_equivalence]
#[test]
fn fork_produces_identical_results_across_multiple_fetches() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open");
    let mut forked = prototype.fork().expect("fork");

    let mut first_pass = Vec::new();
    for &(contig, start, end) in CONTIGS {
        first_pass.push(fetch_record_positions(&mut forked, contig, start, end));
    }

    let mut second_pass = Vec::new();
    for &(contig, start, end) in CONTIGS {
        second_pass.push(fetch_record_positions(&mut forked, contig, start, end));
    }

    assert_eq!(first_pass, second_pass, "repeated fetches on same fork diverged");
}

// r[verify bam.reader.fork]
// r[verify bam.reader.fork_arc_identity]
#[test]
fn fork_shares_arc_with_prototype() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open");
    let forked = prototype.fork().expect("fork");

    assert!(
        Arc::ptr_eq(prototype.shared(), forked.shared()),
        "forked reader should share the same Arc allocation as the prototype"
    );
}

// r[verify bam.reader.fork_arc_identity]
#[test]
fn multiple_forks_share_same_arc() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open");
    let fork1 = prototype.fork().expect("fork1");
    let fork2 = prototype.fork().expect("fork2");
    let fork3 = fork1.fork().expect("fork from fork");

    assert!(Arc::ptr_eq(prototype.shared(), fork1.shared()));
    assert!(Arc::ptr_eq(prototype.shared(), fork2.shared()));
    assert!(Arc::ptr_eq(prototype.shared(), fork3.shared()));
}

// r[verify bam.reader.shared_state]
#[test]
fn fork_header_is_identical() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open");
    let forked = prototype.fork().expect("fork");

    assert_eq!(prototype.header().target_count(), forked.header().target_count());
    for tid in 0..prototype.header().target_count() as u32 {
        assert_eq!(prototype.header().target_name(tid), forked.header().target_name(tid),);
        assert_eq!(prototype.header().target_len(tid), forked.header().target_len(tid),);
    }
}

// r[verify bam.reader.fork_independence]
#[test]
fn fork_fetches_are_independent() {
    let prototype = IndexedBamReader::open(test_bam_path()).expect("open");
    let mut fork_a = prototype.fork().expect("fork_a");
    let mut fork_b = prototype.fork().expect("fork_b");

    let &(contig, start, end) = &CONTIGS[0];
    let tid = fork_a.header().tid(contig).expect("tid");

    let mut store_a = RecordStore::new();
    fork_a
        .fetch_into(
            tid,
            Pos::<Zero>::new(start as u32).unwrap(),
            Pos::<Zero>::new(end as u32).unwrap(),
            &mut store_a,
        )
        .expect("fetch_a");

    // Fetch a different region on fork_b, then the same region
    let &(contig2, start2, end2) = &CONTIGS[1];
    let tid2 = fork_b.header().tid(contig2).expect("tid2");
    let mut store_b = RecordStore::new();
    fork_b
        .fetch_into(
            tid2,
            Pos::<Zero>::new(start2 as u32).unwrap(),
            Pos::<Zero>::new(end2 as u32).unwrap(),
            &mut store_b,
        )
        .expect("fetch_b different region");

    // Now fetch the original region on fork_b
    fork_b
        .fetch_into(
            tid,
            Pos::<Zero>::new(start as u32).unwrap(),
            Pos::<Zero>::new(end as u32).unwrap(),
            &mut store_b,
        )
        .expect("fetch_b same region");

    assert_eq!(store_a.len(), store_b.len(), "independent forks should produce same count");

    for i in 0..store_a.len() as u32 {
        assert_eq!(store_a.record(i).pos, store_b.record(i).pos);
        assert_eq!(store_a.record(i).end_pos, store_b.record(i).end_pos);
        assert_eq!(store_a.record(i).flags, store_b.record(i).flags);
    }
}

// r[verify bam.reader.fork_concurrent]
#[test]
fn forks_work_concurrently_across_threads() {
    let prototype = Arc::new(IndexedBamReader::open(test_bam_path()).expect("open"));

    let handles: Vec<_> = CONTIGS
        .iter()
        .map(|&(contig, start, end)| {
            let proto = Arc::clone(&prototype);
            let contig = contig.to_string();
            thread::spawn(move || {
                let mut forked = proto.fork().expect("fork in thread");
                fetch_record_positions(&mut forked, &contig, start, end)
            })
        })
        .collect();

    // Collect results and compare against sequential reads
    let thread_results: Vec<_> =
        handles.into_iter().map(|h| h.join().expect("thread panicked")).collect();

    let mut sequential = IndexedBamReader::open(test_bam_path()).expect("open sequential");
    for (i, &(contig, start, end)) in CONTIGS.iter().enumerate() {
        let sequential_records = fetch_record_positions(&mut sequential, contig, start, end);
        assert_eq!(
            thread_results[i], sequential_records,
            "{contig}: concurrent fork produced different results than sequential"
        );
    }
}

// r[verify bam.reader.fork_concurrent]
#[test]
fn many_forks_same_region_concurrent() {
    let prototype = Arc::new(IndexedBamReader::open(test_bam_path()).expect("open"));
    let &(contig, start, end) = &CONTIGS[0];

    let handles: Vec<_> = (0..8)
        .map(|_| {
            let proto = Arc::clone(&prototype);
            let contig = contig.to_string();
            thread::spawn(move || {
                let mut forked = proto.fork().expect("fork");
                fetch_record_positions(&mut forked, &contig, start, end)
            })
        })
        .collect();

    let results: Vec<_> = handles.into_iter().map(|h| h.join().expect("join")).collect();

    // All threads should produce identical results
    for (i, result) in results.iter().enumerate().skip(1) {
        assert_eq!(&results[0], result, "thread 0 vs thread {i} diverged");
    }
}
