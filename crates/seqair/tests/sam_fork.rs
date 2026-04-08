//! Thread-safety tests for SAM reader forking.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
use seqair::bam::{Pos, RecordStore, Zero};
use seqair::sam::reader::IndexedSamReader;
use std::path::Path;
use std::process::Command;
use std::sync::Arc;
use std::thread;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

fn create_sam_gz(dir: &Path) -> std::path::PathBuf {
    let sam_gz = dir.join("test.sam.gz");
    let status = Command::new("samtools")
        .args(["view", "-h", "--output-fmt", "SAM,level=6", "-o"])
        .arg(&sam_gz)
        .arg(test_bam_path())
        .status()
        .expect("samtools not found");
    assert!(status.success());
    let status =
        Command::new("tabix").args(["-p", "sam"]).arg(&sam_gz).status().expect("tabix not found");
    assert!(status.success());
    sam_gz
}

fn fetch_record_positions(
    reader: &mut IndexedSamReader,
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

// r[verify sam.reader.fetch_into]
#[test]
fn fork_produces_identical_results() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let prototype = IndexedSamReader::open(&sam_gz).expect("open");
    let mut forked = prototype.fork().expect("fork");
    let mut independent = IndexedSamReader::open(&sam_gz).expect("open2");

    for &(contig, start, end) in CONTIGS {
        let forked_recs = fetch_record_positions(&mut forked, contig, start, end);
        let independent_recs = fetch_record_positions(&mut independent, contig, start, end);

        assert_eq!(
            forked_recs, independent_recs,
            "{contig}: forked and independent SAM readers differ"
        );
        assert!(!forked_recs.is_empty(), "{contig}: no records");
    }
}

#[test]
fn fork_shares_arc() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let prototype = IndexedSamReader::open(&sam_gz).expect("open");
    let forked = prototype.fork().expect("fork");

    assert!(Arc::ptr_eq(prototype.shared(), forked.shared()));
}

#[test]
fn concurrent_forks_across_threads() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let prototype = Arc::new(IndexedSamReader::open(&sam_gz).expect("open"));

    let handles: Vec<_> = CONTIGS
        .iter()
        .map(|&(contig, start, end)| {
            let proto = Arc::clone(&prototype);
            let contig = contig.to_string();
            thread::spawn(move || {
                let mut forked = proto.fork().expect("fork");
                fetch_record_positions(&mut forked, &contig, start, end)
            })
        })
        .collect();

    let thread_results: Vec<_> =
        handles.into_iter().map(|h| h.join().expect("thread panicked")).collect();

    let mut sequential = IndexedSamReader::open(&sam_gz).expect("open seq");
    for (i, &(contig, start, end)) in CONTIGS.iter().enumerate() {
        let seq_recs = fetch_record_positions(&mut sequential, contig, start, end);
        assert_eq!(
            thread_results[i], seq_recs,
            "{contig}: concurrent fork differs from sequential"
        );
    }
}

#[test]
fn many_forks_same_region() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let prototype = Arc::new(IndexedSamReader::open(&sam_gz).expect("open"));
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

    for (i, result) in results.iter().enumerate().skip(1) {
        assert_eq!(&results[0], result, "thread 0 vs thread {i} diverged");
    }
}
