//! Where readers look for an index: `r[index.locate]` and `r[index.explicit]`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use seqair::bam::RecordStore;
use seqair::reader::IndexedReader;
use seqair_types::Pos0;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data").join(name)
}

/// Copy `from` (in `tests/data`) to `name` in `dir`.
fn place(dir: &Path, from: &str, name: &str) -> PathBuf {
    let to = dir.join(name);
    std::fs::copy(data(from), &to).unwrap();
    to
}

/// Something at an index's name that is not an index.
fn decoy(dir: &Path, name: &str) {
    std::fs::write(dir.join(name), b"not an index").unwrap();
}

/// Records the reader returns over the chr19 reads of the test files.
fn chr19_records(reader: &mut IndexedReader) -> usize {
    let tid = reader.header().tid("chr19").unwrap();
    let span = (Pos0::new(6_100_000).unwrap()..=Pos0::new(6_145_000).unwrap()).into();
    reader.fetch_into(tid, span, &mut RecordStore::default()).unwrap()
}

fn open_bam(path: &Path) -> IndexedReader {
    IndexedReader::open(path).unwrap()
}

// r[verify index.locate]
#[test]
fn a_bam_without_the_bam_extension_finds_its_appended_index() {
    let dir = tempfile::tempdir().unwrap();
    let bam = place(dir.path(), "test.bam", "sample");
    place(dir.path(), "test.bam.bai", "sample.bai");
    assert!(chr19_records(&mut open_bam(&bam)) > 0);
}

// r[verify index.locate]
/// `with_extension("bam.bai")` turned `sample.v2` into `sample.bam.bai`, and
/// `with_extension("bai")` into `sample.bai`: another file's name.
#[test]
fn a_dotted_stem_without_the_bam_extension_finds_its_own_index() {
    let dir = tempfile::tempdir().unwrap();
    let bam = place(dir.path(), "test.bam", "sample.v2");
    place(dir.path(), "test.bam.bai", "sample.v2.bai");
    decoy(dir.path(), "sample.bai");
    decoy(dir.path(), "sample.bam.bai");
    assert!(chr19_records(&mut open_bam(&bam)) > 0);
}

// r[verify index.locate]
#[test]
fn a_dotted_stem_finds_an_index_in_place_of_its_extension() {
    let dir = tempfile::tempdir().unwrap();
    let bam = place(dir.path(), "test.bam", "sample.v2.bam");
    place(dir.path(), "test.bam.bai", "sample.v2.bai");
    assert!(chr19_records(&mut open_bam(&bam)) > 0);
}

// r[verify index.locate]
/// Appended names come first, the CSI before the BAI in each form.
#[test]
fn an_appended_index_wins_over_a_replaced_one() {
    let dir = tempfile::tempdir().unwrap();
    let bam = place(dir.path(), "test.bam", "sample.v2.bam");
    place(dir.path(), "test.bam.bai", "sample.v2.bam.bai");
    decoy(dir.path(), "sample.v2.csi");
    decoy(dir.path(), "sample.v2.bai");
    assert!(chr19_records(&mut open_bam(&bam)) > 0);

    // And the CSI is the first one looked for.
    place(dir.path(), "test.bam.samtools.csi", "sample.v2.bam.csi");
    let reader = open_bam(&bam);
    let IndexedReader::Bam(bam_reader) = &reader else { panic!("a BAM") };
    assert!(matches!(bam_reader.shared().index(), seqair::bam::AlignmentIndex::Csi(_)));
}

// r[verify index.explicit]
#[test]
fn an_explicit_index_wins_over_the_one_beside_the_file() {
    let dir = tempfile::tempdir().unwrap();
    let bam = place(dir.path(), "test.bam", "sample.bam");
    decoy(dir.path(), "sample.bam.csi");
    decoy(dir.path(), "sample.bam.bai");
    assert!(IndexedReader::open(&bam).is_err(), "the decoy beside the BAM is not an index");

    let elsewhere = place(dir.path(), "test.bam.bai", "kept-apart.idx");
    let mut reader =
        IndexedReader::open_with_reference_and_index(&bam, &data("test.fasta.gz"), &elsewhere)
            .unwrap();
    assert!(chr19_records(&mut reader) > 0);

    // A CSI is told from a BAI by content, not by name.
    let csi = place(dir.path(), "test.bam.samtools.csi", "kept-apart.bai");
    let mut reader =
        IndexedReader::open_with_reference_and_index(&bam, &data("test.fasta.gz"), &csi).unwrap();
    assert!(chr19_records(&mut reader) > 0);
}

// r[verify index.locate]
// r[verify index.explicit]
#[test]
fn a_cram_looks_up_its_crai_the_same_way() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = data("test.fasta.gz");
    let cram = place(dir.path(), "test_v30.cram", "sample.v2.cram");
    place(dir.path(), "test_v30.cram.crai", "sample.v2.crai");
    let mut reader = IndexedReader::open_with_reference(&cram, &fasta).unwrap();
    let found = chr19_records(&mut reader);
    assert!(found > 0);

    std::fs::remove_file(dir.path().join("sample.v2.crai")).unwrap();
    decoy(dir.path(), "sample.v2.cram.crai");
    let elsewhere = place(dir.path(), "test_v30.cram.crai", "kept-apart.idx");
    let mut reader =
        IndexedReader::open_with_reference_and_index(&cram, &fasta, &elsewhere).unwrap();
    assert_eq!(chr19_records(&mut reader), found);
}
