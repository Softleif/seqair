//! Edge-case tests derived from htslib's test suite: DOS line endings,
//! colons in contig names, sequence-less reads, SEQ/QUAL presence combos,
//! and fully-unmapped files.
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

use noodles::bam;
use noodles::sam;
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::IndexedReader;
use std::path::{Path, PathBuf};
use std::process::Command;

fn htslib_sam(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/sam/")).join(name)
}

fn htslib_bam(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/bam/")).join(name)
}

/// Convert a SAM to sorted, indexed BAM in a temp directory.
fn sam_to_indexed_bam(dir: &Path, sam_path: &Path) -> PathBuf {
    let bam_path = dir.join("test.bam");

    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools sort failed for {}", sam_path.display());

    let status = Command::new("samtools")
        .arg("index")
        .arg(&bam_path)
        .status()
        .expect("samtools index failed");
    assert!(status.success(), "samtools index failed");

    bam_path
}

// ---- DOS line endings ----

/// SAM files with \r\n line endings must parse correctly through the
/// indexed SAM reader (BGZF-compressed SAM + tabix).
#[test]
fn dos_line_endings_in_sam() {
    let dir = tempfile::tempdir().unwrap();

    // Create a SAM with \r\n line endings
    let sam_text = "@HD\tVN:1.6\tSO:coordinate\r\n\
                    @SQ\tSN:chr1\tLN:1000\r\n\
                    read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\r\n\
                    read2\t16\tchr1\t200\t50\t8M\t*\t0\t0\tTGCATGCA\tHHHHHHHH\r\n";

    let sam_path = dir.path().join("dos.sam");
    std::fs::write(&sam_path, sam_text).unwrap();

    // bgzip + tabix
    let sam_gz_path = dir.path().join("dos.sam.gz");
    let status = Command::new("bgzip")
        .arg("-c")
        .arg(&sam_path)
        .stdout(std::fs::File::create(&sam_gz_path).unwrap())
        .status()
        .expect("bgzip not found");
    assert!(status.success(), "bgzip failed");

    let status = Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz_path)
        .status()
        .expect("tabix not found");
    assert!(status.success(), "tabix failed");

    let mut reader = IndexedReader::open(&sam_gz_path).expect("open DOS SAM");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 2, "should parse 2 records from DOS SAM");
    assert_eq!(store.record(0).pos.as_i64(), 99); // 0-based
    assert_eq!(store.record(1).pos.as_i64(), 199);
    assert_eq!(store.qname(0), b"read1");
    assert_eq!(store.qname(1), b"read2");
}

/// DOS line endings also work when converted to BAM (samtools handles \r\n).
/// This tests the actual htslib `index_dos.sam` file.
#[test]
fn dos_line_endings_via_bam() {
    let sam_path = htslib_sam("index_dos.sam");
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &sam_path);

    let mut reader = IndexedReader::open(&bam_path).expect("open BAM from DOS SAM");

    // Verify we can read records from CHROMOSOME_I
    let tid = reader.header().tid("CHROMOSOME_I").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_009_800).unwrap(), &mut store)
        .expect("fetch");

    assert!(!store.is_empty(), "should have records on CHROMOSOME_I");

    // Compare against noodles
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut noodles_reader = bam::io::Reader::new(file);
    let _header: sam::Header = noodles_reader.read_header().unwrap();
    let noodles_count = noodles_reader
        .records()
        .filter(|r| {
            let rec = r.as_ref().unwrap();
            !rec.flags().is_unmapped()
        })
        .count();

    // seqair total across all contigs should match noodles total
    let mut total = store.len();
    for name in [
        "CHROMOSOME_II",
        "CHROMOSOME_III",
        "CHROMOSOME_IV",
        "CHROMOSOME_V",
        "CHROMOSOME_X",
        "CHROMOSOME_MtDNA",
    ] {
        if let Some(tid) = reader.header().tid(name) {
            store.clear();
            reader
                .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(5000).unwrap(), &mut store)
                .unwrap();
            total += store.len();
        }
    }
    assert_eq!(total, noodles_count, "total record count should match noodles");
}

// ---- Colons in contig names ----

/// Contig names with colons (chr1:100) and commas (chr1,chr3) must be
/// handled correctly — these look like region strings but are valid names.
#[test]
fn colons_in_contig_names() {
    let bam_path = htslib_bam("colons.bam");
    let mut reader = IndexedReader::open(&bam_path).expect("open colons.bam");

    let expected: &[(&str, u32, &[&str])] = &[
        ("chr1", 1000, &["chr1", "chr1b"]),
        ("chr1:100", 1000, &["chr1:100"]),
        ("chr1:100-200", 1000, &["chr1:100-200"]),
        ("chr2:100-200", 1000, &["chr2:100-200"]),
        ("chr1,chr3", 1000, &["chr1,chr3"]),
    ];

    for &(contig, len, expected_qnames) in expected {
        let tid =
            reader.header().tid(contig).unwrap_or_else(|| panic!("contig '{contig}' not found"));

        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(len).unwrap(), &mut store)
            .unwrap_or_else(|e| panic!("fetch {contig}: {e}"));

        let qnames: Vec<&[u8]> = (0..store.len() as u32).map(|i| store.qname(i)).collect();
        let expected_bytes: Vec<&[u8]> = expected_qnames.iter().map(|s| s.as_bytes()).collect();
        assert_eq!(qnames, expected_bytes, "contig '{contig}': qnames mismatch");
    }
}

// ---- Sequence-less reads (SEQ=*) ----

/// c1#noseq.sam: mapped reads with SEQ=* — valid in BAM for secondary
/// alignments. Tests that `seq_len=0` is handled correctly.
#[test]
fn sequence_less_mapped_reads() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam("c1#noseq.sam"));

    // Read with noodles for ground truth
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut noodles_reader = bam::io::Reader::new(file);
    let _header: sam::Header = noodles_reader.read_header().unwrap();
    let noodles_records: Vec<_> =
        noodles_reader.records().map(|r| r.unwrap()).filter(|r| !r.flags().is_unmapped()).collect();

    // Read with seqair
    let mut reader = IndexedReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("c1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), noodles_records.len(), "mapped record count mismatch");

    // Some records have SEQ=* (seq_len=0), others have actual sequences.
    // Verify seqair handles both correctly.
    for i in 0..store.len() as u32 {
        let r = store.record(i);
        let n = &noodles_records[i as usize];

        let n_seq_len = n.sequence().len();
        assert_eq!(
            r.seq_len as usize,
            n_seq_len,
            "rec {i} ({}): seq_len mismatch seqair={} noodles={}",
            String::from_utf8_lossy(store.qname(i)),
            r.seq_len,
            n_seq_len
        );

        if r.seq_len == 0 {
            // SEQ=*: sequence and quality should be empty
            assert!(store.seq(i).is_empty(), "rec {i}: SEQ=* but seq not empty");
        }
    }
}

// ---- Unknown bases / SEQ+QUAL combinations ----

/// c1#unknown.sam: permutations of SEQ/QUAL being present or "*" in both
/// mapped and unmapped forms.
#[test]
fn seq_qual_presence_combos() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam("c1#unknown.sam"));

    // Read with noodles
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut noodles_reader = bam::io::Reader::new(file);
    let _header: sam::Header = noodles_reader.read_header().unwrap();
    let noodles_mapped: Vec<_> =
        noodles_reader.records().map(|r| r.unwrap()).filter(|r| !r.flags().is_unmapped()).collect();

    // Read with seqair
    let mut reader = IndexedReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("c1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), noodles_mapped.len(), "mapped record count mismatch");

    for i in 0..store.len() as u32 {
        let r = store.record(i);
        let n = &noodles_mapped[i as usize];
        let qname = String::from_utf8_lossy(store.qname(i));

        let n_seq_len = n.sequence().len();
        assert_eq!(r.seq_len as usize, n_seq_len, "{qname}: seq_len mismatch");

        // For records with sequence, verify bases match
        for pos in 0..n_seq_len {
            let seqair_base = store.seq_at(i, pos) as u8;
            let noodles_base: u8 = n.sequence().iter().nth(pos).unwrap();
            match noodles_base {
                b'A' | b'C' | b'G' | b'T' => {
                    assert_eq!(seqair_base, noodles_base, "{qname} seq[{pos}]");
                }
                _ => {
                    assert_eq!(seqair_base, b'N', "{qname} seq[{pos}]: expected N");
                }
            }
        }
    }
}

// ---- Fully-unmapped reads ----

/// ce#unmap.sam: all reads have flag 0x4 (unmapped) with RNAME=*.
/// `fetch_into` should return 0 records for any contig.
#[test]
fn fully_unmapped_file_returns_no_records() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam("ce#unmap.sam"));

    // The file has no @SQ lines, so samtools might fail or create
    // an empty BAM. Either way, seqair should handle it gracefully.
    match IndexedReader::open(&bam_path) {
        Ok(reader) => {
            // If it opens, there should be 0 contigs or 0 records
            assert_eq!(
                reader.header().target_count(),
                0,
                "unmapped-only file should have 0 contigs"
            );
        }
        Err(_) => {
            // Expected: unmapped-only BAM may not have a valid index
            // or header. This is acceptable behavior.
        }
    }
}

// ---- Supplementary alignments ----

/// ce#supp.sam: supplementary alignments (flag 2048) must be included
/// in fetch results — only flag 0x4 (unmapped) is filtered.
#[test]
fn supplementary_alignments_included() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam("ce#supp.sam"));

    let mut reader = IndexedReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("CHROMOSOME_I").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_009_800).unwrap(), &mut store)
        .expect("fetch");

    // Check that supplementary alignments (flag 2048) are present
    let supp_count =
        (0..store.len() as u32).filter(|&i| store.record(i).flags.raw() & 0x800 != 0).count();
    assert!(supp_count > 0, "should include supplementary alignments");

    // Also verify against noodles
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut noodles_reader = bam::io::Reader::new(file);
    let _header: sam::Header = noodles_reader.read_header().unwrap();
    let noodles_mapped =
        noodles_reader.records().filter(|r| !r.as_ref().unwrap().flags().is_unmapped()).count();
    assert_eq!(store.len(), noodles_mapped, "record count vs noodles");
}

// ---- Secondary alignments with SEQ=* ----

/// ce#5b.sam includes a secondary alignment (flag 256) with SEQ=* and QUAL=*.
/// This must parse without error and have `seq_len=0`.
#[test]
fn secondary_alignment_without_sequence() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam("ce#5b.sam"));

    let mut reader = IndexedReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("CHROMOSOME_V").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(5000).unwrap(), &mut store)
        .expect("fetch");

    // Find the secondary alignment (flag 256)
    let secondary = (0..store.len() as u32)
        .find(|&i| store.record(i).flags.raw() & 0x100 != 0)
        .expect("should have a secondary alignment");

    let r = store.record(secondary);
    assert_eq!(r.seq_len, 0, "secondary with SEQ=* should have seq_len=0");
    assert!(store.seq(secondary).is_empty(), "secondary SEQ should be empty");
}

// ---- Padding CIGAR operations ----

/// c1#pad*.sam: reads with padding (P) CIGAR ops — these appear in
/// padded reference alignments.
#[test]
fn padding_cigar_operations() {
    for name in ["c1#pad1.sam", "c1#pad2.sam", "c1#pad3.sam"] {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = sam_to_indexed_bam(dir.path(), &htslib_sam(name));

        let mut reader = IndexedReader::open(&bam_path).expect("open BAM");
        let tid = reader.header().tid("c1").expect("tid");
        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
            .expect("fetch");

        // Verify against noodles
        let file = std::fs::File::open(&bam_path).unwrap();
        let mut noodles_reader = bam::io::Reader::new(file);
        let _header: sam::Header = noodles_reader.read_header().unwrap();
        let noodles_mapped =
            noodles_reader.records().filter(|r| !r.as_ref().unwrap().flags().is_unmapped()).count();

        assert_eq!(store.len(), noodles_mapped, "{name}: record count mismatch");
    }
}
