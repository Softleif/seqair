//! CRAM version matrix tests: generate CRAM v3.0 and v3.1 from htslib test
//! SAMs via samtools, then read with seqair and compare against noodles.
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

use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record::Sequence as _;
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::Readers;
use std::path::{Path, PathBuf};
use std::process::Command;

fn htslib_sam(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/sam/")).join(name)
}

fn htslib_fasta(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/fasta/")).join(name)
}

/// Convert a SAM to CRAM with a given version, create CRAI index.
/// Uses a two-step process (SAM -> sorted BAM -> CRAM) so that samtools
/// can compute M5 checksums from the reference FASTA.
fn sam_to_cram(
    dir: &Path,
    sam_path: &Path,
    fasta_path: &Path,
    version: &str,
    extra_opts: &[&str],
) -> PathBuf {
    let bam_path = dir.join("sorted.bam");
    let cram_path = dir.join(format!("test_{version}.cram"));

    // Step 1: sort SAM -> BAM
    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools sort failed for {}", sam_path.display());

    // Step 2: BAM -> CRAM with proper reference (gives M5 tags)
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "-C"]).arg("--output-fmt-option").arg(format!("version={version}"));
    for opt in extra_opts {
        cmd.arg("--output-fmt-option").arg(opt);
    }
    cmd.arg("-T").arg(fasta_path).arg("-o").arg(&cram_path).arg(&bam_path);

    let status = cmd.status().expect("samtools view -C failed");
    assert!(status.success(), "samtools view -C failed");

    // Step 3: index
    let status = Command::new("samtools")
        .arg("index")
        .arg(&cram_path)
        .status()
        .expect("samtools index failed");
    assert!(status.success(), "samtools index failed");

    cram_path
}

struct NoodlesRecord {
    pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq: Vec<u8>,
    seq_len: usize,
}

/// Read all mapped records from a CRAM with noodles, grouped by reference index.
fn read_noodles_cram(
    cram_path: &Path,
    fasta_path: &Path,
) -> (sam::Header, Vec<Vec<NoodlesRecord>>) {
    let fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .expect("open indexed FASTA");
    let repository =
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta_reader));

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .expect("open CRAM");

    let header = reader.read_header().expect("read header");
    let n_contigs = header.reference_sequences().len();
    let mut grouped: Vec<Vec<NoodlesRecord>> = (0..n_contigs).map(|_| Vec::new()).collect();

    for result in reader.records(&header) {
        let record = result.expect("noodles: read record");
        let flags = record.flags();

        if flags.is_unmapped() {
            continue;
        }

        let Some(ref_seq_id) = record.reference_sequence_id() else {
            continue;
        };

        let pos = record.alignment_start().map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);

        let mapq = record.mapping_quality().map(u8::from).unwrap_or(255);
        let qname: Vec<u8> = record.name().map(|n| n.to_vec()).unwrap_or_default();
        let seq: Vec<u8> = record.sequence().iter().collect();
        let seq_len = seq.len();

        grouped[ref_seq_id].push(NoodlesRecord {
            pos,
            flags: u16::from(flags),
            mapq,
            qname,
            seq,
            seq_len,
        });
    }

    (header, grouped)
}

/// Core comparison: read CRAM with seqair and noodles, compare record fields.
fn assert_cram_parity(cram_path: &Path, fasta_path: &Path, label: &str) {
    let (header, noodles_grouped) = read_noodles_cram(cram_path, fasta_path);
    let mut readers = Readers::open(cram_path, fasta_path).expect("seqair: open CRAM+FASTA");

    let mut total = 0usize;

    for (contig_idx, noodles_records) in noodles_grouped.iter().enumerate() {
        if noodles_records.is_empty() {
            continue;
        }

        let (contig_name_bytes, contig_map) =
            header.reference_sequences().get_index(contig_idx).unwrap();
        let contig_name = std::str::from_utf8(contig_name_bytes).unwrap();
        let contig_len = usize::from(contig_map.length()) as u32;

        let tid = readers
            .header()
            .tid(contig_name)
            .unwrap_or_else(|| panic!("{label}: seqair missing contig '{contig_name}'"));

        let mut store = RecordStore::new();
        readers
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(contig_len).unwrap(), &mut store)
            .unwrap_or_else(|e| panic!("{label}/{contig_name}: fetch failed: {e:?}"));

        assert_eq!(
            store.len(),
            noodles_records.len(),
            "{label}/{contig_name}: record count seqair={} noodles={}",
            store.len(),
            noodles_records.len()
        );

        for (i, n) in noodles_records.iter().enumerate() {
            let idx = i as u32;
            let r = store.record(idx);
            let ctx = format!("{label}/{contig_name}[{i}]");

            assert_eq!(r.pos.as_i64(), n.pos, "{ctx}: pos");
            assert_eq!(r.flags.raw(), n.flags, "{ctx}: flags");
            assert_eq!(r.mapq, n.mapq, "{ctx}: mapq");
            assert_eq!(store.qname(idx), n.qname.as_slice(), "{ctx}: qname");

            // CRAM may reconstruct the sequence from the reference for
            // secondary alignments that had SEQ=* in the original SAM.
            // Only compare seq when noodles also reports a sequence.
            if n.seq_len > 0 {
                assert_eq!(r.seq_len as usize, n.seq_len, "{ctx}: seq_len");

                for pos in 0..n.seq_len {
                    let seqair_base = store.seq_at(idx, pos) as u8;
                    let noodles_base = n.seq[pos];
                    match noodles_base {
                        b'A' | b'C' | b'G' | b'T' => {
                            assert_eq!(seqair_base, noodles_base, "{ctx} seq[{pos}]");
                        }
                        _ => {
                            assert_eq!(seqair_base, b'N', "{ctx} seq[{pos}]: expected N");
                        }
                    }
                }
            }
        }

        total += noodles_records.len();
    }

    assert!(total > 0, "{label}: no records compared");
}

// --- CRAM v3.0 tests ---

/// ce#5.sam round-trip through CRAM v3.0.
#[test]
fn cram_v30_ce5() {
    let dir = tempfile::tempdir().unwrap();
    let cram = sam_to_cram(dir.path(), &htslib_sam("ce#5.sam"), &htslib_fasta("ce.fa"), "3.0", &[]);
    assert_cram_parity(&cram, &htslib_fasta("ce.fa"), "v3.0/ce#5");
}

/// c1#clip.sam (soft/hard clips, introns) through CRAM v3.0.
#[test]
fn cram_v30_c1_clip() {
    let dir = tempfile::tempdir().unwrap();
    let cram =
        sam_to_cram(dir.path(), &htslib_sam("c1#clip.sam"), &htslib_fasta("c1.fa"), "3.0", &[]);
    assert_cram_parity(&cram, &htslib_fasta("c1.fa"), "v3.0/c1#clip");
}

/// xx#pair.sam (paired-end reads) through CRAM v3.0.
#[test]
fn cram_v30_xx_pair() {
    let dir = tempfile::tempdir().unwrap();
    let cram =
        sam_to_cram(dir.path(), &htslib_sam("xx#pair.sam"), &htslib_fasta("xx.fa"), "3.0", &[]);
    assert_cram_parity(&cram, &htslib_fasta("xx.fa"), "v3.0/xx#pair");
}

// --- CRAM v3.1 tests ---

/// ce#5.sam round-trip through CRAM v3.1.
#[test]
fn cram_v31_ce5() {
    let dir = tempfile::tempdir().unwrap();
    let cram = sam_to_cram(dir.path(), &htslib_sam("ce#5.sam"), &htslib_fasta("ce.fa"), "3.1", &[]);
    assert_cram_parity(&cram, &htslib_fasta("ce.fa"), "v3.1/ce#5");
}

/// c1#clip.sam through CRAM v3.1.
#[test]
fn cram_v31_c1_clip() {
    let dir = tempfile::tempdir().unwrap();
    let cram =
        sam_to_cram(dir.path(), &htslib_sam("c1#clip.sam"), &htslib_fasta("c1.fa"), "3.1", &[]);
    assert_cram_parity(&cram, &htslib_fasta("c1.fa"), "v3.1/c1#clip");
}

/// xx#pair.sam through CRAM v3.1.
#[test]
fn cram_v31_xx_pair() {
    let dir = tempfile::tempdir().unwrap();
    let cram =
        sam_to_cram(dir.path(), &htslib_sam("xx#pair.sam"), &htslib_fasta("xx.fa"), "3.1", &[]);
    assert_cram_parity(&cram, &htslib_fasta("xx.fa"), "v3.1/xx#pair");
}

/// ce#5b.sam (multi-contig + secondary with SEQ=*) through CRAM v3.1.
#[test]
fn cram_v31_ce5b() {
    let dir = tempfile::tempdir().unwrap();
    let cram =
        sam_to_cram(dir.path(), &htslib_sam("ce#5b.sam"), &htslib_fasta("ce.fa"), "3.1", &[]);
    assert_cram_parity(&cram, &htslib_fasta("ce.fa"), "v3.1/ce#5b");
}

// --- Multi-ref CRAM containers ---

/// ce#5b.sam has records across 5 contigs. With `multi_seq_per_slice=1`,
/// samtools puts reads from different references into the same slice
/// (`ref_seq_id` == -2). seqair must handle this.
#[test]
fn cram_multi_ref_container() {
    let dir = tempfile::tempdir().unwrap();
    let cram = sam_to_cram(
        dir.path(),
        &htslib_sam("ce#5b.sam"),
        &htslib_fasta("ce.fa"),
        "3.0",
        &["multi_seq_per_slice=1"],
    );
    assert_cram_parity(&cram, &htslib_fasta("ce.fa"), "multi-ref/ce#5b");
}

/// ce#1000.sam has many records on `CHROMOSOME_I` — `multi_seq_per_slice`
/// shouldn't change behavior for single-contig data, but exercises the
/// container routing.
#[test]
fn cram_multi_ref_single_contig() {
    let dir = tempfile::tempdir().unwrap();
    let cram = sam_to_cram(
        dir.path(),
        &htslib_sam("ce#1000.sam"),
        &htslib_fasta("ce.fa"),
        "3.0",
        &["multi_seq_per_slice=1"],
    );
    assert_cram_parity(&cram, &htslib_fasta("ce.fa"), "multi-ref/ce#1000");
}

// --- Embedded reference ---

/// CRAM with embedded reference: the CRAM file contains the reference
/// sequence, so it can be decoded without an external FASTA file.
/// Verified against samtools (noodles has a bug reading `embed_ref=2` CRAMs).
#[test]
fn cram_embedded_reference() {
    let dir = tempfile::tempdir().unwrap();
    let cram = sam_to_cram(
        dir.path(),
        &htslib_sam("c1#clip.sam"),
        &htslib_fasta("c1.fa"),
        "3.0",
        &["embed_ref=2"],
    );

    // samtools can read the embedded-ref CRAM without -T
    let output =
        Command::new("samtools").args(["view", "-c"]).arg(&cram).output().expect("samtools view");
    assert!(output.status.success(), "samtools should read embedded-ref CRAM");
    let samtools_count: usize = String::from_utf8(output.stdout).unwrap().trim().parse().unwrap();
    assert!(samtools_count > 0, "embedded-ref CRAM should have records");

    // seqair reads the same CRAM with the external FASTA
    let mut readers = Readers::open(&cram, &htslib_fasta("c1.fa")).expect("seqair: open");
    let tid = readers.header().tid("c1").expect("tid");
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(
        store.len(),
        samtools_count,
        "embed-ref: seqair={} samtools={samtools_count}",
        store.len()
    );
}
