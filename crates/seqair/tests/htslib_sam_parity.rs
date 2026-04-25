//! Parity tests: parse htslib's own test SAM files (converted to BAM via
//! samtools) and compare every record field against noodles.
//!
//! These tests import the canonical htslib test data to verify that seqair
//! produces identical results on the same inputs htslib is tested against.
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
use seqair_types::BaseQuality;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn htslib_sam(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/sam/")).join(name)
}

/// Convert a SAM to sorted, indexed BAM in a temp directory.
fn sam_to_indexed_bam(dir: &Path, sam_path: &Path) -> PathBuf {
    let bam_path = dir.join("test.bam");

    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found — required for htslib parity tests");
    assert!(status.success(), "samtools sort failed for {}", sam_path.display());

    let status = Command::new("samtools")
        .arg("index")
        .arg(&bam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index failed");
    assert!(status.success(), "samtools index failed");

    bam_path
}

struct NoodlesRecord {
    pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_len: usize,
    cigar_ops: Vec<(u32, u8)>,
}

/// Read all mapped records from a BAM with noodles, grouped by reference index.
fn read_noodles_grouped(bam_path: &Path) -> (sam::Header, Vec<Vec<NoodlesRecord>>) {
    let file = std::fs::File::open(bam_path).expect("open BAM for noodles");
    let mut reader = bam::io::Reader::new(file);
    let header: sam::Header = reader.read_header().expect("noodles: read header");

    let n_contigs = header.reference_sequences().len();
    let mut grouped: Vec<Vec<NoodlesRecord>> = (0..n_contigs).map(|_| Vec::new()).collect();

    for result in reader.records() {
        let record = result.expect("noodles: read record");
        let flags = record.flags();

        if flags.is_unmapped() {
            continue;
        }

        let Some(ref_seq_id) = record.reference_sequence_id().and_then(|r| r.ok()) else {
            continue;
        };

        let aln_start = match record.alignment_start().and_then(|r| r.ok()) {
            Some(pos) => usize::from(pos) as i64 - 1,
            None => continue,
        };

        let mapq = record.mapping_quality().map(u8::from).unwrap_or(255);
        let qname: Vec<u8> = record.name().map(|n| n.to_vec()).unwrap_or_default();
        let seq: Vec<u8> = record.sequence().iter().collect();
        let seq_len = seq.len();
        let qual: Vec<u8> = record.quality_scores().iter().collect();

        let cigar_ops: Vec<(u32, u8)> = record
            .cigar()
            .iter()
            .map(|res| {
                use noodles::sam::alignment::record::cigar::op::Kind;
                let op = res.expect("cigar op");
                let code = match op.kind() {
                    Kind::Match => 0,
                    Kind::Insertion => 1,
                    Kind::Deletion => 2,
                    Kind::Skip => 3,
                    Kind::SoftClip => 4,
                    Kind::HardClip => 5,
                    Kind::Pad => 6,
                    Kind::SequenceMatch => 7,
                    Kind::SequenceMismatch => 8,
                };
                (op.len() as u32, code)
            })
            .collect();

        grouped[ref_seq_id].push(NoodlesRecord {
            pos: aln_start,
            flags: flags.bits(),
            mapq,
            qname,
            seq,
            qual,
            seq_len,
            cigar_ops,
        });
    }

    (header, grouped)
}

/// Decode seqair's typed CIGAR ops into (length, `op_code`) pairs for comparison.
fn decode_cigar(ops: &[seqair::bam::CigarOp]) -> Vec<(u32, u8)> {
    ops.iter().map(|op| (op.len(), op.op_code())).collect()
}

/// Core comparison: convert SAM to BAM, read with both seqair and noodles,
/// assert every record field matches.
fn assert_bam_parity(sam_name: &str) {
    let sam_path = htslib_sam(sam_name);
    assert!(sam_path.exists(), "htslib test file not found: {}", sam_path.display());

    let dir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(dir.path(), &sam_path);

    let (header, noodles_grouped) = read_noodles_grouped(&bam_path);
    let mut reader = IndexedReader::open(&bam_path).expect("seqair: open BAM");

    let mut total_compared = 0usize;

    for (contig_idx, noodles_records) in noodles_grouped.iter().enumerate() {
        if noodles_records.is_empty() {
            continue;
        }

        let (contig_name_bytes, contig_map) =
            header.reference_sequences().get_index(contig_idx).unwrap();
        let contig_name = std::str::from_utf8(contig_name_bytes).unwrap();
        let contig_len = usize::from(contig_map.length()) as u32;

        let tid = reader
            .header()
            .tid(contig_name)
            .unwrap_or_else(|| panic!("{sam_name}: seqair missing contig '{contig_name}'"));

        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(contig_len).unwrap(), &mut store)
            .unwrap_or_else(|e| panic!("{sam_name}/{contig_name}: fetch failed: {e}"));

        assert_eq!(
            store.len(),
            noodles_records.len(),
            "{sam_name}/{contig_name}: record count seqair={} noodles={}",
            store.len(),
            noodles_records.len()
        );

        for (i, n) in noodles_records.iter().enumerate() {
            let idx = i as u32;
            let r = store.record(idx);
            let ctx = format!("{sam_name}/{contig_name}[{i}]");

            assert_eq!(r.pos.as_i64(), n.pos, "{ctx}: pos");
            assert_eq!(r.flags.raw(), n.flags, "{ctx}: flags");
            assert_eq!(r.mapq, n.mapq, "{ctx}: mapq");
            assert_eq!(store.qname(idx), n.qname.as_slice(), "{ctx}: qname");

            // CIGAR
            let seqair_cigar = decode_cigar(store.cigar(idx));
            assert_eq!(seqair_cigar, n.cigar_ops, "{ctx}: cigar");

            // Sequence (seqair normalises non-ACGT to Unknown/'N')
            assert_eq!(r.seq_len as usize, n.seq_len, "{ctx}: seq_len");
            for pos in 0..n.seq_len {
                let seqair_base = store.seq_at(idx, pos) as u8;
                let noodles_base = n.seq[pos];
                match noodles_base {
                    b'A' | b'C' | b'G' | b'T' => {
                        assert_eq!(
                            seqair_base, noodles_base,
                            "{ctx} seq[{pos}]: seqair='{}' noodles='{}'",
                            seqair_base as char, noodles_base as char
                        );
                    }
                    _ => {
                        assert_eq!(
                            seqair_base, b'N',
                            "{ctx} seq[{pos}]: expected 'N' for IUPAC '{}', got '{}'",
                            noodles_base as char, seqair_base as char
                        );
                    }
                }
            }

            // Quality scores — noodles returns empty when all bytes are 0xFF
            // (QUAL=* in SAM); seqair preserves the raw 0xFF bytes.
            let seqair_qual = store.qual(idx);
            if n.qual.is_empty() {
                assert!(
                    seqair_qual.is_empty() || seqair_qual.iter().all(|q| q.get().is_none()),
                    "{ctx}: noodles reports empty qual but seqair has non-unavailable: {seqair_qual:?}"
                );
            } else {
                assert_eq!(
                    BaseQuality::slice_to_bytes(seqair_qual),
                    n.qual.as_slice(),
                    "{ctx}: qual"
                );
            }
        }

        total_compared += noodles_records.len();
    }

    assert!(total_compared > 0, "{sam_name}: no records compared — file may be empty");
}

// --- Individual test functions ---

/// Single C. elegans record with 27M1D73M CIGAR (deletion).
#[test]
fn htslib_ce1() {
    assert_bam_parity("ce#1.sam");
}

/// 5 records on `CHROMOSOME_I` plus one supplementary alignment (flag 2048)
/// with a large 27M100000D73M CIGAR.
#[test]
fn htslib_ce5() {
    assert_bam_parity("ce#5.sam");
}

/// Records spread across 5 contigs; includes a secondary alignment (flag 256)
/// with SEQ=* and QUAL=*, and a complex CIGAR (7S20M1D23M10I30M10S).
#[test]
fn htslib_ce5b() {
    assert_bam_parity("ce#5b.sam");
}

/// Soft clips (S), hard clips (H), introns/ref-skips (N op), and insertions
/// on a 10bp reference.
#[test]
fn htslib_c1_clip() {
    assert_bam_parity("c1#clip.sam");
}

/// Reads at positions 1, 2, 3 on a 10bp contig — alignment extends past
/// reference end.
#[test]
fn htslib_c1_bounds() {
    assert_bam_parity("c1#bounds.sam");
}

/// Paired-end reads: flags 99 (read1, proper pair, mate reverse) and 147
/// (read2, proper pair, reverse).
#[test]
fn htslib_xx_pair() {
    assert_bam_parity("xx#pair.sam");
}

/// Read groups (RG:Z: aux tags), @RG/@PG/@CO header lines.
#[test]
fn htslib_xx_rg() {
    assert_bam_parity("xx#rg.sam");
}

/// Comprehensive aux tag coverage: all integer widths (c/C/s/S/i/I), floats,
/// strings (Z), hex (H), and typed B-arrays of every integer width.
#[test]
fn htslib_auxf_values() {
    assert_bam_parity("auxf#values.sam");
}
