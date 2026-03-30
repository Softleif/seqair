//! Compares seqair BAM record reading against noodles-bam.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use noodles::bam;
use noodles::sam;
use seqair::bam::{Pos, Zero};
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

/// All contigs with their full covered range (from samtools idxstats + depth).
const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

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

fn read_noodles_records(
    path: &Path,
    contig_name: &str,
    start: u64,
    end: u64,
) -> Vec<NoodlesRecord> {
    let file = std::fs::File::open(path).expect("open bam");
    let mut reader = bam::io::Reader::new(file);
    let header: sam::Header = reader.read_header().expect("read header");

    let contig_idx = header
        .reference_sequences()
        .get_index_of(contig_name.as_bytes())
        .expect("contig not found in header");

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.expect("read record");

        let flags = record.flags();
        if flags.is_unmapped() {
            continue;
        }

        let ref_seq_id = match record.reference_sequence_id().and_then(|r| r.ok()) {
            Some(id) => id,
            None => continue,
        };
        if ref_seq_id != contig_idx {
            continue;
        }

        let aln_start = match record.alignment_start().and_then(|r| r.ok()) {
            Some(pos) => usize::from(pos) as i64 - 1,
            None => continue,
        };
        if aln_start < start as i64 || aln_start >= end as i64 {
            continue;
        }

        let mapq = record.mapping_quality().map(u8::from).unwrap_or(255);

        let qname: Vec<u8> = record.name().map(|n| n.to_vec()).unwrap_or_default();

        let seq: Vec<u8> = record.sequence().iter().collect();
        let seq_len = seq.len();

        let qual: Vec<u8> = record.quality_scores().iter().collect();

        // noodles Cigar iterates Result<Op, io::Error>; each Op has len() -> usize and kind() -> Kind.
        // Convert Kind to BAM op code: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8.
        let cigar_ops: Vec<(u32, u8)> = record
            .cigar()
            .iter()
            .map(|res| {
                use noodles::sam::alignment::record::cigar::op::Kind;
                let op = res.expect("cigar op decode");
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

        records.push(NoodlesRecord {
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

    records
}

#[test]
fn bam_record_count_matches_noodles() {
    for &(contig, start, end) in CONTIGS {
        let noodles = read_noodles_records(test_bam_path(), contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32),
                Pos::<Zero>::new(end as u32),
                &mut store,
            )
            .expect("fetch");

        // seqair's index-based fetch may include records starting before
        // `start` whose alignment span overlaps the region. Filter to match
        // the noodles sequential-read filter (pos >= start).
        let seqair_count = (0..store.len() as u32)
            .filter(|&i| store.record(i).pos.as_i64() >= start as i64)
            .count();
        assert_eq!(
            seqair_count,
            noodles.len(),
            "{contig}: record count mismatch seqair={seqair_count} noodles={}",
            noodles.len(),
        );
    }
}

// r[verify bam.record.decode]
// r[verify bam.record.fields]
#[test]
fn bam_record_fields_match_noodles() {
    for &(contig, start, end) in CONTIGS {
        let noodles = read_noodles_records(test_bam_path(), contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32),
                Pos::<Zero>::new(end as u32),
                &mut store,
            )
            .expect("fetch");

        // Filter seqair records to match noodles' sequential filter (pos >= start).
        let seqair_indices: Vec<u32> = (0..store.len() as u32)
            .filter(|&i| store.record(i).pos.as_i64() >= start as i64)
            .collect();
        assert_eq!(seqair_indices.len(), noodles.len(), "{contig}: record count mismatch");

        for (i, n) in noodles.iter().enumerate() {
            let idx = seqair_indices[i];
            let r = store.record(idx);

            assert_eq!(r.pos.as_i64(), n.pos, "{contig} rec {i}: pos");
            assert_eq!(r.flags, n.flags, "{contig} rec {i}: flags");
            assert_eq!(r.mapq, n.mapq, "{contig} rec {i}: mapq");
            assert_eq!(store.qname(idx), n.qname.as_slice(), "{contig} rec {i}: qname");
        }
    }
}

// r[verify bam.record.seq_4bit]
// r[verify base_decode.decode]
#[test]
fn bam_sequence_matches_noodles() {
    for &(contig, start, end) in CONTIGS {
        let noodles = read_noodles_records(test_bam_path(), contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32),
                Pos::<Zero>::new(end as u32),
                &mut store,
            )
            .expect("fetch");

        let seqair_indices: Vec<u32> = (0..store.len() as u32)
            .filter(|&i| store.record(i).pos.as_i64() >= start as i64)
            .collect();
        assert_eq!(seqair_indices.len(), noodles.len(), "{contig}: record count mismatch");

        for (i, n) in noodles.iter().enumerate() {
            let idx = seqair_indices[i];
            let r = store.record(idx);

            assert_eq!(r.seq_len as usize, n.seq_len, "{contig} rec {i}: seq_len");

            for pos in 0..n.seq_len {
                let seqair_byte = store.seq_at(idx, pos) as u8;
                let noodles_byte = n.seq[pos];

                // Only compare standard bases (A=65, C=67, G=71, T=84).
                // Noodles may emit IUPAC codes (e.g. 'N', 'M') for non-ACGT;
                // seqair normalises all of those to Unknown=b'N'.
                match noodles_byte {
                    b'A' | b'C' | b'G' | b'T' => {
                        assert_eq!(
                            seqair_byte, noodles_byte,
                            "{contig} rec {i} pos {pos}: base mismatch \
                             seqair={seqair_byte} noodles={noodles_byte}",
                        );
                    }
                    _ => {}
                }
            }
        }
    }
}

// r[verify bam.record.fields]
#[test]
fn bam_quality_scores_match_noodles() {
    for &(contig, start, end) in CONTIGS {
        let noodles = read_noodles_records(test_bam_path(), contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32),
                Pos::<Zero>::new(end as u32),
                &mut store,
            )
            .expect("fetch");

        let seqair_indices: Vec<u32> = (0..store.len() as u32)
            .filter(|&i| store.record(i).pos.as_i64() >= start as i64)
            .collect();
        assert_eq!(seqair_indices.len(), noodles.len(), "{contig}: record count mismatch");

        for (i, n) in noodles.iter().enumerate() {
            let idx = seqair_indices[i];
            assert_eq!(
                store.qual(idx),
                n.qual.as_slice(),
                "{contig} rec {i}: quality scores mismatch",
            );
        }
    }
}

// r[verify cigar.operations]
#[test]
fn bam_cigar_matches_noodles() {
    for &(contig, start, end) in CONTIGS {
        let noodles = read_noodles_records(test_bam_path(), contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32),
                Pos::<Zero>::new(end as u32),
                &mut store,
            )
            .expect("fetch");

        let seqair_indices: Vec<u32> = (0..store.len() as u32)
            .filter(|&i| store.record(i).pos.as_i64() >= start as i64)
            .collect();
        assert_eq!(seqair_indices.len(), noodles.len(), "{contig}: record count mismatch");

        for (i, n) in noodles.iter().enumerate() {
            let idx = seqair_indices[i];
            let cigar_bytes = store.cigar(idx);

            // Decode seqair's packed u32 CIGAR: lower 4 bits = op code, upper 28 bits = length.
            let n_ops = cigar_bytes.len() / 4;
            let seqair_ops: Vec<(u32, u8)> = (0..n_ops)
                .map(|j| {
                    let b = cigar_bytes.get(j * 4..j * 4 + 4).expect("cigar slice in bounds");
                    let packed = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
                    let len = packed >> 4;
                    let op = (packed & 0xF) as u8;
                    (len, op)
                })
                .collect();

            assert_eq!(
                seqair_ops.len(),
                n.cigar_ops.len(),
                "{contig} rec {i}: CIGAR op count seqair={} noodles={}",
                seqair_ops.len(),
                n.cigar_ops.len(),
            );

            for (op_idx, (s, noo)) in seqair_ops.iter().zip(&n.cigar_ops).enumerate() {
                assert_eq!(
                    s.1, noo.1,
                    "{contig} rec {i} cigar op {op_idx}: op code seqair={} noodles={}",
                    s.1, noo.1,
                );
                assert_eq!(
                    s.0, noo.0,
                    "{contig} rec {i} cigar op {op_idx}: op length seqair={} noodles={}",
                    s.0, noo.0,
                );
            }
        }
    }
}
