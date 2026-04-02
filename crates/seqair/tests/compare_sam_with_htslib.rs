//! Cross-validates seqair's SAM reader against htslib reading the same .sam.gz file.
//! This catches bugs that SAM-vs-BAM tests would miss (if both our readers had the same bug).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
#![allow(clippy::arithmetic_side_effects)]
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use seqair::bam::{Pos, RecordStore, Zero};
use seqair::sam::reader::IndexedSamReader;
use std::path::Path;
use std::process::Command;

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

struct HtsRecord {
    pos: i64,
    end_pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_len: usize,
    rg: Option<String>,
}

/// Read from the original BAM with htslib as the ground truth.
fn fetch_hts_records(contig: &str, start: u64, end: u64) -> Vec<HtsRecord> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).expect("htslib open");
    let tid = reader.header().tid(contig.as_bytes()).expect("tid");
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).expect("fetch");

    let mut records = Vec::new();
    let mut record = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut record) {
        if record.flags() & 0x4 != 0 {
            continue;
        }
        let seq: Vec<u8> = (0..record.seq_len()).map(|i| record.seq()[i]).collect();
        record.cache_cigar();
        let end_pos = record.cigar_cached().unwrap().end_pos();
        let rg = match record.aux(b"RG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.to_owned()),
            _ => None,
        };

        records.push(HtsRecord {
            pos: record.pos(),
            end_pos,
            flags: record.flags(),
            mapq: record.mapq(),
            qname: record.qname().to_vec(),
            seq,
            qual: record.qual().to_vec(),
            seq_len: record.seq_len(),
            rg,
        });
    }
    records
}

// r[verify sam.reader.fetch_into]
// r[verify sam.record.parse]
// r[verify sam.record.coordinate_conversion]
#[test]
fn sam_record_count_matches_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut reader = IndexedSamReader::open(&sam_gz).expect("rio open");

    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_records(contig, start, end);

        let tid = reader.header().tid(contig).expect("tid");
        let mut store = RecordStore::new();
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store,
            )
            .expect("rio fetch");

        assert_eq!(
            store.len(),
            hts.len(),
            "{contig}: count mismatch rio={} hts={}",
            store.len(),
            hts.len(),
        );
    }
}

// r[verify sam.record.cigar_parse]
// r[verify sam.record.seq_decode]
// r[verify sam.record.qual_decode]
#[test]
fn sam_record_fields_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut reader = IndexedSamReader::open(&sam_gz).expect("rio open");

    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_records(contig, start, end);

        let tid = reader.header().tid(contig).expect("tid");
        let mut store = RecordStore::new();
        reader
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store,
            )
            .expect("rio fetch");

        for (i, h) in hts.iter().enumerate() {
            let idx = i as u32;
            let r = store.record(idx);

            assert_eq!(r.pos.as_i64(), h.pos, "{contig} rec {i}: pos");
            // htslib end_pos is exclusive (past-the-end), ours is inclusive
            assert_eq!(r.end_pos.as_i64(), h.end_pos - 1, "{contig} rec {i}: end_pos");
            assert_eq!(r.flags, h.flags, "{contig} rec {i}: flags");
            assert_eq!(r.mapq, h.mapq, "{contig} rec {i}: mapq");
            assert_eq!(store.qname(idx), h.qname.as_slice(), "{contig} rec {i}: qname");
            assert_eq!(r.seq_len as usize, h.seq_len, "{contig} rec {i}: seq_len");

            // Quality scores
            assert_eq!(store.qual(idx), h.qual.as_slice(), "{contig} rec {i}: qual");

            // Sequence (Base vs ASCII)
            for pos in 0..h.seq_len {
                let base = store.seq_at(idx, pos);
                let hts_base = seqair_types::Base::from(h.seq[pos]);
                assert_eq!(base, hts_base, "{contig} rec {i} pos {pos}: base");
            }
        }
    }
}

// r[verify sam.record.aux_tags]
#[test]
fn sam_aux_tags_present() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut reader = IndexedSamReader::open(&sam_gz).expect("rio open");

    let tid = reader.header().tid("chr19").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(6_105_700).unwrap(),
            Pos::<Zero>::new(6_105_800).unwrap(),
            &mut store,
        )
        .expect("fetch");

    // Every record in the test data should have aux tags (at least RG)
    let mut has_aux = 0;
    for i in 0..store.len() as u32 {
        if !store.aux(i).is_empty() {
            has_aux += 1;
        }
    }
    assert!(
        has_aux > 0,
        "expected at least some records with aux tags, got none out of {}",
        store.len()
    );
}

// r[verify sam.record.aux_tags]
#[test]
fn sam_aux_rg_tag_matches_htslib() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut reader = IndexedSamReader::open(&sam_gz).expect("rio open");

    let hts = fetch_hts_records("chr19", 6_105_700, 6_105_800);

    let tid = reader.header().tid("chr19").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(6_105_700).unwrap(),
            Pos::<Zero>::new(6_105_800).unwrap(),
            &mut store,
        )
        .expect("fetch");

    assert_eq!(store.len(), hts.len());

    for (i, h) in hts.iter().enumerate() {
        let aux = store.aux(i as u32);

        if let Some(hts_rg) = &h.rg {
            let rg_tag_found = find_z_tag(aux, b"RG");
            assert!(
                rg_tag_found.is_some(),
                "rec {i}: htslib has RG={hts_rg} but rio aux has no RG tag (aux len={})",
                aux.len()
            );
            let rg_value = rg_tag_found.unwrap();
            assert_eq!(
                rg_value, hts_rg,
                "rec {i}: RG tag value mismatch: rio={rg_value:?} hts={hts_rg}"
            );
        }
    }
}

// TODO (nice-to-have): Systematically compare all aux tags, not just RG.
// TODO (nice-to-have): Compare with noodles-sam as second oracle.

/// Find a Z-type tag in BAM binary aux data, return its string value.
fn find_z_tag<'a>(aux: &'a [u8], tag: &[u8; 2]) -> Option<&'a str> {
    let mut pos = 0;
    while pos + 3 <= aux.len() {
        let t = aux.get(pos..pos + 2)?;
        let typ = *aux.get(pos + 2)?;
        pos += 3;

        if typ == b'Z' {
            let end = aux.get(pos..)?.iter().position(|&b| b == 0)?;
            let val = std::str::from_utf8(aux.get(pos..pos + end)?).ok()?;
            if t == tag {
                return Some(val);
            }
            pos += end + 1;
        } else if typ == b'A' || typ == b'c' || typ == b'C' {
            pos += 1;
        } else if typ == b's' || typ == b'S' {
            pos += 2;
        } else if typ == b'i' || typ == b'I' || typ == b'f' {
            pos += 4;
        } else if typ == b'H' {
            let end = aux.get(pos..)?.iter().position(|&b| b == 0)?;
            pos += end + 1;
        } else if typ == b'B' {
            let subtype = *aux.get(pos)?;
            pos += 1;
            let count = u32::from_le_bytes([
                *aux.get(pos)?,
                *aux.get(pos + 1)?,
                *aux.get(pos + 2)?,
                *aux.get(pos + 3)?,
            ]) as usize;
            pos += 4;
            let elem_size = match subtype {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            pos += count * elem_size;
        } else {
            return None;
        }
    }
    None
}

// r[verify sam.header.parse]
// r[verify sam.header.sq_required]
#[test]
fn sam_header_matches_htslib_header() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let reader = IndexedSamReader::open(&sam_gz).expect("rio open");
    let hts_reader = bam::IndexedReader::from_path(test_bam_path()).expect("htslib open");

    let h = reader.header();
    let hts_h = hts_reader.header();

    // Same number of targets
    let hts_count = hts_h.target_count();
    assert_eq!(h.target_count(), hts_count as usize);

    // Same names and lengths
    for tid in 0..hts_count {
        let hts_name = std::str::from_utf8(hts_h.tid2name(tid).unwrap()).unwrap();
        let name = h.target_name(tid).unwrap();
        assert_eq!(name, hts_name, "tid {tid}: name");

        let hts_len = hts_h.target_len(tid).unwrap();
        let len = h.target_len(tid).unwrap();
        assert_eq!(len, hts_len, "tid {tid}: length");
    }
}
