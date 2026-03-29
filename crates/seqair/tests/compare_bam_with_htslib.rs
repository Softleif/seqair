//! Compares seqair against htslib across the entire test BAM file
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use rust_htslib::bam::{self, FetchDefinition, Read as _, record::Aux};
use seqair::bam::aux::{AuxValue, find_tag as find_aux_tag};
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

// ---- Record-level comparison across all contigs ----

struct HtsRecord {
    pos: i64,
    end_pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_len: usize,
}

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

        records.push(HtsRecord {
            pos: record.pos(),
            end_pos,
            flags: record.flags(),
            mapq: record.mapq(),
            qname: record.qname().to_vec(),
            seq,
            qual: record.qual().to_vec(),
            seq_len: record.seq_len(),
        });
    }
    records
}

// r[verify bam.reader.fetch_into+2]
// r[verify bam.record.decode]
// r[verify bam.record.fields]
#[test]
fn all_contigs_record_count_matches() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_records(contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        assert_eq!(
            store.len(),
            hts.len(),
            "{contig}: record count mismatch rio={} hts={}",
            store.len(),
            hts.len(),
        );
    }
}

// r[verify bam.record.fields]
// r[verify bam.record.end_pos]
// r[verify bam.record.seq_4bit]
// r[verify record_store.field_access]
// r[verify base_decode.decode]
#[test]
fn all_contigs_record_fields_match() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_records(contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        for (i, h) in hts.iter().enumerate() {
            let idx = i as u32;
            let r = store.record(idx);

            assert_eq!(r.pos, h.pos, "{contig} rec {i}: pos");
            assert_eq!(r.end_pos, h.end_pos - 1, "{contig} rec {i}: end_pos");
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

// ---- Pileup-level comparison across all contigs ----

struct HtsPileupCol {
    pos: i64,
    qpos_flags: Vec<(usize, u16)>,
    bases: Vec<u8>,
}

fn fetch_hts_pileup(contig: &str, start: u64, end: u64) -> Vec<HtsPileupCol> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).expect("htslib open");
    let tid = reader.header().tid(contig.as_bytes()).expect("tid");
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).expect("fetch");

    let mut columns = Vec::new();
    for p in reader.pileup() {
        let p = p.expect("pileup");
        let pos = p.pos() as i64;
        if (pos as u64) < start || (pos as u64) > end {
            continue;
        }
        let mut qpos_flags = Vec::new();
        let mut bases = Vec::new();
        for a in p.alignments() {
            if let Some(qpos) = a.qpos() {
                qpos_flags.push((qpos, a.record().flags()));
                let record = a.record();
                bases.push(record.seq()[qpos]);
            }
        }
        qpos_flags.sort();
        bases.sort();
        columns.push(HtsPileupCol { pos, qpos_flags, bases });
    }
    columns
}

// r[verify pileup.htslib_compat]
// r[verify pileup.position_iteration]
// r[verify pileup.column_contents]
// r[verify pileup.qpos]
#[test]
fn all_contigs_pileup_positions_and_depth_match() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_pileup(contig, start, end);
        if hts.is_empty() {
            continue;
        }

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        let engine = seqair::bam::PileupEngine::new(store, start as i64, end as i64);
        let rio: Vec<_> = engine.collect();

        assert_eq!(
            rio.len(),
            hts.len(),
            "{contig}: pileup column count mismatch rio={} hts={}",
            rio.len(),
            hts.len(),
        );

        for (col_idx, (r, h)) in rio.iter().zip(&hts).enumerate() {
            assert_eq!(r.pos(), h.pos, "{contig} col {col_idx}: position mismatch");
            // Compare match-depth only (htslib's qpos_flags excludes deletions/refskips)
            let match_depth = r.alignments().filter(|a| a.qpos().is_some()).count();
            assert_eq!(
                match_depth,
                h.qpos_flags.len(),
                "{contig} col {col_idx} (pos {}): match-depth mismatch",
                h.pos,
            );
        }
    }
}

// r[verify pileup.qpos]
// r[verify cigar.qpos_accuracy]
#[test]
fn all_contigs_pileup_qpos_and_flags_match() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_pileup(contig, start, end);
        if hts.is_empty() {
            continue;
        }

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        let engine = seqair::bam::PileupEngine::new(store, start as i64, end as i64);
        let rio: Vec<_> = engine.collect();

        for (col_idx, (r, h)) in rio.iter().zip(&hts).enumerate() {
            let mut qf: Vec<(usize, u16)> =
                r.alignments().filter_map(|a| a.qpos().map(|q| (q, a.flags))).collect();
            qf.sort();

            assert_eq!(
                qf.len(),
                h.qpos_flags.len(),
                "{contig} col {col_idx} (pos {}): alignment count",
                h.pos,
            );

            for (aln_idx, (ra, ha)) in qf.iter().zip(&h.qpos_flags).enumerate() {
                assert_eq!(
                    ra.0, ha.0,
                    "{contig} col {col_idx} (pos {}) aln {aln_idx}: qpos",
                    h.pos,
                );
                assert_eq!(
                    ra.1, ha.1,
                    "{contig} col {col_idx} (pos {}) aln {aln_idx}: flags",
                    h.pos,
                );
            }
        }
    }
}

// r[verify base_decode.decode]
// r[verify pileup.htslib_compat]
#[test]
fn all_contigs_pileup_bases_match() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_pileup(contig, start, end);
        if hts.is_empty() {
            continue;
        }

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        let engine = seqair::bam::PileupEngine::new(store, start as i64, end as i64);
        let rio: Vec<_> = engine.collect();

        for (col_idx, (r, h)) in rio.iter().zip(&hts).enumerate() {
            let mut bases: Vec<u8> =
                r.alignments().filter_map(|a| a.base().map(|b| b as u8)).collect();
            bases.sort();

            assert_eq!(bases, h.bases, "{contig} col {col_idx} (pos {}): bases mismatch", h.pos,);
        }
    }
}

// ---- Aux tag comparison ----

struct HtsAuxRecord {
    rg: Option<String>,
    nm: Option<i64>,
    md: Option<String>,
}

fn fetch_hts_aux_records(contig: &str, start: u64, end: u64) -> Vec<HtsAuxRecord> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).expect("htslib open");
    let tid = reader.header().tid(contig.as_bytes()).expect("tid");
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).expect("fetch");

    let mut records = Vec::new();
    let mut record = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut record) {
        if record.flags() & 0x4 != 0 {
            continue;
        }

        let rg = match record.aux(b"RG") {
            Ok(Aux::String(s)) => Some(s.to_owned()),
            _ => None,
        };
        let nm = match record.aux(b"NM") {
            Ok(Aux::I8(v)) => Some(v as i64),
            Ok(Aux::U8(v)) => Some(v as i64),
            Ok(Aux::I16(v)) => Some(v as i64),
            Ok(Aux::U16(v)) => Some(v as i64),
            Ok(Aux::I32(v)) => Some(v as i64),
            Ok(Aux::U32(v)) => Some(v as i64),
            _ => None,
        };
        let md = match record.aux(b"MD") {
            Ok(Aux::String(s)) => Some(s.to_owned()),
            _ => None,
        };

        records.push(HtsAuxRecord { rg, nm, md });
    }
    records
}

/// Extract a Z-type (string) aux tag from raw BAM aux bytes.
fn aux_z_tag<'a>(aux: &'a [u8], tag: &[u8; 2]) -> Option<&'a str> {
    match find_aux_tag(aux, *tag)? {
        AuxValue::String(bytes) => std::str::from_utf8(bytes).ok(),
        _ => None,
    }
}

/// Extract an integer aux tag (any integer sub-type) from raw BAM aux bytes.
fn aux_int_tag(aux: &[u8], tag: &[u8; 2]) -> Option<i64> {
    match find_aux_tag(aux, *tag)? {
        AuxValue::I8(v) => Some(v as i64),
        AuxValue::U8(v) => Some(v as i64),
        AuxValue::I16(v) => Some(v as i64),
        AuxValue::U16(v) => Some(v as i64),
        AuxValue::I32(v) => Some(v as i64),
        AuxValue::U32(v) => Some(v as i64),
        _ => None,
    }
}

// r[verify bam.record.fields]
#[test]
fn all_contigs_aux_tags_match() {
    for &(contig, start, end) in CONTIGS {
        let hts = fetch_hts_aux_records(contig, start, end);

        let mut reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("open");
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(contig).expect("tid");
        reader.fetch_into(tid, start, end, &mut store).expect("fetch");

        assert_eq!(
            store.len(),
            hts.len(),
            "{contig}: record count mismatch store={} hts={}",
            store.len(),
            hts.len(),
        );

        for (i, h) in hts.iter().enumerate() {
            let idx = i as u32;
            let aux = store.aux(idx);

            // RG (read group) — Z-type string
            let rio_rg = aux_z_tag(aux, b"RG");
            assert_eq!(
                rio_rg,
                h.rg.as_deref(),
                "{contig} rec {i}: RG mismatch rio={rio_rg:?} hts={:?}",
                h.rg,
            );

            // NM (edit distance) — integer
            let rio_nm = aux_int_tag(aux, b"NM");
            assert_eq!(rio_nm, h.nm, "{contig} rec {i}: NM mismatch rio={rio_nm:?} hts={:?}", h.nm,);

            // MD (mismatching positions) — Z-type string
            let rio_md = aux_z_tag(aux, b"MD");
            assert_eq!(
                rio_md,
                h.md.as_deref(),
                "{contig} rec {i}: MD mismatch rio={rio_md:?} hts={:?}",
                h.md,
            );
        }
    }
}

// TODO (nice-to-have): Compare @PG and @RG header lines, not just target names/lengths.
// TODO (nice-to-have): Compare additional aux tag types beyond RG (MD, NM, XS, etc.).
