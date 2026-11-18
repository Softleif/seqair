#![allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]

use super::*;
use crate::bam::aux_data::AuxData;
use crate::bam::cigar::{CigarOp, CigarOpType};
use crate::cram::stream::CramStreamReader;
use seqair_types::{BamFlags, IupacBase, MappingQuality};
use std::path::PathBuf;

const BASES: &[u8] = b"ACGTACGTTTGACCAGTNNACGTAcgtacgtACGGGCATTACG";

/// A one-contig FASTA (`c0`) and its index.
fn fasta(dir: &Path) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut text = b">c0\n".to_vec();
    text.extend_from_slice(BASES);
    text.push(b'\n');
    std::fs::write(&path, text).unwrap();
    let n = BASES.len();
    std::fs::write(dir.join("ref.fa.fai"), format!("c0\t{n}\t4\t{n}\t{}\n", n + 1)).unwrap();
    path
}

fn header(extra: &str) -> BamHeader {
    BamHeader::from_sam_text(&format!("@SQ\tSN:c0\tLN:{}{extra}\n", BASES.len())).unwrap()
}

fn mapped(pos: u32, name: &[u8], seq: &[u8]) -> OwnedBamRecord {
    let n = u32::try_from(seq.len()).unwrap();
    OwnedBamRecord::builder(Some(crate::bam::test_util::tid(0)), Pos0::new(pos), name.to_vec())
        .mapq(MappingQuality::from_byte(60))
        .cigar(vec![CigarOp::new(CigarOpType::Match, n)])
        .seq(IupacBase::from_ascii_slice(seq))
        .build()
        .unwrap()
}

fn writer(dir: &Path, header: &BamHeader) -> CramWriter<Vec<u8>> {
    let reference = IndexedFastaReader::open(&fasta(dir)).unwrap();
    CramWriterBuilder::to_writer(Vec::new(), header, reference).write_index(true).build().unwrap()
}

/// The records of `cram` as seqair decodes them, against the fixture FASTA.
fn read_back(dir: &Path, cram: &[u8]) -> Vec<crate::bam::RawBamRecord> {
    let mut reader = CramStreamReader::new(cram).unwrap();
    let mut decoder = reader.decoder(&dir.join("ref.fa")).unwrap();
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    records
}

// r[verify cram_writer.validate_before_write]
// r[verify cram_writer.sort_order]
// r[verify cram_writer.record.name]
// r[verify cram_writer.record.unmapped]
// r[verify cram_writer.tags]
/// Each refusal is typed, and a refused record leaves no trace: the records
/// around it are written as if it had not been offered.
#[test]
fn refused_records_leave_the_writer_usable() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    w.write(&mapped(5, b"first", b"ACGTTT")).unwrap();

    let mut unsorted = mapped(2, b"r", b"GT");
    let mut unmapped_with_cigar = mapped(5, b"r", b"AC");
    unmapped_with_cigar.flags = BamFlags::from(0x4);
    let mut tab = mapped(5, b"r", b"AC");
    // `set_string` refuses a TAB; a BAM read from disk can hold one.
    tab.aux = AuxData::from_bytes(b"ZZZa\tb\0".to_vec());
    let mut reserved = mapped(5, b"r", b"AC");
    reserved.aux.set_int(*b"cF", 1).unwrap();
    let mut bad_tag_name = mapped(5, b"r", b"AC");
    bad_tag_name.aux = AuxData::from_bytes(b"1XC\x05".to_vec());
    let mut double = mapped(5, b"r", b"AC");
    double.aux.set_double(*b"XD", 1.0);
    let mut short_seq = mapped(5, b"r", b"AC");
    short_seq.cigar = vec![CigarOp::new(CigarOpType::Match, 3)];
    let mut no_reference = mapped(5, b"r", b"AC");
    no_reference.tid = None;
    let past_end = mapped(u32::try_from(BASES.len()).unwrap(), b"r", b"AC");
    unsorted.pos = Pos0::new(2);

    type IsExpected = fn(&CramWriteError) -> bool;
    let cases: Vec<(OwnedBamRecord, IsExpected)> = vec![
        (unsorted, |e| matches!(e, CramWriteError::Unsorted { .. })),
        (mapped(5, b"", b"AC"), |e| matches!(e, CramWriteError::EmptyReadName)),
        (unmapped_with_cigar, |e| {
            matches!(e, CramWriteError::UnmappedWithAlignment { cigar_ops: 1, mapq: 60 })
        }),
        (mapped(5, b"a\0b", b"AC"), |e| matches!(e, CramWriteError::NulInReadName)),
        (bad_tag_name, |e| matches!(e, CramWriteError::InvalidTagName { .. })),
        (tab, |e| matches!(e, CramWriteError::TabInStringTag { .. })),
        (reserved, |e| matches!(e, CramWriteError::ReservedTag { .. })),
        (double, |e| matches!(e, CramWriteError::UnsupportedTagType { ty: 'd', .. })),
        (short_seq, |e| matches!(e, CramWriteError::CigarSeqLengthMismatch { cigar: 3, seq: 2 })),
        (no_reference, |e| matches!(e, CramWriteError::InvalidPlacement { .. })),
        (past_end, |e| matches!(e, CramWriteError::PositionPastReferenceEnd { .. })),
    ];
    for (record, expected) in cases {
        let err = w.write(&record).unwrap_err();
        assert!(expected(&err), "{err:?}");
    }
    w.write(&mapped(6, b"second", b"CGTTTG")).unwrap();
    let (cram, index) = w.finish().unwrap();

    let records = read_back(dir.path(), &cram);
    let names: Vec<&[u8]> = records.iter().map(|r| r.qname()).collect();
    assert_eq!(names, [&b"first"[..], b"second"]);
    assert_eq!(index.unwrap().entries().len(), 1);
}

// r[verify cram_writer.header_m5]
/// A missing `M5` is added from the FASTA's uppercased bases; a wrong `LN`
/// fails the build, a wrong `M5` the first slice on the reference.
#[test]
fn the_header_names_the_reference_it_was_written_against() {
    let dir = tempfile::tempdir().unwrap();
    let (cram, _) = writer(dir.path(), &header("")).finish().unwrap();
    let text = CramStreamReader::new(&cram[..]).unwrap().into_header().header_text().to_owned();
    let m5 = format!("{:x}", md5::compute(BASES.to_ascii_uppercase()));
    assert!(text.contains(&format!("M5:{m5}")), "{text}");

    let reference = IndexedFastaReader::open(&fasta(dir.path())).unwrap();
    let long = BamHeader::from_sam_text("@SQ\tSN:c0\tLN:9999\n").unwrap();
    let err = CramWriterBuilder::to_writer(Vec::new(), &long, reference).build().err();
    assert!(matches!(err, Some(CramWriteError::ReferenceLengthMismatch { fasta: 43, .. })));

    let mut w = writer(dir.path(), &header("\tM5:00000000000000000000000000000000"));
    let err = w.write(&mapped(0, b"r", b"AC")).unwrap_err();
    assert!(matches!(err, CramWriteError::ReferenceMd5Mismatch { .. }), "{err:?}");
}

// r[verify cram_writer.md_nm]
/// A correct `MD`/`NM` is dropped and regenerated; a wrong one is kept.
#[test]
fn md_nm_are_dropped_only_where_the_decode_regenerates_them() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    // c0 at 0 is ACGTACGT; one substitution at offset 2.
    let mut right = mapped(0, b"right", b"ACTTACGT");
    right.aux = AuxData::new();
    right.aux.set_string(*b"MD", b"2G5").unwrap();
    right.aux.set_int(*b"NM", 1).unwrap();
    let mut wrong = mapped(0, b"wrong", b"ACTTACGT");
    wrong.aux.set_string(*b"MD", b"8").unwrap();
    wrong.aux.set_int(*b"NM", 0).unwrap();
    w.write(&right).unwrap();
    w.write(&wrong).unwrap();
    let (cram, _) = w.finish().unwrap();
    let records = read_back(dir.path(), &cram);
    let aux = |i: usize| records[i].aux().as_bytes().to_vec();
    assert_eq!(aux(0), b"MDZ2G5\0NMC\x01");
    assert_eq!(aux(1), b"MDZ8\0NMC\x00");
}
