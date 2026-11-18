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
    let refused = w.finish().unwrap();

    let mut clean = writer(dir.path(), &header(""));
    clean.write(&mapped(5, b"first", b"ACGTTT")).unwrap();
    clean.write(&mapped(6, b"second", b"CGTTTG")).unwrap();
    let clean = clean.finish().unwrap();
    assert!(refused.0 == clean.0, "a refused record left a trace in the file");
    assert_eq!(refused.1.unwrap().entries(), clean.1.unwrap().entries());
}

/// A sink that fails once it has taken `left` bytes.
struct Failing {
    left: usize,
}

impl std::io::Write for Failing {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if buf.len() > self.left {
            return Err(std::io::Error::other("full"));
        }
        self.left -= buf.len();
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// r[verify cram_writer.validate_before_write]
/// A sink error writing a slice poisons the writer, threaded or not.
#[test]
fn a_sink_error_poisons_the_writer() {
    let dir = tempfile::tempdir().unwrap();
    for threads in [0, 2] {
        let reference = IndexedFastaReader::open(&fasta(dir.path())).unwrap();
        let mut w = CramWriterBuilder::to_writer(Failing { left: 1000 }, &header(""), reference)
            .compression_threads(threads)
            .records_per_slice(1)
            .build()
            .unwrap();
        let mut failed = None;
        for i in 0..20 {
            if let Err(e) = w.write(&mapped(i, b"r", b"ACGT")) {
                failed = Some(e);
                break;
            }
        }
        assert!(matches!(failed, Some(CramWriteError::Io { .. })), "{failed:?}");
        assert!(matches!(w.write(&mapped(30, b"r", b"AC")), Err(CramWriteError::Poisoned)));
        assert!(matches!(w.finish(), Err(CramWriteError::Poisoned)));
    }
}

// r[verify cram_writer.tags]
// r[verify cram_writer.blocks]
/// Hundreds of tag lines in one slice (a `TL` of more than one ITF8 byte),
/// and zero-length CIGAR operations, which read back as nothing.
#[test]
fn many_tag_lines_and_zero_length_operations_read_back() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    let mut expected = Vec::new();
    for i in 0..300u32 {
        let mut r = mapped(1, format!("r{i}").as_bytes(), b"CGTA");
        for bit in 0..9 {
            if i >> bit & 1 == 1 {
                r.aux.set_int([b'X', b'0' + u8::try_from(bit).unwrap()], i64::from(i)).unwrap();
            }
        }
        expected.push(r.aux.as_bytes().to_vec());
        w.write(&r).unwrap();
    }
    let mut zero = mapped(2, b"zero", b"GTAC");
    zero.cigar = vec![
        CigarOp::new(CigarOpType::Match, 2),
        CigarOp::new(CigarOpType::Insertion, 0),
        CigarOp::new(CigarOpType::Match, 0),
        CigarOp::new(CigarOpType::Deletion, 0),
        CigarOp::new(CigarOpType::Match, 2),
    ];
    w.write(&zero).unwrap();
    let (cram, _) = w.finish().unwrap();
    let records = read_back(dir.path(), &cram);
    // Matches only: the decode generates MD/NM after the stored tags.
    for (record, aux) in records.iter().zip(&expected) {
        assert!(record.aux().as_bytes().starts_with(aux));
    }
    let last = records.last().unwrap();
    assert_eq!(last.qname(), b"zero");
    assert_eq!(last.cigar().iter().collect::<Vec<_>>(), [CigarOp::new(CigarOpType::Match, 4)]);
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

// r[verify cram_writer.header_m5.reuse]
/// The bases `build()` hashed serve the first slice: a FASTA changed on disk
/// after `build()` (same layout, other bases) does not reach the file, which
/// decodes, MD5s checked, against the original.
#[test]
fn the_first_slice_uses_the_bases_build_hashed() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    let path = dir.path().join("ref.fa");
    let original = std::fs::read(&path).unwrap();
    let mut changed = original.clone();
    for b in changed.iter_mut().skip(4).filter(|b| b.is_ascii_alphabetic()) {
        *b = b'C';
    }
    // In place, so the writer's open handle sees it.
    std::fs::OpenOptions::new().write(true).open(&path).unwrap().write_all(&changed).unwrap();
    let records = [mapped(0, b"a", b"ACGTACGTTT"), mapped(24, b"b", b"ACGTACGA")];
    for r in &records {
        w.write(r).unwrap();
    }
    let (cram, _) = w.finish().unwrap();
    std::fs::write(&path, original).unwrap();
    let back = read_back(dir.path(), &cram);
    assert_eq!(back.len(), records.len());
    for (r, b) in records.iter().zip(&back) {
        let mut seq = Vec::new();
        b.iupac_bases_into(&mut seq);
        assert_eq!(seq, r.seq);
    }
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

// r[verify cram_writer.record.unmapped]
/// An unmapped record's CIGAR and MAPQ are dropped, as htslib drops them.
#[test]
fn an_unmapped_record_loses_its_cigar_and_mapq() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    let mut r = mapped(5, b"u", b"ACGT");
    r.flags = BamFlags::from(0x4);
    w.write(&r).unwrap();
    let (cram, _) = w.finish().unwrap();
    let records = read_back(dir.path(), &cram);
    assert_eq!(records[0].mapq().as_byte(), 0);
    assert!(records[0].cigar().is_empty());
    assert_eq!(records[0].pos(), Pos0::new(5));
    assert_eq!(records[0].seq_len(), 4);
}

// r[verify cram_writer.header_m5]
/// Every reference's `M5` is the MD5 of its uppercased bases, the first
/// one's hashed whole and the others' a chunk at a time.
#[test]
fn every_reference_gets_its_md5() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("two.fa");
    let second = b"ttGACCAGTnacgTAGGGCATTACGAAAC";
    let mut text = format!(">c0\n{}\n>c1\n", std::str::from_utf8(BASES).unwrap()).into_bytes();
    text.extend_from_slice(second);
    text.push(b'\n');
    std::fs::write(&path, &text).unwrap();
    let (n0, n1) = (BASES.len(), second.len());
    let fai =
        format!("c0\t{n0}\t4\t{n0}\t{}\nc1\t{n1}\t{}\t{n1}\t{}\n", n0 + 1, 4 + n0 + 1 + 4, n1 + 1);
    std::fs::write(dir.path().join("two.fa.fai"), fai).unwrap();
    let header =
        BamHeader::from_sam_text(&format!("@SQ\tSN:c0\tLN:{n0}\n@SQ\tSN:c1\tLN:{n1}\n")).unwrap();
    let reference = IndexedFastaReader::open(&path).unwrap();
    let (cram, _) = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
        .build()
        .unwrap()
        .finish()
        .unwrap();
    let text = CramStreamReader::new(&cram[..]).unwrap().into_header().header_text().to_owned();
    for bases in [BASES, &second[..]] {
        let m5 = format!("{:x}", md5::compute(bases.to_ascii_uppercase()));
        assert!(text.contains(&format!("M5:{m5}")), "{text}");
    }
}
