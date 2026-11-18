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
// r[verify cram_writer.blocks+2]
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

// r[verify cram_writer.md_nm+2]
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

// r[verify cram_writer.substitution+2]
/// The slice's `SM`, read by the decoder's own parser, turns every coded
/// `BS` back into the read base, and codes each row's bases from the most
/// substituted down.
#[hegel::test]
fn the_substitution_matrix_orders_bases_by_frequency(tc: hegel::TestCase) {
    use crate::cram::compression_header::SubstitutionMatrix;
    use hegel::generators as gs;
    let n = tc.draw(gs::integers::<usize>().max_value(400));
    // Skewed draws, so some rows have a clear favourite and some ties.
    let pairs: Vec<(usize, usize)> = (0..n)
        .map(|_| {
            let row = tc.draw(gs::integers::<usize>().max_value(4));
            let alt = tc
                .draw(gs::integers::<usize>().max_value(3))
                .min(tc.draw(gs::integers::<usize>().max_value(3)));
            (row, alt)
        })
        .collect();
    let refs = *b"ACGTN";
    let others = |row: usize| refs.iter().copied().filter(move |&b| b != refs[row]).take(4);
    let encoded: Vec<u8> =
        pairs.iter().map(|&(row, alt)| u8::try_from(row * 4 + alt).unwrap()).collect();
    let (sm, coded) = container::substitution_matrix(&encoded);
    let matrix = SubstitutionMatrix::parse(&sm);
    for (&(row, alt), &code) in pairs.iter().zip(&coded) {
        assert_eq!(matrix.substitute(refs[row], code), others(row).nth(alt).unwrap());
    }
    for (row, codes) in matrix.matrix.iter().enumerate() {
        let count = |base: u8| {
            pairs.iter().filter(|&&(r, a)| r == row && others(row).nth(a) == Some(base)).count()
        };
        let counts: Vec<usize> = codes.iter().map(|&b| count(b)).collect();
        assert!(counts.windows(2).all(|w| w[0] >= w[1]), "row {row}: {codes:?} {counts:?}");
    }
}

// r[verify cram_writer.data_series+2]
// r[verify cram_writer.blocks+2]
/// Two reads, alike but for their names and positions, leave most series
/// constant: those take no block, and the slice reads back.
#[test]
fn constant_series_take_no_block() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(dir.path(), &header(""));
    let records = [mapped(0, b"a", b"ACGT"), mapped(3, b"bb", b"TACG")];
    for r in &records {
        w.write(r).unwrap();
    }
    let (cram, _) = w.finish().unwrap();
    let mut reader = CramStreamReader::new(&cram[..]).unwrap();
    let container = reader.next_container().unwrap().unwrap();
    // Both match the reference, have no qualities, and are detached: only
    // `RN` and `AP` vary. Compression header, slice header, core, and those two.
    assert_eq!(container.header.num_blocks, 5);
    let back = read_back(dir.path(), &cram);
    assert_eq!(back.len(), 2);
    for (r, b) in records.iter().zip(&back) {
        assert_eq!(b.qname(), &r.qname[..]);
        assert_eq!(b.pos(), r.pos);
    }
}

#[test]
fn a_constant_is_one_value_or_none() {
    use crate::cram::varint::write_itf8;
    let itf8 = |values: &[i32]| {
        let mut out = Vec::new();
        for &v in values {
            write_itf8(&mut out, v);
        }
        out
    };
    let c = |s: container::Series, data: &[u8]| s.constant(data);
    assert_eq!(c(container::Series::Mq, &itf8(&[60, 60, 60])), Some(60));
    assert_eq!(c(container::Series::Rg, &itf8(&[-1, -1])), Some(-1));
    assert_eq!(c(container::Series::Ap, &itf8(&[300_000, 300_000])), Some(300_000));
    assert_eq!(c(container::Series::Ap, &itf8(&[1, 300_000])), None);
    assert_eq!(c(container::Series::Tl, &[]), Some(0));
    assert_eq!(c(container::Series::Qs, &[0xFF; 9]), Some(0xFF));
    assert_eq!(c(container::Series::Qs, &[30, 31]), None);
    assert_eq!(c(container::Series::In, &[]), None);
}

/// One slice's worth of `tests/data/test.bam`, 10 000 records as htslib cuts
/// slices: its `QS` block (`0xFF` for a record without qualities, as the
/// writer stores it) and records, its names NUL-terminated, its aux bytes.
#[derive(Default)]
struct Slice {
    quals: Vec<u8>,
    records: Vec<fqzcomp::QualRecord>,
    names: Vec<u8>,
    aux: Vec<u8>,
}

fn test_bam_slices() -> Vec<Slice> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");
    let mut reader =
        crate::bam::stream::BamStreamReader::new(std::fs::File::open(path).unwrap()).unwrap();
    let mut record = crate::bam::RawBamRecord::default();
    let mut slices = Vec::new();
    let mut s = Slice::default();
    while reader.read_record(&mut record).unwrap() {
        let start = s.quals.len();
        match record.qual() {
            Some(q) => s.quals.extend(q.iter().map(|q| q.as_byte())),
            None => s.quals.resize(start + record.seq_len(), 0xFF),
        }
        let flags = record.flags();
        s.records.push(fqzcomp::QualRecord {
            len: u32::try_from(s.quals.len() - start).unwrap(),
            reverse: flags.is_reverse(),
            read2: flags.is_second_in_template(),
        });
        s.names.extend_from_slice(record.qname());
        s.names.push(0);
        s.aux.extend_from_slice(record.aux().as_bytes());
        if s.records.len() == 10_000 {
            slices.push(std::mem::take(&mut s));
        }
    }
    if !s.records.is_empty() {
        slices.push(s);
    }
    slices
}

fn htscodecs_records(s: &Slice) -> crate::cram::fqzcomp_reference::htscodecs_fqz::Records {
    use crate::cram::fqzcomp_reference::htscodecs_fqz::{FREAD2, FREVERSE, Records};
    let flag = |on: bool, f: u32| if on { f } else { 0 };
    Records {
        quals: s.quals.clone(),
        lens: s.records.iter().map(|r| r.len).collect(),
        flags: s
            .records
            .iter()
            .map(|r| flag(r.reverse, FREVERSE) | flag(r.read2, FREAD2))
            .collect(),
        alphabet: Vec::new(),
    }
}

// r[verify cram_writer.fqzcomp]
// r[verify cram_writer.arith]
/// Real qualities: every fqzcomp strategy and every arith combination of
/// the profile round-trips and is byte for byte htscodecs' stream.
#[test]
fn test_bam_qualities_match_htscodecs() {
    use crate::cram::arith::htscodecs_oracle::htscodecs_compress as arith_htscodecs;
    use crate::cram::fqzcomp_reference::htscodecs_fqz::htscodecs_compress as fqz_htscodecs;
    let slices = test_bam_slices();
    assert!(!slices.is_empty());
    for s in &slices {
        for (strat, strategy) in (0..).zip(fqzcomp::Strategy::ALL) {
            let ours = fqzcomp::encode(&s.quals, &s.records, strategy).unwrap();
            assert!(crate::cram::fqzcomp::decode(&ours).unwrap() == s.quals, "{strategy:?}");
            let theirs = fqz_htscodecs(&htscodecs_records(s), 3, strat, None).unwrap();
            assert!(ours == theirs, "{strategy:?}: {} vs htscodecs' {}", ours.len(), theirs.len());
        }
        for flags in arith::PROFILE {
            let ours = arith::encode(&s.quals, flags).unwrap();
            assert!(crate::cram::arith::decode(&ours, 0).unwrap() == s.quals, "{flags:?}");
            assert!(ours == arith_htscodecs(&s.quals, flags.byte().into()), "{flags:?}");
        }
    }
}

/// Sizes and speeds of the archive codecs on `test.bam`'s slices, against
/// rANS Nx16's normal profile and htscodecs. Run with `cargo nextest run
/// --release --run-ignored only -E 'test(archive_codec_report)' --no-capture`.
#[test]
#[ignore = "a report, not a check"]
#[allow(clippy::cast_precision_loss, clippy::print_stdout, reason = "a report")]
fn archive_codec_report() {
    use crate::cram::arith::htscodecs_oracle::htscodecs_compress as arith_htscodecs;
    use crate::cram::fqzcomp_reference::htscodecs_fqz::htscodecs_compress as fqz_htscodecs;
    use std::time::{Duration, Instant};

    fn time(mut f: impl FnMut() -> usize) -> (usize, Duration) {
        let start = Instant::now();
        let mut size = 0;
        let mut runs = 0u32;
        while runs < 3 || start.elapsed() < Duration::from_millis(500) {
            size = f();
            runs += 1;
        }
        (size, start.elapsed() / runs)
    }

    let slices = test_bam_slices();
    let mb = |bytes: usize, t: Duration| bytes as f64 / 1e6 / t.as_secs_f64();
    type Block = fn(&Slice) -> &[u8];
    let blocks: [(&str, Block); 3] =
        [("QS", |s| &s.quals), ("RN", |s| &s.names), ("aux", |s| &s.aux)];
    for (name, block) in blocks {
        let raw: usize = slices.iter().map(|s| block(s).len()).sum();
        println!("{name}: {raw} bytes in {} slices", slices.len());
        let (nx16, t) = time(|| {
            slices
                .iter()
                .map(|s| {
                    let d = block(s);
                    rans_nx16::normal_profile(d.len())
                        .iter()
                        .map(|&f| rans_nx16::encode(d, f).unwrap().len())
                        .min()
                        .unwrap()
                })
                .sum()
        });
        println!("  rANS Nx16 normal profile best {nx16:>9}  {:6.1} MB/s (all six)", mb(raw, t));
        let (best, t) = time(|| {
            slices
                .iter()
                .map(|s| {
                    arith::PROFILE.iter().map(|&f| arith::encode(block(s), f).unwrap().len()).min()
                })
                .map(Option::unwrap)
                .sum()
        });
        println!("  arith profile best            {best:>9}  {:6.1} MB/s (all eight)", mb(raw, t));
        for flags in arith::PROFILE {
            let (ours, t) =
                time(|| slices.iter().map(|s| arith::encode(block(s), flags).unwrap().len()).sum());
            let (theirs, th) = time(|| {
                slices.iter().map(|s| arith_htscodecs(block(s), flags.byte().into()).len()).sum()
            });
            println!(
                "  arith {:>3}  {ours:>9}  {:6.1} MB/s | htscodecs {theirs:>9}  {:6.1} MB/s",
                flags.byte(),
                mb(raw, t),
                mb(raw, th)
            );
        }
        if name == "QS" {
            let records: Vec<_> = slices.iter().map(htscodecs_records).collect();
            for (strat, strategy) in (0..).zip(fqzcomp::Strategy::ALL) {
                let (ours, t) = time(|| {
                    slices
                        .iter()
                        .map(|s| fqzcomp::encode(&s.quals, &s.records, strategy).unwrap().len())
                        .sum()
                });
                let (theirs, th) = time(|| {
                    records.iter().map(|r| fqz_htscodecs(r, 3, strat, None).unwrap().len()).sum()
                });
                println!(
                    "  fqzcomp {strategy:?}  {ours:>9}  {:6.1} MB/s | htscodecs {theirs:>9}  {:6.1} MB/s",
                    mb(raw, t),
                    mb(raw, th)
                );
            }
        }
    }
}
