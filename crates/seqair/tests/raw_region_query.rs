//! `query_raw` on the fixture BAM and a bgzipped SAM of it, against two
//! oracles: what `samtools view <file> <region>` prints, and the records the
//! streaming readers yield, filtered to the region by a CIGAR walk of their
//! own.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use core::range::RangeInclusive;
use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::stream::BamStreamReader;
use seqair::bam::{BgzfStream, RawBamRecord};
use seqair::reader::{IndexedReader, ReaderError, ResolveTid};
use seqair::sam::SamStreamReader;
use seqair::sam::format::format_sam_record;
use seqair_types::Pos0;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data").join(name)
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().unwrap();
    assert!(out.status.success(), "{cmd:?}: {}", out.stderr.escape_ascii());
    out.stdout
}

struct Fixture {
    bam: PathBuf,
    sam: PathBuf,
    bam_records: Vec<RawBamRecord>,
    sam_records: Vec<RawBamRecord>,
    _dir: tempfile::TempDir,
}

/// The fixture BAM, the same records as bgzipped, tabix-indexed SAM, and
/// every record of each, streamed.
fn fixture() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let bam = data("test.bam");
        let sam = dir.path().join("test.sam.gz");
        run(Command::new("samtools")
            .args(["view", "--no-PG", "-h", "--output-fmt", "SAM,level=6", "-o"])
            .arg(&sam)
            .arg(&bam));
        run(Command::new("tabix").args(["-p", "sam"]).arg(&sam));
        let open = |path: &Path| BgzfStream::builder(std::fs::File::open(path).unwrap()).build();
        let bam_records = BamStreamReader::from_bgzf(open(&bam))
            .unwrap()
            .records()
            .collect::<Result<_, _>>()
            .unwrap();
        let sam_records =
            SamStreamReader::new(open(&sam)).unwrap().records().collect::<Result<_, _>>().unwrap();
        Fixture { bam, sam, bam_records, sam_records, _dir: dir }
    })
}

/// Last reference base a record covers, from its decoded CIGAR: its POS when
/// it is unmapped or spans no reference (htslib's `bam_endpos` minus one).
fn last_base(record: &RawBamRecord, pos: u32) -> u32 {
    if record.flags().is_unmapped() {
        return pos;
    }
    let rlen: u32 = record
        .cigar()
        .map(|op| op.to_bam_u32())
        .filter(|op| matches!(op & 0xF, 0 | 2 | 3 | 7 | 8))
        .map(|op| op >> 4)
        .sum();
    pos + rlen.saturating_sub(1)
}

/// The streamed records on `tid` that overlap `span`, in file order.
fn streamed(records: &[RawBamRecord], tid: u32, span: RangeInclusive<Pos0>) -> Vec<&[u8]> {
    let (start, last) = (span.start.as_u32(), span.last.as_u32());
    records
        .iter()
        .filter(|r| {
            let (Some(t), Some(pos)) = (r.tid(), r.pos()) else { return false };
            let pos = pos.as_u32();
            t.as_u32() == tid && pos <= last && last_base(r, pos) >= start
        })
        .map(RawBamRecord::as_payload)
        .collect()
}

/// A region of a contig `len` long: anchored a base or two off a streamed
/// record's first or last base, or anywhere, a few widths wide.
#[hegel::composite]
fn region(tc: &TestCase, len: u32, ends: Vec<(u32, u32)>) -> RangeInclusive<Pos0> {
    let last = len - 1;
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, 5_000]));
    let start = if !ends.is_empty() && tc.draw(gs::weighted_booleans(0.7)) {
        let (pos, end) = ends[tc.draw(gs::integers::<usize>().max_value(ends.len() - 1))];
        let anchor = if tc.draw(gs::booleans()) { pos } else { end };
        anchor.saturating_add_signed(tc.draw(gs::integers::<i32>().min_value(-2).max_value(2)))
    } else {
        tc.draw(gs::integers::<u32>().max_value(last))
    }
    .min(last);
    let start = Pos0::new(start).unwrap();
    RangeInclusive {
        start,
        last: Pos0::new(start.as_u32().saturating_add(width).min(last)).unwrap(),
    }
}

/// `query_raw`'s records, read into one reused record, checked against the
/// same query's iterator.
fn query(reader: &mut IndexedReader, tid: u32, span: RangeInclusive<Pos0>) -> Vec<RawBamRecord> {
    let tid = tid.resolve_tid(reader.header()).unwrap();
    let mut q = reader.query_raw(tid, span).unwrap();
    let mut got = Vec::new();
    let mut record = RawBamRecord::default();
    while q.read_record(&mut record).unwrap() {
        got.push(record.clone());
    }
    let last = record.clone();
    assert!(!q.read_record(&mut record).unwrap(), "a spent query stays spent");
    assert!(record == last, "the end leaves the record as it was");
    drop(q);
    let iterated: Vec<RawBamRecord> =
        reader.query_raw(tid, span).unwrap().collect::<Result<_, _>>().unwrap();
    assert!(iterated == got, "read_record and the iterator differ");
    got
}

fn payloads(records: &[RawBamRecord]) -> Vec<&[u8]> {
    records.iter().map(RawBamRecord::as_payload).collect()
}

// r[verify bam.reader.query_raw]
// r[verify bam.reader.query_raw.errors]
// r[verify sam.reader.query_raw]
// r[verify unified.query_raw]
// r[verify bam.reader.resume]
/// A few regions in a row on one BAM and one SAM reader (so the resume
/// shortcut is in play): each query gives exactly the streamed records that
/// overlap it, byte for byte and in file order, and, formatted as SAM, the
/// lines `samtools view` prints for it.
#[hegel::test(test_cases = pinned::cases(12))]
fn raw_regions_are_the_streamed_records_and_what_samtools_prints(tc: TestCase) {
    let f = fixture();
    let mut bam = IndexedReader::open(&f.bam).unwrap();
    let mut sam = IndexedReader::open(&f.sam).unwrap();
    let header = bam.header().clone();
    let n_targets = header.target_count();

    for _ in 0..tc.draw(gs::integers::<usize>().min_value(1).max_value(4)) {
        let tid = tc.draw(gs::integers::<u32>().max_value(n_targets as u32 - 1));
        let resolved = tid.resolve_tid(&header).unwrap();
        let name = header.target_name(resolved).unwrap().to_owned();
        let len = header.target_len(resolved).unwrap() as u32;
        let ends: Vec<(u32, u32)> = f
            .bam_records
            .iter()
            .filter(|r| r.tid().is_some_and(|t| t.as_u32() == tid))
            .filter_map(|r| r.pos().map(|p| (p.as_u32(), last_base(r, p.as_u32()))))
            .collect();
        let span = tc.draw(region(len, ends));

        let from_bam = query(&mut bam, tid, span);
        assert!(
            payloads(&from_bam) == streamed(&f.bam_records, tid, span),
            "BAM {name}:{}-{} (0-based) differs from the stream",
            span.start.as_u32(),
            span.last.as_u32()
        );
        let from_sam = query(&mut sam, tid, span);
        assert!(
            payloads(&from_sam) == streamed(&f.sam_records, tid, span),
            "SAM {name}:{}-{} (0-based) differs from the stream",
            span.start.as_u32(),
            span.last.as_u32()
        );

        let region = format!("{name}:{}-{}", span.start.as_u32() + 1, span.last.as_u32() + 1);
        let want = run(Command::new("samtools").arg("view").arg(&f.bam).arg(&region));
        for (format, records) in [("BAM", &from_bam), ("SAM", &from_sam)] {
            let mut text = Vec::new();
            for record in records.iter() {
                format_sam_record(record, &header, &mut text).unwrap();
                text.push(b'\n');
            }
            assert!(
                text == want,
                "{format} {region} differs from samtools:\n{}\nwant\n{}",
                text.escape_ascii(),
                want.escape_ascii()
            );
        }
        tc.event_value("records", from_bam.len() as f64);
    }
}

// r[verify unified.query_raw]
#[test]
fn a_cram_raw_query_is_unsupported() {
    let mut reader =
        IndexedReader::open_with_reference(&data("test_v30.cram"), &data("test.fasta.gz")).unwrap();
    let tid = reader.header().tid("chr19").unwrap();
    let span = (Pos0::new(6_105_700).unwrap()..=Pos0::new(6_105_800).unwrap()).into();
    let err = reader.query_raw(tid, span).unwrap_err();
    assert!(matches!(err, ReaderError::RawQueryUnsupported), "{err:?}");
}

/// What a generated record has wrong with it, if anything.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Flaw {
    None,
    /// One base more than its CIGAR covers: htslib rejects it, and neither
    /// `RecordStore::push_raw` nor the indexed SAM reader checks.
    CigarSeqMismatch,
    /// RNEXT (in the SAM) names a contig the header lacks: htslib keeps the
    /// line with no mate reference.
    UnknownRnext,
}

/// A record on `chr1` of the generated files: its 0-based start, the bases
/// its `M` covers, and its flaw.
#[derive(Debug, Clone, Copy)]
struct Generated {
    pos: u32,
    len: u32,
    flaw: Flaw,
}

const GENERATED_LEN: u32 = 4_000;

/// Records sorted by start, some of them flawed.
#[hegel::composite]
fn generated_records(tc: &TestCase) -> Vec<Generated> {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(60));
    let mut records: Vec<Generated> = (0..n)
        .map(|_| Generated {
            pos: tc.draw(gs::integers::<u32>().max_value(GENERATED_LEN - 200)),
            len: tc.draw(gs::integers::<u32>().min_value(1).max_value(150)),
            flaw: tc.draw(gs::sampled_from(vec![
                Flaw::None,
                Flaw::None,
                Flaw::CigarSeqMismatch,
                Flaw::UnknownRnext,
            ])),
        })
        .collect();
    records.sort_by_key(|r| r.pos);
    records
}

struct GeneratedFiles {
    bam: PathBuf,
    sam: PathBuf,
    _dir: tempfile::TempDir,
}

/// The records as an indexed BAM (written by seqair, which does not check
/// CIGAR against SEQ) and as a bgzipped, tabix-indexed SAM. Read `i` is
/// named `r<i>`; an unknown RNEXT is in the BAM as htslib converts it.
fn write_generated(records: &[Generated]) -> GeneratedFiles {
    use seqair::bam::MateFields;
    use seqair::bam::cigar::{CigarOp, CigarOpType};
    use seqair::bam::header::BamHeader;
    use seqair::bam::owned_record::OwnedBamRecord;
    use seqair::bam::writer::BamWriterBuilder;
    use seqair_types::{BamFlags, Base, BaseQuality, MappingQuality};

    let dir = tempfile::tempdir().unwrap();
    let header_text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{GENERATED_LEN}\n");
    let header = BamHeader::from_sam_text(&header_text).unwrap();
    let tid = "chr1".resolve_tid(&header).unwrap();
    let seq_len = |r: &Generated| r.len + u32::from(r.flaw == Flaw::CigarSeqMismatch);

    let bam = dir.path().join("generated.bam");
    let mut writer = BamWriterBuilder::to_path(&bam, &header).write_index(true).build().unwrap();
    for (i, r) in records.iter().enumerate() {
        let n = seq_len(r) as usize;
        let record = OwnedBamRecord::builder(
            Some(tid),
            Some(Pos0::new(r.pos).unwrap()),
            format!("r{i}").into_bytes(),
        )
        .flags(BamFlags::empty())
        .mapq(MappingQuality::from_byte(60))
        .cigar(vec![CigarOp::new(CigarOpType::Match, r.len)])
        .seq(vec![Base::A; n])
        .qual(vec![BaseQuality::from_byte(30); n])
        .mate_fields(MateFields {
            tid: None,
            pos: (r.flaw == Flaw::UnknownRnext).then(|| Pos0::new(6).unwrap()),
            template_len: None,
        })
        .build()
        .unwrap();
        writer.write(&record).unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let bai = std::fs::File::create(bam.with_extension("bam.bai")).unwrap();
    index.unwrap().write_bai(bai, 1).unwrap();

    let mut text = header_text;
    for (i, r) in records.iter().enumerate() {
        let (rnext, pnext) = if r.flaw == Flaw::UnknownRnext { ("chrZ", 7) } else { ("*", 0) };
        let seq = "A".repeat(seq_len(r) as usize);
        text.push_str(&format!(
            "r{i}\t1\tchr1\t{}\t60\t{}M\t{rnext}\t{pnext}\t0\t{seq}\t*\n",
            r.pos + 1,
            r.len
        ));
    }
    let plain = dir.path().join("generated.sam");
    std::fs::write(&plain, text).unwrap();
    run(Command::new("bgzip").arg(&plain));
    let sam = dir.path().join("generated.sam.gz");
    run(Command::new("tabix").args(["-p", "sam"]).arg(&sam));
    GeneratedFiles { bam, sam, _dir: dir }
}

/// A span of the generated contig, a few widths wide.
#[hegel::composite]
fn generated_span(tc: &TestCase) -> RangeInclusive<Pos0> {
    let start = tc.draw(gs::integers::<u32>().max_value(GENERATED_LEN - 1));
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, GENERATED_LEN]));
    RangeInclusive {
        start: Pos0::new(start).unwrap(),
        last: Pos0::new(start.saturating_add(width).min(GENERATED_LEN - 1)).unwrap(),
    }
}

/// The read index a generated name `r<i>` carries.
fn generated_index(qname: &[u8]) -> usize {
    std::str::from_utf8(qname).unwrap().trim_start_matches('r').parse().unwrap()
}

/// The read name of the BAM record at virtual offset `offset`, read with
/// nothing but a BGZF reader.
fn qname_at(bam: &Path, offset: u64) -> Vec<u8> {
    use seqair::bam::bgzf::{BgzfReader, VirtualOffset};
    let mut bgzf = BgzfReader::open(bam).unwrap();
    bgzf.seek_virtual(VirtualOffset::new(offset >> 16, (offset & 0xFFFF) as u16)).unwrap();
    let mut record = vec![0; bgzf.read_u32().unwrap() as usize];
    bgzf.read_exact_into(&mut record).unwrap();
    let l_read_name = usize::from(record[8]);
    record[32..32 + l_read_name - 1].to_vec()
}

/// Whether `e` is the error a record that does not validate is reported as,
/// which the query continues past.
fn is_invalid_record(e: &ReaderError) -> bool {
    matches!(
        e,
        ReaderError::Bam { source: seqair::bam::BamError::InvalidRecord { .. } }
            | ReaderError::Sam { source: seqair::sam::SamError::MalformedRecord { .. } }
    )
}

// r[verify bam.reader.query_raw]
// r[verify bam.reader.query_raw.errors]
// r[verify sam.reader.query_raw]
// r[verify unified.query_raw]
/// Every record `fetch_into` keeps comes out of `query_raw` in its place:
/// as itself, or, where it does not validate (a CIGAR that does not cover
/// SEQ), as an error the query continues past. A SAM line naming an
/// unknown mate reference is kept, without one, as `fetch_into` keeps it.
#[hegel::test(test_cases = pinned::cases(40))]
fn a_record_that_does_not_validate_is_an_error_in_its_place(tc: TestCase) {
    let records = tc.draw(generated_records());
    let files = write_generated(&records);
    let span = tc.draw(generated_span());
    for path in [&files.bam, &files.sam] {
        let mut reader = IndexedReader::open(path).unwrap();
        let tid = "chr1".resolve_tid(reader.header()).unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, span, &mut store).unwrap();
        let want: Vec<Option<usize>> = store
            .records()
            .map(|r| generated_index(r.qname(&store).unwrap()))
            .map(|i| (records[i].flaw != Flaw::CigarSeqMismatch).then_some(i))
            .collect();

        let mut query = reader.query_raw(tid, span).unwrap();
        let mut record = RawBamRecord::default();
        let mut got = Vec::new();
        loop {
            match query.read_record(&mut record) {
                Ok(true) => {
                    let i = generated_index(record.qname());
                    if records[i].flaw == Flaw::UnknownRnext {
                        let mate = record.mate_fields();
                        assert_eq!((mate.tid, mate.pos), (None, Pos0::new(6)), "{path:?} r{i}");
                    }
                    got.push(Some(i));
                }
                Ok(false) => break,
                Err(e) => {
                    assert!(is_invalid_record(&e), "{path:?}: {e:?}");
                    if let ReaderError::Bam {
                        source: seqair::bam::BamError::InvalidRecord { offset, .. },
                    } = e
                    {
                        let i = generated_index(&qname_at(path, offset));
                        assert_eq!(records[i].flaw, Flaw::CigarSeqMismatch, "r{i} at {offset:#x}");
                    }
                    got.push(None);
                }
            }
        }
        assert_eq!(got, want, "{path:?} {}..={}", span.start.as_u32(), span.last.as_u32());
        assert!(!query.read_record(&mut record).unwrap(), "{path:?}: spent, not poisoned");
        drop(query);

        let iterated: Vec<Option<usize>> = reader
            .query_raw(tid, span)
            .unwrap()
            .map(|r| match r {
                Ok(r) => Some(generated_index(r.qname())),
                Err(e) => {
                    assert!(is_invalid_record(&e), "{path:?}: {e:?}");
                    None
                }
            })
            .collect();
        assert_eq!(iterated, want, "{path:?}: the iterator");
        tc.event_value("errors", want.iter().filter(|w| w.is_none()).count() as f64);
    }
}

/// Start and length of each BGZF block in `bytes`, by the `BSIZE` of the
/// `BC` subfield htslib and seqair write first.
fn bgzf_blocks(bytes: &[u8]) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let len = usize::from(u16::from_le_bytes([bytes[at + 16], bytes[at + 17]])) + 1;
        blocks.push((at, len));
        at += len;
    }
    blocks
}

/// Enough unflawed records for several BGZF blocks.
fn many_records() -> Vec<Generated> {
    (0..3_000).map(|i| Generated { pos: i, len: 150, flaw: Flaw::None }).collect()
}

fn whole_contig() -> RangeInclusive<Pos0> {
    RangeInclusive { start: Pos0::new(0).unwrap(), last: Pos0::new(GENERATED_LEN - 1).unwrap() }
}

// r[verify bam.reader.query_raw.errors]
/// A BAM cut off inside its last block: the query yields records up to
/// there, then an error that is not about one record, and from then on only
/// `Poisoned`; the iterator stops after the error.
#[test]
fn a_truncated_bam_poisons_its_raw_query() {
    use seqair::bam::{BamError, IndexedBamReader};
    let files = write_generated(&many_records());
    let bytes = std::fs::read(&files.bam).unwrap();
    let blocks = bgzf_blocks(&bytes);
    assert!(blocks.len() > 3, "{} blocks", blocks.len());
    let (last_data, len) = blocks[blocks.len() - 2];
    std::fs::write(&files.bam, &bytes[..last_data + len / 2]).unwrap();

    let mut reader = IndexedBamReader::open(&files.bam).unwrap();
    let tid = "chr1".resolve_tid(reader.header()).unwrap();
    let mut query = reader.query_raw(tid, whole_contig()).unwrap();
    let mut record = RawBamRecord::default();
    let mut read = 0;
    let err = loop {
        match query.read_record(&mut record) {
            Ok(true) => read += 1,
            Ok(false) => panic!("the truncated query ended cleanly after {read} records"),
            Err(e) => break e,
        }
    };
    assert!(read > 0 && read < 3_000, "{read}");
    assert!(!matches!(err, BamError::InvalidRecord { .. }), "{err:?}");
    assert!(matches!(query.read_record(&mut record), Err(BamError::Poisoned)));
    assert!(matches!(query.read_record(&mut record), Err(BamError::Poisoned)));
    drop(query);

    let mut items = reader.query_raw(tid, whole_contig()).unwrap();
    let oks = items.by_ref().take_while(Result::is_ok).count();
    assert_eq!(oks, read);
    assert!(items.next().is_none(), "the iterator ends at the error");
}

// r[verify sam.reader.query_raw]
/// A bgzipped SAM with a corrupt block: as for BAM, the error ends the query
/// and later reads are `Poisoned`.
#[test]
fn a_corrupt_sam_block_poisons_its_raw_query() {
    use seqair::sam::SamError;
    use seqair::sam::reader::IndexedSamReader;
    let files = write_generated(&many_records());
    let mut bytes = std::fs::read(&files.sam).unwrap();
    let blocks = bgzf_blocks(&bytes);
    assert!(blocks.len() > 3, "{} blocks", blocks.len());
    let (middle, len) = blocks[blocks.len() / 2];
    bytes[middle + len - 8] ^= 0xFF; // the block's CRC32
    std::fs::write(&files.sam, &bytes).unwrap();

    let mut reader = IndexedSamReader::open(&files.sam).unwrap();
    let tid = "chr1".resolve_tid(reader.header()).unwrap();
    let mut query = reader.query_raw(tid, whole_contig()).unwrap();
    let mut record = RawBamRecord::default();
    let mut read = 0;
    let err = loop {
        match query.read_record(&mut record) {
            Ok(true) => read += 1,
            Ok(false) => panic!("the corrupt query ended cleanly after {read} records"),
            Err(e) => break e,
        }
    };
    assert!(read > 0 && read < 3_000, "{read}");
    assert!(!matches!(err, SamError::MalformedRecord { .. }), "{err:?}");
    assert!(matches!(query.read_record(&mut record), Err(SamError::Poisoned)));
    assert!(matches!(query.read_record(&mut record), Err(SamError::Poisoned)));
    drop(query);

    let mut items = reader.query_raw(tid, whole_contig()).unwrap();
    let oks = items.by_ref().take_while(Result::is_ok).count();
    assert_eq!(oks, read);
    assert!(items.next().is_none(), "the iterator ends at the error");
}

/// Every item of a raw query: a record's bytes, or `None` for a record that
/// does not validate.
fn raw_items(reader: &mut IndexedReader, span: RangeInclusive<Pos0>) -> Vec<Option<Vec<u8>>> {
    let tid = "chr1".resolve_tid(reader.header()).unwrap();
    reader
        .query_raw(tid, span)
        .unwrap()
        .map(|r| match r {
            Ok(r) => Some(r.into_payload()),
            Err(e) => {
                assert!(is_invalid_record(&e), "{e:?}");
                None
            }
        })
        .collect()
}

// r[verify bam.reader.resume]
// r[verify bam.reader.query_raw]
// r[verify sam.reader.query_raw]
/// A raw query read only part way — past a record that does not validate,
/// perhaps — and dropped leaves the reader as a full one would: a query
/// after it starting no earlier gives what it gives on a fresh reader.
#[hegel::test(test_cases = pinned::cases(40))]
fn a_query_after_a_partly_read_one_is_as_on_a_fresh_reader(tc: TestCase) {
    let records = tc.draw(generated_records());
    let files = write_generated(&records);
    let first = tc.draw(generated_span());
    let start =
        tc.draw(gs::integers::<u32>().min_value(first.start.as_u32()).max_value(GENERATED_LEN - 1));
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, GENERATED_LEN]));
    let second = RangeInclusive {
        start: Pos0::new(start).unwrap(),
        last: Pos0::new(start.saturating_add(width).min(GENERATED_LEN - 1)).unwrap(),
    };
    let k = tc.draw(gs::integers::<usize>().max_value(records.len()));
    for path in [&files.bam, &files.sam] {
        let mut reader = IndexedReader::open(path).unwrap();
        let tid = "chr1".resolve_tid(reader.header()).unwrap();
        let mut query = reader.query_raw(tid, first).unwrap();
        let mut record = RawBamRecord::default();
        let mut errors = 0;
        for _ in 0..k {
            match query.read_record(&mut record) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => {
                    assert!(is_invalid_record(&e), "{path:?}: {e:?}");
                    errors += 1;
                }
            }
        }
        drop(query);
        let after = raw_items(&mut reader, second);
        let fresh = raw_items(&mut IndexedReader::open(path).unwrap(), second);
        assert!(after == fresh, "{path:?}: {} items after, {} fresh", after.len(), fresh.len());
        tc.event_value("errors read in the first query", f64::from(errors));
    }
}
