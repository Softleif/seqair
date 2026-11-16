//! CIGARs on either side of BAM's 16-bit op count, through every BAM writer
//! and back: our readers and samtools must see the real CIGAR, and the tags
//! without a `CG` tag. A stale `CG` among the tags must not shadow the real
//! one.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::cigar::{CigarOpType, CigarStr};
use seqair::bam::{AuxData, BamHeader, BamWriterBuilder, CigarOp, OwnedBamRecord, RecordStore};
use seqair::reader::AlignmentStream;
use seqair_types::{BaseQuality, IupacBase, Pos0};
use std::process::Command;

const CONTIG_LEN: u32 = 1_000_000;

fn header() -> BamHeader {
    BamHeader::from_sam_text(&format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{CONTIG_LEN}\n"))
        .unwrap()
}

/// A CIGAR: `n` ops alternating `kM` and `1I` (or `1D`), on either side
/// of BAM's 16-bit op count; or a clip of the whole read and a deletion,
/// which starts as BAM's placeholder does.
#[hegel::composite]
fn cigar(tc: &TestCase) -> Vec<CigarOp> {
    if tc.draw(gs::weighted_booleans(0.2)) {
        let clip = tc.draw(gs::integers::<u32>().min_value(1).max_value(9));
        return vec![
            CigarOp::new(CigarOpType::SoftClip, clip),
            CigarOp::new(CigarOpType::Deletion, 3),
        ];
    }
    let n = if tc.draw(gs::booleans()) {
        tc.draw(gs::integers::<usize>().min_value(65_533).max_value(65_538))
    } else {
        tc.draw(gs::integers::<usize>().min_value(1).max_value(6))
    };
    let gap = tc.draw(gs::sampled_from(vec![CigarOpType::Insertion, CigarOpType::Deletion]));
    // One draw for the match lengths: a draw per op would be 65k round trips.
    let seed = tc.draw(gs::integers::<u32>());
    (0..n)
        .map(|i| {
            if i % 2 == 0 {
                let len = 1 + (seed.rotate_left(i as u32 % 32) % 3);
                CigarOp::new(CigarOpType::Match, len)
            } else {
                CigarOp::new(gap, 1)
            }
        })
        .collect()
}

/// Whether a reader takes a `CG` tag for this CIGAR: it starts with a clip
/// of the whole read, on a placed record.
fn starts_like_placeholder(ops: &[CigarOp], placed: bool) -> bool {
    let qlen = seqair::bam::cigar::query_len(ops);
    placed
        && ops
            .first()
            .is_some_and(|op| op.op_type() == CigarOpType::SoftClip && u64::from(op.len()) == qlen)
}

/// Tags, with a stale `CG` tag among them when `stale`: two ops, so that a
/// reader would take it for a two-op CIGAR that starts as the placeholder.
fn aux(stale: bool) -> AuxData {
    let mut aux = AuxData::new();
    aux.set_int(*b"NM", 3).unwrap();
    if stale {
        let ops = [CigarOp::new(CigarOpType::Match, 7), CigarOp::new(CigarOpType::Deletion, 1)];
        aux.set_array_u32(*b"CG", &ops.map(CigarOp::to_bam_u32)).unwrap();
    }
    aux.set_string(*b"RG", b"g1").unwrap();
    aux
}

fn record(cigar: Vec<CigarOp>, aux: AuxData, placed: bool) -> OwnedBamRecord {
    let qlen = seqair::bam::cigar::query_len(&cigar) as usize;
    let (tid, pos, flags) =
        if placed { (header().tid("c"), Pos0::new(99), 0) } else { (None, None, 0x4) };
    OwnedBamRecord::builder(tid, pos, b"r".to_vec())
        .flags(seqair_types::BamFlags::from(flags))
        .cigar(cigar)
        .seq(vec![IupacBase::A; qlen])
        .qual(vec![BaseQuality::from_byte(30); qlen])
        .aux(aux)
        .build()
        .unwrap()
}

/// The tags as our reader shows them: `aux` without `CG` when a reader
/// would take one for the CIGAR.
fn expected_tags(aux: &AuxData, drop_cg: bool) -> Vec<u8> {
    let mut want = aux.clone();
    if drop_cg {
        want.remove(*b"CG");
    }
    want.as_bytes().to_vec()
}

fn samtools_fields(bam: &[u8]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("r.bam");
    std::fs::write(&path, bam).unwrap();
    let out = Command::new("samtools").arg("view").arg(&path).output().expect("samtools");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let line = String::from_utf8(out.stdout).unwrap();
    line.trim_end().split('\t').map(str::to_owned).collect()
}

// r[verify bam.owned_record.long_cigar]
// r[verify bam_writer.store_long_cigar]
// r[verify sam.encode.long_cigar+1]
// r[verify record_store.long_cigar+1]
/// An `OwnedBamRecord` and the same record copied through a `RecordStore`
/// write the same bytes, which read back — by seqair and by samtools — as
/// the real CIGAR and the record's tags, a stale `CG` dropped exactly when a
/// reader would take one for the CIGAR. An unplaced record cannot carry a
/// long CIGAR, as no reader would take it out of the tag.
#[hegel::test(test_cases = pinned::cases(16))]
fn long_cigars_survive_every_writer(tc: TestCase) {
    let ops = tc.draw(cigar());
    let stale = tc.draw(gs::booleans());
    let placed = tc.draw(gs::weighted_booleans(0.8));
    let aux = aux(stale);
    let long = ops.len() > usize::from(u16::MAX);
    let drop_cg = long || starts_like_placeholder(&ops, placed);
    let owned = record(ops.clone(), aux.clone(), placed);
    let header = header();

    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    match writer.write(&owned) {
        Err(seqair::bam::BamWriteError::Record {
            source:
                seqair::bam::OwnedRecordError::LongCigar(seqair::bam::LongCigarError::Unplaced {
                    n_ops,
                }),
        }) => {
            assert!(long && !placed && n_ops == ops.len());
            return;
        }
        other => other.unwrap(),
    }
    let (from_owned, _) = writer.finish().unwrap();

    // Read it into a store and write the store's record.
    let mut stream =
        AlignmentStream::builder(std::io::Cursor::new(from_owned.clone())).build().unwrap();
    let raw = stream.records().next().unwrap().unwrap();
    assert_eq!(raw.has_cg_placeholder(), long);
    assert_eq!(raw.cigar().to_vec(), ops);
    assert_eq!(raw.aux().as_bytes(), expected_tags(&aux, drop_cg));
    let mut written = vec![(from_owned, "owned")];
    // The store holds placed records only.
    if placed {
        let mut store = RecordStore::new();
        let idx = store.push_raw(raw.as_payload(), &header, &mut ()).unwrap().unwrap();
        assert_eq!(store.record(idx).unwrap().cigar(), ops.as_slice());
        assert_eq!(store.record(idx).unwrap().aux(), expected_tags(&aux, drop_cg));
        let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
        writer.write_store_record(&store, idx).unwrap();
        written.push((writer.finish().unwrap().0, "store"));
    }

    for (bam, what) in &written {
        let fields = samtools_fields(bam);
        assert_eq!(fields[5], CigarStr(&ops).to_string(), "{what}: samtools CIGAR");
        let cg_tags = fields[11..].iter().filter(|f| f.starts_with("CG:")).count();
        assert_eq!(cg_tags, usize::from(stale && !drop_cg), "{what}: CG tags samtools shows");
    }
}

/// One BAM record payload on `c` at `pos`, mapped, placed.
fn payload(pos: i32, cigar: &[CigarOp], seq_len: usize, aux: &[u8]) -> Vec<u8> {
    let mut p = Vec::new();
    p.extend_from_slice(&0i32.to_le_bytes());
    p.extend_from_slice(&pos.to_le_bytes());
    p.extend_from_slice(&[2, 60]);
    p.extend_from_slice(&4680u16.to_le_bytes());
    p.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    p.extend_from_slice(&0u16.to_le_bytes());
    p.extend_from_slice(&(seq_len as u32).to_le_bytes());
    p.extend_from_slice(&(-1i32).to_le_bytes());
    p.extend_from_slice(&(-1i32).to_le_bytes());
    p.extend_from_slice(&0i32.to_le_bytes());
    p.extend_from_slice(b"r\0");
    p.extend_from_slice(CigarOp::ops_as_bytes(cigar));
    p.extend(std::iter::repeat_n(0x11, seq_len.div_ceil(2)));
    p.extend(std::iter::repeat_n(30, seq_len));
    p.extend_from_slice(aux);
    p
}

fn cg_tag(subtype: u8, ops: &[CigarOp]) -> Vec<u8> {
    let mut tag = vec![b'C', b'G', b'B', subtype];
    tag.extend_from_slice(&(ops.len() as u32).to_le_bytes());
    tag.extend_from_slice(CigarOp::ops_as_bytes(ops));
    tag
}

// r[verify record_store.long_cigar+1]
// r[verify bam.raw_record.cg_hidden]
/// `RecordStore::push_raw` (the indexed BAM path, which never builds a
/// `RawBamRecord`) and `RawBamRecord` read the same record the same way —
/// CIGAR, end and tags — wherever the `CG` tag sits, whatever other tags and
/// `CG` tags surround it, and however wrong the placeholder's `N` is; and the
/// reader's pre-push overlap test agrees with the store's end.
#[hegel::test(test_cases = pinned::cases(60))]
fn the_store_reads_a_long_cigar_as_the_record_does(tc: TestCase) {
    let real = tc.draw(cigar());
    let seq_len = seqair::bam::cigar::query_len(&real) as usize;
    let span = seqair::bam::cigar::ref_len(&real) as u32;
    // A placeholder whose span may disagree with the real CIGAR's.
    let n = tc.draw(gs::integers::<u32>().min_value(1).max_value(span + 3));
    let placeholder = [
        CigarOp::new(CigarOpType::SoftClip, seq_len as u32),
        CigarOp::new(CigarOpType::RefSkip, n),
    ];
    let mut aux = Vec::new();
    for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        aux.extend_from_slice(b"NMC\x02");
    }
    // A `CG` tag in front that cannot be the CIGAR: the first decides.
    let blocked = tc.draw(gs::booleans());
    if blocked {
        aux.extend(cg_tag(b'S', &[]));
    }
    aux.extend(cg_tag(tc.draw(gs::sampled_from(vec![b'I', b'i'])), &real));
    for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        match tc.draw(gs::booleans()) {
            true => aux.extend(cg_tag(b'I', &[CigarOp::new(CigarOpType::Match, 7)])),
            false => aux.extend_from_slice(b"RGZg\0"),
        }
    }
    let raw = payload(99, &placeholder, seq_len, &aux);
    let header = header();

    let record = seqair::bam::RawBamRecord::from_payload(raw.clone(), &header).unwrap();
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &header, &mut ()).unwrap().unwrap();
    let stored = store.record(idx).unwrap();
    assert_eq!(stored.cigar(), record.cigar().to_vec().as_slice());
    assert_eq!(stored.aux(), record.aux().as_bytes());
    assert_eq!(Some(stored.end_pos), record.end_pos());
    // `bam_tag2cigar` wants at least as many ops as the placeholder's two.
    if blocked || real.len() < placeholder.len() {
        assert_eq!(stored.cigar(), placeholder.as_slice());
    } else {
        assert_eq!(stored.cigar(), real.as_slice());
        assert_eq!(stored.end_pos.as_u32(), 99 + span - 1);
        if real.len() > 65_535 {
            assert!(!seqair::bam::aux::iter_tags(stored.aux()).any(|(t, _)| t == *b"CG"));
        }
    }
    // The reader's region test, from the raw bytes alone.
    let mut reader_store = RecordStore::new();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("r.bam");
    let mut w = seqair::io::BgzfWriter::new(Vec::new());
    let mut bam = b"BAM\x01".to_vec();
    let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{CONTIG_LEN}\n");
    bam.extend_from_slice(&(text.len() as i32).to_le_bytes());
    bam.extend_from_slice(text.as_bytes());
    bam.extend_from_slice(&1i32.to_le_bytes());
    bam.extend_from_slice(&2i32.to_le_bytes());
    bam.extend_from_slice(b"c\0");
    bam.extend_from_slice(&(CONTIG_LEN as i32).to_le_bytes());
    bam.extend_from_slice(&(raw.len() as i32).to_le_bytes());
    bam.extend_from_slice(&raw);
    w.write_all(&bam).unwrap();
    std::fs::write(&path, w.finish().unwrap()).unwrap();
    let status = Command::new("samtools").arg("index").arg(&path).status().unwrap();
    assert!(status.success());
    let mut reader = seqair::bam::IndexedBamReader::open(&path).unwrap();
    let tid = reader.header().tid("c").unwrap();
    // A region starting at the record's last base finds it; one past does not.
    let end = stored.end_pos.as_u32();
    for (start, found) in [(end, 1), (end + 1, 0)] {
        reader_store.clear();
        let span = (Pos0::new(start).unwrap()..=Pos0::new(start + 10).unwrap()).into();
        reader.fetch_into(tid, span, &mut reader_store).unwrap();
        assert_eq!(reader_store.len(), found, "region from {start}");
    }
}

/// One read's tags as a BAM aux block, and the readers' fields for them.
struct Fixture {
    _dir: tempfile::TempDir,
    text: String,
    sam: std::path::PathBuf,
    sam_gz: std::path::PathBuf,
    /// `None` when samtools cannot index it (it refuses the record).
    bam: Option<std::path::PathBuf>,
    /// `None` when htslib refuses the SAM (no CRAM can be made of it).
    cram: Option<std::path::PathBuf>,
    fasta: std::path::PathBuf,
}

fn run(cmd: &mut Command) -> bool {
    cmd.output().unwrap().status.success()
}

/// The fixture's contig: long enough for its reads, short enough for a
/// FASTA per case.
const FIXTURE_LEN: u32 = 1_000;

/// The same record as SAM text, bgzipped SAM, BAM and CRAM, the first three
/// keeping the placeholder and the tags exactly; the CRAM is samtools',
/// which resolves a `CG` tag it can as it reads the SAM.
fn fixture(line: &str, payload: &[u8]) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{FIXTURE_LEN}\n{line}\n");
    let sam = dir.path().join("r.sam");
    std::fs::write(&sam, &text).unwrap();
    let sam_gz = dir.path().join("r.sam.gz");
    let gz = Command::new("bgzip").arg("-c").arg(&sam).output().unwrap();
    assert!(gz.status.success());
    std::fs::write(&sam_gz, gz.stdout).unwrap();
    assert!(run(Command::new("tabix").args(["-p", "sam"]).arg(&sam_gz)));

    let bam = dir.path().join("r.bam");
    let mut bytes = b"BAM\x01".to_vec();
    let header_text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{FIXTURE_LEN}\n");
    bytes.extend_from_slice(&(header_text.len() as i32).to_le_bytes());
    bytes.extend_from_slice(header_text.as_bytes());
    bytes.extend_from_slice(&1i32.to_le_bytes());
    bytes.extend_from_slice(&2i32.to_le_bytes());
    bytes.extend_from_slice(b"c\0");
    bytes.extend_from_slice(&(FIXTURE_LEN as i32).to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as i32).to_le_bytes());
    bytes.extend_from_slice(payload);
    let mut w = seqair::io::BgzfWriter::new(Vec::new());
    w.write_all(&bytes).unwrap();
    std::fs::write(&bam, w.finish().unwrap()).unwrap();
    let indexed = run(Command::new("samtools").arg("index").arg(&bam));

    let fasta = dir.path().join("ref.fa");
    let reference: String =
        (0..FIXTURE_LEN as usize).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
    std::fs::write(&fasta, format!(">c\n{reference}\n")).unwrap();
    assert!(run(Command::new("samtools").arg("faidx").arg(&fasta)));
    let cram = dir.path().join("r.cram");
    let made = run(Command::new("samtools")
        .args(["view", "-C", "-T"])
        .arg(&fasta)
        .arg("-o")
        .arg(&cram)
        .arg(&sam))
        && run(Command::new("samtools").arg("index").arg(&cram));
    Fixture {
        _dir: dir,
        text,
        sam,
        sam_gz,
        bam: indexed.then_some(bam),
        cram: made.then_some(cram),
        fasta,
    }
}

/// What samtools makes of the SAM line: its CIGAR and tag names, or `None`
/// when it refuses the record. The oracle.
fn htslib_reads(sam: &std::path::Path) -> Option<(String, Vec<String>)> {
    let out = Command::new("samtools").arg("view").arg(sam).output().unwrap();
    let line = String::from_utf8(out.stdout).unwrap();
    let fields: Vec<&str> = line.trim_end().split('\t').collect();
    if !out.status.success() || fields.len() < 11 {
        return None;
    }
    let tags = fields[11..].iter().map(|f| f[..2].to_owned()).collect();
    Some((fields[5].to_owned(), tags))
}

// r[verify sam.reader.long_cigar]
// r[verify cram.record.long_cigar]
// r[verify record_store.long_cigar+1]
// r[verify bam.raw_record.cigar_seq]
/// The same record — a placeholder CIGAR, maybe a wrong one, its `CG` tag,
/// maybe one that does not cover SEQ, maybe another `CG` tag after it — read
/// by every reader into a store, with every tag, without `CG`, or with `CG`
/// among the tags asked for: they agree with samtools on the CIGAR and on
/// whether the record is read at all, find it in a region its real span
/// reaches, and keep a `CG` tag only when it is not the CIGAR and was asked
/// for.
#[hegel::test(test_cases = pinned::cases(30))]
fn every_reader_takes_a_cg_tag_alike(tc: TestCase) {
    use seqair::reader::{AuxTags, IndexedReader, RecordFields};
    let real = [
        CigarOp::new(CigarOpType::Match, tc.draw(gs::integers::<u32>().min_value(1).max_value(5))),
        CigarOp::new(
            CigarOpType::Deletion,
            tc.draw(gs::integers::<u32>().min_value(1).max_value(5)),
        ),
        CigarOp::new(CigarOpType::Match, 2),
    ];
    let qlen = seqair::bam::cigar::query_len(real) as u32;
    let span = seqair::bam::cigar::ref_len(real) as u32;
    // A CG that does not cover SEQ; or one too short to be taken.
    let covers = tc.draw(gs::weighted_booleans(0.8));
    let seq_len = if covers { qlen } else { qlen + 1 };
    let cg_ops: Vec<CigarOp> =
        if tc.draw(gs::weighted_booleans(0.15)) { real[..1].to_vec() } else { real.to_vec() };
    let n = tc.draw(gs::integers::<u32>().min_value(1).max_value(span + 2));
    let stale_after = tc.draw(gs::booleans());
    let tags = tc.draw(gs::sampled_from(vec![
        AuxTags::All,
        AuxTags::only([*b"NM"]),
        AuxTags::only([*b"NM", *b"CG"]),
    ]));
    let cg_text = |ops: &[CigarOp]| {
        let ops: Vec<String> = ops.iter().map(|o| o.to_bam_u32().to_string()).collect();
        format!("CG:B:I,{}", ops.join(","))
    };
    let stale = [CigarOp::new(CigarOpType::Match, 7)];
    let mut tag_text = format!("NM:i:0\t{}", cg_text(&cg_ops));
    let mut aux = b"NMC\x00".to_vec();
    aux.extend(cg_tag(b'I', &cg_ops));
    if stale_after {
        tag_text.push('\t');
        tag_text.push_str(&cg_text(&stale));
        aux.extend(cg_tag(b'I', &stale));
    }
    let seq: String = "A".repeat(seq_len as usize);
    let line = format!("r\t0\tc\t100\t30\t{seq_len}S{n}N\t*\t0\t0\t{seq}\t*\t{tag_text}");
    let placeholder =
        [CigarOp::new(CigarOpType::SoftClip, seq_len), CigarOp::new(CigarOpType::RefSkip, n)];
    let f = fixture(&line, &payload(99, &placeholder, seq_len as usize, &aux));

    let oracle = htslib_reads(&f.sam);
    let taken = cg_ops.len() >= 2;
    match &oracle {
        None => assert!(taken && !covers, "samtools refused {line}"),
        Some((cigar, _)) => {
            let want = if taken { CigarStr(&real).to_string() } else { format!("{seq_len}S{n}N") };
            assert_eq!(cigar, &want, "samtools on {line}");
        }
    }
    // A region from the real CIGAR's last base: only its span reaches it.
    let last = 99 + span - 1;
    let region = (Pos0::new(last).unwrap()..=Pos0::new(last + 5).unwrap()).into();
    let fields = RecordFields::all().with_aux(tags.clone());

    // Without an index for the BAM, its record still goes through the store.
    let raw_payload = payload(99, &placeholder, seq_len as usize, &aux);
    let pushed = RecordStore::new().push_raw(&raw_payload, &header(), &mut ()).map(|_| ());
    assert_eq!(pushed.is_ok(), oracle.is_some(), "push_raw: {pushed:?}");
    let mut readers: Vec<(&str, seqair::reader::PartialIndexedReader)> =
        vec![("sam", IndexedReader::builder(&f.sam_gz).fields(fields.clone()).open().unwrap())];
    if let Some(bam) = &f.bam {
        readers.push(("bam", IndexedReader::builder(bam).fields(fields.clone()).open().unwrap()));
    }
    if let Some(cram) = &f.cram {
        readers.push((
            "cram",
            IndexedReader::builder(cram).reference(&f.fasta).fields(fields.clone()).open().unwrap(),
        ));
    }
    for (what, mut reader) in readers {
        let tid = reader.header().tid("c").unwrap();
        let mut store = RecordStore::new();
        let fetched = reader.fetch_into(tid, region, &mut store);
        let Some((_, htslib_tags)) = &oracle else {
            assert!(fetched.is_err(), "{what}: a CG tag that does not cover SEQ");
            continue;
        };
        fetched.unwrap();
        if !taken {
            // The placeholder's own span decides, and may fall short.
            assert!(store.len() <= 1, "{what}");
            if store.is_empty() {
                continue;
            }
        }
        assert_eq!(store.len(), 1, "{what}: the real span reaches {last}");
        let rec = store.record(store.indices().next().unwrap()).unwrap();
        let want = if taken { real.as_slice() } else { placeholder.as_slice() };
        assert_eq!(rec.cigar(), want, "{what}: cigar");
        let stored_tags: Vec<[u8; 2]> =
            seqair::bam::aux::iter_tags(rec.aux()).map(|(t, _)| t).collect();
        let cg_kept = stored_tags.iter().filter(|t| **t == *b"CG").count();
        // Not the CIGAR, asked for: kept as htslib keeps it. Once one is the
        // CIGAR, none stays — but samtools wrote the CRAM from its own
        // reading, where a `CG` tag after the one it took is just a tag.
        let htslib_cg = htslib_tags.iter().filter(|t| t.as_str() == "CG").count();
        let want_cg = match () {
            () if !tags.contains(*b"CG") => 0,
            () if what == "cram" => htslib_cg,
            () if taken => 0,
            () => htslib_cg,
        };
        assert_eq!(cg_kept, want_cg, "{what}: CG tags kept");
        // (samtools' CRAM leaves out an `NM` it can recompute.)
        assert!(what == "cram" || stored_tags.contains(b"NM"), "{what}: NM");
    }

    // The SAM stream reads the text through the encoder.
    let mut stream = AlignmentStream::builder(std::io::Cursor::new(f.text.clone().into_bytes()))
        .build()
        .unwrap();
    match (stream.records().next().unwrap(), &oracle) {
        (Err(_), None) => {}
        (Ok(raw), Some(_)) => {
            let want = if taken { real.as_slice() } else { placeholder.as_slice() };
            assert_eq!(raw.cigar().to_vec(), want, "sam stream");
        }
        (got, want) => panic!("sam stream: {got:?} where samtools gave {want:?}"),
    }
}
