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

/// `n` ops alternating `kM` and `1I` (or `1D`), starting with a match, so the
/// first op is never the placeholder's `<l_seq>S`.
#[hegel::composite]
fn cigar(tc: &TestCase) -> Vec<CigarOp> {
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

/// Tags, with a stale `CG` tag among them when `stale`.
fn aux(stale: bool) -> AuxData {
    let mut aux = AuxData::new();
    aux.set_int(*b"NM", 3).unwrap();
    if stale {
        aux.set_array_u32(*b"CG", &[CigarOp::new(CigarOpType::Match, 7).to_bam_u32()]).unwrap();
    }
    aux.set_string(*b"RG", b"g1").unwrap();
    aux
}

fn record(cigar: Vec<CigarOp>, aux: AuxData) -> OwnedBamRecord {
    let qlen = seqair::bam::cigar::query_len(&cigar) as usize;
    OwnedBamRecord::builder(header().tid("c"), Some(Pos0::new(99).unwrap()), b"r".to_vec())
        .cigar(cigar)
        .seq(vec![IupacBase::A; qlen])
        .qual(vec![BaseQuality::from_byte(30); qlen])
        .aux(aux)
        .build()
        .unwrap()
}

/// The tags as our reader shows them: `aux` without `CG` when the record's
/// CIGAR needed one.
fn expected_tags(aux: &AuxData, long: bool) -> Vec<u8> {
    let mut want = aux.clone();
    if long {
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
/// the real CIGAR and the record's tags, a stale `CG` dropped exactly when
/// the CIGAR needed the tag.
#[hegel::test(test_cases = pinned::cases(12))]
fn long_cigars_survive_every_writer(tc: TestCase) {
    let ops = tc.draw(cigar());
    let stale = tc.draw(gs::booleans());
    let aux = aux(stale);
    let long = ops.len() > usize::from(u16::MAX);
    let owned = record(ops.clone(), aux.clone());
    let header = header();

    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    writer.write(&owned).unwrap();
    let (from_owned, _) = writer.finish().unwrap();

    // Read it into a store and write the store's record.
    let mut stream =
        AlignmentStream::builder(std::io::Cursor::new(from_owned.clone())).build().unwrap();
    let raw = stream.records().next().unwrap().unwrap();
    assert_eq!(raw.has_cg_placeholder(), long);
    assert_eq!(raw.cigar().to_vec(), ops);
    assert_eq!(raw.aux().as_bytes(), expected_tags(&aux, long));
    let mut store = RecordStore::new();
    let idx = store.push_raw(raw.as_payload(), &header, &mut ()).unwrap().unwrap();
    assert_eq!(store.record(idx).unwrap().cigar(), ops.as_slice());
    assert_eq!(store.record(idx).unwrap().aux(), expected_tags(&aux, long));
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    writer.write_store_record(&store, idx).unwrap();
    let (from_store, _) = writer.finish().unwrap();

    for (bam, what) in [(&from_owned, "owned"), (&from_store, "store")] {
        let fields = samtools_fields(bam);
        assert_eq!(fields[5], CigarStr(&ops).to_string(), "{what}: samtools CIGAR");
        let cg_tags = fields[11..].iter().filter(|f| f.starts_with("CG:")).count();
        assert_eq!(cg_tags, usize::from(stale && !long), "{what}: CG tags samtools shows");
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
