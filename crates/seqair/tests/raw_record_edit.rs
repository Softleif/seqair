//! `RawBamRecord::set_seq`, `set_base` and `edit_aux`: against a model of the
//! expected bytes, and against `samtools view` on the written file.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/raw_record_gen.rs"]
mod raw_record_gen;

use hegel::TestCase;
use hegel::generators as gs;
use raw_record_gen::{GenRecord, GenTag, HTSLIB, arb_record, array_tag, header};
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{AuxDataError, RawBamRecord, RawRecordEditError};
use seqair_types::{IupacBase, QPos};
use std::process::Command;

const EDITED_KEYS: [[u8; 2]; 6] = [*b"MM", *b"ML", *b"XR", *b"XG", *b"XM", *b"CG"];

fn record_of(g: &GenRecord) -> RawBamRecord {
    RawBamRecord::from_payload(g.payload(), &header()).expect("generated records are valid")
}

fn iupac(tc: &TestCase, len: usize) -> Vec<IupacBase> {
    (0..len).map(|_| tc.draw_silent(gs::sampled_from(IupacBase::ALL.to_vec()))).collect()
}

fn decoded(record: &RawBamRecord) -> Vec<IupacBase> {
    record.iupac_bases().collect()
}

// r[verify bam.raw_record.edit_seq]
#[hegel::test(test_cases = pinned::cases(200))]
fn set_seq_and_set_base_change_only_the_bases(tc: TestCase) {
    let g: GenRecord = tc.draw(arb_record(HTSLIB));
    let mut record = record_of(&g);
    let len = record.seq_len();
    let original = record.as_payload().to_vec();
    let packed_start = original.len()
        - g.tags.iter().map(|t| 2 + t.value.len()).sum::<usize>()
        - g.qual.len()
        - len.div_ceil(2);
    let packed_end = packed_start + len.div_ceil(2);

    let new = iupac(&tc, len);
    record.set_seq(&new).unwrap();
    assert_eq!(decoded(&record), new);
    // Everything but the packed bases, the odd pad nibble included, is as it was.
    assert_eq!(record.as_payload()[..packed_start], original[..packed_start]);
    assert_eq!(record.as_payload()[packed_end..], original[packed_end..]);
    if len % 2 == 1 {
        assert_eq!(record.as_payload()[packed_end - 1] & 15, original[packed_end - 1] & 15);
    }
    let bytes = record.as_payload().to_vec();
    RawBamRecord::from_payload(bytes, &header()).unwrap();

    // A different length is refused and changes nothing.
    let snapshot = record.as_payload().to_vec();
    let mut wrong = new.clone();
    wrong.push(IupacBase::A);
    assert!(matches!(record.set_seq(&wrong), Err(RawRecordEditError::SeqLenMismatch { .. })));
    assert_eq!(record.as_payload(), snapshot);

    // One base: only it moves.
    if len > 0 {
        let i = tc.draw(gs::integers::<usize>().max_value(len - 1));
        let base = tc.draw(gs::sampled_from(IupacBase::ALL.to_vec()));
        let before_base = record.as_payload().to_vec();
        record.set_base(QPos::new(i as u32), base).unwrap();
        // Only the edited base's byte may differ, and within it only its nibble.
        let byte_at = packed_start + i / 2;
        let after_base = record.as_payload();
        assert_eq!(after_base[..byte_at], before_base[..byte_at]);
        assert_eq!(after_base[byte_at + 1..], before_base[byte_at + 1..]);
        let kept_nibble = if i % 2 == 0 { 0x0F } else { 0xF0 };
        assert_eq!(after_base[byte_at] & kept_nibble, before_base[byte_at] & kept_nibble);
        if len % 2 == 1 && i == len - 1 {
            assert_eq!(after_base[byte_at] & 0x0F, before_base[byte_at] & 0x0F, "pad nibble");
        }
        let mut expected = new;
        expected[i] = base;
        assert_eq!(decoded(&record), expected);
        assert_eq!(record.as_payload().len(), snapshot.len());
    }
    assert!(matches!(
        record.set_base(QPos::new(len as u32), IupacBase::A),
        Err(RawRecordEditError::BaseOutOfRange { .. })
    ));
}

/// One aux edit the model can replay.
#[derive(Debug, Clone)]
enum Edit {
    Remove([u8; 2]),
    SetString([u8; 2], Vec<u8>),
    SetBytes([u8; 2], Vec<u8>),
}

fn draw_edits(tc: &TestCase, g: &GenRecord) -> Vec<Edit> {
    let mut keys: Vec<[u8; 2]> = EDITED_KEYS.to_vec();
    keys.extend(g.tags.iter().map(|t| t.key));
    let count = tc.draw(gs::integers::<usize>().max_value(5));
    (0..count)
        .map(|_| {
            let key = tc.draw_silent(gs::sampled_from(keys.clone()));
            match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
                0 => Edit::Remove(key),
                1 => {
                    let text = tc.draw_silent(
                        gs::vecs(gs::sampled_from(b"CT+?,0123".to_vec())).max_size(12),
                    );
                    Edit::SetString(key, text)
                }
                _ => {
                    Edit::SetBytes(key, tc.draw_silent(gs::vecs(gs::integers::<u8>()).max_size(12)))
                }
            }
        })
        .collect()
}

/// The tags after `edits`: a removal drops the first tag of that name, a set
/// drops it and appends the new one.
fn model_tags(tags: &[GenTag], edits: &[Edit]) -> Vec<GenTag> {
    let mut tags = tags.to_vec();
    for edit in edits {
        let (Edit::Remove(key) | Edit::SetString(key, _) | Edit::SetBytes(key, _)) = edit;
        if let Some(at) = tags.iter().position(|t| t.key == *key) {
            tags.remove(at);
        }
        match edit {
            Edit::Remove(_) => {}
            Edit::SetString(key, text) => {
                let mut value = vec![b'Z'];
                value.extend_from_slice(text);
                value.push(0);
                tags.push(GenTag { key: *key, value });
            }
            Edit::SetBytes(key, bytes) => {
                tags.push(array_tag(*key, b'C', bytes, bytes.len() as u32));
            }
        }
    }
    tags
}

fn apply(record: &mut RawBamRecord, edits: &[Edit]) -> Result<(), RawRecordEditError> {
    record.edit_aux(|aux| {
        for edit in edits {
            match edit {
                Edit::Remove(key) => aux.remove(*key),
                Edit::SetString(key, text) => aux.set_string(*key, text),
                Edit::SetBytes(key, bytes) => aux.set_array_u8(*key, bytes)?,
            }
        }
        Ok::<(), AuxDataError>(())
    })
}

// r[verify bam.raw_record.edit_aux]
#[hegel::test(test_cases = pinned::cases(200))]
fn edit_aux_keeps_untouched_tags_byte_for_byte(tc: TestCase) {
    let g: GenRecord = tc.draw(arb_record(HTSLIB));
    let mut record = record_of(&g);
    let original = record.as_payload().to_vec();
    let edits = draw_edits(&tc, &g);
    let expected = GenRecord { tags: model_tags(&g.tags, &edits), ..g.clone() };
    let header = header();

    // The edit is refused exactly when the real CIGAR would differ.
    let expected_record = RawBamRecord::from_payload(expected.payload(), &header);
    let same_cigar = expected_record.as_ref().is_ok_and(|e| {
        e.has_cg_placeholder() == record.has_cg_placeholder()
            && e.cigar().collect::<Vec<_>>() == record.cigar().collect::<Vec<_>>()
    });
    match apply(&mut record, &edits) {
        Ok(()) => {
            assert!(same_cigar);
            assert_eq!(record.as_payload(), expected.payload());
            let reread = RawBamRecord::from_payload(record.into_payload(), &header).unwrap();
            assert_eq!(reread.as_payload(), expected.payload());
        }
        Err(RawRecordEditError::LongCigarChanged) => {
            assert!(!same_cigar);
            assert_eq!(record.as_payload(), original);
        }
        Err(e) => panic!("{e}"),
    }
}

// r[verify bam.raw_record.edit_aux]
#[test]
fn a_failing_aux_edit_leaves_the_record_alone() {
    let mut record = RawBamRecord::default();
    let before = record.as_payload().to_vec();
    let err = record
        .edit_aux(|aux| {
            aux.set_string(*b"XM", b"abc");
            aux.set_char(*b"XR", 0x01)
        })
        .unwrap_err();
    assert!(matches!(err, RawRecordEditError::Aux(AuxDataError::InvalidCharByte { .. })));
    assert_eq!(record.as_payload(), before);
}

// r[verify bam.raw_record.edit_aux]
#[test]
fn the_cg_tag_of_a_long_cigar_cannot_be_removed_or_changed() {
    let (m, s, n) = (raw_record_gen::op(2, 0), raw_record_gen::op(4, 4), raw_record_gen::op(5, 3));
    let real: Vec<u8> =
        [m, raw_record_gen::op(1, 2), m].iter().flat_map(|o| o.to_le_bytes()).collect();
    let g = GenRecord {
        refid: 0,
        pos: 10,
        mapq: 1,
        bin: 0,
        flags: 0,
        next_refid: -1,
        next_pos: -1,
        tlen: 0,
        qname: b"r".to_vec(),
        cigar: vec![s, n],
        seq: vec![1, 2, 4, 8],
        seq_pad: 0,
        qual: vec![30; 4],
        tags: vec![array_tag(*b"CG", b'I', &real, 3)],
    };
    let mut record = record_of(&g);
    assert!(record.has_cg_placeholder());
    let before = record.as_payload().to_vec();
    let removed = record.edit_aux(|aux| {
        aux.remove(*b"CG");
        Ok::<(), AuxDataError>(())
    });
    assert!(matches!(removed, Err(RawRecordEditError::LongCigarChanged)));
    let changed = record.edit_aux(|aux| aux.set_array_u32(*b"CG", &[m, m]));
    assert!(matches!(changed, Err(RawRecordEditError::LongCigarChanged)));
    assert_eq!(record.as_payload(), before);
    // Same ops under another subtype are the same CIGAR.
    let ops: Vec<i32> = [m, raw_record_gen::op(1, 2), m].iter().map(|&o| o as i32).collect();
    record.edit_aux(|aux| aux.set_array_i32(*b"CG", &ops)).unwrap();
    assert_eq!(record.cigar().count(), 3);
}

fn samtools_line(path: &std::path::Path) -> Vec<Vec<u8>> {
    let out = Command::new("samtools").arg("view").arg(path).output().expect("samtools");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let line = out.stdout.strip_suffix(b"\n").expect("one record");
    line.split(|&b| b == b'\t').map(<[u8]>::to_vec).collect()
}

fn write_bam(path: &std::path::Path, record: &RawBamRecord) {
    let header = header();
    let mut writer = BamWriterBuilder::to_path(path, &header).build().unwrap();
    writer.write_raw(record).unwrap();
    writer.finish().unwrap();
}

// r[verify bam.raw_record.edit_seq]
// r[verify bam.raw_record.edit_aux]
/// samtools 1.25 reads the edited record: SEQ is the new bases, the other
/// fixed columns and QUAL are as before, and the tags are the old ones minus
/// the replaced, plus the new ones, as samtools prints them.
#[hegel::test(test_cases = pinned::cases(40))]
fn samtools_reads_the_edited_record(tc: TestCase) {
    let g: GenRecord = tc.draw(arb_record(HTSLIB));
    // CG edits are the model test's business; samtools would also move a CG.
    // Tabs and newlines in names or text tags would break the line split below.
    let clean = |b: &u8| if matches!(b, 9 | 10) { b'x' } else { *b };
    let g = GenRecord {
        tags: g.tags.into_iter().filter(|t| t.key != *b"CG").collect(),
        qname: g.qname.iter().map(clean).collect(),
        ..g
    };
    let mut record = record_of(&g);
    let dir = tempfile::tempdir().unwrap();
    let (before_path, after_path) = (dir.path().join("a.bam"), dir.path().join("b.bam"));
    write_bam(&before_path, &record);
    let before = samtools_line(&before_path);

    let new_seq = iupac(&tc, record.seq_len());
    record.set_seq(&new_seq).unwrap();
    let ml: Vec<u8> = tc.draw(gs::vecs(gs::integers::<u8>()).max_size(8));
    record
        .edit_aux(|aux| {
            for key in [*b"XM", *b"XR", *b"MM"] {
                aux.remove(key);
            }
            aux.set_string(*b"MM", b"C+m,0;");
            aux.set_array_u8(*b"ML", &ml)?;
            aux.set_string(*b"XG", b"CT");
            Ok::<(), AuxDataError>(())
        })
        .unwrap();
    write_bam(&after_path, &record);
    let after = samtools_line(&after_path);

    assert_eq!(after[..9], before[..9]);
    assert_eq!(after[10], before[10], "QUAL");
    let new_ascii = IupacBase::to_ascii_slice(&new_seq);
    let want_seq: &[u8] = if new_ascii.is_empty() { b"*" } else { &new_ascii };
    assert_eq!(after[9], want_seq);

    // Each removed or set key loses its first occurrence in samtools' own
    // rendering of the original.
    let mut tags = before[11..].to_vec();
    for key in [*b"XM", *b"XR", *b"MM", *b"ML", *b"XG"] {
        if let Some(at) = tags.iter().position(|t| t.starts_with(&key)) {
            tags.remove(at);
        }
    }
    let mut ml_text = b"ML:B:C".to_vec();
    for b in &ml {
        ml_text.extend_from_slice(format!(",{b}").as_bytes());
    }
    tags.extend([b"MM:Z:C+m,0;".to_vec(), ml_text, b"XG:Z:CT".to_vec()]);
    assert_eq!(after[11..], tags);
}

// r[verify types.iupac_base.from_ascii]
// r[verify bam.raw_record.iupac_bases]
/// samtools encodes every printable SEQ character (all IUPAC letters in both
/// cases, `U`, `.`, digits, and the rest) to the base `IupacBase::from_ascii`
/// names.
#[test]
fn from_ascii_agrees_with_htslib_on_every_sam_seq_character() {
    // SAM's SEQ grammar is `\*|[A-Za-z=.]+`; the other printable bytes are
    // outside it, but htslib encodes them too, as N.
    let chars: Vec<u8> = (33u8..127).filter(|&b| b != b'*').collect();
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("a.sam");
    let line = format!(
        "r\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\n",
        String::from_utf8(chars.clone()).unwrap(),
        "I".repeat(chars.len())
    );
    std::fs::write(&sam, line).unwrap();
    let bam = dir.path().join("a.bam");
    let out =
        Command::new("samtools").args(["view", "-b", "-o"]).arg(&bam).arg(&sam).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mut reader =
        seqair::bam::stream::BamStreamReader::new(std::fs::File::open(&bam).unwrap()).unwrap();
    let record = reader.records().next().expect("one record").unwrap();
    assert_eq!(record.iupac_bases().collect::<Vec<_>>(), IupacBase::from_ascii_slice(&chars));
}
