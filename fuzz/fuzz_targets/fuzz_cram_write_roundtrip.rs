//! Reads drawn against a drawn reference, written as CRAM by seqair and
//! read back by seqair's stream decoder: every record the writer accepts
//! reads back as it was, but for the normalizations `r[cram_writer.fidelity]`
//! lists, and writing what was read gives the same records again.
#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::RawBamRecord;
use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::mate::MateFields;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair_fuzz::cram_roundtrip;
use seqair_types::{BamFlags, BaseQuality, IupacBase, MappingQuality, Pos0};
use std::num::NonZeroI32;

/// Reference letters, soft-masked and ambiguous ones included.
const REF_LETTERS: &[u8] = b"ACGTACGTACGTNacgtnRYKMU";
/// Read bases: the reference's mostly, else any IUPAC code.
const READ_LETTERS: &[u8] = b"=ACMGRSVTWYHKDBN";
const OPS: [CigarOpType; 9] = [
    CigarOpType::Match,
    CigarOpType::Insertion,
    CigarOpType::Deletion,
    CigarOpType::RefSkip,
    CigarOpType::SoftClip,
    CigarOpType::HardClip,
    CigarOpType::Padding,
    CigarOpType::SeqMatch,
    CigarOpType::SeqMismatch,
];

#[derive(Arbitrary, Debug)]
enum Value {
    Int(i32),
    Char(u8),
    Float(f32),
    Str(Vec<u8>),
    Hex(Vec<u8>),
    Shorts(Vec<i16>),
    Bytes(Vec<u8>),
}

#[derive(Arbitrary, Debug)]
struct Tag {
    name: [u8; 2],
    value: Value,
}

#[derive(Arbitrary, Debug)]
struct Read {
    /// 0, 1: a contig in the FASTA; 2: one it lacks; else unplaced.
    contig: u8,
    pos: u16,
    unmapped: bool,
    flags: u16,
    mapq: u8,
    ops: Vec<(u8, u8)>,
    /// Per base: under 200 copies the reference, else picks a read letter.
    bases: Vec<u8>,
    quals: Option<Vec<u8>>,
    name: Vec<u8>,
    tags: Vec<Tag>,
    mate: Option<(u8, u16, i32)>,
}

#[derive(Arbitrary, Debug)]
struct Input {
    contigs: [Vec<u8>; 2],
    reads: Vec<Read>,
    threads: u8,
    records_per_slice: u8,
}

fn tag_into(aux: &mut AuxData, tag: &Tag) {
    // Valid names only (the writer refuses others; the BAM side checks too).
    let name = [b'X', tag.name[1] % 26 + b'A'];
    let printable = |v: &[u8]| v.iter().map(|b| b % 94 + b'!').collect::<Vec<u8>>();
    match &tag.value {
        Value::Int(v) => drop(aux.set_int(name, i64::from(*v))),
        Value::Char(c) => drop(aux.set_char(name, c % 94 + b'!')),
        Value::Float(f) => aux.set_float(name, *f),
        Value::Str(s) => drop(aux.set_string(name, &printable(s))),
        Value::Hex(h) => {
            let hex: Vec<u8> = h.iter().map(|b| b"0123456789ABCDEF"[usize::from(b % 16)]).collect();
            let even = &hex[..hex.len() & !1];
            drop(aux.set_hex(name, even));
        }
        Value::Shorts(v) => drop(aux.set_array_i16(name, v)),
        Value::Bytes(v) => drop(aux.set_array_u8(name, v)),
    }
}

fn build(read: &Read, refs: &[Vec<u8>; 2], header: &BamHeader, i: usize) -> Option<RawBamRecord> {
    let tid_of = |c: u8| ["c0", "c1", "missing"].get(usize::from(c)).and_then(|n| header.tid(n));
    let tid = tid_of(read.contig);
    let len = match read.contig {
        0 | 1 => refs[usize::from(read.contig)].len(),
        2 => 500,
        _ => 0,
    };
    let pos = tid.map(|_| u32::from(read.pos) % len.max(1) as u32).and_then(Pos0::new);
    let unmapped = read.unmapped || tid.is_none();
    let mut flags = read.flags & 0xfff & !0x4;
    if unmapped {
        flags |= 0x4;
    }

    let reference = refs.get(usize::from(read.contig));
    let mut cigar = Vec::new();
    let mut seq = Vec::new();
    let mut base_stream = read.bases.iter().copied().cycle();
    let mut refpos = pos.map_or(0, |p| p.as_usize());
    if !unmapped {
        for &(op, n) in read.ops.iter().take(12) {
            let op = OPS[usize::from(op) % OPS.len()];
            let n = u32::from(n % 40);
            cigar.push(CigarOp::new(op, n));
            let consumes_ref = matches!(
                op,
                CigarOpType::Match
                    | CigarOpType::SeqMatch
                    | CigarOpType::SeqMismatch
                    | CigarOpType::Deletion
                    | CigarOpType::RefSkip
            );
            let consumes_read = matches!(
                op,
                CigarOpType::Match
                    | CigarOpType::SeqMatch
                    | CigarOpType::SeqMismatch
                    | CigarOpType::Insertion
                    | CigarOpType::SoftClip
            );
            for _ in 0..n {
                if consumes_read {
                    let pick = base_stream.next().unwrap_or(0);
                    let r = reference.and_then(|r| r.get(refpos)).filter(|_| consumes_ref);
                    let b = match r {
                        Some(&r) if pick < 200 => IupacBase::from_ascii(r),
                        _ => IupacBase::from_ascii(READ_LETTERS[usize::from(pick) % 16]),
                    };
                    seq.push(b);
                }
                if consumes_ref {
                    refpos += 1;
                }
            }
        }
    } else {
        let n = read.bases.len().min(60);
        seq.extend(read.bases.iter().take(n).map(|&b| IupacBase::from_ascii(READ_LETTERS[usize::from(b) % 16])));
    }
    if read.bases.first() == Some(&255) {
        seq.clear(); // SEQ `*`
    }
    let qual = match &read.quals {
        Some(q) if !seq.is_empty() && !q.is_empty() => {
            (0..seq.len()).map(|j| BaseQuality::from_byte(q[j % q.len()] % 94)).collect()
        }
        _ => Vec::new(),
    };
    let mut aux = AuxData::new();
    for tag in read.tags.iter().take(6) {
        tag_into(&mut aux, tag);
    }
    let mate = read.mate.map_or_else(MateFields::default, |(c, p, t)| MateFields {
        tid: tid_of(c),
        pos: Pos0::new(u32::from(p)),
        template_len: NonZeroI32::new(t),
    });
    let mut name = format!("r{i}_").into_bytes();
    name.extend(read.name.iter().take(20).map(|b| b % 94 + b'!').filter(|&b| b != b'@'));
    let record = OwnedBamRecord {
        tid,
        pos,
        mapq: MappingQuality::from_byte(read.mapq),
        flags: BamFlags::from(flags),
        mate_fields: mate,
        qname: name,
        cigar,
        seq,
        qual,
        aux,
    };
    let mut payload = Vec::new();
    record.to_bam_bytes(&mut payload).ok()?;
    RawBamRecord::from_payload(payload, header).ok()
}

fuzz_target!(|input: Input| {
    let refs: [Vec<u8>; 2] = input.contigs.map(|c| {
        let mut bases: Vec<u8> =
            c.iter().take(2000).map(|b| REF_LETTERS[usize::from(*b) % REF_LETTERS.len()]).collect();
        if bases.is_empty() {
            bases.push(b'A');
        }
        bases
    });
    let header_text = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c0\tLN:{}\n@SQ\tSN:c1\tLN:{}\n@SQ\tSN:missing\tLN:500\n@RG\tID:a\n@RG\tID:b\n",
        refs[0].len(),
        refs[1].len()
    );
    let header = BamHeader::from_sam_text(&header_text).unwrap();
    let mut records: Vec<RawBamRecord> = input
        .reads
        .iter()
        .take(64)
        .enumerate()
        .filter_map(|(i, r)| build(r, &refs, &header, i))
        .collect();
    records.sort_by_key(|r| (r.tid().is_none(), r.tid(), r.tid().and(r.pos())));

    let dir = tempfile::tempdir().unwrap();
    let contigs = [("c0".to_string(), refs[0].clone()), ("c1".to_string(), refs[1].clone())];
    let fasta = cram_roundtrip::write_fasta(dir.path(), &contigs);
    let rps = u32::from(input.records_per_slice % 8).max(1) * if input.records_per_slice > 128 { 10_000 } else { 1 };
    cram_roundtrip::check(&header, &fasta, &records, usize::from(input.threads % 3), rps);
});
