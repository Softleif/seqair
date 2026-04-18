//! Oracle tests: compare [`seqair::bam::base_mod::BaseModState`] against
//! htslib's `hts_base_mod_state` via rust-htslib. Round-trips a BAM record
//! with known MM/ML tags (including a reverse-strand case) and confirms that
//! both libraries resolve the same (qpos, mod-code, probability) tuples.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test code with known small values"
)]

use rust_htslib::bam::{self, Read as _};
use seqair::bam::aux_data::AuxData;
// Use the top-level `bam` re-export rather than the inner module path so the
// public-API surface from `seqair::bam` is exercised. The inner path
// `seqair::bam::base_mod::*` remains valid; both must compile.
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriter;
use seqair::bam::{BaseModState, ModType};
use seqair_types::{BamFlags, Base, BaseQuality};
const FLAG_REVERSE: u16 = 0x10;
use std::path::Path;

/// A single call, normalized for comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Call {
    qpos: u32,
    /// Ascii code for single-char modifications, negative for `-ChEBI` id.
    mod_code: i32,
    /// htslib uses `-1` for missing; seqair always has a `u8` probability.
    probability: i32,
}

fn make_header() -> BamHeader {
    BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n").unwrap()
}

fn write_bam(dir: &Path, records: &[OwnedBamRecord]) -> std::path::PathBuf {
    let header = make_header();
    let bam_path = dir.join("mods.bam");
    let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();
    for rec in records {
        writer.write(rec).unwrap();
    }
    let (_inner, index_builder) = writer.finish().unwrap();
    if let Some(ib) = index_builder {
        let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
        ib.write_bai(bai_file, header.target_count()).unwrap();
    }
    bam_path
}

fn reverse_complement(seq: &[Base]) -> Vec<Base> {
    seq.iter().rev().map(|b| b.inverse()).collect()
}

fn build_record(
    pos: i64,
    qname: &[u8],
    seq: Vec<Base>,
    cigar: Vec<CigarOp>,
    mm: &[u8],
    ml: &[u8],
    is_reverse: bool,
) -> OwnedBamRecord {
    let mut aux = AuxData::new();
    aux.set_string(*b"MM", mm);
    aux.set_array_u8(*b"ML", ml).unwrap();
    let flags = if is_reverse { BamFlags::from(FLAG_REVERSE) } else { BamFlags::empty() };
    let seq_len = seq.len();
    OwnedBamRecord::builder(0, pos, qname.to_vec())
        .flags(flags)
        .mapq(60)
        .cigar(cigar)
        .seq(seq)
        .qual(vec![BaseQuality::from_byte(30); seq_len])
        .aux(aux)
        .build()
        .unwrap()
}

fn htslib_calls(bam_path: &Path, qname: &[u8]) -> Vec<Call> {
    let mut reader = bam::Reader::from_path(bam_path).expect("htslib open");
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.expect("htslib record");
        if rec.qname() != qname {
            continue;
        }
        let iter = rec.basemods_position_iter().expect("basemods iter");
        for item in iter {
            let (qpos, mods) = item.expect("basemods item");
            for m in mods {
                out.push(Call {
                    qpos: qpos as u32,
                    mod_code: m.modified_base,
                    probability: m.qual,
                });
            }
        }
    }
    out.sort();
    out
}

/// Convert a [`ModType`] into the signed integer code htslib reports.
///
/// htslib uses negative integers to represent `ChEBI` ids (negated) and the
/// ASCII byte for single-character codes. Returns `None` if the `ChEBI` id
/// cannot be represented as a positive `i32` (i.e. `id > i32::MAX`), which
/// would otherwise wrap silently and cause a wrong oracle comparison.
fn mod_code_for(mt: ModType) -> Option<i32> {
    match mt {
        ModType::Code(c) => Some(i32::from(c)),
        ModType::ChEBI(n) => i32::try_from(n).ok().map(|v| -v),
    }
}

fn seqair_calls(state: &BaseModState, seq_len: usize) -> Vec<Call> {
    let mut out = Vec::new();
    for qp in 0..seq_len {
        if let Some(mods) = state.mod_at_qpos(qp) {
            for m in mods {
                let mod_code = mod_code_for(m.mod_type)
                    .unwrap_or_else(|| panic!("`ChEBI` id does not fit in i32: {:?}", m.mod_type));
                // htslib reports qual=-1 when ML entry is missing. We always
                // have a probability (u8) so report it as-is.
                out.push(Call { qpos: qp as u32, mod_code, probability: i32::from(m.probability) });
            }
        }
    }
    out.sort();
    out
}

#[test]
fn mod_code_for_rejects_chebi_overflow() {
    // `ChEBI` id > i32::MAX would silently wrap with `-(n as i32)` and produce
    // a wrong oracle comparison. The new helper returns None instead.
    let huge = u32::MAX;
    assert!(mod_code_for(ModType::ChEBI(huge)).is_none(), "huge `ChEBI` id must not overflow");
    let just_over = (i32::MAX as u32).saturating_add(1);
    assert!(mod_code_for(ModType::ChEBI(just_over)).is_none(), "id > i32::MAX must be rejected");
    // Boundary: i32::MAX is representable.
    assert_eq!(mod_code_for(ModType::ChEBI(i32::MAX as u32)), Some(-i32::MAX));
    // Single-char code passes through.
    assert_eq!(mod_code_for(ModType::Code(b'm')), Some(i32::from(b'm')));
}

// r[verify base_mod.parse_mm]
// r[verify base_mod.resolve_positions]
/// Forward-strand read: seqair's resolved calls match htslib's.
#[test]
fn forward_strand_matches_htslib() {
    let dir = tempfile::tempdir().unwrap();
    // Sequence contains 4 C's at qpos 1, 5, 7, 13.
    //   A C G T A C G C G A C G T C G T
    //   0 1 2 3 4 5 6 7 8 9 ...
    let seq = vec![
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::C,
        Base::G,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::C,
        Base::G,
        Base::T,
    ];
    let rec = build_record(
        100,
        b"fwd",
        seq.clone(),
        vec![CigarOp::new(CigarOpType::Match, seq.len() as u32)],
        b"C+m,0,1,0;",
        &[200, 180, 230],
        false,
    );
    let bam_path = write_bam(dir.path(), &[rec]);

    let oracle = htslib_calls(&bam_path, b"fwd");
    assert!(!oracle.is_empty(), "htslib produced no calls");

    let state = BaseModState::parse(b"C+m,0,1,0;", &[200, 180, 230], &seq, false).unwrap();
    let ours = seqair_calls(&state, seq.len());
    assert_eq!(ours, oracle, "seqair vs htslib disagree on forward strand");
}

// r[verify base_mod.reverse_complement]
/// Reverse-strand read: MM positions are relative to the original
/// (unreversed) sequence. seqair must handle the reverse-complement mapping
/// and produce the same per-stored-qpos calls as htslib.
#[test]
fn reverse_strand_matches_htslib() {
    let dir = tempfile::tempdir().unwrap();
    // Original (unreversed) 5'→3' sequence with C's at original pos 1, 5, 7, 13.
    let original = vec![
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::C,
        Base::G,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::C,
        Base::G,
        Base::T,
    ];
    // BAM stores the sequence reverse-complemented for reverse-strand alignments.
    let stored = reverse_complement(&original);
    let rec = build_record(
        100,
        b"rev",
        stored.clone(),
        vec![CigarOp::new(CigarOpType::Match, stored.len() as u32)],
        b"C+m,0,1,0;",
        &[200, 180, 230],
        true, // is_reverse
    );
    let bam_path = write_bam(dir.path(), &[rec]);

    let oracle = htslib_calls(&bam_path, b"rev");
    assert!(!oracle.is_empty(), "htslib produced no calls");

    let state = BaseModState::parse(b"C+m,0,1,0;", &[200, 180, 230], &stored, true).unwrap();
    let ours = seqair_calls(&state, stored.len());
    assert_eq!(ours, oracle, "seqair vs htslib disagree on reverse strand");
}

// r[verify base_mod.parse_mm]
/// Combined modification codes (`C+mh,…`): htslib agrees with seqair on the
/// interleaved ML assignment.
#[test]
fn combined_codes_match_htslib() {
    let dir = tempfile::tempdir().unwrap();
    // 5 C's at qpos 1, 3, 5, 7, 9.
    let seq = vec![
        Base::A,
        Base::C,
        Base::A,
        Base::C,
        Base::A,
        Base::C,
        Base::A,
        Base::C,
        Base::A,
        Base::C,
    ];
    let rec = build_record(
        100,
        b"combo",
        seq.clone(),
        vec![CigarOp::new(CigarOpType::Match, seq.len() as u32)],
        b"C+mh,0,1;",
        &[200, 100, 220, 120],
        false,
    );
    let bam_path = write_bam(dir.path(), &[rec]);
    let oracle = htslib_calls(&bam_path, b"combo");
    assert_eq!(oracle.len(), 4, "expect 4 mod calls (2 positions × 2 codes), got {oracle:?}");

    let state = BaseModState::parse(b"C+mh,0,1;", &[200, 100, 220, 120], &seq, false).unwrap();
    let ours = seqair_calls(&state, seq.len());
    assert_eq!(ours, oracle);
}
