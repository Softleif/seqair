//! NM and MD tag computation by walking a record's CIGAR against a reference.
//!
//! Per [SAM1] §1.5 (auxiliary fields), the `NM` tag is the *edit distance*
//! between the read and the reference: number of mismatches in match
//! operations + total inserted bases + total deleted bases. The `MD` tag is a
//! string encoding the reference bases at mismatch positions, used to
//! reconstruct the reference without re-fetching it.
//!
//! Both are derivable from [`AlignedPairsWithRef`](super::aligned_pairs_view::AlignedPairsWithRef).
//! The methods on `AlignedPairsWithRef` consume the iterator and walk it to
//! completion.
//!
//! ## NM semantics
//!
//! - `M` (ambiguous): contribute 1 per position where `query != ref_base`
//!   (when the reference is loaded). Positions where the reference is outside
//!   the loaded window are skipped — `nm()` is best-effort.
//! - `=` (sequence match): always contribute 0.
//! - `X` (sequence mismatch): always contribute 1 per position, regardless of
//!   `ref_base` (the aligner already committed).
//! - `I` (insertion): contribute the insertion length.
//! - `D` (deletion): contribute the deletion length.
//! - `N` (reference skip): contribute 0 (introns are not edits).
//! - `S` / `H` / `P`: contribute 0.
//!
//! Comparing to `Base::Unknown`: a query `N` aligned against a reference
//! `A` counts as a mismatch (htslib behavior). Two `N`s also count as a
//! mismatch (they're not provably equal).
//!
//! ## MD semantics
//!
//! Per [SAM1] §1.5 (MD tag):
//! - Numbers: count of consecutive matching reference positions.
//! - Letters: a single ref base at a mismatch position.
//! - `^` followed by ref bases: a deletion run (the deleted ref bases).
//!
//! Insertions and soft clips are NOT in MD. Reference skips (N) are
//! intentionally NOT in MD either — different aligners handle this
//! inconsistently and we follow the conservative path of treating N as a
//! reference-coordinate skip (same as a deletion in span but no `^` prefix).
//! This matches what the htslib `samtools calmd` tool emits for typical
//! RNA-seq reads.
//!
//! Because MD encodes reference bases verbatim, missing reference data
//! yields an incorrect output — `md()` returns [`NmMdError::MissingReference`]
//! rather than silently producing a wrong tag.

use super::aligned_pairs::MatchKind;
use super::aligned_pairs_view::{AlignedPairWithRef, AlignedPairsWithRef};
use seqair_types::Base;
use thiserror::Error;

/// Errors from `md()` (and any future strict NM-style computations).
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum NmMdError {
    /// The reference window did not cover a position required by the CIGAR.
    /// MD encodes ref bases verbatim, so a missing position means the result
    /// would be wrong; we surface this rather than silently corrupt output.
    #[error(
        "reference base missing for MD computation at rpos={rpos}; either the loaded RefSeq \
         doesn't cover this position or the CIGAR walks past the end"
    )]
    MissingReference { rpos: u32 },
}

impl<'cigar, 'read, 'ref_seq> AlignedPairsWithRef<'cigar, 'read, 'ref_seq> {
    // r[impl cigar.aligned_pairs.nm.compute]
    /// Compute the `NM` tag (edit distance) by consuming this iterator.
    ///
    /// Best-effort: positions where the reference is outside the loaded
    /// `RefSeq` window are silently skipped (count 0), since `NM` is a
    /// downstream metric and partial coverage is informative even if not
    /// authoritative. For strict computation, ensure your `RefSeq` covers the
    /// record's full reference span before calling.
    ///
    /// # Example
    /// ```ignore
    /// let nm = slim
    ///     .aligned_pairs_with_read(&store)?
    ///     .with_reference(&ref_seq)
    ///     .nm();
    /// ```
    pub fn nm(self) -> u32 {
        let mut nm: u32 = 0;
        for ev in self {
            match ev {
                AlignedPairWithRef::Match { kind, query, ref_base, .. } => {
                    nm = nm.saturating_add(match_contribution(kind, query, ref_base));
                }
                AlignedPairWithRef::Insertion { query, .. } => {
                    // I contributes its length to edit distance.
                    let len = u32::try_from(query.len()).unwrap_or(u32::MAX);
                    nm = nm.saturating_add(len);
                }
                AlignedPairWithRef::Deletion { del_len, .. } => {
                    nm = nm.saturating_add(del_len);
                }
                AlignedPairWithRef::RefSkip { .. }
                | AlignedPairWithRef::SoftClip { .. }
                | AlignedPairWithRef::Padding { .. }
                | AlignedPairWithRef::Unknown { .. } => {
                    // Don't contribute to edit distance.
                }
            }
        }
        nm
    }

    // r[impl cigar.aligned_pairs.md.compute]
    /// Compute the `MD` tag bytes by consuming this iterator.
    ///
    /// Returns [`NmMdError::MissingReference`] if any required reference base
    /// (M-op position or D-op span) is outside the loaded `RefSeq` window —
    /// MD encodes ref bases verbatim so a missing position would corrupt the
    /// output. Insertions and soft clips don't appear in MD; reference skips
    /// (N) are excluded by convention to match `samtools calmd` behavior.
    ///
    /// # Example
    /// ```ignore
    /// let md = slim
    ///     .aligned_pairs_with_read(&store)?
    ///     .with_reference(&ref_seq)
    ///     .md()?;
    /// // md is e.g. b"5G2^AC3"
    /// ```
    pub fn md(self) -> Result<Vec<u8>, NmMdError> {
        let mut out: Vec<u8> = Vec::new();
        let mut run: u32 = 0; // current matching-base counter
        let mut last_was_deletion = false;

        for ev in self {
            match ev {
                AlignedPairWithRef::Match { kind, query, ref_base, rpos, .. } => {
                    let rb = ref_base.ok_or(NmMdError::MissingReference { rpos: *rpos })?;
                    let is_match = match kind {
                        MatchKind::SeqMatch => true,
                        MatchKind::SeqMismatch => false,
                        MatchKind::Match => bases_equal(query, rb),
                    };
                    if is_match {
                        run = run.saturating_add(1);
                        last_was_deletion = false;
                    } else {
                        // Emit run, then the ref base. MD spec requires a
                        // number BEFORE every letter — even 0.
                        push_number(&mut out, run);
                        out.push(base_to_md_byte(rb));
                        run = 0;
                        last_was_deletion = false;
                    }
                }
                AlignedPairWithRef::Deletion { del_len, ref_bases, rpos } => {
                    let bases = ref_bases.ok_or(NmMdError::MissingReference { rpos: *rpos })?;
                    debug_assert_eq!(bases.len(), del_len as usize);
                    // MD requires a number before `^` (even 0). However, two
                    // adjacent deletions need a `0` between them — `^A^C`
                    // would be ambiguous; canonical form is `^A0^C`.
                    if last_was_deletion {
                        out.push(b'0');
                    }
                    push_number(&mut out, run);
                    out.push(b'^');
                    for &b in bases {
                        out.push(base_to_md_byte(b));
                    }
                    run = 0;
                    last_was_deletion = true;
                }
                AlignedPairWithRef::Insertion { .. }
                | AlignedPairWithRef::RefSkip { .. }
                | AlignedPairWithRef::SoftClip { .. }
                | AlignedPairWithRef::Padding { .. }
                | AlignedPairWithRef::Unknown { .. } => {
                    // None of these contribute to MD.
                    last_was_deletion = false;
                }
            }
        }
        // MD must always end with a number (possibly 0).
        push_number(&mut out, run);
        Ok(out)
    }
}

/// Per-position NM contribution for a `Match` event.
fn match_contribution(kind: MatchKind, query: Base, ref_base: Option<Base>) -> u32 {
    match kind {
        MatchKind::SeqMatch => 0,
        MatchKind::SeqMismatch => 1,
        MatchKind::Match => {
            // Ambiguous: compare query against ref_base. If ref isn't loaded,
            // skip (count 0) — best-effort semantics.
            match ref_base {
                Some(rb) if !bases_equal(query, rb) => 1,
                _ => 0,
            }
        }
    }
}

/// Two `Base` values are "equal" for NM/MD purposes if they're identical AND
/// both are concrete (A/C/G/T). Two `N`s are NOT equal — htslib treats them as
/// a mismatch since neither is provably equal to the other.
fn bases_equal(a: Base, b: Base) -> bool {
    a == b && a != Base::Unknown
}

/// Convert a `Base` to its single-byte ASCII letter for MD encoding.
fn base_to_md_byte(b: Base) -> u8 {
    // Base is `#[repr(u8)]` with ASCII discriminants; this is the canonical
    // mapping used everywhere in seqair (see seqair-types/src/base.rs).
    b as u8
}

/// Append the decimal representation of `n` to `out`. We do not use `format!`
/// to avoid an allocation per run.
fn push_number(out: &mut Vec<u8>, n: u32) {
    if n == 0 {
        out.push(b'0');
        return;
    }
    // u32 max is 10 digits; stack buffer is plenty.
    let mut buf = [0u8; 10];
    let mut i = buf.len();
    let mut n = n;
    while n > 0 {
        i = i.saturating_sub(1);
        let digit = (n % 10) as u8;
        #[allow(clippy::indexing_slicing, reason = "i bounded by buf.len() above")]
        {
            buf[i] = b'0'.saturating_add(digit);
        }
        n /= 10;
    }
    #[allow(clippy::indexing_slicing, reason = "i is a valid start index into buf")]
    out.extend_from_slice(&buf[i..]);
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on small values")]
mod tests {
    use super::super::cigar::{CigarOp, CigarOpType};
    use super::super::owned_record::OwnedBamRecord;
    use super::super::pileup::RefSeq;
    use super::super::record_store::RecordStore;
    use super::*;
    use seqair_types::{BamFlags, Base, BaseQuality, Pos0};
    use std::rc::Rc;

    fn op(t: CigarOpType, len: u32) -> CigarOp {
        CigarOp::new(t, len)
    }

    fn p0(v: u32) -> Pos0 {
        Pos0::new(v).unwrap()
    }

    fn store_with_record(pos: Pos0, cigar: Vec<CigarOp>, seq: Vec<Base>) -> RecordStore<()> {
        let qual = vec![BaseQuality::from_byte(30); seq.len()];
        let owned = OwnedBamRecord::builder(0, Some(pos), b"r".to_vec())
            .flags(BamFlags::empty())
            .cigar(cigar)
            .seq(seq)
            .qual(qual)
            .build()
            .unwrap();
        let mut buf = Vec::new();
        owned.to_bam_bytes(&mut buf).unwrap();
        let mut store = RecordStore::<()>::new();
        let _ = store.push_raw(&buf, &mut ()).unwrap();
        store
    }

    fn ref_window(start: u32, ascii: &[u8]) -> RefSeq {
        let bases: Vec<Base> = ascii.iter().map(|&b| Base::from(b)).collect();
        RefSeq::new(Rc::from(bases), p0(start))
    }

    // ── NM ───────────────────────────────────────────────────────────────

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_match_only_no_mismatches_is_zero() {
        // 5M, query == ref everywhere
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"ACGTA");
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 0);
    }

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_counts_mismatches_in_match_op() {
        // 5M, 2 mismatches: query ACGTA vs ref ACGAA → only position 3 differs
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"ACGAA"); // pos 3 is A vs query T
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 1);
    }

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_counts_insertions_and_deletions() {
        // 2M + 3I + 2M + 1D + 1M
        // NM = 0 mismatches + 3 (I) + 1 (D) = 4
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Insertion, 3),
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 1),
            op(CigarOpType::Match, 1),
        ];
        let seq: Vec<Base> = b"AAGGGCCAA".iter().map(|&b| Base::from(b)).collect(); // 2+3+2+0+1=8 query
        // Wait: 2M(AA) + 3I(GGG) + 2M(CC) + 1D + 1M(A) = 2+3+2+0+1 = 8 query
        // But seq has 9 bytes. Adjust:
        let _ = seq.len(); // shut up
        let seq: Vec<Base> = b"AAGGGCCAA"[..8].iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        // Ref covers pos 100..106 (M+M+M+D+M = 6 ref positions)
        let ref_seq = ref_window(100, b"AACCXA"); // 6 bases
        // Query M positions: 0,1 → AA vs AA (match), 5,6 → CC vs CC (match),
        //                    7 → A vs A (match). 0 mismatches.
        // NM = 0 + 3 (insert) + 1 (delete) = 4
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 4);
    }

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_seqmatch_and_seqmismatch_dont_consult_ref() {
        // 2= + 3X + 2= → NM = 0 + 3 + 0 = 3 (without consulting ref).
        let cigar = vec![
            op(CigarOpType::SeqMatch, 2),
            op(CigarOpType::SeqMismatch, 3),
            op(CigarOpType::SeqMatch, 2),
        ];
        let seq: Vec<Base> = b"AAACCCC".iter().map(|&b| Base::from(b)).collect(); // 7 query
        let store = store_with_record(p0(100), cigar, seq);
        // Use a ref that AGREES with the read everywhere — = and X must
        // override, so NM stays 3 regardless.
        let ref_seq = ref_window(100, b"AAACCCC");
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 3, "X always contributes regardless of ref agreement");
    }

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_skips_match_positions_outside_loaded_ref() {
        // Read at pos 100 with 5M, but ref window starts at 200 — entire
        // record's match positions are outside the loaded window. NM is
        // best-effort, so it returns 0 (skipped, not errored).
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(200, b"AAAAA");
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 0);
    }

    // r[verify cigar.aligned_pairs.nm.compute]
    #[test]
    fn nm_treats_n_as_mismatch() {
        // Query N vs ref A — htslib counts this as a mismatch (1).
        let cigar = vec![op(CigarOpType::Match, 1)];
        let seq = vec![Base::Unknown];
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"A");
        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        assert_eq!(nm, 1);
    }

    // ── MD ───────────────────────────────────────────────────────────────

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_perfect_match_is_just_length_number() {
        // 5M with all matches — MD = "5"
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"ACGTA");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(md, b"5");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_single_mismatch_emits_ref_base() {
        // 5M, query ACGTA vs ref ACTTA → mismatch at pos 2: G vs T
        // MD = "2T2"
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"ACTTA");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(md, b"2T2");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_consecutive_mismatches() {
        // 5M, query ACGTA vs ref TCGGT → mismatches at 0, 3, 4
        // MD = "0T2G0T0"
        // Reading: 0 matches, ref T, 2 matches, ref G, 0 matches, ref T, 0 matches end
        let cigar = vec![op(CigarOpType::Match, 5)];
        let seq: Vec<Base> = b"ACGTA".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"TCGGT");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        // pos 0: A vs T mismatch → "0T"
        // pos 1: C vs C match → run=1
        // pos 2: G vs G match → run=2
        // pos 3: T vs G mismatch → emit "2G", run=0
        // pos 4: A vs T mismatch → emit "0T", run=0
        // end: emit "0"
        // Total: "0T2G0T0"
        assert_eq!(std::str::from_utf8(&md).unwrap(), "0T2G0T0");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_deletion_emits_caret_and_ref_bases() {
        // 2M + 3D + 2M, all matches in M ops, deletion of "ACG"
        // MD = "2^ACG2"
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 3),
            op(CigarOpType::Match, 2),
        ];
        let seq: Vec<Base> = b"AATT".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"AAACGTT");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(std::str::from_utf8(&md).unwrap(), "2^ACG2");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_skips_insertions_and_soft_clips() {
        // 2S + 2M + 3I + 2M, no mismatches
        // MD ignores S and I → "4"
        let cigar = vec![
            op(CigarOpType::SoftClip, 2),
            op(CigarOpType::Match, 2),
            op(CigarOpType::Insertion, 3),
            op(CigarOpType::Match, 2),
        ];
        let seq: Vec<Base> = b"NNAAGGGCC".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"AACC");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(md, b"4");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_starts_with_zero_when_first_position_mismatches() {
        // 3M, query TGC vs ref AGC → mismatch at pos 0: T vs A
        // MD must start with a number, even 0: "0A2"
        let cigar = vec![op(CigarOpType::Match, 3)];
        let seq: Vec<Base> = b"TGC".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"AGC");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(md, b"0A2");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_seqmatch_and_seqmismatch_use_ref_base_for_mismatch_letter() {
        // 2= + 1X + 2= — for X, MD must emit the REF base, not the query.
        // ref ACGTT, query AC*TT where * is whatever (ignored — X says it's a
        // mismatch).
        let cigar = vec![
            op(CigarOpType::SeqMatch, 2),
            op(CigarOpType::SeqMismatch, 1),
            op(CigarOpType::SeqMatch, 2),
        ];
        // We choose query = "ACATT" (X position has 'A'); ref has 'G' there.
        let seq: Vec<Base> = b"ACATT".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(100, b"ACGTT");
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(md, b"2G2");
    }

    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn md_errors_on_missing_reference() {
        // Read at pos 100 with 3M, ref window starts at 200 — first M position
        // has no ref. md() returns Err.
        let cigar = vec![op(CigarOpType::Match, 3)];
        let seq: Vec<Base> = b"ACG".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        let ref_seq = ref_window(200, b"AAA");
        let result =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).md();
        assert!(matches!(result, Err(NmMdError::MissingReference { rpos: 100 })));
    }

    // ── End-to-end: NM + MD on a realistic record ─────────────────────────

    // r[verify cigar.aligned_pairs.nm.compute]
    // r[verify cigar.aligned_pairs.md.compute]
    #[test]
    fn nm_and_md_realistic_record() {
        // 2S + 3M + 1I + 2M + 1D + 2M
        // ref window covers M+M+D+M positions:
        //   M(3): ACG vs ref ACT (1 mismatch at pos 2: G vs T)
        //   M(2): GG vs ref GG (matches)
        //   D(1): ref T (deleted)
        //   M(2): TT vs ref TT (matches)
        // Insertion of 1 query base, soft clip of 2 query bases.
        // NM = 1 (M mismatch) + 1 (insert) + 1 (delete) = 3
        // MD: 2T2 ... 2 ... ^T2 → walk:
        //   pos 0 (M): A vs A match → run=1
        //   pos 1 (M): C vs C match → run=2
        //   pos 2 (M): G vs T mismatch → "2T", run=0
        //   I: skip
        //   pos 3 (M): G vs G match → run=1
        //   pos 4 (M): G vs G match → run=2
        //   D: emit "2^T", run=0
        //   pos 5 (M): T vs T match → run=1
        //   pos 6 (M): T vs T match → run=2
        //   end: emit "2"
        // MD = "2T2^T2"
        let cigar = vec![
            op(CigarOpType::SoftClip, 2),
            op(CigarOpType::Match, 3),
            op(CigarOpType::Insertion, 1),
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 1),
            op(CigarOpType::Match, 2),
        ];
        // Query (10): NN ACG X GG TT → 2S + 3M(ACG) + 1I(X) + 2M(GG) + 0D + 2M(TT)
        let seq: Vec<Base> = b"NNACGXGGTT".iter().map(|&b| Base::from(b)).collect();
        let store = store_with_record(p0(100), cigar, seq);
        // Ref bases at positions 100..108 (8 ref bases for 3M+2M+1D+2M = 8)
        let ref_seq = ref_window(100, b"ACTGGTTT");

        let nm =
            store.record(0).aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).nm();
        let md = store
            .record(0)
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .md()
            .unwrap();
        assert_eq!(nm, 3);
        assert_eq!(md, b"2T2^T2");
    }

    // ── Cross-validation: NM and MD must agree ────────────────────────────

    /// Parse the count of mismatch letters and deleted bases from an MD tag.
    /// Returns `(mismatch_count, deleted_count)`.
    ///
    /// MD grammar: `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`. Digits = match-run lengths
    /// (skipped); letters = mismatches; `^[A-Z]+` = deletion (each letter
    /// counts toward `deleted_count`).
    fn parse_md_for_counts(md: &[u8]) -> (u32, u32) {
        let mut mismatches = 0u32;
        let mut deletions = 0u32;
        let mut i = 0;
        while i < md.len() {
            #[allow(clippy::indexing_slicing, reason = "i bounded by md.len()")]
            let b = md[i];
            match b {
                b'0'..=b'9' => {
                    i = i.saturating_add(1);
                }
                b'^' => {
                    // Consume all letters following ^
                    i = i.saturating_add(1);
                    while i < md.len() {
                        #[allow(clippy::indexing_slicing, reason = "i bounded by md.len()")]
                        let c = md[i];
                        if c.is_ascii_alphabetic() {
                            deletions = deletions.saturating_add(1);
                            i = i.saturating_add(1);
                        } else {
                            break;
                        }
                    }
                }
                b'A'..=b'Z' | b'a'..=b'z' => {
                    mismatches = mismatches.saturating_add(1);
                    i = i.saturating_add(1);
                }
                _ => {
                    i = i.saturating_add(1);
                }
            }
        }
        (mismatches, deletions)
    }

    #[test]
    fn parse_md_for_counts_matches_known_strings() {
        assert_eq!(parse_md_for_counts(b"5"), (0, 0));
        assert_eq!(parse_md_for_counts(b"2T2"), (1, 0));
        assert_eq!(parse_md_for_counts(b"2^ACG2"), (0, 3));
        assert_eq!(parse_md_for_counts(b"0T2G0T0"), (3, 0));
        assert_eq!(parse_md_for_counts(b"2T2^T2"), (1, 1));
    }

    mod proptests {
        use super::super::super::cigar::CigarOp;
        use super::*;
        use proptest::prelude::*;

        /// Generate a CIGAR using only ops that work cleanly with NM/MD.
        /// Avoids `H` (zero contribution but complicates query length),
        /// `P`/Unknown (filtered out by default), `B` (back).
        fn arb_cigar_op() -> impl Strategy<Value = CigarOp> {
            // Op types that consume query and/or ref. Lengths 1..=8 to keep
            // cigars tractable.
            let ops = prop_oneof![
                Just(CigarOpType::Match),
                Just(CigarOpType::Insertion),
                Just(CigarOpType::Deletion),
                Just(CigarOpType::SoftClip),
                Just(CigarOpType::SeqMatch),
                Just(CigarOpType::SeqMismatch),
            ];
            (ops, 1u32..=8u32).prop_map(|(t, len)| CigarOp::new(t, len))
        }

        proptest! {
            /// **NM/MD consistency.** For any record where the reference covers
            /// the full alignment, computed NM must equal the sum of:
            /// - mismatches encoded in MD (letter count, excluding `^` deletions)
            /// - deletions encoded in MD (bases after `^`)
            /// - insertions in the CIGAR
            ///
            /// This is a cross-validation: if NM and MD ever disagree, one of
            /// them is wrong. Both are implemented in this file by walking the
            /// same iterator, so this catches divergent bugs that fixture tests
            /// might miss.
            #[test]
            fn nm_equals_md_mismatches_plus_md_deletions_plus_cigar_insertions(
                ops in proptest::collection::vec(arb_cigar_op(), 1..=8),
                read_pos in 0u32..=1_000u32,
            ) {
                let qlen: u32 = ops
                    .iter()
                    .filter(|o| o.op_type().consumes_query())
                    .map(|o| o.len())
                    .sum();
                let rlen: u32 = ops
                    .iter()
                    .filter(|o| o.op_type().consumes_ref())
                    .map(|o| o.len())
                    .sum();
                if qlen == 0 {
                    return Ok(()); // need at least 1 query base for OwnedBamRecord
                }

                // Synthetic seq: rotate through ACGT so we get realistic mismatches.
                let bases = [Base::A, Base::C, Base::G, Base::T];
                let seq: Vec<Base> =
                    (0..qlen).map(|i| bases[(i % 4) as usize]).collect();
                let qual: Vec<BaseQuality> =
                    (0..qlen).map(|_| BaseQuality::from_byte(30)).collect();

                let owned = OwnedBamRecord::builder(0, Some(p0(read_pos)), b"r".to_vec())
                    .flags(BamFlags::empty())
                    .cigar(ops.clone())
                    .seq(seq)
                    .qual(qual)
                    .build()
                    .unwrap();
                let mut buf = Vec::new();
                owned.to_bam_bytes(&mut buf).unwrap();
                let mut store = RecordStore::<()>::new();
                let _ = store.push_raw(&buf, &mut ()).unwrap();

                // Reference window covering the full ref span. Different
                // pattern (CGTA rotation) so M-op mismatches actually happen.
                let ref_pat = [Base::C, Base::G, Base::T, Base::A];
                let ref_buf: Vec<Base> =
                    (0..rlen).map(|i| ref_pat[(i % 4) as usize]).collect();
                let ref_seq = RefSeq::new(Rc::from(ref_buf), p0(read_pos));

                let nm = store
                    .record(0)
                    .aligned_pairs_with_read(&store)
                    .unwrap()
                    .with_reference(&ref_seq)
                    .nm();
                let md = store
                    .record(0)
                    .aligned_pairs_with_read(&store)
                    .unwrap()
                    .with_reference(&ref_seq)
                    .md()
                    .unwrap();

                let (md_mismatches, md_deletions) = parse_md_for_counts(&md);
                let cigar_insertions: u32 = ops
                    .iter()
                    .filter(|o| matches!(o.op_type(), CigarOpType::Insertion))
                    .map(|o| o.len())
                    .sum();

                let derived_nm = md_mismatches
                    .saturating_add(md_deletions)
                    .saturating_add(cigar_insertions);
                prop_assert_eq!(
                    nm, derived_nm,
                    "NM ({}) != mismatches({}) + deletions({}) + insertions({}); MD={:?}",
                    nm, md_mismatches, md_deletions, cigar_insertions, std::str::from_utf8(&md)
                );

                // Also: MD must always start with a digit and end with a digit.
                prop_assert!(!md.is_empty(), "MD must never be empty");
                prop_assert!(
                    md[0].is_ascii_digit(),
                    "MD must start with a digit: {:?}",
                    std::str::from_utf8(&md),
                );
                prop_assert!(
                    md[md.len() - 1].is_ascii_digit(),
                    "MD must end with a digit: {:?}",
                    std::str::from_utf8(&md),
                );

            }
        }
    }

    // ── push_number unit tests (internal helper) ─────────────────────────

    #[test]
    fn push_number_handles_zero_and_multidigit() {
        let mut buf = Vec::new();
        push_number(&mut buf, 0);
        assert_eq!(buf, b"0");

        let mut buf = Vec::new();
        push_number(&mut buf, 1);
        assert_eq!(buf, b"1");

        let mut buf = Vec::new();
        push_number(&mut buf, 9);
        assert_eq!(buf, b"9");

        let mut buf = Vec::new();
        push_number(&mut buf, 10);
        assert_eq!(buf, b"10");

        let mut buf = Vec::new();
        push_number(&mut buf, 12345);
        assert_eq!(buf, b"12345");

        let mut buf = Vec::new();
        push_number(&mut buf, u32::MAX);
        assert_eq!(buf, b"4294967295");
    }
}
