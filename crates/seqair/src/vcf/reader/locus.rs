//! [`Locus`]: where a record is (`r[vcf_fields.locus+2]`).

use super::error::{RecordError, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::raw::RawRecord;
use crate::vcf::header::{HeaderId, VcfHeader};
use crate::vcf::record_encoder::ContigId;
use seqair_types::Pos1;

/// A VCF position: POS 0 (BCF −1), the telomere before a contig's first
/// base, or a base. POS N+1, the telomere after the last base, is an
/// ordinary `At`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VcfPos {
    Telomere,
    At(Pos1),
}

impl VcfPos {
    /// The 1-based position `pos` (0 for the telomere).
    pub(crate) fn from_one_based(pos: i64) -> Option<Self> {
        match pos {
            0 => Some(Self::Telomere),
            p => Pos1::new(u32::try_from(p).ok()?).map(Self::At),
        }
    }
}

// r[impl vcf_fields.locus+2]
/// A record's contig, start and inclusive end. Cheap to clone; reading it
/// skips everything past a record's fixed fields.
///
/// The end is the start plus the reference length minus one, as the index
/// was built: for BCF the stored `rlen`, which accounts for INFO END,
/// symbolic alleles' SVLEN and gVCF blocks.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Locus {
    contig: ContigId,
    start: VcfPos,
    end: VcfPos,
}

impl Locus {
    /// A record's locus from its 0-based BCF position and reference length.
    pub(crate) fn from_bcf(contig: ContigId, pos: i32, rlen: i32) -> Result<Self, RecordError> {
        Self::from_zero_based(contig, i64::from(pos), i64::from(rlen))
    }

    // r[impl vcf_values.rlen]
    /// A record's locus from its 0-based position (−1 for the telomere) and
    /// reference length: the end is `start + rlen − 1`.
    pub(crate) fn from_zero_based(
        contig: ContigId,
        pos: i64,
        rlen: i64,
    ) -> Result<Self, RecordError> {
        let start = pos.saturating_add(1);
        let end = start.saturating_add(rlen).saturating_sub(1);
        let start = VcfPos::from_one_based(start).ok_or(RecordError::InvalidPosition { pos })?;
        let end = VcfPos::from_one_based(end).ok_or(RecordError::InvalidSpan { pos, rlen })?;
        Ok(Self { contig, start, end })
    }

    pub fn contig(&self) -> &ContigId {
        &self.contig
    }

    pub fn start(&self) -> VcfPos {
        self.start
    }

    /// The last reference base the record spans.
    pub fn end(&self) -> VcfPos {
        self.end
    }
}

impl Fields for Locus {
    type Plan = HeaderId;

    fn plan(header: &VcfHeader, _: &ReadOptions) -> Result<HeaderId, PlanError> {
        Ok(header.id())
    }

    fn needs(_: &HeaderId) -> Sections {
        Sections::FIXED
    }

    fn new(_: &HeaderId, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        raw.locus()
    }

    fn fill(&mut self, _: &HeaderId, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        // `locus` clones the contig: a `SmolStr` of a name is inline.
        *self = raw.locus()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_fields.locus+2]
    #[test]
    fn bcf_positions_and_spans() {
        let header = VcfHeader::builder().build().unwrap();
        let chr = ContigId { tid: 0, name: "chr1".into(), header: header.id() };
        let at = |p| VcfPos::At(Pos1::new(p).unwrap());
        let locus = |pos, rlen| Locus::from_bcf(chr.clone(), pos, rlen);
        let span = |l: Locus| (l.start(), l.end());
        assert_eq!(span(locus(99, 1).unwrap()), (at(100), at(100)));
        assert_eq!(span(locus(99, 4).unwrap()), (at(100), at(103)));
        assert_eq!(span(locus(-1, 1).unwrap()), (VcfPos::Telomere, VcfPos::Telomere));
        assert_eq!(span(locus(-1, 2).unwrap()), (VcfPos::Telomere, at(1)));
        assert!(matches!(locus(-2, 1), Err(RecordError::InvalidPosition { pos: -2 })));
        assert!(matches!(locus(0, -5), Err(RecordError::InvalidSpan { .. })));
        assert!(matches!(locus(i32::MAX, 2), Err(RecordError::InvalidPosition { .. })));
        assert!(VcfPos::Telomere < at(1));
    }
}
