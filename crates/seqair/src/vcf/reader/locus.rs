//! [`Locus`]: where a record is (`r[vcf_fields.locus+3]`).

use super::error::{RecordError, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::raw::RawRecord;
use crate::vcf::header::{HeaderId, VcfHeader};
use crate::vcf::record_encoder::ContigId;
use seqair_types::{Pos0, Pos1};
use std::ops::RangeInclusive;

/// A VCF position: POS 0 (BCF −1), the telomere before a contig's first
/// base, or a base. POS N+1, the telomere after the last base, is an
/// ordinary `At`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VcfPos {
    Telomere,
    At(Pos1),
}

impl std::fmt::Display for VcfPos {
    /// The VCF POS: 0 for the telomere.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Telomere => f.write_str("0"),
            Self::At(pos) => write!(f, "{}", pos.as_u32()),
        }
    }
}

impl VcfPos {
    // r[impl vcf_fields.locus+3]
    /// The base, 1-based; `None` for the telomere.
    pub fn pos1(self) -> Option<Pos1> {
        match self {
            Self::Telomere => None,
            Self::At(pos) => Some(pos),
        }
    }

    /// The base, 0-based; `None` for the telomere.
    pub fn pos0(self) -> Option<Pos0> {
        self.pos1().map(Pos1::to_zero_based)
    }

    /// The 1-based position `pos` (0 for the telomere).
    pub(crate) fn from_one_based(pos: i64) -> Option<Self> {
        match pos {
            0 => Some(Self::Telomere),
            p => Pos1::new(u32::try_from(p).ok()?).map(Self::At),
        }
    }
}

// r[impl vcf_values.rlen]
/// A record's start and inclusive end from its 0-based position (−1 for the
/// telomere) and reference length: the end is `start + rlen − 1`.
pub(crate) fn span_of(pos: i64, rlen: i64) -> Result<(VcfPos, VcfPos), RecordError> {
    let start = pos.saturating_add(1);
    let end = start.saturating_add(rlen).saturating_sub(1);
    let start = VcfPos::from_one_based(start).ok_or(RecordError::InvalidPosition { pos })?;
    let end = VcfPos::from_one_based(end).ok_or(RecordError::InvalidSpan { pos, rlen })?;
    Ok((start, end))
}

// r[impl vcf_fields.locus+3]
/// A record's contig, start and inclusive end. Cheap to clone; reading it
/// skips everything past a record's fixed fields.
///
/// The end is the start plus the reference length minus one, as the index
/// was built: for BCF the stored `rlen`, for VCF text the length htslib
/// computes; both account for INFO END, symbolic alleles' SVLEN and gVCF
/// blocks' FORMAT LEN.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Locus {
    contig: ContigId,
    start: VcfPos,
    end: VcfPos,
}

// r[impl vcf_reader.debug]
/// `chr1:100`, or `chr1:100-103` for a record spanning more than POS.
impl std::fmt::Display for Locus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.contig.name(), self.start)?;
        if self.end != self.start {
            write!(f, "-{}", self.end)?;
        }
        Ok(())
    }
}

impl Locus {
    /// A locus whose span [`span_of`] checked.
    pub(crate) fn new_unchecked(contig: ContigId, start: VcfPos, end: VcfPos) -> Self {
        Self { contig, start, end }
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

    // r[impl vcf_fields.locus+3]
    /// The reference bases the record spans, from base 1 when it starts at
    /// the telomere before it; `None` when it spans that telomere only.
    pub fn bases(&self) -> Option<RangeInclusive<Pos1>> {
        let end = self.end.pos1()?;
        Some(self.start.pos1().unwrap_or(Pos1::MIN)..=end)
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
        let (start, end) = raw.span()?;
        // Sorted files change contig rarely: compare before cloning.
        if self.contig != *raw.contig() {
            self.contig = raw.contig().clone();
        }
        self.start = start;
        self.end = end;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_fields.locus+3]
    #[test]
    fn positions_and_spans() {
        let at = |p| VcfPos::At(Pos1::new(p).unwrap());
        let locus = span_of;
        assert_eq!(locus(99, 1).unwrap(), (at(100), at(100)));
        assert_eq!(locus(99, 4).unwrap(), (at(100), at(103)));
        assert_eq!(locus(-1, 1).unwrap(), (VcfPos::Telomere, VcfPos::Telomere));
        assert_eq!(locus(-1, 2).unwrap(), (VcfPos::Telomere, at(1)));
        assert!(matches!(locus(-2, 1), Err(RecordError::InvalidPosition { pos: -2 })));
        assert!(matches!(locus(0, -5), Err(RecordError::InvalidSpan { .. })));
        assert!(matches!(locus(i64::from(i32::MAX), 2), Err(RecordError::InvalidPosition { .. })));
        assert!(VcfPos::Telomere < at(1));
        assert_eq!((at(7).pos1(), at(7).pos0()), (Pos1::new(7), Some(Pos0::new(6).unwrap())));
        assert_eq!((VcfPos::Telomere.pos1(), VcfPos::Telomere.pos0()), (None, None));
        let contig = ContigId {
            tid: 0,
            name: "c".into(),
            header: VcfHeader::builder().build().unwrap().id(),
        };
        let bases = |s, e| Locus::new_unchecked(contig.clone(), s, e).bases();
        assert_eq!(bases(at(5), at(7)), Some(Pos1::new(5).unwrap()..=Pos1::new(7).unwrap()));
        assert_eq!(bases(VcfPos::Telomere, at(1)), Some(Pos1::MIN..=Pos1::MIN));
        assert_eq!(bases(VcfPos::Telomere, VcfPos::Telomere), None);
    }
}
