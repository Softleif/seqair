//! [`Site`]: ID, REF, ALT, QUAL and FILTER (`r[vcf_fields.site+2]`).

use super::bcf_value::{BcfValueError, TypedValue, skip_typed_value};
use super::error::{KeyKind, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::file::FileShared;
use super::index::DictId;
use super::raw::{Col, RawRecord, TextCols};
use super::text_value::parse_f32;
use crate::vcf::alleles::Alleles;
use crate::vcf::bcf_encoding::FLOAT_MISSING;
use crate::vcf::error::AllelesError;
use crate::vcf::header::VcfHeader;
use crate::vcf::record_encoder::FilterId;
use seqair_types::{Base, SmallVec, SmolStr};
use std::ops::Range;
use std::sync::Arc;

/// One allele of a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Allele<'a> {
    /// Bases (`C`, `CTT`), as stored; REF is always `Bases`.
    Bases(&'a [u8]),
    /// A symbolic allele, brackets included: `<DEL>`, `<*>`, `<NON_REF>`.
    Symbolic(&'a str),
    /// A breakend: `G]17:198982]`, `.A`.
    Breakend(&'a str),
    /// `*`: the allele is deleted by an overlapping variant.
    Overlapping,
    /// `.` inside a list of ALTs (`T,.`): malformed, but kept so the
    /// alleles after it keep their GT index. A record without ALT has no
    /// alleles past REF instead.
    Missing,
}

impl<'a> Allele<'a> {
    /// The allele as written: `b"*"` for [`Overlapping`](Self::Overlapping),
    /// `b"."` for [`Missing`](Self::Missing).
    pub fn as_bytes(&self) -> &'a [u8] {
        match *self {
            Self::Bases(b) => b,
            Self::Symbolic(s) | Self::Breakend(s) => s.as_bytes(),
            Self::Overlapping => b"*",
            Self::Missing => b".",
        }
    }

    fn classify(allele: &'a str) -> Self {
        let bytes = allele.as_bytes();
        if bytes == b"*" {
            Self::Overlapping
        } else if bytes == b"." {
            Self::Missing
        } else if bytes.first() == Some(&b'<') && bytes.last() == Some(&b'>') {
            Self::Symbolic(allele)
        } else if bytes.iter().any(|&b| b == b'[' || b == b']')
            || (bytes.len() > 1 && (bytes.first() == Some(&b'.') || bytes.last() == Some(&b'.')))
        {
            Self::Breakend(allele)
        } else {
            Self::Bases(bytes)
        }
    }
}

/// A record's site columns, decoded once by the reader, whatever the
/// format: where the ID and each allele lie in the record's bytes (checked
/// to be UTF-8), the FILTER ids (checked to name FILTERs), and QUAL.
#[derive(Debug, Clone, Default)]
pub(crate) struct SiteParts {
    pub(crate) id: Range<u32>,
    /// REF, then each ALT.
    pub(crate) alleles: Vec<Range<u32>>,
    pub(crate) filter: FilterColumn,
    /// The FILTER ids in record order (for `Failed`, and `PASS` for `Pass`).
    pub(crate) filter_ids: Vec<DictId>,
    pub(crate) qual: Option<f32>,
}

/// What a FILTER column says, before its ids are listed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) enum FilterColumn {
    #[default]
    Missing,
    Pass,
    Failed,
}

impl SiteParts {
    /// Clear for the next record, keeping the buffers.
    pub(crate) fn begin(&mut self) {
        self.alleles.clear();
        self.filter_ids.clear();
    }

    /// The bytes the ID and alleles span: one run in either format.
    fn span(&self) -> Range<usize> {
        let end = self.alleles.last().map_or(self.id.end, |a| a.end);
        self.id.start as usize..end as usize
    }

    /// `other`, its spans moved `shift` bytes towards the start.
    fn rebased_from(&mut self, other: &Self, shift: u32) {
        self.id = shifted(&other.id, shift);
        self.alleles.clear();
        self.alleles.extend(other.alleles.iter().map(|r| shifted(r, shift)));
        self.filter = other.filter;
        self.filter_ids.clone_from(&other.filter_ids);
        self.qual = other.qual;
    }

    // r[impl vcf_fields.site+2]
    /// Decode the site columns of a BCF shared section (after its 24 fixed
    /// bytes): IDs and alleles UTF-8, every FILTER index a FILTER of the
    /// header. Returns where INFO starts.
    pub(crate) fn decode_bcf(
        &mut self,
        shared: &[u8],
        n_allele: u16,
        qual: u32,
        file: &FileShared,
    ) -> Result<usize, RecordError> {
        self.begin();
        let mut at = 24usize;
        let string = |at: &mut usize, part: RecordPart| -> Result<Range<u32>, RecordError> {
            let rest = shared.get(*at..).unwrap_or_default();
            let (value, len) = TypedValue::parse(rest).map_err(|e| value_error(part, e))?;
            let s = value.str().map_err(|e| value_error(part, e))?;
            // The payload, as `str` read it: after the descriptor, cut at NUL.
            let start = at.saturating_add(len).saturating_sub(value.stored_len());
            *at = at.saturating_add(len);
            span(start, s.len())
        };
        self.id = string(&mut at, RecordPart::Id)?;
        for _ in 0..n_allele {
            let allele = string(&mut at, RecordPart::Alleles)?;
            self.alleles.push(allele);
        }
        // r[impl vcf_fields.site+2]
        // htslib writes a text ALT `.` as one `.` allele after REF.
        if let [_, alt] = self.alleles.as_slice()
            && shared.get(alt.start as usize..alt.end as usize) == Some(b".")
        {
            self.alleles.truncate(1);
        }
        let rest = shared.get(at..).unwrap_or_default();
        let (value, len) =
            TypedValue::parse(rest).map_err(|e| value_error(RecordPart::Filter, e))?;
        let ids = value.ints().map_err(|e| value_error(RecordPart::Filter, e))?;
        for index in ids.flatten() {
            let declared = usize::try_from(index).ok().and_then(|i| file.filters.get(i)?.as_ref());
            let id = declared
                .and_then(|f| DictId::new(f.dict_idx))
                .ok_or(RecordError::UnknownKey { index, kind: KeyKind::Filter })?;
            self.filter_ids.push(id);
        }
        self.filter = match self.filter_ids.as_slice() {
            [] => FilterColumn::Missing,
            [id] if id.get() == 0 => FilterColumn::Pass,
            _ => FilterColumn::Failed,
        };
        self.qual = (qual != FLOAT_MISSING).then(|| f32::from_bits(qual));
        Ok(at.saturating_add(len))
    }
}

impl SiteParts {
    // r[impl vcf_fields.site+2]
    // r[impl vcf_reader.undeclared.filter]
    /// Decode the site columns of a VCF text line: IDs and alleles UTF-8,
    /// QUAL a float or `.`, FILTER names looked up in the header. A FILTER
    /// name the header does not declare still fails the record (htslib
    /// adds it to its header); like an undeclared INFO key it is logged
    /// once per file, and no `FilterId` can name it.
    pub(crate) fn decode_text(
        &mut self,
        line: &[u8],
        cols: &TextCols,
        file: &FileShared,
    ) -> Result<(), RecordError> {
        self.begin();
        let utf8 = |col: Col, part: RecordPart| -> Result<Range<u32>, RecordError> {
            std::str::from_utf8(cols.col(line, col))
                .map_err(|source| RecordError::Utf8 { part, source })?;
            Ok(cols.range(col))
        };
        self.id = utf8(Col::Id, RecordPart::Id)?;
        self.alleles.push(utf8(Col::Ref, RecordPart::Alleles)?);
        let alt = cols.range(Col::Alt);
        let alt_bytes = cols.col(line, Col::Alt);
        if alt_bytes != b"." {
            std::str::from_utf8(alt_bytes)
                .map_err(|source| RecordError::Utf8 { part: RecordPart::Alleles, source })?;
            let mut start = alt.start;
            for allele in alt_bytes.split(|&b| b == b',') {
                // r[impl vcf_reader.limits]
                if self.alleles.len() >= usize::from(u16::MAX) {
                    let n = alt_bytes.iter().filter(|&&b| b == b',').count().saturating_add(2);
                    return Err(RecordError::TooManyAlleles { n });
                }
                // An allele is at most the line, which fits u32.
                let end = start.saturating_add(u32::try_from(allele.len()).unwrap_or(u32::MAX));
                self.alleles.push(start..end);
                start = end.saturating_add(1);
            }
        }
        self.qual = parse_f32(cols.col(line, Col::Qual))
            .map_err(|source| RecordError::Text { part: RecordPart::Qual, source })?;
        let filter = cols.col(line, Col::Filter);
        self.filter = match filter {
            b"." | b"" => FilterColumn::Missing,
            b"PASS" => FilterColumn::Pass,
            _ => FilterColumn::Failed,
        };
        if self.filter != FilterColumn::Missing {
            for name in filter.split(|&b| b == b';') {
                let declared = file
                    .keys
                    .get(name)
                    .filter(|id| file.filters.get(id.get() as usize).is_some_and(Option::is_some));
                match declared {
                    Some(id) => self.filter_ids.push(id),
                    None => file.warn_undeclared(KeyKind::Filter, name),
                }
            }
        }
        Ok(())
    }
}

fn shifted(r: &Range<u32>, shift: u32) -> Range<u32> {
    r.start.saturating_sub(shift)..r.end.saturating_sub(shift)
}

fn value_error(part: RecordPart, e: BcfValueError) -> RecordError {
    match e {
        BcfValueError::Utf8(source) => RecordError::Utf8 { part, source },
        source => RecordError::Value { part, source },
    }
}

/// `start..start + len` as `u32`s: a record is at most 256 MiB.
pub(crate) fn span(start: usize, len: usize) -> Result<Range<u32>, RecordError> {
    let end = start.saturating_add(len);
    match (u32::try_from(start), u32::try_from(end)) {
        (Ok(s), Ok(e)) => Ok(s..e),
        _ => Err(RecordError::TooLarge { size: end as u64, max: u64::from(u32::MAX) }),
    }
}

/// Where INFO starts in a BCF shared section, its site columns skipped
/// unchecked (for a type that reads INFO but not the site).
pub(crate) fn bcf_info_start(shared: &[u8], n_allele: u16) -> Result<usize, RecordError> {
    let skip = |at: usize, part: RecordPart| -> Result<usize, RecordError> {
        let rest = shared.get(at..).unwrap_or_default();
        let len = skip_typed_value(rest).map_err(|source| RecordError::Value { part, source })?;
        Ok(at.saturating_add(len))
    };
    let mut at = skip(24, RecordPart::Id)?;
    for _ in 0..n_allele {
        at = skip(at, RecordPart::Alleles)?;
    }
    skip(at, RecordPart::Filter)
}

/// The text of `range` in `bytes`; decoding checked it is UTF-8.
fn text<'a>(bytes: &'a [u8], range: &Range<u32>) -> &'a str {
    let b = bytes.get(range.start as usize..range.end as usize).unwrap_or_default();
    std::str::from_utf8(b).unwrap_or_default()
}

/// A record's alleles: REF, then the ALTs.
#[derive(Debug, Clone, Copy)]
pub struct AllelesRef<'a> {
    bytes: &'a [u8],
    spans: &'a [Range<u32>],
}

impl<'a> AllelesRef<'a> {
    /// The REF bases as stored.
    pub fn reference(&self) -> &'a [u8] {
        self.spans.first().map_or(&[], |r| text(self.bytes, r).as_bytes())
    }

    // r[impl vcf_fields.site+2]
    /// Allele `i` as a GT names it: 0 is REF, 1 the first ALT; `None` past
    /// the last.
    pub fn get(&self, i: u16) -> Option<Allele<'a>> {
        let span = self.spans.get(usize::from(i))?;
        let allele = text(self.bytes, span);
        Some(if i == 0 { Allele::Bases(allele.as_bytes()) } else { Allele::classify(allele) })
    }

    /// REF, then the ALTs, in GT index order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Allele<'a>> + Clone + use<'a> {
        let this = *self;
        // At most 65 535 alleles (`r[vcf_reader.limits]`).
        let n = u16::try_from(self.spans.len()).unwrap_or(u16::MAX);
        // `get` finds every index below the count.
        (0..n).map(move |i| this.get(i).unwrap_or(Allele::Missing))
    }

    /// How many alleles there are, REF included.
    pub fn len(&self) -> usize {
        self.spans.len()
    }

    /// Never true for a decoded record, which has a REF.
    pub fn is_empty(&self) -> bool {
        self.spans.is_empty()
    }

    fn alt_strs(&self) -> impl ExactSizeIterator<Item = &'a str> + Clone + use<'a> {
        let bytes = self.bytes;
        self.spans.get(1..).unwrap_or_default().iter().map(move |r| text(bytes, r))
    }

    /// The ALT alleles; none for a record without ALT (`.`).
    pub fn alts(&self) -> impl ExactSizeIterator<Item = Allele<'a>> + Clone + use<'a> {
        self.alt_strs().map(Allele::classify)
    }

    /// The writer's typed alleles: an SNV, insertion, deletion or
    /// reference-only site where the bases fit one, else `Complex`.
    pub fn to_alleles(&self) -> Result<Alleles, AllelesError> {
        let reference = self.reference();
        let alts: SmallVec<&str, 2> = self.alt_strs().collect();
        let bases = |s: &[u8]| s.iter().map(|&b| base(b)).collect::<Option<SmallVec<Base, 4>>>();
        let complex = || {
            let reference = std::str::from_utf8(reference).unwrap_or_default();
            Alleles::complex(
                SmolStr::from(reference),
                alts.iter().map(|&a| SmolStr::from(a)).collect(),
            )
        };
        let Some(r) = bases(reference) else { return Ok(complex()) };
        let Some(&anchor) = r.first() else { return Ok(complex()) };
        let alt_bases: Option<SmallVec<SmallVec<Base, 4>, 2>> =
            alts.iter().map(|a| bases(a.as_bytes())).collect();
        let Some(alt_bases) = alt_bases else { return Ok(complex()) };
        match (&*r, &*alt_bases) {
            ([_], []) => Ok(Alleles::reference(anchor)),
            ([_], alts) if alts.iter().all(|a| a.len() == 1) => {
                let alts: SmallVec<Base, 2> =
                    alts.iter().filter_map(|a| a.first().copied()).collect();
                Alleles::snv_multi(anchor, &alts).or_else(|_| Ok(complex()))
            }
            ([_], [alt]) if alt.len() > 1 && alt.first() == Some(&anchor) => {
                Alleles::insertion(anchor, alt.get(1..).unwrap_or_default())
            }
            ([_, deleted @ ..], [alt]) if **alt == [anchor] => Alleles::deletion(anchor, deleted),
            _ => Ok(complex()),
        }
    }
}

fn base(b: u8) -> Option<Base> {
    match b {
        b'A' => Some(Base::A),
        b'C' => Some(Base::C),
        b'G' => Some(Base::G),
        b'T' => Some(Base::T),
        _ => None,
    }
}

/// A record's FILTER column.
#[derive(Debug, Clone, Copy)]
pub struct FiltersRef<'a> {
    column: FilterColumn,
    ids: &'a [DictId],
    file: &'a FileShared,
}

impl<'a> FiltersRef<'a> {
    /// `Pass` (exactly `PASS`), `Missing` (`.`) or `Failed` with the
    /// filters that failed. htslib's `bcf_has_filter` counts `.` as PASS;
    /// here each caller decides.
    ///
    /// A VCF text FILTER name the header does not declare still fails the
    /// record, but no [`FilterId`] can name it, so `Failed` does not list it.
    pub fn status(&self) -> FilterStatus<'a> {
        match self.column {
            FilterColumn::Missing => FilterStatus::Missing,
            FilterColumn::Pass => FilterStatus::Pass,
            FilterColumn::Failed => FilterStatus::Failed(FailedFilters {
                ids: self.ids.iter(),
                table: &self.file.filters,
            }),
        }
    }

    // r[impl vcf_values.header_brand]
    /// Whether the record lists `filter` (`PASS` only when it passed).
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a filter of another header.
    pub fn contains(&self, filter: &FilterId) -> Result<bool, VcfReadError> {
        if let Some(key) = filter.header
            && key != self.file.header.id()
        {
            return Err(VcfReadError::ForeignKey { name: filter.name.clone() });
        }
        Ok(self.ids.iter().any(|id| id.get() == filter.dict_idx))
    }
}

/// Whether a record passed its filters.
#[derive(Debug, Clone, PartialEq)]
pub enum FilterStatus<'a> {
    /// FILTER is exactly `PASS`.
    Pass,
    /// FILTER is `.` (no filters applied).
    Missing,
    /// The filters the record failed.
    Failed(FailedFilters<'a>),
}

impl Eq for FilterStatus<'_> {}

/// The filters a record failed, in record order.
#[derive(Debug, Clone)]
pub struct FailedFilters<'a> {
    ids: std::slice::Iter<'a, DictId>,
    table: &'a [Option<FilterId>],
}

impl<'a> Iterator for FailedFilters<'a> {
    type Item = &'a FilterId;

    fn next(&mut self) -> Option<&'a FilterId> {
        // Decoding checked every id names a FILTER.
        self.ids.find_map(|id| self.table.get(id.get() as usize)?.as_ref())
    }
}

impl PartialEq for FailedFilters<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.clone().eq(other.clone())
    }
}

// r[impl vcf_fields.site+2]
/// A record's ID, REF, ALT, QUAL and FILTER, borrowed.
#[derive(Debug, Clone, Copy)]
pub struct SiteRef<'a> {
    bytes: &'a [u8],
    parts: &'a SiteParts,
    file: &'a FileShared,
}

impl<'a> SiteRef<'a> {
    /// The site of a record whose `bytes` `parts` decoded.
    pub(crate) fn new(bytes: &'a [u8], parts: &'a SiteParts, file: &'a FileShared) -> Self {
        Self { bytes, parts, file }
    }

    /// The IDs (`;`-separated in VCF), none for `.`.
    pub fn ids(&self) -> impl Iterator<Item = &'a str> + Clone + use<'a> {
        let id = text(self.bytes, &self.parts.id);
        let id = if id == "." { "" } else { id };
        id.split(';').filter(|s| !s.is_empty())
    }

    pub fn alleles(&self) -> AllelesRef<'a> {
        AllelesRef { bytes: self.bytes, spans: &self.parts.alleles }
    }

    /// The writer's typed alleles (see [`AllelesRef::to_alleles`]).
    pub fn to_alleles(&self) -> Result<Alleles, AllelesError> {
        self.alleles().to_alleles()
    }

    // r[impl vcf_values.qual]
    /// QUAL as stored, `None` for `.`. Not a `Phred`: QUAL may be negative,
    /// huge or NaN, and reaches the thousands.
    pub fn qual(&self) -> Option<f32> {
        self.parts.qual
    }

    pub fn filters(&self) -> FiltersRef<'a> {
        FiltersRef { column: self.parts.filter, ids: &self.parts.filter_ids, file: self.file }
    }
}

/// A record's ID, REF, ALT, QUAL and FILTER, owned.
#[derive(Debug, Clone)]
pub struct Site {
    /// The ID and allele columns, as stored.
    bytes: Vec<u8>,
    parts: SiteParts,
    file: Arc<FileShared>,
}

impl Site {
    fn view(&self) -> SiteRef<'_> {
        SiteRef::new(&self.bytes, &self.parts, &self.file)
    }

    /// The IDs (`;`-separated in VCF), none for `.`.
    pub fn ids(&self) -> impl Iterator<Item = &str> + Clone + '_ {
        self.view().ids()
    }

    pub fn alleles(&self) -> AllelesRef<'_> {
        self.view().alleles()
    }

    /// The writer's typed alleles (see [`AllelesRef::to_alleles`]).
    pub fn to_alleles(&self) -> Result<Alleles, AllelesError> {
        self.view().to_alleles()
    }

    /// QUAL as stored, `None` for `.` (`r[vcf_values.qual]`).
    pub fn qual(&self) -> Option<f32> {
        self.parts.qual
    }

    pub fn filters(&self) -> FiltersRef<'_> {
        self.view().filters()
    }
}

impl Fields for Site {
    type Plan = ();

    fn plan(_: &VcfHeader, _: &ReadOptions) -> Result<(), PlanError> {
        Ok(())
    }

    fn needs(_: &()) -> Sections {
        Sections::SITE
    }

    fn new(plan: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut site =
            Self { bytes: Vec::new(), parts: SiteParts::default(), file: Arc::clone(raw.file) };
        site.fill(plan, raw)?;
        Ok(site)
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        let (bytes, parts) = raw.site_parts()?;
        let span = parts.span();
        self.bytes.clear();
        self.bytes.extend_from_slice(bytes.get(span.clone()).unwrap_or_default());
        // A record is at most 256 MiB: its offsets fit u32.
        self.parts.rebased_from(parts, u32::try_from(span.start).unwrap_or(u32::MAX));
        raw.keep_file(&mut self.file);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alleles_are_classified() {
        assert_eq!(Allele::classify("ACGT"), Allele::Bases(b"ACGT"));
        assert_eq!(Allele::classify("*"), Allele::Overlapping);
        assert_eq!(Allele::classify("<DEL>"), Allele::Symbolic("<DEL>"));
        assert_eq!(Allele::classify("<*>"), Allele::Symbolic("<*>"));
        assert_eq!(Allele::classify("G]17:198982]"), Allele::Breakend("G]17:198982]"));
        assert_eq!(Allele::classify(".A"), Allele::Breakend(".A"));
        assert_eq!(Allele::classify("A."), Allele::Breakend("A."));
        assert_eq!(Allele::classify("."), Allele::Missing);
        for a in ["ACGT", "*", "<DEL>", "G]17:198982]", "."] {
            assert_eq!(Allele::classify(a).as_bytes(), a.as_bytes());
        }
    }

    // r[verify vcf_fields.site+2]
    #[test]
    fn to_alleles_finds_the_writers_shapes() {
        let typed = |alleles: &[&str]| {
            let bytes = alleles.concat();
            let mut spans = Vec::new();
            let mut at = 0u32;
            for a in alleles {
                let len = u32::try_from(a.len()).unwrap();
                spans.push(at..at + len);
                at += len;
            }
            AllelesRef { bytes: bytes.as_bytes(), spans: &spans }.to_alleles().unwrap()
        };
        assert_eq!(typed(&["A"]), Alleles::reference(Base::A));
        assert_eq!(
            typed(&["A", "T", "C"]),
            Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap()
        );
        assert_eq!(typed(&["A", "ACG"]), Alleles::insertion(Base::A, &[Base::C, Base::G]).unwrap());
        assert_eq!(typed(&["ACG", "A"]), Alleles::deletion(Base::A, &[Base::C, Base::G]).unwrap());
        let complex = |r: &str, alts: &[&str]| {
            Alleles::complex(r.into(), alts.iter().map(|&a| SmolStr::from(a)).collect())
        };
        assert_eq!(typed(&["AC", "GT"]), complex("AC", &["GT"]));
        assert_eq!(typed(&["A", "<DEL>"]), complex("A", &["<DEL>"]));
        assert_eq!(typed(&["A", "A"]), complex("A", &["A"]));
        assert_eq!(typed(&["a", "c"]), complex("a", &["c"]));
        assert_eq!(typed(&["A", "C", "CT"]), complex("A", &["C", "CT"]));
    }
}
