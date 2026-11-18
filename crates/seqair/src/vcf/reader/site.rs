//! [`Site`]: ID, REF, ALT, QUAL and FILTER (`r[vcf_fields.site]`).

use super::bcf_value::{BcfValueError, BcfVec, TypedValue, skip_typed_value};
use super::error::{KeyKind, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::file::FileShared;
use super::raw::RawRecord;
use crate::vcf::alleles::Alleles;
use crate::vcf::bcf_encoding::FLOAT_MISSING;
use crate::vcf::error::AllelesError;
use crate::vcf::header::VcfHeader;
use crate::vcf::record_encoder::FilterId;
use seqair_types::{Base, SmallVec, SmolStr};
use std::ops::Range;
use std::sync::Arc;

/// One ALT allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Allele<'a> {
    /// Bases (`C`, `CTT`), as stored.
    Bases(&'a [u8]),
    /// A symbolic allele, brackets included: `<DEL>`, `<*>`, `<NON_REF>`.
    Symbolic(&'a str),
    /// A breakend: `G]17:198982]`, `.A`.
    Breakend(&'a str),
    /// `*`: the allele is deleted by an overlapping variant.
    Overlapping,
}

impl<'a> Allele<'a> {
    fn classify(allele: &'a str) -> Self {
        let bytes = allele.as_bytes();
        if bytes == b"*" {
            Self::Overlapping
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

/// The typed strings of a record's alleles, checked to parse as UTF-8.
#[derive(Debug, Clone)]
struct TypedStrings<'a> {
    rest: &'a [u8],
    left: usize,
}

impl<'a> Iterator for TypedStrings<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<&'a str> {
        self.left = self.left.checked_sub(1)?;
        let (value, len) = TypedValue::parse(self.rest).ok()?;
        self.rest = self.rest.get(len..)?;
        value.str().ok()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.left, Some(self.left))
    }
}

impl ExactSizeIterator for TypedStrings<'_> {}

/// A record's alleles: REF, then the ALTs.
#[derive(Debug, Clone, Copy)]
pub struct AllelesRef<'a> {
    bytes: &'a [u8],
    n: u16,
}

impl<'a> AllelesRef<'a> {
    fn strings(&self) -> TypedStrings<'a> {
        TypedStrings { rest: self.bytes, left: usize::from(self.n) }
    }

    /// The REF bases as stored.
    pub fn reference(&self) -> &'a [u8] {
        self.strings().next().map_or(&[], str::as_bytes)
    }

    /// The ALT alleles; none for a record without ALT (`.`).
    pub fn alts(&self) -> impl ExactSizeIterator<Item = Allele<'a>> + Clone + use<'a> {
        self.strings().skip(1).map(Allele::classify)
    }

    /// The writer's typed alleles: an SNV, insertion, deletion or
    /// reference-only site where the bases fit one, else `Complex`.
    pub fn to_alleles(&self) -> Result<Alleles, AllelesError> {
        let reference = self.reference();
        let alts: SmallVec<&str, 2> = self.strings().skip(1).collect();
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
    /// The typed int vector of dictionary indices, checked to name FILTERs.
    bytes: &'a [u8],
    file: Option<&'a FileShared>,
}

impl<'a> FiltersRef<'a> {
    /// The dictionary indices; checked when the record was read, and an
    /// unfilled site has no bytes, which is no filters.
    fn ids(&self) -> Option<BcfVec<'a, i32>> {
        TypedValue::parse(self.bytes).and_then(|(v, _)| v.ints()).ok()
    }

    /// `Pass` (exactly `PASS`), `Missing` (`.`) or `Failed` with the
    /// filters that failed. htslib's `bcf_has_filter` counts `.` as PASS;
    /// here each caller decides.
    pub fn status(&self) -> FilterStatus<'a> {
        let Some(ids) = self.ids() else { return FilterStatus::Missing };
        match (ids.len(), ids.get(0)) {
            (0, _) => FilterStatus::Missing,
            (1, Some(Some(0))) => FilterStatus::Pass,
            _ => FilterStatus::Failed(FailedFilters {
                ids,
                table: self.file.map_or(&[], |f| f.filters.as_slice()),
            }),
        }
    }

    // r[impl vcf_values.header_brand]
    /// Whether the record lists `filter` (`PASS` only when it passed).
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a filter of another header.
    pub fn contains(&self, filter: &FilterId) -> Result<bool, VcfReadError> {
        if let (Some(key), Some(file)) = (filter.header, self.file)
            && key != file.header.id()
        {
            return Err(VcfReadError::ForeignKey { name: filter.name.clone() });
        }
        let wanted = i32::try_from(filter.dict_idx).ok();
        Ok(self.ids().into_iter().flatten().any(|id| id.is_some() && id == wanted))
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
    ids: BcfVec<'a, i32>,
    table: &'a [Option<FilterId>],
}

impl<'a> Iterator for FailedFilters<'a> {
    type Item = &'a FilterId;

    fn next(&mut self) -> Option<&'a FilterId> {
        loop {
            let id = self.ids.next()?;
            let found = id.and_then(|i| self.table.get(usize::try_from(i).ok()?)?.as_ref());
            if found.is_some() {
                return found;
            }
        }
    }
}

impl PartialEq for FailedFilters<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.clone().eq(other.clone())
    }
}

/// Where the site columns sit in a BCF record's shared section.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct SiteSpans {
    pub(crate) id: Range<usize>,
    pub(crate) alleles: Range<usize>,
    pub(crate) filter: Range<usize>,
    pub(crate) info: Range<usize>,
}

impl SiteSpans {
    /// Find the columns in `shared`, the whole shared section, whose site
    /// columns start after the 24 fixed bytes.
    pub(crate) fn find(shared: &[u8], n_allele: u16) -> Result<Self, RecordError> {
        let value = |at: usize, part: RecordPart| -> Result<usize, RecordError> {
            let rest = shared.get(at..).unwrap_or_default();
            let len =
                skip_typed_value(rest).map_err(|source| RecordError::Value { part, source })?;
            Ok(at.saturating_add(len))
        };
        let id_end = value(24, RecordPart::Id)?;
        let mut alleles_end = id_end;
        for _ in 0..n_allele {
            alleles_end = value(alleles_end, RecordPart::Alleles)?;
        }
        let filter_end = value(alleles_end, RecordPart::Filter)?;
        Ok(Self {
            id: 24..id_end,
            alleles: id_end..alleles_end,
            filter: alleles_end..filter_end,
            info: filter_end..shared.len(),
        })
    }
}

// r[impl vcf_fields.site]
/// A record's ID, REF, ALT, QUAL and FILTER, borrowed.
#[derive(Debug, Clone, Copy)]
pub struct SiteRef<'a> {
    id: &'a str,
    alleles: AllelesRef<'a>,
    filters: FiltersRef<'a>,
    qual: u32,
}

impl<'a> SiteRef<'a> {
    /// The site of a BCF record, every part checked: IDs and alleles are
    /// UTF-8, every FILTER index names a FILTER of the header.
    pub(crate) fn check(
        shared: &'a [u8],
        spans: &SiteSpans,
        n_allele: u16,
        qual: u32,
        file: &'a FileShared,
    ) -> Result<Self, RecordError> {
        let part = |r: &Range<usize>| shared.get(r.clone()).unwrap_or_default();
        let value = |part: RecordPart, e: BcfValueError| match e {
            BcfValueError::Utf8(source) => RecordError::Utf8 { part, source },
            source => RecordError::Value { part, source },
        };
        let id = TypedValue::parse(part(&spans.id))
            .and_then(|(v, _)| v.str())
            .map_err(|e| value(RecordPart::Id, e))?;
        let alleles = AllelesRef { bytes: part(&spans.alleles), n: n_allele };
        let mut rest = alleles.bytes;
        for _ in 0..n_allele {
            let (v, len) = TypedValue::parse(rest).map_err(|e| value(RecordPart::Alleles, e))?;
            v.str().map_err(|e| value(RecordPart::Alleles, e))?;
            rest = rest.get(len..).unwrap_or_default();
        }
        let filter_bytes = part(&spans.filter);
        let ids = TypedValue::parse(filter_bytes)
            .and_then(|(v, _)| v.ints())
            .map_err(|e| value(RecordPart::Filter, e))?;
        for index in ids.flatten() {
            let declared = usize::try_from(index).ok().and_then(|i| file.filters.get(i)?.as_ref());
            if declared.is_none() {
                return Err(RecordError::UnknownKey { index, kind: KeyKind::Filter });
            }
        }
        let filters = FiltersRef { bytes: filter_bytes, file: Some(file) };
        Ok(Self { id, alleles, filters, qual })
    }

    /// The IDs (`;`-separated in VCF), none for `.`.
    pub fn ids(&self) -> impl Iterator<Item = &'a str> + Clone + use<'a> {
        let id = if self.id == "." { "" } else { self.id };
        id.split(';').filter(|s| !s.is_empty())
    }

    pub fn alleles(&self) -> AllelesRef<'a> {
        self.alleles
    }

    /// The writer's typed alleles (see [`AllelesRef::to_alleles`]).
    pub fn to_alleles(&self) -> Result<Alleles, AllelesError> {
        self.alleles.to_alleles()
    }

    // r[impl vcf_values.qual]
    /// QUAL as stored, `None` for `.`. Not a `Phred`: QUAL may be negative,
    /// huge or NaN, and reaches the thousands.
    pub fn qual(&self) -> Option<f32> {
        (self.qual != FLOAT_MISSING).then(|| f32::from_bits(self.qual))
    }

    pub fn filters(&self) -> FiltersRef<'a> {
        self.filters
    }
}

/// A record's ID, REF, ALT, QUAL and FILTER, owned.
#[derive(Debug, Clone)]
pub struct Site {
    /// The shared section's site columns, as stored.
    bytes: Vec<u8>,
    spans: SiteSpans,
    n_allele: u16,
    qual: u32,
    file: Option<Arc<FileShared>>,
}

impl Site {
    fn view(&self) -> SiteRef<'_> {
        let part = |r: &Range<usize>| self.bytes.get(r.clone()).unwrap_or_default();
        let id =
            TypedValue::parse(part(&self.spans.id)).and_then(|(v, _)| v.str()).unwrap_or_default();
        SiteRef {
            id,
            alleles: AllelesRef { bytes: part(&self.spans.alleles), n: self.n_allele },
            filters: FiltersRef { bytes: part(&self.spans.filter), file: self.file.as_deref() },
            qual: self.qual,
        }
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
        self.view().qual()
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
        let mut site = Self {
            bytes: Vec::new(),
            spans: SiteSpans::default(),
            n_allele: 0,
            qual: FLOAT_MISSING,
            file: None,
        };
        site.fill(plan, raw)?;
        Ok(site)
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        raw.site()?;
        let (shared, spans, n_allele, qual) = raw.bcf_site()?;
        let start = spans.id.start;
        let end = spans.filter.end;
        self.bytes.clear();
        self.bytes.extend_from_slice(shared.get(start..end).unwrap_or_default());
        let shift = |r: &Range<usize>| r.start.saturating_sub(start)..r.end.saturating_sub(start);
        self.spans = SiteSpans {
            id: shift(&spans.id),
            alleles: shift(&spans.alleles),
            filter: shift(&spans.filter),
            info: 0..0,
        };
        self.n_allele = n_allele;
        self.qual = qual;
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
    }

    fn typed_strings(alleles: &[&str]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for a in alleles {
            crate::vcf::bcf_encoding::encode_typed_string(&mut bytes, a.as_bytes());
        }
        bytes
    }

    // r[verify vcf_fields.site]
    #[test]
    fn to_alleles_finds_the_writers_shapes() {
        let typed = |alleles: &[&str]| {
            let bytes = typed_strings(alleles);
            let n = u16::try_from(alleles.len()).unwrap();
            AllelesRef { bytes: &bytes, n }.to_alleles().unwrap()
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
