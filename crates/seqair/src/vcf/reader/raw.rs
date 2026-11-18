//! [`RawRecord`]: one record as the reader found it, and the borrowed views
//! a [`Fields`](super::Fields) type decodes it through
//! (`r[vcf_fields.user_types]`).

use super::error::{RecordAt, RecordError, RecordSection, VcfReadError};
use super::fields::{SampleSelection, Sections};
use super::file::FileShared;
use super::info::InfoRef;
use super::locus::{Locus, VcfPos};
use super::samples::{SampleRef, SamplesRef};
use super::site::{SiteParts, SiteRef};
use crate::vcf::header::{SampleIdx, VcfHeader};
use crate::vcf::record_encoder::ContigId;
use seqair_types::SmolStr;
use std::sync::Arc;

/// The 24 fixed bytes at the start of a BCF record's shared section.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Fixed {
    pub(crate) rid: i32,
    pub(crate) pos: i32,
    pub(crate) rlen: i32,
    pub(crate) qual: u32,
    pub(crate) n_info: u16,
    pub(crate) n_allele: u16,
    pub(crate) n_sample: u32,
    pub(crate) n_fmt: u8,
}

impl Fixed {
    /// The fixed fields of `shared`; `None` if it is shorter than 24 bytes.
    pub(crate) fn parse(shared: &[u8]) -> Option<Self> {
        let &[
            r0,
            r1,
            r2,
            r3,
            p0,
            p1,
            p2,
            p3,
            l0,
            l1,
            l2,
            l3,
            q0,
            q1,
            q2,
            q3,
            i0,
            i1,
            a0,
            a1,
            s0,
            s1,
            s2,
            n_fmt,
        ] = shared.first_chunk::<24>()?;
        Some(Self {
            rid: i32::from_le_bytes([r0, r1, r2, r3]),
            pos: i32::from_le_bytes([p0, p1, p2, p3]),
            rlen: i32::from_le_bytes([l0, l1, l2, l3]),
            qual: u32::from_le_bytes([q0, q1, q2, q3]),
            n_info: u16::from_le_bytes([i0, i1]),
            n_allele: u16::from_le_bytes([a0, a1]),
            n_sample: u32::from_le_bytes([s0, s1, s2, 0]),
            n_fmt,
        })
    }
}

/// A BCF record's bytes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BcfRecord<'a> {
    /// The whole record as stored: `l_shared`, `l_indiv`, the shared
    /// section and the individual section (empty when no sample is read).
    pub(crate) bytes: &'a [u8],
    pub(crate) shared: &'a [u8],
    pub(crate) indiv: &'a [u8],
    pub(crate) fixed: Fixed,
    /// Where INFO starts in `shared`; found only when it is read.
    pub(crate) info_start: usize,
}

/// How a record is stored. VCF text records are the next phase: a variant
/// holding the line and its column offsets.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Repr<'a> {
    Bcf(BcfRecord<'a>),
}

// r[impl vcf_fields.user_types]
/// One record as the reader found it, not yet decoded: what a
/// [`Fields::new`](super::Fields::new) or `fill` reads, through the views it hands
/// out. Only the [`Sections`] the type's `needs` named can be viewed; the
/// others fail with [`VcfReadError::NotRequested`], as their bytes may have
/// been skipped.
#[derive(Debug, Clone, Copy)]
pub struct RawRecord<'a> {
    pub(crate) file: &'a Arc<FileShared>,
    pub(crate) contig: &'a ContigId,
    pub(crate) at: RecordAt,
    pub(crate) sections: Sections,
    /// The site columns, decoded when [`Sections::SITE`] is read.
    pub(crate) site: Option<&'a SiteParts>,
    pub(crate) repr: Repr<'a>,
}

impl<'a> RawRecord<'a> {
    fn bcf(&self) -> &BcfRecord<'a> {
        match &self.repr {
            Repr::Bcf(record) => record,
        }
    }

    fn error(&self, source: RecordError) -> VcfReadError {
        self.at.error(source)
    }

    fn require(&self, ok: bool, section: RecordSection) -> Result<(), VcfReadError> {
        if ok { Ok(()) } else { Err(VcfReadError::NotRequested { section }) }
    }

    /// The header of the file the record is from.
    pub fn header(&self) -> &'a VcfHeader {
        &self.file.header
    }

    /// The record's contig.
    pub fn contig(&self) -> &'a ContigId {
        self.contig
    }

    /// The record's start and inclusive end.
    pub(crate) fn span(&self) -> Result<(VcfPos, VcfPos), VcfReadError> {
        let l = self.locus()?;
        Ok((l.start(), l.end()))
    }

    /// The record's contig, start and end.
    ///
    /// # Errors
    /// [`VcfReadError::Record`] for a position or span outside VCF's range.
    pub fn locus(&self) -> Result<Locus, VcfReadError> {
        let f = self.bcf().fixed;
        Locus::from_bcf(self.contig.clone(), f.pos, f.rlen).map_err(|e| self.error(e))
    }

    /// The ID, REF, ALT, QUAL and FILTER columns; needs [`Sections::SITE`].
    ///
    /// # Errors
    /// [`VcfReadError::NotRequested`], or [`VcfReadError::Record`] if a
    /// column is malformed.
    pub fn site(&self) -> Result<SiteRef<'a>, VcfReadError> {
        let (bytes, parts) = self.site_parts()?;
        Ok(SiteRef::new(bytes, parts, self.file))
    }

    /// The INFO column; needs [`Sections::INFO`]. Each `get` walks the
    /// column: a type that reads many fields should keep an [`Info`](super::Info).
    ///
    /// # Errors
    /// [`VcfReadError::NotRequested`].
    pub fn info(&self) -> Result<InfoRef<'a>, VcfReadError> {
        let (bytes, n) = self.bcf_info()?;
        Ok(InfoRef::scan(bytes, n, self.file, self.at))
    }

    /// One sample's FORMAT values; needs that sample or all of them.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a sample of another header,
    /// [`VcfReadError::NotRequested`], or [`VcfReadError::Record`] if the
    /// record's sample count is not the header's.
    pub fn sample(&self, sample: SampleIdx) -> Result<SampleRef<'a>, VcfReadError> {
        if sample.header != self.file.header.id() {
            let name = SmolStr::from(format!("sample {}", sample.index()));
            return Err(VcfReadError::ForeignKey { name });
        }
        let (bytes, n_fmt, n_samples) = self.bcf_samples(Some(sample))?;
        Ok(SamplesRef::scan(bytes, n_fmt, n_samples, self.file, self.at).sample(sample))
    }

    /// Every sample's FORMAT values; needs [`Sections::SAMPLES`].
    ///
    /// # Errors
    /// As [`Self::sample`].
    pub fn samples(&self) -> Result<SamplesRef<'a>, VcfReadError> {
        let (bytes, n_fmt, n_samples) = self.bcf_samples(None)?;
        Ok(SamplesRef::scan(bytes, n_fmt, n_samples, self.file, self.at))
    }

    pub(crate) fn at(&self) -> RecordAt {
        self.at
    }

    /// Point `slot` at the record's file, cloning the `Arc` only when it
    /// changes, so filling a component does not touch the refcount.
    pub(crate) fn keep_file(&self, slot: &mut Arc<FileShared>) {
        if !Arc::ptr_eq(slot, self.file) {
            *slot = Arc::clone(self.file);
        }
    }

    /// The whole record as stored.
    pub(crate) fn bytes(&self) -> &'a [u8] {
        self.bcf().bytes
    }

    /// The bytes the site columns were decoded from, and their parts.
    pub(crate) fn site_parts(&self) -> Result<(&'a [u8], &'a SiteParts), VcfReadError> {
        let parts = self.site.ok_or(VcfReadError::NotRequested { section: RecordSection::Site })?;
        Ok((self.bcf().shared, parts))
    }

    /// The INFO entries and their count.
    pub(crate) fn bcf_info(&self) -> Result<(&'a [u8], usize), VcfReadError> {
        self.require(self.sections.info(), RecordSection::Info)?;
        let r = self.bcf();
        Ok((r.shared.get(r.info_start..).unwrap_or_default(), usize::from(r.fixed.n_info)))
    }

    /// The individual section, `n_fmt` and the sample count, for `sample`
    /// (`None`: all samples).
    pub(crate) fn bcf_samples(
        &self,
        sample: Option<SampleIdx>,
    ) -> Result<(&'a [u8], usize, usize), VcfReadError> {
        match sample {
            Some(s) => self.require(
                self.sections.has_sample(s.index()),
                RecordSection::Sample { index: s.index() },
            )?,
            None => self
                .require(self.sections.samples() == SampleSelection::All, RecordSection::Samples)?,
        }
        let r = self.bcf();
        let expected = self.file.header.samples().len();
        if r.fixed.n_fmt > 0 && r.fixed.n_sample as usize != expected {
            return Err(self.error(RecordError::SampleCount { found: r.fixed.n_sample, expected }));
        }
        Ok((r.indiv, usize::from(r.fixed.n_fmt), expected))
    }
}
