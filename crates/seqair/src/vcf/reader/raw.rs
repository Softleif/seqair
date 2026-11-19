//! [`RawRecord`]: one record as the reader found it, and the borrowed views
//! a [`Fields`](super::Fields) type decodes it through
//! (`r[vcf_fields.user_types]`).

use super::entries::{FormatIndex, LazyIndex};
use super::error::{RecordAt, RecordError, RecordPart, RecordSection, VcfReadError};
use super::fields::{SampleSelection, Sections};
use super::file::FileShared;
use super::info::InfoRef;
use super::locus::{Locus, VcfPos, span_of};
use super::rlen::{RlenColumns, text_rlen};
use super::samples::{FormatSource, FormatView, SampleRef, SamplesRef};
use super::site::{SiteParts, SiteRef, span};
use super::text_value::TextValueError;
use crate::io::text_scan::{find_byte, for_each_byte};
use crate::vcf::header::{SampleIdx, VcfHeader};
use crate::vcf::record_encoder::ContigId;
use seqair_types::SmolStr;
use std::ops::{ControlFlow, Range};
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

/// Where the columns of a VCF text line are.
#[derive(Debug, Clone, Default)]
pub(crate) struct TextCols {
    /// CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
    fixed: [Range<u32>; 8],
    /// Where FORMAT starts, if the line has more than 8 columns.
    format: Option<u32>,
}

/// The fixed columns of a VCF text line, by position.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Col {
    Chrom = 0,
    Pos = 1,
    Id = 2,
    Ref = 3,
    Alt = 4,
    Qual = 5,
    Filter = 6,
    Info = 7,
}

impl TextCols {
    // r[impl vcf_reader.text_scan]
    /// Find the columns of `line` with one vectorised TAB scan that stops
    /// at FORMAT.
    pub(crate) fn find(&mut self, line: &[u8]) -> Result<(), RecordError> {
        let mut n = 0usize;
        let mut start = 0usize;
        self.format = None;
        let fixed = &mut self.fixed;
        let flow = for_each_byte(line, b'\t', |tab| {
            match (fixed.get_mut(n), span(start, tab.saturating_sub(start))) {
                (Some(col), Ok(r)) => *col = r,
                (_, Err(e)) => return ControlFlow::Break(Err(e)),
                (None, Ok(_)) => return ControlFlow::Break(Ok(())),
            }
            n = n.saturating_add(1);
            start = tab.saturating_add(1);
            if n == 8 { ControlFlow::Break(Ok(())) } else { ControlFlow::Continue(()) }
        });
        if let ControlFlow::Break(Err(e)) = flow {
            return Err(e);
        }
        if n < 7 {
            return Err(RecordError::TooFewColumns { found: n.saturating_add(1) });
        }
        if n == 7 {
            let info = span(start, line.len().saturating_sub(start))?;
            if let Some(col) = self.fixed.get_mut(7) {
                *col = info;
            }
        } else {
            self.format = Some(span(start, 0)?.start);
        }
        debug_assert!(n <= 8, "the scan stops at FORMAT");
        Ok(())
    }

    pub(crate) fn col<'a>(&self, line: &'a [u8], col: Col) -> &'a [u8] {
        self.fixed
            .get(col as usize)
            .and_then(|r| line.get(r.start as usize..r.end as usize))
            .unwrap_or_default()
    }

    pub(crate) fn range(&self, col: Col) -> Range<u32> {
        self.fixed.get(col as usize).cloned().unwrap_or_default()
    }

    /// FORMAT and the sample columns after it; empty without them.
    pub(crate) fn rest<'a>(&self, line: &'a [u8]) -> &'a [u8] {
        self.format.and_then(|f| line.get(f as usize..)).unwrap_or_default()
    }
}

/// A VCF text record: its line (without the newline) and columns.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TextRecord<'a> {
    pub(crate) line: &'a [u8],
    pub(crate) cols: &'a TextCols,
}

/// How a record is stored.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Repr<'a> {
    Bcf(BcfRecord<'a>),
    Text(TextRecord<'a>),
}

// r[impl vcf_fields.user_types+2]
/// One record as the reader found it, not yet decoded: what a
/// [`Fields::new`](super::Fields::new) or `fill` reads, through the views
/// it hands out. Only the [`Sections`] the type's `needs` named can be
/// viewed; the others fail with [`VcfReadError::NotRequested`], as their
/// bytes may have been skipped. The INFO and FORMAT views index their
/// section on first access, once per record.
#[derive(Debug, Clone, Copy)]
pub struct RawRecord<'a> {
    pub(crate) file: &'a Arc<FileShared>,
    pub(crate) contig: &'a ContigId,
    pub(crate) at: RecordAt,
    pub(crate) sections: Sections,
    /// The site columns, decoded when [`Sections::SITE`] is read.
    pub(crate) site: Option<&'a SiteParts>,
    /// The reader's indexes of this record's INFO and FORMAT data, for the
    /// views; reset for every record.
    pub(crate) info_index: &'a LazyIndex,
    pub(crate) format_index: &'a LazyIndex<FormatIndex>,
    pub(crate) repr: Repr<'a>,
}

impl<'a> RawRecord<'a> {
    #[cold]
    fn error(&self, source: RecordError) -> VcfReadError {
        self.at.error(self.file, source)
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

    /// The record's contig, start and end.
    ///
    /// # Errors
    /// [`VcfReadError::Record`] for a position or span outside VCF's range,
    /// or a VCF text POS that is not a number.
    pub fn locus(&self) -> Result<Locus, VcfReadError> {
        let (start, end) = self.span()?;
        Ok(Locus::new_unchecked(self.contig.clone(), start, end))
    }

    /// The record's start and inclusive end.
    #[inline]
    pub(crate) fn span(&self) -> Result<(VcfPos, VcfPos), VcfReadError> {
        match self.repr {
            Repr::Bcf(r) => {
                span_of(i64::from(r.fixed.pos), i64::from(r.fixed.rlen)).map_err(|e| self.error(e))
            }
            Repr::Text(t) => self.text_span(t),
        }
    }

    /// A VCF text record's span: POS, which the reader parsed into `at`
    /// already, and the length htslib computes.
    #[inline(never)]
    fn text_span(
        &self,
        TextRecord { line, cols }: TextRecord<'a>,
    ) -> Result<(VcfPos, VcfPos), VcfReadError> {
        self.at
            .pos
            .ok_or_else(|| RecordError::Text {
                part: RecordPart::Pos,
                source: TextValueError::invalid_integer(cols.col(line, Col::Pos)),
            })
            .and_then(|pos0| {
                let rlen = RlenColumns {
                    reference: cols.col(line, Col::Ref),
                    alt: cols.col(line, Col::Alt),
                    info: cols.col(line, Col::Info),
                    format: cols.rest(line),
                };
                span_of(pos0, text_rlen(&rlen, self.file, self.contig, pos0))
            })
            .map_err(|e| self.error(e))
    }

    /// The ID, REF, ALT, QUAL and FILTER columns; needs [`Sections::SITE`].
    ///
    /// # Errors
    /// [`VcfReadError::NotRequested`].
    pub fn site(&self) -> Result<SiteRef<'a>, VcfReadError> {
        let (bytes, parts) = self.site_parts()?;
        Ok(SiteRef::new(bytes, parts, self.file))
    }

    /// The INFO column; needs [`Sections::INFO`].
    ///
    /// # Errors
    /// [`VcfReadError::NotRequested`].
    pub fn info(&self) -> Result<InfoRef<'a>, VcfReadError> {
        let (bytes, n) = self.info_bytes()?;
        Ok(InfoRef::new(bytes, n, self.info_index, self.file, self.at))
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
        let all = self.sections.samples() == SampleSelection::All;
        // The view holds what the reader read: every sample, or (VCF text)
        // the one column `needs` asked for.
        let source = self.format_source(if all { None } else { Some(sample) })?;
        let column = match source {
            FormatSource::Text { .. } if !all => 0,
            _ => sample.index(),
        };
        let view = FormatView::new(source, self.format_index, self.file, self.at);
        Ok(SampleRef::new(view, column, Some(sample)))
    }

    /// Every sample's FORMAT values; needs [`Sections::SAMPLES`].
    ///
    /// # Errors
    /// As [`Self::sample`].
    pub fn samples(&self) -> Result<SamplesRef<'a>, VcfReadError> {
        let source = self.format_source(None)?;
        Ok(SamplesRef::new(FormatView::new(source, self.format_index, self.file, self.at)))
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

    /// The whole record as stored: a BCF record's bytes, or a VCF line.
    pub(crate) fn bytes(&self) -> &'a [u8] {
        match self.repr {
            Repr::Bcf(r) => r.bytes,
            Repr::Text(t) => t.line,
        }
    }

    /// The bytes the site columns were decoded from, and their parts.
    pub(crate) fn site_parts(&self) -> Result<(&'a [u8], &'a SiteParts), VcfReadError> {
        let parts = self.site.ok_or(VcfReadError::NotRequested { section: RecordSection::Site })?;
        let bytes = match self.repr {
            Repr::Bcf(r) => r.shared,
            Repr::Text(t) => t.line,
        };
        Ok((bytes, parts))
    }

    /// The INFO bytes and (BCF) their entry count.
    pub(crate) fn info_bytes(&self) -> Result<(&'a [u8], usize), VcfReadError> {
        self.require(self.sections.info(), RecordSection::Info)?;
        Ok(match self.repr {
            Repr::Bcf(r) => {
                (r.shared.get(r.info_start..).unwrap_or_default(), usize::from(r.fixed.n_info))
            }
            Repr::Text(t) => (t.cols.col(t.line, Col::Info), 0),
        })
    }

    /// The FORMAT data of `sample` (`None`: of all samples). BCF gives the
    /// whole individual section either way; VCF text only that sample's
    /// column, found by scanning no further.
    pub(crate) fn format_source(
        &self,
        sample: Option<SampleIdx>,
    ) -> Result<FormatSource<'a>, VcfReadError> {
        match sample {
            Some(s) => self.require(
                self.sections.has_sample(s.index()),
                RecordSection::Sample { index: s.index() },
            )?,
            None => self
                .require(self.sections.samples() == SampleSelection::All, RecordSection::Samples)?,
        }
        let expected = self.file.header.samples().len();
        match self.repr {
            Repr::Bcf(r) => {
                if r.fixed.n_fmt > 0 && r.fixed.n_sample as usize != expected {
                    let found = r.fixed.n_sample;
                    return Err(self.error(RecordError::SampleCount { found, expected }));
                }
                let n_fmt = usize::from(r.fixed.n_fmt);
                Ok(FormatSource::Bcf { indiv: r.indiv, n_fmt, n_samples: expected })
            }
            Repr::Text(t) => {
                let rest = t.cols.rest(t.line);
                let (keys, samples) = match find_byte(rest, b'\t') {
                    Some(tab) => (
                        rest.get(..tab).unwrap_or_default(),
                        rest.get(tab.saturating_add(1)..).unwrap_or_default(),
                    ),
                    None => (rest, &[][..]),
                };
                let Some(s) = sample else {
                    return Ok(FormatSource::Text { keys, samples, n_samples: expected });
                };
                if keys.is_empty() || keys == b"." {
                    return Ok(FormatSource::Text { keys: b"", samples: b"", n_samples: 1 });
                }
                let column = nth_column(samples, s.index()).ok_or_else(|| {
                    let found = samples.iter().filter(|&&b| b == b'\t').count().saturating_add(1);
                    let found = u32::try_from(found).unwrap_or(u32::MAX);
                    self.error(RecordError::SampleCount { found, expected })
                })?;
                Ok(FormatSource::Text { keys, samples: column, n_samples: 1 })
            }
        }
    }
}

/// Column `n` of the TAB-separated `columns`, scanning no further.
fn nth_column(columns: &[u8], n: usize) -> Option<&[u8]> {
    let mut start = 0usize;
    let mut seen = 0usize;
    let flow = for_each_byte(columns, b'\t', |tab| {
        if seen == n {
            return ControlFlow::Break(start..tab);
        }
        seen = seen.saturating_add(1);
        start = tab.saturating_add(1);
        ControlFlow::Continue(())
    });
    match flow {
        ControlFlow::Break(r) => columns.get(r),
        ControlFlow::Continue(()) if seen == n => columns.get(start..),
        ControlFlow::Continue(()) => None,
    }
}

/// A VCF text POS: decimal digits only.
pub(crate) fn parse_pos(pos: &[u8]) -> Option<i64> {
    if pos.is_empty() || pos.len() > 18 || !pos.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(pos).ok()?.parse().ok()
}
