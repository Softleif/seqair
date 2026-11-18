//! [`FullRecord`]: every part of a record (`r[vcf_fields.full_record]`).

use super::error::VcfReadError;
use super::fields::{Fields, PlanContext, PlanError, Sections};
use super::info::Info;
use super::locus::Locus;
use super::raw::RawRecord;
use super::samples::Samples;
use super::site::Site;
use crate::vcf::header::HeaderId;

// r[impl vcf_fields.full_record]
/// Every part of a record, plus its bytes as stored and its header's id, so
/// a writer with the same header can copy it unchanged.
#[derive(Clone)]
pub struct FullRecord {
    locus: Locus,
    site: Site,
    info: Info,
    samples: Samples,
    bytes: Vec<u8>,
    header: HeaderId,
}

// r[impl vcf_reader.debug]
impl std::fmt::Debug for FullRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FullRecord")
            .field("locus", &format_args!("{}", self.locus))
            .field("site", &self.site)
            .field("info", &self.info)
            .field("samples", &self.samples)
            .finish()
    }
}

impl FullRecord {
    pub fn locus(&self) -> &Locus {
        &self.locus
    }

    pub fn site(&self) -> &Site {
        &self.site
    }

    pub fn info(&self) -> &Info {
        &self.info
    }

    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// The record as stored: for BCF, `l_shared`, `l_indiv` and both
    /// sections.
    pub fn raw_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// The header the record was read with.
    pub fn header_id(&self) -> HeaderId {
        self.header
    }
}

impl Fields for FullRecord {
    type Plan = HeaderId;

    fn plan(cx: &PlanContext<'_>) -> Result<HeaderId, PlanError> {
        Locus::plan(cx)
    }

    fn needs(_: &HeaderId) -> Sections {
        Sections::SITE | Sections::INFO | Sections::SAMPLES
    }

    fn new(header: &HeaderId, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        Ok(Self {
            locus: Locus::new(header, raw)?,
            site: Site::new(&(), raw)?,
            info: Info::new(&(), raw)?,
            samples: Samples::new(&(), raw)?,
            bytes: raw.bytes().to_vec(),
            header: *header,
        })
    }

    fn fill(&mut self, header: &HeaderId, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.locus.fill(header, raw)?;
        self.site.fill(&(), raw)?;
        self.info.fill(&(), raw)?;
        self.samples.fill(&(), raw)?;
        self.bytes.clear();
        self.bytes.extend_from_slice(raw.bytes());
        self.header = *header;
        Ok(())
    }
}
