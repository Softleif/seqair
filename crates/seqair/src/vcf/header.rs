//! VCF header model with builder, field definitions, and BCF string map.

use super::error::VcfHeaderError;
use crate::bam::header::BamHeader;
use indexmap::IndexMap;
use seqair_types::SmolStr;
use std::marker::PhantomData;

// r[impl vcf_header.info_def]
/// Number of values for an INFO or FORMAT field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Number {
    /// Fixed count.
    Count(u32),
    /// One per ALT allele (A).
    AlternateBases,
    /// One per allele including REF (R).
    ReferenceAlternateBases,
    /// One per genotype combination (G).
    Genotypes,
    /// One per possible base modification (M). VCF 4.2+ extension.
    BaseModification,
    /// Unknown/variable count (.).
    Unknown,
}

impl Number {
    /// VCF text representation.
    pub fn as_str(&self) -> SmolStr {
        match self {
            Self::Count(n) => SmolStr::from(n.to_string()),
            Self::AlternateBases => SmolStr::from("A"),
            Self::ReferenceAlternateBases => SmolStr::from("R"),
            Self::Genotypes => SmolStr::from("G"),
            Self::BaseModification => SmolStr::from("M"),
            Self::Unknown => SmolStr::from("."),
        }
    }
}

// r[impl vcf_header.info_def]
/// Value type for an INFO or FORMAT field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl ValueType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Integer => "Integer",
            Self::Float => "Float",
            Self::Flag => "Flag",
            Self::Character => "Character",
            Self::String => "String",
        }
    }
}

/// INFO field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InfoDef {
    pub number: Number,
    pub typ: ValueType,
    pub description: SmolStr,
}

/// FORMAT field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatDef {
    pub number: Number,
    pub typ: ValueType,
    pub description: SmolStr,
}

/// FILTER definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterDef {
    pub description: SmolStr,
}

/// Contig definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigDef {
    pub length: Option<u64>,
}

// r[impl vcf_header.string_map]
/// Ordered mapping of header string IDs to BCF dictionary indices.
///
/// The dictionary is small (typically 10–30 entries), so a `Vec` with linear
/// scan is faster than hashing and uses less memory.
#[derive(Debug, Clone)]
pub struct StringMap {
    entries: Vec<SmolStr>,
}

impl Default for StringMap {
    fn default() -> Self {
        Self::new()
    }
}

impl StringMap {
    fn new() -> Self {
        Self { entries: Vec::new() }
    }

    fn insert(&mut self, id: SmolStr) -> usize {
        if let Some(idx) = self.get(&id) {
            return idx;
        }
        let idx = self.entries.len();
        self.entries.push(id);
        idx
    }

    /// Look up the BCF dictionary index for an ID.
    pub fn get(&self, id: &str) -> Option<usize> {
        self.entries.iter().position(|s| s == id)
    }
}

// r[impl vcf_header.file_format]
// r[impl vcf_header.ordered_maps]
// r[impl vcf_header.smolstr]
/// VCF header with typed field definitions in insertion order.
#[derive(Debug, Clone)]
pub struct VcfHeader {
    file_format: SmolStr,
    infos: IndexMap<SmolStr, InfoDef>,
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
    string_map: StringMap,
}

impl VcfHeader {
    /// Start building a new header.
    pub fn builder() -> VcfHeaderBuilder {
        VcfHeaderBuilder::new()
    }

    pub fn file_format(&self) -> &str {
        &self.file_format
    }

    pub fn infos(&self) -> &IndexMap<SmolStr, InfoDef> {
        &self.infos
    }

    pub fn formats(&self) -> &IndexMap<SmolStr, FormatDef> {
        &self.formats
    }

    pub fn filters(&self) -> &IndexMap<SmolStr, FilterDef> {
        &self.filters
    }

    pub fn contigs(&self) -> &IndexMap<SmolStr, ContigDef> {
        &self.contigs
    }

    pub fn samples(&self) -> &[SmolStr] {
        &self.samples
    }

    /// Look up the integer index for a contig name (for BCF encoding).
    pub fn contig_id(&self, name: &str) -> Result<usize, VcfHeaderError> {
        self.contigs
            .get_index_of(name)
            .ok_or_else(|| VcfHeaderError::MissingContig { name: SmolStr::from(name) })
    }

    /// BCF string dictionary.
    pub fn string_map(&self) -> &StringMap {
        &self.string_map
    }

    // r[impl vcf_header.serialization]
    /// Serialize as VCF header text (all ## lines + #CHROM line).
    pub fn to_vcf_text(&self) -> String {
        let mut out = String::new();

        // fileformat first
        out.push_str("##fileformat=");
        out.push_str(&self.file_format);
        out.push('\n');

        // PASS filter — always emitted immediately after fileformat, before metadata.
        // This matches htslib's convention and ensures PASS gets BCF dict index 0.
        if let Some(pass_def) = self.filters.get("PASS") {
            out.push_str("##FILTER=<ID=PASS,Description=\"");
            out.push_str(&pass_def.description);
            out.push_str("\">\n");
        }

        // Other lines (metadata such as ##rastairVersion, ##rastairCommand) come next.
        for line in &self.other_lines {
            out.push_str("##");
            out.push_str(line);
            out.push('\n');
        }

        // Remaining FILTER definitions (skip PASS, already emitted above).
        // r[impl vcf_header.serialization]
        for (id, def) in &self.filters {
            if id == "PASS" {
                continue;
            }
            out.push_str("##FILTER=<ID=");
            out.push_str(id);
            out.push_str(",Description=\"");
            out.push_str(&def.description);
            out.push_str("\">\n");
        }

        // INFO
        for (id, def) in &self.infos {
            out.push_str("##INFO=<ID=");
            out.push_str(id);
            out.push_str(",Number=");
            out.push_str(&def.number.as_str());
            out.push_str(",Type=");
            out.push_str(def.typ.as_str());
            out.push_str(",Description=\"");
            out.push_str(&def.description);
            out.push_str("\">\n");
        }

        // FORMAT
        for (id, def) in &self.formats {
            out.push_str("##FORMAT=<ID=");
            out.push_str(id);
            out.push_str(",Number=");
            out.push_str(&def.number.as_str());
            out.push_str(",Type=");
            out.push_str(def.typ.as_str());
            out.push_str(",Description=\"");
            out.push_str(&def.description);
            out.push_str("\">\n");
        }

        // contig
        for (name, def) in &self.contigs {
            out.push_str("##contig=<ID=");
            out.push_str(name);
            if let Some(len) = def.length {
                out.push_str(",length=");
                out.push_str(&len.to_string());
            }
            out.push_str(">\n");
        }

        // #CHROM line
        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.push_str("\tFORMAT");
            for sample in &self.samples {
                out.push('\t');
                out.push_str(sample);
            }
        }
        out.push('\n');

        out
    }
}

// ── Builder typestate phases ──────────────────────────────────────────
//
// The builder progresses through phases: Contigs → Filters → Infos →
// Formats → Samples.  Each phase only exposes methods for that stage.
// The BCF string dictionary is built incrementally in the guaranteed-
// correct order (FILTER → INFO → FORMAT) — no runtime order check needed.

/// Result of [`VcfHeaderBuilder::from_bam_header`]: a builder pre-populated
/// with contigs and the resolved [`ContigId`](super::record_encoder::ContigId)
/// handles needed by the record encoder.
#[derive(Debug, Clone)]
pub struct FromBamHeader {
    /// The builder, ready for further field registration.
    pub builder: VcfHeaderBuilder<Contigs>,
    /// One [`ContigId`](super::record_encoder::ContigId) per `@SQ` line, in BAM tid order.
    pub contigs: Vec<super::record_encoder::ContigId>,
}

impl FromBamHeader {
    /// Look up the [`ContigId`](super::record_encoder::ContigId) registered
    /// for `tid`.
    ///
    /// `tid` is a validated index into the BAM header that produced this
    /// `FromBamHeader`, so the lookup succeeds whenever both come from the
    /// same header. Returns `None` if the `tid` was resolved against a
    /// different header (mismatched `Readers` foot-gun).
    #[must_use]
    pub fn contig(&self, tid: crate::reader::Tid) -> Option<&super::record_encoder::ContigId> {
        // tid.as_u32() is bounded by the header's target count (Tid invariant).
        // u32 -> usize widens on 64-bit and is equal-width on 32-bit; never
        // truncates on supported platforms.
        self.contigs.get(tid.as_u32() as usize)
    }
}

/// Initial phase: set file format, contigs, and metadata.
#[derive(Debug, Clone, Copy)]
pub struct Contigs;
/// Register FILTER definitions (after contigs, before INFO).
#[derive(Debug, Clone, Copy)]
pub struct Filters;
/// Register INFO field definitions (after filters, before FORMAT).
#[derive(Debug, Clone, Copy)]
pub struct Infos;
/// Register FORMAT field definitions (after INFO, before samples).
#[derive(Debug, Clone, Copy)]
pub struct Formats;
/// Add sample names (final phase before build).
#[derive(Debug, Clone, Copy)]
pub struct Samples;

// r[impl vcf_header.builder]
/// Typestate builder for [`VcfHeader`].
///
/// Fields must be registered in BCF string dictionary order:
/// contigs → filters → INFO → FORMAT → samples.  The `Phase` type
/// parameter tracks the current stage and prevents out-of-order
/// registration at compile time.
///
/// Transition methods ([`filters`](VcfHeaderBuilder::filters),
/// [`infos`](VcfHeaderBuilder::infos), etc.) consume `self` and return
/// the next phase.  Phases can be skipped (e.g. go from `Contigs`
/// straight to `Infos`).  [`build`](VcfHeaderBuilder::build) is
/// available from every phase.
///
/// Out-of-order registration is a compile error:
///
/// ```compile_fail
/// use seqair::vcf::{VcfHeader, Number, ValueType};
/// use seqair::vcf::record_encoder::{InfoFieldDef, Scalar};
///
/// let mut builder = VcfHeader::builder();
/// // ERROR: register_info is not available in the Contigs phase
/// builder.register_info(&InfoFieldDef::<Scalar<i32>>::new(
///     "DP", Number::Count(1), ValueType::Integer, "Depth",
/// )).unwrap();
/// ```
pub struct VcfHeaderBuilder<Phase = Contigs> {
    file_format: SmolStr,
    infos: IndexMap<SmolStr, InfoDef>,
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
    string_map: StringMap,
    _phase: PhantomData<Phase>,
}

// Manual impls to avoid Phase: Debug/Clone bounds from derive macros.
impl<P> std::fmt::Debug for VcfHeaderBuilder<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VcfHeaderBuilder")
            .field("file_format", &self.file_format)
            .field("infos", &self.infos)
            .field("formats", &self.formats)
            .field("filters", &self.filters)
            .field("contigs", &self.contigs)
            .field("samples", &self.samples)
            .field("other_lines", &self.other_lines)
            .finish()
    }
}

impl<P> Clone for VcfHeaderBuilder<P> {
    fn clone(&self) -> Self {
        Self {
            file_format: self.file_format.clone(),
            infos: self.infos.clone(),
            formats: self.formats.clone(),
            filters: self.filters.clone(),
            contigs: self.contigs.clone(),
            samples: self.samples.clone(),
            other_lines: self.other_lines.clone(),
            string_map: self.string_map.clone(),
            _phase: PhantomData,
        }
    }
}

// ── Methods available on every phase ─────────────────────────────────

impl<P> VcfHeaderBuilder<P> {
    /// Transition to a different phase (private — callers use named methods).
    fn into_phase<Q>(self) -> VcfHeaderBuilder<Q> {
        VcfHeaderBuilder {
            file_format: self.file_format,
            infos: self.infos,
            formats: self.formats,
            filters: self.filters,
            contigs: self.contigs,
            samples: self.samples,
            other_lines: self.other_lines,
            string_map: self.string_map,
            _phase: PhantomData,
        }
    }

    /// Insert an ID into the string map, returning the dict index.
    fn insert_string_map_entry(&mut self, id: &SmolStr) -> Result<u32, VcfHeaderError> {
        u32::try_from(self.string_map.insert(id.clone())).map_err(|_| VcfHeaderError::TooManyFields)
    }

    // r[impl vcf_header.file_format]
    /// Set the VCF file format version (default: VCFv4.3).
    ///
    /// Available from any phase — metadata does not affect BCF string
    /// dictionary ordering.
    pub fn file_format(&mut self, version: impl Into<SmolStr>) {
        self.file_format = version.into();
    }

    /// Append a `##key=value` metadata line (e.g. `##source=myapp`).
    ///
    /// Available from any phase — metadata does not affect BCF string
    /// dictionary ordering.
    pub fn add_other_line(&mut self, line: impl Into<SmolStr>) {
        self.other_lines.push(line.into());
    }

    // r[impl vcf_header.builder]
    // r[impl vcf_header.string_map]
    /// Finalize the builder and produce a [`VcfHeader`].
    ///
    /// Available from any phase — all fields after contigs are optional.
    #[must_use = "build() returns the header; ignoring it discards all configuration"]
    pub fn build(self) -> Result<VcfHeader, VcfHeaderError> {
        Ok(VcfHeader {
            file_format: self.file_format,
            infos: self.infos,
            formats: self.formats,
            filters: self.filters,
            contigs: self.contigs,
            samples: self.samples,
            other_lines: self.other_lines,
            string_map: self.string_map,
        })
    }
}

// ── Contigs phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Contigs> {
    fn new() -> Self {
        // PASS is always the first filter and dict index 0.
        let mut filters = IndexMap::new();
        filters.insert(
            SmolStr::from("PASS"),
            FilterDef { description: SmolStr::from("All filters passed") },
        );
        let mut string_map = StringMap::new();
        string_map.insert(SmolStr::from("PASS"));

        Self {
            file_format: SmolStr::from("VCFv4.3"),
            infos: IndexMap::new(),
            formats: IndexMap::new(),
            filters,
            contigs: IndexMap::new(),
            samples: Vec::new(),
            other_lines: Vec::new(),
            string_map,
            _phase: PhantomData,
        }
    }

    // r[impl vcf_header.contig_required]
    pub fn add_contig(
        &mut self,
        name: impl Into<SmolStr>,
        def: ContigDef,
    ) -> Result<(), VcfHeaderError> {
        let name = name.into();
        if self.contigs.contains_key(&name) {
            return Err(VcfHeaderError::DuplicateContig { name });
        }
        self.contigs.insert(name, def);
        Ok(())
    }

    // r[impl record_encoder.register_contig]
    /// Register a contig: adds the header entry and returns a resolved [`ContigId`](super::record_encoder::ContigId).
    pub fn register_contig(
        &mut self,
        name: impl Into<SmolStr>,
        def: ContigDef,
    ) -> Result<super::record_encoder::ContigId, VcfHeaderError> {
        let name = name.into();
        if self.contigs.contains_key(&name) {
            return Err(VcfHeaderError::DuplicateContig { name });
        }
        let tid = u32::try_from(self.contigs.len()).map_err(|_| VcfHeaderError::TooManyContigs)?;
        self.contigs.insert(name.clone(), def);
        Ok(super::record_encoder::ContigId { tid, name })
    }

    // r[impl vcf_header.from_bam_header]
    /// Build a `VcfHeaderBuilder` from a BAM header, copying contig names and lengths.
    ///
    /// Returns a [`FromBamHeader`] containing the builder and all resolved
    /// [`ContigId`](super::record_encoder::ContigId) handles, so callers
    /// can use them with the record encoder without re-registering.
    pub fn from_bam_header(header: &BamHeader) -> Result<FromBamHeader, VcfHeaderError> {
        let mut builder = Self::new();
        let mut contigs = Vec::new();
        for (tid, name) in header.target_names().enumerate() {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "BAM header n_ref ≤ 1M (enforced at parse), so tid fits in u32"
            )]
            let length = header.target_len(tid as u32);
            let def = ContigDef { length };
            contigs.push(builder.register_contig(name, def)?);
        }
        Ok(FromBamHeader { builder, contigs })
    }

    /// Advance to the [`Filters`] phase.
    pub fn filters(self) -> VcfHeaderBuilder<Filters> {
        self.into_phase()
    }

    /// Skip filters, advance to the [`Infos`] phase.
    pub fn infos(self) -> VcfHeaderBuilder<Infos> {
        self.into_phase()
    }

    /// Skip filters and infos, advance to the [`Formats`] phase.
    pub fn formats(self) -> VcfHeaderBuilder<Formats> {
        self.into_phase()
    }

    /// Skip filters, infos, and formats, advance to the [`Samples`] phase.
    pub fn samples(self) -> VcfHeaderBuilder<Samples> {
        self.into_phase()
    }
}

// ── Filters phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Filters> {
    // r[impl vcf_header.filter_def]
    /// Insert a FILTER definition, returning its BCF dict index.
    fn insert_filter(&mut self, id: SmolStr, def: FilterDef) -> Result<u32, VcfHeaderError> {
        if self.filters.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFilter { id });
        }
        self.filters.insert(id.clone(), def);
        self.insert_string_map_entry(&id)
    }

    pub fn add_filter(
        &mut self,
        id: impl Into<SmolStr>,
        def: FilterDef,
    ) -> Result<(), VcfHeaderError> {
        self.insert_filter(id.into(), def)?;
        Ok(())
    }

    // r[impl record_encoder.register]
    /// Register a FILTER: adds the header entry and returns a resolved [`FilterId`](super::record_encoder::FilterId).
    pub fn register_filter(
        &mut self,
        def: &super::record_encoder::FilterFieldDef,
    ) -> Result<super::record_encoder::FilterId, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        let dict_idx = self
            .insert_filter(id.clone(), FilterDef { description: SmolStr::from(def.description) })?;
        Ok(super::record_encoder::FilterId { dict_idx, name: id })
    }

    /// Advance to the [`Infos`] phase.
    pub fn infos(self) -> VcfHeaderBuilder<Infos> {
        self.into_phase()
    }

    /// Skip infos, advance to the [`Formats`] phase.
    pub fn formats(self) -> VcfHeaderBuilder<Formats> {
        self.into_phase()
    }

    /// Skip infos and formats, advance to the [`Samples`] phase.
    pub fn samples(self) -> VcfHeaderBuilder<Samples> {
        self.into_phase()
    }
}

// ── Infos phase ──────────────────────────────────────────────────────

impl VcfHeaderBuilder<Infos> {
    // r[impl vcf_header.info_def]
    /// Insert an INFO definition, returning its BCF dict index.
    fn insert_info(&mut self, id: SmolStr, def: InfoDef) -> Result<u32, VcfHeaderError> {
        if self.infos.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateInfo { id });
        }
        if def.typ == ValueType::Flag && def.number != Number::Count(0) {
            return Err(VcfHeaderError::FlagNumberMismatch { id });
        }
        self.infos.insert(id.clone(), def);
        self.insert_string_map_entry(&id)
    }

    pub fn add_info(&mut self, id: impl Into<SmolStr>, def: InfoDef) -> Result<(), VcfHeaderError> {
        self.insert_info(id.into(), def)?;
        Ok(())
    }

    // r[impl record_encoder.register]
    /// Register an INFO field: adds the header entry and returns a resolved typed key.
    pub fn register_info<V>(
        &mut self,
        def: &super::record_encoder::InfoFieldDef<V>,
    ) -> Result<super::record_encoder::InfoKey<V>, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        let dict_idx = self.insert_info(
            id.clone(),
            InfoDef {
                number: def.number,
                typ: def.value_type,
                description: SmolStr::from(def.description),
            },
        )?;
        Ok(super::record_encoder::InfoKey(
            super::record_encoder::FieldId { dict_idx, name: id },
            PhantomData,
        ))
    }

    /// Advance to the [`Formats`] phase.
    pub fn formats(self) -> VcfHeaderBuilder<Formats> {
        self.into_phase()
    }

    /// Skip formats, advance to the [`Samples`] phase.
    pub fn samples(self) -> VcfHeaderBuilder<Samples> {
        self.into_phase()
    }
}

// ── Formats phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Formats> {
    // r[impl vcf_header.format_def]
    /// Insert a FORMAT definition, returning its BCF dict index.
    fn insert_format(&mut self, id: SmolStr, def: FormatDef) -> Result<u32, VcfHeaderError> {
        if self.formats.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFormat { id });
        }
        if def.typ == ValueType::Flag {
            return Err(VcfHeaderError::FormatFlagNotAllowed { id });
        }
        self.formats.insert(id.clone(), def);
        self.insert_string_map_entry(&id)
    }

    pub fn add_format(
        &mut self,
        id: impl Into<SmolStr>,
        def: FormatDef,
    ) -> Result<(), VcfHeaderError> {
        self.insert_format(id.into(), def)?;
        Ok(())
    }

    // r[impl record_encoder.register]
    /// Register a FORMAT field: adds the header entry and returns a resolved typed key.
    pub fn register_format<V>(
        &mut self,
        def: &super::record_encoder::FormatFieldDef<V>,
    ) -> Result<super::record_encoder::FormatKey<V>, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        let dict_idx = self.insert_format(
            id.clone(),
            FormatDef {
                number: def.number,
                typ: def.value_type,
                description: SmolStr::from(def.description),
            },
        )?;
        Ok(super::record_encoder::FormatKey(
            super::record_encoder::FieldId { dict_idx, name: id },
            PhantomData,
        ))
    }

    /// Advance to the [`Samples`] phase.
    pub fn samples(self) -> VcfHeaderBuilder<Samples> {
        self.into_phase()
    }
}

// ── Samples phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Samples> {
    // r[impl vcf_header.sample_names]
    pub fn add_sample(&mut self, name: impl Into<SmolStr>) -> Result<(), VcfHeaderError> {
        let name = name.into();
        if self.samples.contains(&name) {
            return Err(VcfHeaderError::DuplicateSample { name });
        }
        self.samples.push(name);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_header.from_bam_header]
    /// `FromBamHeader::contig(tid)` must return the same `ContigId` as
    /// indexing `from_bam.contigs[tid as usize]`, and `None` for an
    /// out-of-range tid (the cross-header foot-gun).
    #[test]
    fn from_bam_header_contig_lookup() {
        let sam =
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chrM\tLN:16000\n";
        let header = BamHeader::from_sam_text(sam).unwrap();
        let from_bam = VcfHeaderBuilder::from_bam_header(&header).unwrap();

        for tid_u32 in 0..3u32 {
            let tid = crate::reader::ResolveTid::resolve_tid(&tid_u32, &header).unwrap();
            let via_helper = from_bam.contig(tid).expect("tid in range");
            let via_index = &from_bam.contigs[tid_u32 as usize];
            // ContigId doesn't impl PartialEq; compare its public fields.
            assert_eq!(via_helper.tid(), via_index.tid());
            assert_eq!(via_helper.name(), via_index.name());
        }

        // Out-of-range tid (e.g. resolved against a different header). Build
        // a Tid with a numeric value past the contig list so the lookup
        // misses, but stay within u32.
        let bigger =
            "@HD\tVN:1.6\n@SQ\tSN:a\tLN:1\n@SQ\tSN:b\tLN:1\n@SQ\tSN:c\tLN:1\n@SQ\tSN:d\tLN:1\n";
        let bigger_header = BamHeader::from_sam_text(bigger).unwrap();
        let stale_tid = crate::reader::ResolveTid::resolve_tid(&3u32, &bigger_header).unwrap();
        assert!(from_bam.contig(stale_tid).is_none());
    }

    // r[verify vcf_header.builder]
    // r[verify vcf_header.file_format]
    #[test]
    fn basic_header_build() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef { length: Some(248_956_422) }).unwrap();
        let mut builder = builder.infos();
        builder
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Total Depth"),
                },
            )
            .unwrap();
        let mut builder = builder.formats();
        builder
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("sample0").unwrap();
        let header = builder.build().unwrap();

        assert_eq!(header.file_format(), "VCFv4.3");
        assert!(header.contigs().contains_key("chr1"));
        assert!(header.infos().contains_key("DP"));
        assert!(header.formats().contains_key("GT"));
        assert_eq!(header.samples(), &[SmolStr::from("sample0")]);
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_contig() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef { length: None }).unwrap();
        let result = builder.add_contig("chr1", ContigDef { length: None });
        assert!(result.is_err());
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_info() {
        let mut builder = VcfHeader::builder().infos();
        builder
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("x"),
                },
            )
            .unwrap();
        let result = builder.add_info(
            "DP",
            InfoDef {
                number: Number::Count(1),
                typ: ValueType::Integer,
                description: SmolStr::from("y"),
            },
        );
        assert!(result.is_err());
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_sample() {
        let mut builder = VcfHeader::builder().samples();
        builder.add_sample("s1").unwrap();
        let result = builder.add_sample("s1");
        assert!(result.is_err());
    }

    // r[verify vcf_header.filter_def]
    #[test]
    fn pass_filter_auto_inserted() {
        let header = VcfHeader::builder().build().unwrap();
        assert!(header.filters().contains_key("PASS"));
        // PASS should be first (index 0)
        let first = header.filters().keys().next().unwrap();
        assert_eq!(first.as_str(), "PASS");
    }

    // r[verify vcf_header.string_map]
    #[test]
    fn string_map_pass_is_zero() {
        let mut builder = VcfHeader::builder().infos();
        builder
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap();
        let header = builder.build().unwrap();

        assert_eq!(header.string_map().get("PASS"), Some(0));
        assert_eq!(header.string_map().get("DP"), Some(1));
    }

    // r[verify vcf_header.info_def]
    #[test]
    fn flag_requires_number_zero() {
        let mut builder = VcfHeader::builder().infos();
        let result = builder.add_info(
            "DB",
            InfoDef {
                number: Number::Count(1),
                typ: ValueType::Flag,
                description: SmolStr::from("dbSNP"),
            },
        );
        assert!(result.is_err());

        // Number=0 + Flag should work
        let result = builder.add_info(
            "DB",
            InfoDef {
                number: Number::Count(0),
                typ: ValueType::Flag,
                description: SmolStr::from("dbSNP"),
            },
        );
        assert!(result.is_ok());
    }

    // r[verify vcf_header.format_def]
    #[test]
    fn format_rejects_flag() {
        let mut builder = VcfHeader::builder().formats();
        let result = builder.add_format(
            "X",
            FormatDef {
                number: Number::Count(0),
                typ: ValueType::Flag,
                description: SmolStr::from("bad"),
            },
        );
        assert!(result.is_err());
    }

    // r[verify vcf_header.serialization]
    #[test]
    fn serialization_format() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let mut builder = builder.infos();
        builder
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();
        let header = builder.build().unwrap();

        let text = header.to_vcf_text();
        assert!(text.starts_with("##fileformat=VCFv4.3\n"));
        assert!(text.contains("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">"));
        assert!(text.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">"));
        assert!(text.contains("##contig=<ID=chr1,length=1000>"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"));
    }

    // r[verify vcf_header.serialization]
    #[test]
    fn serialization_no_samples_omits_format_column() {
        let header = VcfHeader::builder().build().unwrap();
        let text = header.to_vcf_text();
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
        assert!(!text.contains("FORMAT"));
    }

    // r[verify vcf_header.string_map]
    // Typestate now prevents out-of-order registration at compile time —
    // the old `register_out_of_order_detected` runtime test is replaced
    // by this compile-time guarantee test.
    #[test]
    fn register_canonical_order_accepted() {
        use crate::vcf::record_encoder::{
            FilterFieldDef, FormatFieldDef, Gt, InfoFieldDef, Scalar,
        };

        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef { length: Some(1000) }).unwrap();

        // Typestate enforces: Contigs → Filters → Infos → Formats → Samples
        let mut builder = builder.filters();
        builder.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let mut builder = builder.infos();
        builder
            .register_info(&InfoFieldDef::<Scalar<i32>>::new(
                "DP",
                Number::Count(1),
                ValueType::Integer,
                "Depth",
            ))
            .unwrap();
        let mut builder = builder.formats();
        builder
            .register_format(&FormatFieldDef::<Gt>::new(
                "GT",
                Number::Count(1),
                ValueType::String,
                "Genotype",
            ))
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();

        let header = builder.build();
        assert!(header.is_ok(), "canonical order should succeed");
    }
}
