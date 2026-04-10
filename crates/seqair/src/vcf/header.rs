//! VCF header model with builder, field definitions, and BCF string map.

use super::error::VcfHeaderError;
use crate::bam::header::BamHeader;
use indexmap::IndexMap;
use seqair_types::SmolStr;

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

// r[impl vcf_header.builder]
/// Builder for `VcfHeader` with validation at build time.
#[derive(Debug, Clone)]
pub struct VcfHeaderBuilder {
    file_format: SmolStr,
    infos: IndexMap<SmolStr, InfoDef>,
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
    /// Incremental string map, populated by `register_*` methods.
    /// `None` when only `add_*` methods are used (built fresh in `build()`).
    string_map: Option<StringMap>,
}

impl VcfHeaderBuilder {
    fn new() -> Self {
        Self {
            file_format: SmolStr::from("VCFv4.3"),
            infos: IndexMap::new(),
            formats: IndexMap::new(),
            filters: IndexMap::new(),
            contigs: IndexMap::new(),
            samples: Vec::new(),
            other_lines: Vec::new(),
            string_map: None,
        }
    }

    pub fn file_format(mut self, version: impl Into<SmolStr>) -> Self {
        self.file_format = version.into();
        self
    }

    // r[impl vcf_header.contig_required]
    pub fn add_contig(
        mut self,
        name: impl Into<SmolStr>,
        def: ContigDef,
    ) -> Result<Self, VcfHeaderError> {
        let name = name.into();
        if self.contigs.contains_key(&name) {
            return Err(VcfHeaderError::DuplicateContig { name });
        }
        self.contigs.insert(name, def);
        Ok(self)
    }

    // r[impl vcf_header.info_def]
    pub fn add_info(
        mut self,
        id: impl Into<SmolStr>,
        def: InfoDef,
    ) -> Result<Self, VcfHeaderError> {
        let id = id.into();
        if self.infos.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateInfo { id });
        }
        if def.typ == ValueType::Flag && def.number != Number::Count(0) {
            return Err(VcfHeaderError::FlagNumberMismatch { id });
        }
        self.infos.insert(id, def);
        Ok(self)
    }

    // r[impl vcf_header.format_def]
    pub fn add_format(
        mut self,
        id: impl Into<SmolStr>,
        def: FormatDef,
    ) -> Result<Self, VcfHeaderError> {
        let id = id.into();
        if self.formats.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFormat { id });
        }
        if def.typ == ValueType::Flag {
            return Err(VcfHeaderError::FormatFlagNotAllowed { id });
        }
        self.formats.insert(id, def);
        Ok(self)
    }

    // r[impl vcf_header.filter_def]
    pub fn add_filter(
        mut self,
        id: impl Into<SmolStr>,
        def: FilterDef,
    ) -> Result<Self, VcfHeaderError> {
        let id = id.into();
        if self.filters.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFilter { id });
        }
        self.filters.insert(id, def);
        Ok(self)
    }

    // r[impl vcf_header.sample_names]
    pub fn add_sample(mut self, name: impl Into<SmolStr>) -> Result<Self, VcfHeaderError> {
        let name = name.into();
        if self.samples.contains(&name) {
            return Err(VcfHeaderError::DuplicateSample { name });
        }
        self.samples.push(name);
        Ok(self)
    }

    pub fn add_other_line(mut self, line: impl Into<SmolStr>) -> Self {
        self.other_lines.push(line.into());
        self
    }

    // ── Registration methods (return typed keys) ───────────────────────

    // r[impl record_encoder.register]

    /// Register an INFO field: adds the header entry and returns a resolved typed key.
    pub fn register_info<V>(
        &mut self,
        def: &super::record_encoder::InfoFieldDef<V>,
    ) -> Result<super::record_encoder::InfoKey<V>, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        if self.infos.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateInfo { id });
        }
        if def.value_type == ValueType::Flag && def.number != Number::Count(0) {
            return Err(VcfHeaderError::FlagNumberMismatch { id });
        }
        self.infos.insert(
            id.clone(),
            InfoDef {
                number: def.number,
                typ: def.value_type,
                description: SmolStr::from(def.description),
            },
        );
        let dict_idx = self.ensure_string_map_entry(&id)?;
        Ok(super::record_encoder::InfoKey(
            super::record_encoder::FieldId { dict_idx, name: id },
            std::marker::PhantomData,
        ))
    }

    /// Register a FORMAT field: adds the header entry and returns a resolved typed key.
    pub fn register_format<V>(
        &mut self,
        def: &super::record_encoder::FormatFieldDef<V>,
    ) -> Result<super::record_encoder::FormatKey<V>, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        if self.formats.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFormat { id });
        }
        if def.value_type == ValueType::Flag {
            return Err(VcfHeaderError::FormatFlagNotAllowed { id });
        }
        self.formats.insert(
            id.clone(),
            FormatDef {
                number: def.number,
                typ: def.value_type,
                description: SmolStr::from(def.description),
            },
        );
        let dict_idx = self.ensure_string_map_entry(&id)?;
        Ok(super::record_encoder::FormatKey(
            super::record_encoder::FieldId { dict_idx, name: id },
            std::marker::PhantomData,
        ))
    }

    /// Register a FILTER: adds the header entry and returns a resolved [`FilterId`](super::record_encoder::FilterId).
    pub fn register_filter(
        &mut self,
        def: &super::record_encoder::FilterFieldDef,
    ) -> Result<super::record_encoder::FilterId, VcfHeaderError> {
        let id = SmolStr::from(def.name);
        if self.filters.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFilter { id });
        }
        self.filters.insert(id.clone(), FilterDef { description: SmolStr::from(def.description) });
        let dict_idx = self.ensure_string_map_entry(&id)?;
        Ok(super::record_encoder::FilterId { dict_idx, name: id })
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

    /// Ensure a string map entry exists, returning the dict index.
    /// Lazily initializes the string map with PASS at index 0.
    fn ensure_string_map_entry(&mut self, id: &SmolStr) -> Result<u32, VcfHeaderError> {
        let map = self.string_map.get_or_insert_with(|| {
            let mut m = StringMap::new();
            // PASS must always be at index 0
            m.insert(SmolStr::from("PASS"));
            m
        });
        u32::try_from(map.insert(id.clone())).map_err(|_| VcfHeaderError::TooManyFields)
    }

    // r[impl vcf_header.builder]
    // r[impl vcf_header.string_map]
    #[must_use = "build() returns the header; ignoring it discards all configuration"]
    pub fn build(mut self) -> Result<VcfHeader, VcfHeaderError> {
        // Ensure PASS filter is present and first
        if !self.filters.contains_key("PASS") {
            // Insert PASS at position 0
            let mut new_filters = IndexMap::new();
            new_filters.insert(
                SmolStr::from("PASS"),
                FilterDef { description: SmolStr::from("All filters passed") },
            );
            for (id, def) in self.filters {
                new_filters.insert(id, def);
            }
            self.filters = new_filters;
        }

        // Build BCF string map in the same order as to_vcf_text() emits header lines.
        // This ensures dict indices match what noodles/htslib compute from the header text.
        // Order: FILTER (PASS first), then INFO, then FORMAT.
        // r[impl vcf_header.string_map]
        //
        // If register_* methods were used, the string map was built incrementally
        // and already has correct entries. We still need to ensure all add_*-only
        // fields are included.
        let mut string_map = self.string_map.unwrap_or_default();
        for id in self.filters.keys() {
            string_map.insert(id.clone());
        }
        for id in self.infos.keys() {
            string_map.insert(id.clone());
        }
        for id in self.formats.keys() {
            string_map.insert(id.clone());
        }

        Ok(VcfHeader {
            file_format: self.file_format,
            infos: self.infos,
            formats: self.formats,
            filters: self.filters,
            contigs: self.contigs,
            samples: self.samples,
            other_lines: self.other_lines,
            string_map,
        })
    }

    // r[impl vcf_header.from_bam_header]
    /// Build a `VcfHeaderBuilder` from a BAM header, copying contig names and lengths.
    pub fn from_bam_header(header: &BamHeader) -> Result<Self, VcfHeaderError> {
        let mut builder = Self::new();
        for (tid, name) in header.target_names().enumerate() {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "BAM header n_ref ≤ 1M (enforced at parse), so tid fits in u32"
            )]
            let length = header.target_len(tid as u32);
            let def = ContigDef { length };
            builder = builder.add_contig(name, def)?;
        }
        Ok(builder)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_header.builder]
    // r[verify vcf_header.file_format]
    #[test]
    fn basic_header_build() {
        let header = VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(248_956_422) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Total Depth"),
                },
            )
            .unwrap()
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap()
            .add_sample("sample0")
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(header.file_format(), "VCFv4.3");
        assert!(header.contigs().contains_key("chr1"));
        assert!(header.infos().contains_key("DP"));
        assert!(header.formats().contains_key("GT"));
        assert_eq!(header.samples(), &[SmolStr::from("sample0")]);
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_contig() {
        let result = VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: None })
            .unwrap()
            .add_contig("chr1", ContigDef { length: None });
        assert!(result.is_err());
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_info() {
        let result = VcfHeader::builder()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("x"),
                },
            )
            .unwrap()
            .add_info(
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
        let result = VcfHeader::builder().add_sample("s1").unwrap().add_sample("s1");
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
        let header = VcfHeader::builder()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(header.string_map().get("PASS"), Some(0));
        assert_eq!(header.string_map().get("DP"), Some(1));
    }

    // r[verify vcf_header.info_def]
    #[test]
    fn flag_requires_number_zero() {
        let result = VcfHeader::builder().add_info(
            "DB",
            InfoDef {
                number: Number::Count(1),
                typ: ValueType::Flag,
                description: SmolStr::from("dbSNP"),
            },
        );
        assert!(result.is_err());

        // Number=0 + Flag should work
        let result = VcfHeader::builder().add_info(
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
        let result = VcfHeader::builder().add_format(
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
        let header = VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(1000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .add_sample("S1")
            .unwrap()
            .build()
            .unwrap();

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
}
