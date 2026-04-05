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
/// Bidirectional mapping of header string IDs to BCF dictionary indices.
#[derive(Debug, Clone)]
pub struct StringMap {
    /// ID → index (for encoding).
    to_idx: IndexMap<SmolStr, usize>,
}

impl StringMap {
    fn new() -> Self {
        Self { to_idx: IndexMap::new() }
    }

    fn insert(&mut self, id: SmolStr) -> usize {
        let next = self.to_idx.len();
        *self.to_idx.entry(id).or_insert(next)
    }

    /// Look up the BCF dictionary index for an ID.
    pub fn get(&self, id: &str) -> Option<usize> {
        self.to_idx.get_index_of(id)
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

        // FILTER first — PASS is always first in the IndexMap, giving it dict index 0.
        // This ordering ensures the BCF string dictionary matches the header text order.
        // r[impl vcf_header.serialization]
        for (id, def) in &self.filters {
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

        // other lines
        for line in &self.other_lines {
            out.push_str("##");
            out.push_str(line);
            out.push('\n');
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
/// Builder for VcfHeader with validation at build time.
#[derive(Debug, Clone)]
pub struct VcfHeaderBuilder {
    file_format: SmolStr,
    infos: IndexMap<SmolStr, InfoDef>,
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
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

    // r[impl vcf_header.builder]
    // r[impl vcf_header.string_map]
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
        let mut string_map = StringMap::new();
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
    /// Build a VcfHeaderBuilder from a BAM header, copying contig names and lengths.
    pub fn from_bam_header(header: &BamHeader) -> Result<Self, VcfHeaderError> {
        let mut builder = Self::new();
        for (tid, name) in header.target_names().enumerate() {
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
