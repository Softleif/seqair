//! VCF header model with builder, field definitions, and BCF string map.

use super::error::{KeyError, KeyKind, VcfHeaderError};
use super::record_encoder::{
    ContigId, FieldId, FieldKind, FilterId, FormatKey, InfoKey, Section, check_genotype, check_kind,
};
use crate::bam::header::BamHeader;
use indexmap::IndexMap;
use seqair_types::SmolStr;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU64, Ordering};

mod parse;

// r[impl vcf_header.info_def]
/// Number of values for an INFO or FORMAT field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Number {
    /// Fixed count.
    Count(u32),
    /// One per ALT allele (A).
    AlternateBases,
    /// One per allele including REF (R).
    ReferenceAlternateBases,
    /// One per genotype combination (G).
    Genotypes,
    /// FORMAT only: one per allele in the sample's GT (P). VCF 4.4+.
    Ploidy,
    /// FORMAT only: like `A`, over the alternate alleles in `LAA` (LA). VCF 4.5.
    LocalAlternateBases,
    /// FORMAT only: like `R`, over the alleles in `LAA` (LR). VCF 4.5.
    LocalReferenceAlternateBases,
    /// FORMAT only: like `G`, over the alleles in `LAA` (LG). VCF 4.5.
    LocalGenotypes,
    /// One per possible base modification (M). **VCF 4.5**: "one value for each
    /// possible base modification for the corresponding `ChEBI` ID" — not the
    /// "4.2+ extension" this comment used to claim. Built headers declare 4.5
    /// ([`VcfHeader::WRITE_VERSION`]), so it is covered.
    ///
    /// Note that a strict reader may still refuse it: noodles, and so every
    /// tool built on it, rejects a header carrying `Number=M` regardless of the
    /// declared version. [`Number::Unknown`] states the same cardinality in a
    /// form every reader accepts.
    BaseModification,
    /// Unknown/variable count (.).
    Unknown,
}

impl Number {
    // r[impl vcf_header.file_format.single_version+4]
    /// The first VCF version that defines this cardinality: 4.4 for `P`, 4.5
    /// for `M`, `LA`, `LR` and `LG`, 4.0 for the rest.
    pub const fn min_version(self) -> VcfVersion {
        match self {
            Self::Ploidy => VcfVersion::V4_4,
            Self::BaseModification
            | Self::LocalAlternateBases
            | Self::LocalReferenceAlternateBases
            | Self::LocalGenotypes => VcfVersion::V4_5,
            Self::Count(_)
            | Self::AlternateBases
            | Self::ReferenceAlternateBases
            | Self::Genotypes
            | Self::Unknown => VcfVersion::V4_0,
        }
    }

    /// VCF text representation.
    pub fn as_str(&self) -> SmolStr {
        match self {
            Self::Count(n) => SmolStr::from(n.to_string()),
            Self::AlternateBases => SmolStr::from("A"),
            Self::ReferenceAlternateBases => SmolStr::from("R"),
            Self::Genotypes => SmolStr::from("G"),
            Self::Ploidy => SmolStr::from("P"),
            Self::LocalAlternateBases => SmolStr::from("LA"),
            Self::LocalReferenceAlternateBases => SmolStr::from("LR"),
            Self::LocalGenotypes => SmolStr::from("LG"),
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

impl std::fmt::Display for Number {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Count(n) => write!(f, "{n}"),
            _ => f.write_str(&self.as_str()),
        }
    }
}

impl std::fmt::Display for ValueType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
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

/// An attribute of a structured header line that the model does not
/// interpret — `Source`/`Version` on INFO, `md5`/`assembly` on contig — kept
/// in order so a parsed header writes it back out.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attribute {
    pub key: SmolStr,
    /// Unescaped: a quoted value is written back with `"` and `\\` escaped.
    pub value: SmolStr,
    pub quoted: bool,
}

/// INFO field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InfoDef {
    pub number: Number,
    pub typ: ValueType,
    pub description: SmolStr,
    pub attributes: Vec<Attribute>,
}

impl InfoDef {
    pub fn new(number: Number, typ: ValueType, description: impl Into<SmolStr>) -> Self {
        Self { number, typ, description: description.into(), attributes: Vec::new() }
    }
}

/// FORMAT field definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatDef {
    pub number: Number,
    pub typ: ValueType,
    pub description: SmolStr,
    pub attributes: Vec<Attribute>,
}

impl FormatDef {
    pub fn new(number: Number, typ: ValueType, description: impl Into<SmolStr>) -> Self {
        Self { number, typ, description: description.into(), attributes: Vec::new() }
    }
}

/// FILTER definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterDef {
    pub description: SmolStr,
    pub attributes: Vec<Attribute>,
}

impl FilterDef {
    pub fn new(description: impl Into<SmolStr>) -> Self {
        Self { description: description.into(), attributes: Vec::new() }
    }
}

/// Contig definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigDef {
    pub length: Option<u64>,
    pub attributes: Vec<Attribute>,
}

impl ContigDef {
    pub fn new(length: Option<u64>) -> Self {
        Self { length, attributes: Vec::new() }
    }
}

// r[impl vcf_header.dict_capacity+2]
/// Maximum number of FILTER/INFO/FORMAT entries in the BCF string dictionary,
/// so every dictionary index fits a `u16`. Exceeding it returns
/// [`VcfHeaderError::TooManyFields`].
pub const MAX_DICT_ENTRIES: usize = 65_536;

// r[impl vcf_header.parse.limits]
/// Longest header text [`VcfHeader::parse`] accepts, as for a BAM header.
pub const MAX_HEADER_TEXT: usize = 256 << 20;

// r[impl vcf_header.parse.limits]
/// Most contigs a header may declare, for the builder and the parser alike.
pub const MAX_CONTIGS: usize = 1_048_576;

// r[impl vcf_header.parse.limits]
/// Most samples a header may declare.
pub const MAX_SAMPLES: usize = 10_000_000;

// r[impl vcf_header.string_map]
// r[impl vcf_header.parse.dictionary]
/// A dictionary of header ids to BCF indices: FILTER/INFO/FORMAT ids in the
/// string dictionary ([`VcfHeader::string_map`]), contig names in the contig
/// dictionary ([`VcfHeader::contig_map`]).
///
/// A parsed header may leave holes (`IDX=` values with gaps, as htslib writes
/// after removing a field), so an index names at most one id, and some none.
/// Hashed, not scanned: a header may declare up to [`MAX_DICT_ENTRIES`] ids.
///
/// Two maps are equal when they give every id the same index: a hole past
/// the highest index, or one a removed id keeps, names nothing.
#[derive(Debug, Clone, Default)]
pub struct StringMap {
    ids: HashMap<SmolStr, u32>,
    slots: Vec<Option<SmolStr>>,
    /// Removed ids and the slot each keeps for itself (`r[vcf_header.remove+2]`).
    retired: HashMap<SmolStr, u32>,
}

impl PartialEq for StringMap {
    fn eq(&self, other: &Self) -> bool {
        self.ids == other.ids
    }
}

impl Eq for StringMap {}

impl StringMap {
    /// Look up the BCF dictionary index for an ID.
    pub fn get(&self, id: &str) -> Option<usize> {
        self.index(id).map(|idx| idx as usize)
    }

    fn index(&self, id: &str) -> Option<u32> {
        self.ids.get(id).copied()
    }

    /// The id at `idx`, or `None` past the end or in a hole.
    pub fn name(&self, idx: usize) -> Option<&str> {
        self.slots.get(idx)?.as_deref()
    }

    /// One past the highest index, holes included: the size of a table with
    /// one entry per index.
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// The ids in index order, skipping holes.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &str)> {
        self.slots.iter().enumerate().filter_map(|(idx, id)| Some((idx, id.as_deref()?)))
    }

    /// The index of `id`, appending it past the highest index (as htslib does,
    /// never filling a hole) if it is new. A removed id gets its own slot back.
    /// `None`, with the map unchanged, when the new index would reach `max`.
    pub(crate) fn push(&mut self, id: &SmolStr, max: usize) -> Option<u32> {
        if let Some(&idx) = self.ids.get(id) {
            return Some(idx);
        }
        // r[impl vcf_header.remove+2]
        if let Some(idx) = self.retired.remove(id) {
            if let Some(slot) = self.slots.get_mut(idx as usize) {
                *slot = Some(id.clone());
            }
            self.ids.insert(id.clone(), idx);
            return Some(idx);
        }
        let idx = self.slots.len();
        if idx >= max {
            return None;
        }
        let idx32 = u32::try_from(idx).ok()?;
        self.ids.insert(id.clone(), idx32);
        self.slots.push(Some(id.clone()));
        Some(idx32)
    }

    // r[impl vcf_header.remove+2]
    /// Make `id`'s slot a hole that stays `id`'s, as htslib's `bcf_hdr_remove`
    /// keeps the id in its dictionary: every other id keeps its index, and a
    /// `push` of `id` takes the slot back.
    pub(crate) fn remove(&mut self, id: &str) {
        if let Some((id, idx)) = self.ids.remove_entry(id) {
            if let Some(slot) = self.slots.get_mut(idx as usize) {
                *slot = None;
            }
            self.retired.insert(id, idx);
        }
    }

    /// Put a new `id` at `idx`, growing the map with holes up to it. On a
    /// conflict, returns the id already there; the map is unchanged.
    pub(crate) fn insert_at(&mut self, id: &SmolStr, idx: u32) -> Result<(), SmolStr> {
        let at = idx as usize;
        if let Some(Some(taken)) = self.slots.get(at) {
            return Err(taken.clone());
        }
        if at >= self.slots.len() {
            self.slots.resize(at.saturating_add(1), None);
        }
        if let Some(slot) = self.slots.get_mut(at) {
            *slot = Some(id.clone());
        }
        self.ids.insert(id.clone(), idx);
        Ok(())
    }
}

// r[impl vcf_fields.sample_idx]
/// One sample of one header: obtainable only from that header
/// ([`VcfHeader::sample`], [`VcfHeader::sample_indices`]) and carrying its
/// [`HeaderId`], so a record of another file refuses it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SampleIdx {
    pub(crate) index: u32,
    pub(crate) header: HeaderId,
}

impl SampleIdx {
    /// The sample's column, 0 for the first sample.
    pub fn index(&self) -> usize {
        self.index as usize
    }
}

fn unknown_field(name: &str, kind: KeyKind) -> KeyError {
    KeyError::UnknownField { name: SmolStr::from(name), kind }
}

// r[impl vcf_header.brand+3]
/// Identifies one header within the process. Every header built or parsed
/// gets a fresh one; a clone keeps it. Keys resolved from a header carry it,
/// so a key cannot silently read a field of another file that happens to
/// share its dictionary index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HeaderId(u64);

impl HeaderId {
    pub(crate) fn next() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self(NEXT.fetch_add(1, Ordering::Relaxed))
    }
}

// r[impl vcf_header.file_format+2]
/// A VCF version, as a header's `##fileformat=VCFv{major}.{minor}` line declares it.
///
/// Ordered, so "is this header new enough" is a comparison:
/// `header.version() >= VcfVersion::V4_4`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VcfVersion {
    pub major: u8,
    pub minor: u8,
}

impl VcfVersion {
    pub const V4_0: Self = Self::new(4, 0);
    pub const V4_1: Self = Self::new(4, 1);
    pub const V4_2: Self = Self::new(4, 2);
    pub const V4_3: Self = Self::new(4, 3);
    pub const V4_4: Self = Self::new(4, 4);
    pub const V4_5: Self = Self::new(4, 5);

    pub const fn new(major: u8, minor: u8) -> Self {
        Self { major, minor }
    }
}

impl std::fmt::Display for VcfVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "VCFv{}.{}", self.major, self.minor)
    }
}

// r[impl vcf_header.ordered_maps]
// r[impl vcf_header.smolstr]
/// VCF header with typed field definitions in insertion order.
#[derive(Debug, Clone)]
pub struct VcfHeader {
    id: HeaderId,
    version: VcfVersion,
    infos: IndexMap<SmolStr, InfoDef>,
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
    string_map: StringMap,
    contig_map: StringMap,
    /// The header `to_builder` derived this one from (`r[vcf_header.lineage]`).
    parent: Option<HeaderId>,
    /// Names removed from the builder since, by section.
    removed: Vec<(KeyKind, SmolStr)>,
}

impl VcfHeader {
    // r[impl vcf_header.file_format.single_version+4]
    /// The VCF version a built header declares. There is no setter: a
    /// cardinality is versioned (`Number=M` is 4.5 and nothing earlier), so a
    /// header that could declare an older version could promise a grammar it
    /// then violates. Emitting one version removes the question.
    pub const WRITE_VERSION: VcfVersion = VcfVersion::V4_5;

    /// Start building a new header.
    pub fn builder() -> VcfHeaderBuilder {
        VcfHeaderBuilder::new()
    }

    // r[impl vcf_header.to_builder]
    /// A builder seeded with every line of this header — version, other
    /// lines, definitions, contigs, samples — and its dictionaries as they
    /// are, under a fresh [`HeaderId`]: `h.to_builder().build()` is `h` but
    /// for the id, so keys resolved from `h` must be resolved again from the
    /// new header.
    ///
    /// The builder starts in the [`Contigs`] phase, so every phase is
    /// reachable. A registration appends to the dictionary past every index
    /// already taken (a new FILTER may land after the copied INFO ids), but an
    /// id removed and registered again gets its old index back;
    /// [`to_vcf_text`](Self::to_vcf_text) then writes `IDX=` so a reader
    /// assigns the indices the writer encodes with.
    ///
    /// ```
    /// use seqair::vcf::record_encoder::{FormatFieldDef, Scalar, Str};
    /// use seqair::vcf::{Number, VcfHeader};
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// // An input file's header, which declares PS as a String.
    /// let mut input = VcfHeader::builder().formats();
    /// input.register_format(&FormatFieldDef::<Str>::new("PS", Number::Count(1), "Old phase set"))?;
    /// let mut input = input.samples();
    /// input.add_sample("S1")?;
    /// let input = input.build()?;
    ///
    /// let mut builder = input.to_builder();
    /// assert!(builder.remove_format("PS"));
    /// builder.add_other_line("phaser_version=1.0");
    /// let mut builder = builder.formats();
    /// let ps = builder.register_format(&FormatFieldDef::<Scalar<i32>>::new("PS", "Phase set"))?;
    /// let header = builder.build()?;
    /// assert_eq!(header.samples(), ["S1"]);
    /// assert_eq!(header.version(), input.version());
    /// // As htslib: PS gets its index back.
    /// assert_eq!(header.string_map().get("PS"), input.string_map().get("PS"));
    /// # let _ = ps;
    /// # Ok(())
    /// # }
    /// ```
    pub fn to_builder(&self) -> VcfHeaderBuilder<Contigs> {
        let mut header = self.clone();
        header.id = HeaderId::next();
        // r[impl vcf_header.lineage]
        header.parent = Some(self.id);
        header.removed.clear();
        VcfHeaderBuilder { header, _phase: PhantomData }
    }

    // r[impl vcf_header.file_format+2]
    // r[impl vcf_header.lineage]
    /// The header this one was derived from with [`to_builder`](Self::to_builder).
    pub fn parent(&self) -> Option<HeaderId> {
        self.parent
    }

    // r[impl vcf_header.lineage]
    /// Whether `name` was removed from `kind`'s definitions since this
    /// header was derived, whether or not it was declared again.
    pub(crate) fn was_removed(&self, kind: KeyKind, name: &str) -> bool {
        self.removed.iter().any(|(k, n)| *k == kind && n == name)
    }

    /// The VCF version this header declares: [`Self::WRITE_VERSION`] for a
    /// built header, the file's own for a parsed one.
    pub fn version(&self) -> VcfVersion {
        self.version
    }

    /// The INFO definitions by ID, in header line order (not
    /// dictionary order: see [`Self::string_map`]).
    pub fn infos(&self) -> &IndexMap<SmolStr, InfoDef> {
        &self.infos
    }

    /// The FORMAT definitions by ID, in header line order (not
    /// dictionary order: see [`Self::string_map`]).
    pub fn formats(&self) -> &IndexMap<SmolStr, FormatDef> {
        &self.formats
    }

    /// The FILTER definitions by ID, in header line order (not
    /// dictionary order: see [`Self::string_map`]).
    pub fn filters(&self) -> &IndexMap<SmolStr, FilterDef> {
        &self.filters
    }

    // r[impl vcf_header.contigs]
    /// Every contig with its definition, in contig dictionary order: the
    /// BCF `rid` order, which [`ContigId::index`] gives and htslib numbers
    /// contigs by.
    pub fn contigs(&self) -> impl Iterator<Item = (ContigId, &ContigDef)> + '_ {
        self.contig_map.iter().filter_map(|(rid, name)| {
            let (name, def) = self.contigs.get_key_value(name)?;
            let id =
                ContigId { tid: u32::try_from(rid).ok()?, name: name.clone(), header: self.id };
            Some((id, def))
        })
    }

    /// The contig definitions by name, in header line order.
    pub(crate) fn contig_defs(&self) -> &IndexMap<SmolStr, ContigDef> {
        &self.contigs
    }

    pub fn samples(&self) -> &[SmolStr] {
        &self.samples
    }

    /// The `##key=value` lines the model does not interpret, without `##`, in order.
    pub fn other_lines(&self) -> &[SmolStr] {
        &self.other_lines
    }

    /// Look up the integer index for a contig name (for BCF encoding).
    pub fn contig_id(&self, name: &str) -> Result<usize, VcfHeaderError> {
        self.contig_map
            .get(name)
            .ok_or_else(|| VcfHeaderError::MissingContig { name: SmolStr::from(name) })
    }

    /// BCF string dictionary (FILTER, INFO and FORMAT ids).
    pub fn string_map(&self) -> &StringMap {
        &self.string_map
    }

    /// BCF contig dictionary: contig name to `rid`.
    pub(crate) fn contig_map(&self) -> &StringMap {
        &self.contig_map
    }

    // r[impl vcf_values.keys+2]
    /// Resolve an INFO key, checking marker `V` against the declaration
    /// ([`FieldKind`]). The marker is usually inferred from the binding:
    /// `let dp: InfoInt = header.info_key("DP")?;`.
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        let def = self.infos.get(name).ok_or_else(|| unknown_field(name, KeyKind::Info))?;
        check_kind::<V>(Section::Info, name, def.typ, def.number)?;
        Ok(InfoKey(self.field_id(name, def.number, KeyKind::Info)?, PhantomData))
    }

    // r[impl vcf_values.keys+2]
    /// Resolve a FORMAT key, checking marker `V` against the declaration.
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        let def = self.formats.get(name).ok_or_else(|| unknown_field(name, KeyKind::Format))?;
        check_kind::<V>(Section::Format, name, def.typ, def.number)?;
        Ok(FormatKey(self.field_id(name, def.number, KeyKind::Format)?, PhantomData))
    }

    fn field_id(&self, name: &str, number: Number, kind: KeyKind) -> Result<FieldId, KeyError> {
        let dict_idx = self.string_map.index(name).ok_or_else(|| unknown_field(name, kind))?;
        Ok(FieldId { dict_idx, name: SmolStr::from(name), number, header: self.id })
    }

    // r[impl vcf_values.keys+2]
    pub fn filter(&self, name: &str) -> Result<FilterId, KeyError> {
        let unknown = || KeyError::UnknownFilter { name: SmolStr::from(name) };
        if !self.filters.contains_key(name) {
            return Err(unknown());
        }
        let dict_idx = self.string_map.index(name).ok_or_else(unknown)?;
        Ok(FilterId { dict_idx, name: SmolStr::from(name), header: Some(self.id) })
    }

    // r[impl vcf_values.keys+2]
    pub fn contig(&self, name: &str) -> Result<ContigId, KeyError> {
        let tid = self
            .contig_map
            .index(name)
            .ok_or_else(|| KeyError::UnknownContig { name: SmolStr::from(name) })?;
        Ok(ContigId { tid, name: SmolStr::from(name), header: self.id })
    }

    // r[impl vcf_fields.sample_idx]
    /// The sample named `name`, as a handle valid for this header's records only.
    pub fn sample(&self, name: &str) -> Result<SampleIdx, KeyError> {
        // A scan, not a map: it runs once per reader, and a map would double a
        // 10-million-sample header's memory.
        let index = self
            .samples
            .iter()
            .position(|s| s == name)
            .and_then(|i| u32::try_from(i).ok())
            .ok_or_else(|| KeyError::UnknownSample { name: SmolStr::from(name) })?;
        Ok(SampleIdx { index, header: self.id })
    }

    // r[impl vcf_fields.sample_idx]
    /// Every sample, in column order.
    pub fn sample_indices(&self) -> impl ExactSizeIterator<Item = SampleIdx> + '_ {
        // At most MAX_SAMPLES (builder and parser both check), so this never clamps.
        let n = u32::try_from(self.samples.len()).unwrap_or(u32::MAX);
        (0..n).map(|index| SampleIdx { index, header: self.id })
    }

    // r[impl vcf_header.brand+3]
    /// This header's identity: distinct from every other header's in the
    /// process, kept by a clone.
    pub fn id(&self) -> HeaderId {
        self.id
    }

    /// A header with nothing but PASS, version `version` and a fresh id.
    fn empty(version: VcfVersion) -> Self {
        let mut header = Self {
            id: HeaderId::next(),
            version,
            infos: IndexMap::new(),
            formats: IndexMap::new(),
            filters: IndexMap::new(),
            contigs: IndexMap::new(),
            samples: Vec::new(),
            other_lines: Vec::new(),
            parent: None,
            removed: Vec::new(),
            string_map: StringMap::default(),
            contig_map: StringMap::default(),
        };
        // PASS is always the first filter and dict index 0.
        let pass = SmolStr::new_static("PASS");
        header.filters.insert(pass.clone(), FilterDef::new("All filters passed"));
        header.string_map.push(&pass, MAX_DICT_ENTRIES);
        header
    }

    // r[impl vcf_header.serialization]
    // r[impl vcf_header.parse.preserve]
    /// Serialize as VCF header text (all ## lines + #CHROM line).
    ///
    /// A parsed header whose dictionary a reader of this text would not
    /// reproduce (its file declared ids in another order than the emission
    /// order below, or left `IDX` gaps) gets `IDX=` on every definition line,
    /// as htslib writes BCF headers, so a re-read assigns the same indices.
    pub fn to_vcf_text(&self) -> String {
        let (field_idx, contig_idx) = self.needs_idx();
        let field = |id: &str| if field_idx { self.string_map.get(id) } else { None };
        let mut out = String::new();

        out.push_str("##fileformat=");
        out.push_str(&self.version.to_string());
        out.push('\n');

        // PASS right after fileformat, before metadata, as htslib does: it is
        // dictionary index 0.
        if let Some(def) = self.filters.get("PASS") {
            let line = DefLine { key: "FILTER", id: "PASS", kind: None, idx: field("PASS") };
            line.push(&mut out, &def.description, &def.attributes);
        }

        for line in &self.other_lines {
            out.push_str("##");
            out.push_str(line);
            out.push('\n');
        }

        for (id, def) in self.filters.iter().filter(|(id, _)| *id != "PASS") {
            let line = DefLine { key: "FILTER", id, kind: None, idx: field(id) };
            line.push(&mut out, &def.description, &def.attributes);
        }
        for (id, def) in &self.infos {
            let kind = Some((def.number, def.typ));
            let line = DefLine { key: "INFO", id, kind, idx: field(id) };
            line.push(&mut out, &def.description, &def.attributes);
        }
        for (id, def) in &self.formats {
            let kind = Some((def.number, def.typ));
            let line = DefLine { key: "FORMAT", id, kind, idx: field(id) };
            line.push(&mut out, &def.description, &def.attributes);
        }

        for (name, def) in &self.contigs {
            out.push_str("##contig=<ID=");
            push_value(&mut out, name);
            if let Some(len) = def.length {
                out.push_str(",length=");
                out.push_str(&len.to_string());
            }
            push_attributes(&mut out, &def.attributes);
            if let Some(idx) = self.contig_map.get(name).filter(|_| contig_idx) {
                out.push_str(",IDX=");
                out.push_str(&idx.to_string());
            }
            out.push_str(">\n");
        }

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

    /// Whether the string and the contig dictionary differ from what a reader
    /// assigns to `to_vcf_text`'s lines without `IDX` (in order of appearance).
    fn needs_idx(&self) -> (bool, bool) {
        let mut fields = StringMap::default();
        let ids = std::iter::once("PASS")
            .chain(self.filters.keys().map(SmolStr::as_str).filter(|id| *id != "PASS"))
            .chain(self.infos.keys().map(SmolStr::as_str))
            .chain(self.formats.keys().map(SmolStr::as_str));
        for id in ids {
            fields.push(&SmolStr::from(id), usize::MAX);
        }
        let mut contigs = StringMap::default();
        for name in self.contigs.keys() {
            contigs.push(name, usize::MAX);
        }
        (fields != self.string_map, contigs != self.contig_map)
    }
}

/// `##KEY=<ID=id[,Number=n,Type=t],Description="d"[,attributes][,IDX=i]>`.
struct DefLine<'a> {
    key: &'a str,
    id: &'a str,
    kind: Option<(Number, ValueType)>,
    idx: Option<usize>,
}

impl DefLine<'_> {
    fn push(&self, out: &mut String, description: &str, attributes: &[Attribute]) {
        out.push_str("##");
        out.push_str(self.key);
        out.push_str("=<ID=");
        push_value(out, self.id);
        if let Some((number, typ)) = self.kind {
            out.push_str(",Number=");
            out.push_str(&number.as_str());
            out.push_str(",Type=");
            out.push_str(typ.as_str());
        }
        out.push_str(",Description=");
        push_quoted(out, description);
        push_attributes(out, attributes);
        if let Some(idx) = self.idx {
            out.push_str(",IDX=");
            out.push_str(&idx.to_string());
        }
        out.push_str(">\n");
    }
}

fn push_attributes(out: &mut String, attributes: &[Attribute]) {
    for attr in attributes {
        out.push(',');
        out.push_str(&attr.key);
        out.push('=');
        if attr.quoted {
            push_quoted(out, &attr.value);
        } else {
            out.push_str(&attr.value);
        }
    }
}

/// `value` as it was read: bare, unless a reader would take a bare one
/// to end early (at `,` or `>`) or to be quoted or bracketed.
fn push_value(out: &mut String, value: &str) {
    if value.contains([',', '>']) || value.starts_with(['"', '[']) {
        push_quoted(out, value);
    } else {
        out.push_str(value);
    }
}

/// A double-quoted value with `"` and `\\` escaped, as the VCF spec requires.
fn push_quoted(out: &mut String, value: &str) {
    out.push('"');
    for c in value.chars() {
        if matches!(c, '"' | '\\') {
            out.push('\\');
        }
        out.push(c);
    }
    out.push('"');
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
#[derive(Debug)]
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
/// In order, it compiles — and this example is the guard for the rejection
/// below it, which uses the same names. Without it, renaming
/// `register_info` or `InfoFieldDef` would leave a `compile_fail` block that
/// passes for the wrong reason:
///
/// ```
/// use seqair::vcf::VcfHeader;
/// use seqair::vcf::record_encoder::{InfoFieldDef, Scalar};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let builder = VcfHeader::builder();
/// let mut builder = builder.infos();
/// let _dp: seqair::vcf::InfoInt = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new(
///     "DP", "Depth",
/// ))?;
/// # Ok(())
/// # }
/// ```
///
/// Out-of-order registration is a compile error:
///
/// ```compile_fail,E0599
/// use seqair::vcf::VcfHeader;
/// use seqair::vcf::record_encoder::{InfoFieldDef, Scalar};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut builder = VcfHeader::builder();
/// // ERROR: register_info is not available in the Contigs phase
/// builder.register_info(&InfoFieldDef::<Scalar<i32>>::new(
///     "DP", "Depth",
/// ))?;
/// # Ok(())
/// # }
/// ```
pub struct VcfHeaderBuilder<Phase = Contigs> {
    header: VcfHeader,
    _phase: PhantomData<Phase>,
}

// Manual impl to avoid a Phase: Debug bound from the derive. No Clone: the
// header's id is drawn when the builder starts (keys registered on it carry
// it), and two headers built from clones would share it.
impl<P> std::fmt::Debug for VcfHeaderBuilder<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VcfHeaderBuilder").field("header", &self.header).finish()
    }
}

// ── Methods available on every phase ─────────────────────────────────

impl<P> VcfHeaderBuilder<P> {
    /// Transition to a different phase (private — callers use named methods).
    fn into_phase<Q>(self) -> VcfHeaderBuilder<Q> {
        VcfHeaderBuilder { header: self.header, _phase: PhantomData }
    }

    /// The dictionary index `id` gets, without inserting it. The dictionary is
    /// capped at [`MAX_DICT_ENTRIES`], so every index fits a `u16`.
    fn push_dict_entry(&mut self, id: &SmolStr) -> Result<u32, VcfHeaderError> {
        self.header
            .string_map
            .push(id, MAX_DICT_ENTRIES)
            .ok_or(VcfHeaderError::TooManyFields { max: MAX_DICT_ENTRIES })
    }

    // r[impl vcf_header.set_samples]
    /// Append a `##key=value` metadata line, given without the `##` (e.g.
    /// `source=myapp`).
    ///
    /// Available from any phase — metadata does not affect BCF string
    /// dictionary ordering.
    pub fn add_other_line(&mut self, line: impl Into<SmolStr>) {
        self.header.other_lines.push(line.into());
    }

    // r[impl vcf_header.remove+2]
    /// Remove the INFO definition `id`; returns whether there was one. Its
    /// dictionary slot becomes a hole unless a FILTER or FORMAT `id` keeps it.
    pub fn remove_info(&mut self, id: &str) -> bool {
        let removed = self.header.infos.shift_remove(id).is_some();
        self.release(id);
        self.note_removed(KeyKind::Info, id, removed);
        removed
    }

    // r[impl vcf_header.remove+2]
    /// Remove the FORMAT definition `id`; returns whether there was one. Its
    /// dictionary slot becomes a hole unless a FILTER or INFO `id` keeps it.
    pub fn remove_format(&mut self, id: &str) -> bool {
        let removed = self.header.formats.shift_remove(id).is_some();
        self.release(id);
        self.note_removed(KeyKind::Format, id, removed);
        removed
    }

    // r[impl vcf_header.remove+2]
    /// Remove the FILTER definition `id`; returns whether there was one.
    /// `PASS` is index 0 in every header and is not removed: `false`.
    pub fn remove_filter(&mut self, id: &str) -> bool {
        if id == "PASS" {
            return false;
        }
        let removed = self.header.filters.shift_remove(id).is_some();
        self.release(id);
        self.note_removed(KeyKind::Filter, id, removed);
        removed
    }

    // r[impl vcf_header.lineage]
    fn note_removed(&mut self, kind: KeyKind, id: &str, removed: bool) {
        if removed && self.header.parent.is_some() && !self.header.was_removed(kind, id) {
            self.header.removed.push((kind, SmolStr::from(id)));
        }
    }

    // r[impl vcf_header.file_format.single_version+4]
    /// Refuse a `Number` newer than the header's version, which a derived
    /// header keeps from its source.
    fn check_number_version(&self, id: &SmolStr, number: Number) -> Result<(), VcfHeaderError> {
        let (version, needs) = (self.header.version, number.min_version());
        if version < needs {
            return Err(VcfHeaderError::NumberNeedsNewerVersion {
                id: id.clone(),
                number,
                version,
                needs,
            });
        }
        Ok(())
    }

    /// Drop `id` from the dictionary if no category declares it any more.
    fn release(&mut self, id: &str) {
        let h = &self.header;
        if !(h.filters.contains_key(id) || h.infos.contains_key(id) || h.formats.contains_key(id)) {
            self.header.string_map.remove(id);
        }
    }

    // r[impl vcf_header.builder]
    // r[impl vcf_header.string_map]
    /// Finalize the builder and produce a [`VcfHeader`].
    ///
    /// Available from any phase — all fields after contigs are optional.
    #[must_use = "build() returns the header; ignoring it discards all configuration"]
    pub fn build(self) -> Result<VcfHeader, VcfHeaderError> {
        Ok(self.header)
    }
}

// ── Contigs phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Contigs> {
    fn new() -> Self {
        Self { header: VcfHeader::empty(VcfHeader::WRITE_VERSION), _phase: PhantomData }
    }

    // r[impl vcf_header.contig_required]
    pub fn add_contig(
        &mut self,
        name: impl Into<SmolStr>,
        def: ContigDef,
    ) -> Result<(), VcfHeaderError> {
        self.register_contig(name, def)?;
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
        if self.header.contigs.contains_key(&name) {
            return Err(VcfHeaderError::DuplicateContig { name });
        }
        let tid = self
            .header
            .contig_map
            .push(&name, MAX_CONTIGS)
            .ok_or(VcfHeaderError::TooManyContigs { max: MAX_CONTIGS })?;
        self.header.contigs.insert(name.clone(), def);
        Ok(super::record_encoder::ContigId { tid, name, header: self.header.id })
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
        for target in header.targets() {
            let def = ContigDef::new(Some(target.length));
            contigs.push(builder.register_contig(target.name, def)?);
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
    /// Insert a FILTER definition, returning its BCF dict index. A refused
    /// definition leaves the builder unchanged.
    fn insert_filter(&mut self, id: SmolStr, def: FilterDef) -> Result<u32, VcfHeaderError> {
        if self.header.filters.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFilter { id });
        }
        let dict_idx = self.push_dict_entry(&id)?;
        self.header.filters.insert(id, def);
        Ok(dict_idx)
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
        let dict_idx = self.insert_filter(id.clone(), FilterDef::new(def.description))?;
        Ok(super::record_encoder::FilterId { dict_idx, name: id, header: Some(self.header.id) })
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
    /// Insert an INFO definition, returning its BCF dict index. A refused
    /// definition leaves the builder unchanged.
    fn insert_info(&mut self, id: SmolStr, def: InfoDef) -> Result<u32, VcfHeaderError> {
        if self.header.infos.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateInfo { id });
        }
        if def.typ == ValueType::Flag && def.number != Number::Count(0) {
            return Err(VcfHeaderError::FlagNumberMismatch { id });
        }
        self.check_number_version(&id, def.number)?;
        let dict_idx = self.push_dict_entry(&id)?;
        self.header.infos.insert(id, def);
        Ok(dict_idx)
    }

    pub fn add_info(&mut self, id: impl Into<SmolStr>, def: InfoDef) -> Result<(), VcfHeaderError> {
        self.insert_info(id.into(), def)?;
        Ok(())
    }

    // r[impl record_encoder.register]
    /// Register an INFO field: adds the header entry and returns a resolved typed key.
    pub fn register_info<V: FieldKind>(
        &mut self,
        def: &super::record_encoder::InfoFieldDef<V>,
    ) -> Result<super::record_encoder::InfoKey<V>, VcfHeaderError> {
        check_genotype::<V>(Section::Info, def.name)?;
        let id = SmolStr::from(def.name);
        let dict_idx =
            self.insert_info(id.clone(), InfoDef::new(def.number, V::VALUE_TYPE, def.description))?;
        let field = super::record_encoder::FieldId {
            dict_idx,
            name: id,
            number: def.number,
            header: self.header.id,
        };
        Ok(super::record_encoder::InfoKey(field, PhantomData))
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
    /// Insert a FORMAT definition, returning its BCF dict index. A refused
    /// definition leaves the builder unchanged.
    fn insert_format(&mut self, id: SmolStr, def: FormatDef) -> Result<u32, VcfHeaderError> {
        if self.header.formats.contains_key(&id) {
            return Err(VcfHeaderError::DuplicateFormat { id });
        }
        if def.typ == ValueType::Flag {
            return Err(VcfHeaderError::FormatFlagNotAllowed { id });
        }
        self.check_number_version(&id, def.number)?;
        let dict_idx = self.push_dict_entry(&id)?;
        self.header.formats.insert(id, def);
        Ok(dict_idx)
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
    pub fn register_format<V: FieldKind>(
        &mut self,
        def: &super::record_encoder::FormatFieldDef<V>,
    ) -> Result<super::record_encoder::FormatKey<V>, VcfHeaderError> {
        check_genotype::<V>(Section::Format, def.name)?;
        let id = SmolStr::from(def.name);
        let dict_idx = self.insert_format(
            id.clone(),
            FormatDef::new(def.number, V::VALUE_TYPE, def.description),
        )?;
        let field = super::record_encoder::FieldId {
            dict_idx,
            name: id,
            number: def.number,
            header: self.header.id,
        };
        Ok(super::record_encoder::FormatKey(field, PhantomData))
    }

    /// Advance to the [`Samples`] phase.
    pub fn samples(self) -> VcfHeaderBuilder<Samples> {
        self.into_phase()
    }
}

// ── Samples phase ────────────────────────────────────────────────────

impl VcfHeaderBuilder<Samples> {
    // r[impl vcf_header.set_samples]
    /// Replace the samples with `names`, in order. On an error (a name given
    /// twice, more than [`MAX_SAMPLES`]) the samples are left as they were.
    pub fn set_samples<S: Into<SmolStr>>(
        &mut self,
        names: impl IntoIterator<Item = S>,
    ) -> Result<(), VcfHeaderError> {
        let mut samples: Vec<SmolStr> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for name in names {
            let name = name.into();
            if samples.len() >= MAX_SAMPLES {
                return Err(VcfHeaderError::TooManySamples { max: MAX_SAMPLES });
            }
            if !seen.insert(name.clone()) {
                return Err(VcfHeaderError::DuplicateSample { name });
            }
            samples.push(name);
        }
        self.header.samples = samples;
        Ok(())
    }

    // r[impl vcf_header.sample_names]
    pub fn add_sample(&mut self, name: impl Into<SmolStr>) -> Result<(), VcfHeaderError> {
        let name = name.into();
        if self.header.samples.contains(&name) {
            return Err(VcfHeaderError::DuplicateSample { name });
        }
        if self.header.samples.len() >= MAX_SAMPLES {
            return Err(VcfHeaderError::TooManySamples { max: MAX_SAMPLES });
        }
        self.header.samples.push(name);
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
    // r[verify vcf_header.file_format+2]
    #[test]
    fn basic_header_build() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef::new(Some(248_956_422))).unwrap();
        let mut builder = builder.infos();
        builder
            .add_info("DP", InfoDef::new(Number::Count(1), ValueType::Integer, "Total Depth"))
            .unwrap();
        let mut builder = builder.formats();
        builder
            .add_format("GT", FormatDef::new(Number::Count(1), ValueType::String, "Genotype"))
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("sample0").unwrap();
        let header = builder.build().unwrap();

        assert_eq!(header.version(), VcfVersion::V4_5);
        assert!(header.contig_defs().contains_key("chr1"));
        assert!(header.infos().contains_key("DP"));
        assert!(header.formats().contains_key("GT"));
        assert_eq!(header.samples(), &[SmolStr::from("sample0")]);
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_contig() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef::new(None)).unwrap();
        let result = builder.add_contig("chr1", ContigDef::new(None));
        assert!(result.is_err());
    }

    // r[verify vcf_header.no_duplicates]
    #[test]
    fn rejects_duplicate_info() {
        let mut builder = VcfHeader::builder().infos();
        builder.add_info("DP", InfoDef::new(Number::Count(1), ValueType::Integer, "x")).unwrap();
        let result =
            builder.add_info("DP", InfoDef::new(Number::Count(1), ValueType::Integer, "y"));
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

    /// An INFO builder whose dictionary is full: PASS plus `MAX_DICT_ENTRIES - 1` INFO ids.
    fn full_dictionary() -> VcfHeaderBuilder<Infos> {
        let mut builder = VcfHeader::builder().infos();
        for i in 1..MAX_DICT_ENTRIES {
            let def = InfoDef::new(Number::Count(1), ValueType::Integer, "x");
            builder.add_info(format!("I{i}"), def).unwrap();
        }
        builder
    }

    // r[verify vcf_header.dict_capacity+2]
    #[test]
    fn dict_capacity_rejects_an_entry_past_the_limit() {
        let mut builder = full_dictionary();
        let def = InfoDef::new(Number::Count(1), ValueType::Integer, "x");
        let err = builder.add_info("overflow", def).unwrap_err();
        assert!(matches!(err, VcfHeaderError::TooManyFields { max: 65_536 }), "{err:?}");
    }

    // r[verify vcf_header.dict_capacity+2]
    #[test]
    fn dict_capacity_allows_exactly_the_limit() {
        let header = full_dictionary().build().unwrap();
        let last = header.string_map().get(&format!("I{}", MAX_DICT_ENTRIES - 1));
        assert_eq!(last, Some(usize::from(u16::MAX)), "the highest index fits a u16");
    }

    // r[verify vcf_header.dict_capacity+2]
    /// A registration the full dictionary refuses leaves the builder as it was:
    /// no header entry without a dictionary index.
    #[test]
    fn a_refused_registration_leaves_the_builder_unchanged() {
        let mut builder = full_dictionary();
        let before = (builder.header.infos.len(), builder.header.string_map.clone());
        let def = InfoDef::new(Number::Count(1), ValueType::Integer, "x");
        builder.add_info("overflow", def).unwrap_err();
        assert_eq!((builder.header.infos.len(), builder.header.string_map.clone()), before);
        assert!(!builder.header.infos.contains_key("overflow"));

        // An id already in the dictionary still fits: FORMAT shares INFO's index.
        let mut builder = builder.formats();
        let def = FormatDef::new(Number::Count(1), ValueType::Integer, "x");
        builder.add_format("I1", def).unwrap();
        let def = FormatDef::new(Number::Count(1), ValueType::Integer, "x");
        builder.add_format("overflow", def).unwrap_err();
        assert!(!builder.header.formats.contains_key("overflow"));
    }

    // r[verify vcf_header.brand+3]
    #[test]
    fn every_header_has_its_own_id_and_a_clone_keeps_it() {
        use crate::vcf::record_encoder::{InfoFieldDef, Scalar};
        let ids: Vec<HeaderId> =
            (0..64).map(|_| VcfHeader::builder().build().unwrap().id()).collect();
        let distinct: std::collections::HashSet<_> = ids.iter().collect();
        assert_eq!(distinct.len(), ids.len());

        let mut builder = VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef::new(None)).unwrap();
        let mut builder = builder.infos();
        let dp = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "")).unwrap();
        let header = builder.build().unwrap();
        assert_eq!(header.clone().id(), header.id());
        assert_eq!(contig.header, header.id(), "a registered contig carries the header's id");
        assert_eq!(dp.0.header, header.id(), "a registered key carries the header's id");
    }

    // r[verify vcf_header.brand+3]
    // r[verify vcf_values.header_brand+2]
    #[test]
    fn the_writer_refuses_a_contig_of_another_header() {
        use crate::vcf::{Alleles, OutputFormat, VcfError, Writer};
        use seqair_types::{Base, Pos1};
        let make = || {
            let mut builder = VcfHeader::builder();
            let contig = builder.register_contig("chr1", ContigDef::new(None)).unwrap();
            (builder.build().unwrap(), contig)
        };
        let (header, _) = make();
        let (_, foreign) = make();
        for format in [OutputFormat::Vcf, OutputFormat::Bcf] {
            let mut writer = Writer::new(Vec::new(), format).write_header(&header).unwrap();
            let alleles = Alleles::snv(Base::A, Base::T).unwrap();
            let err = writer.begin_record(&foreign, Pos1::new(1).unwrap(), &alleles, None).err();
            assert!(
                matches!(err, Some(VcfError::ForeignKey { ref name }) if name == "chr1"),
                "{err:?}"
            );
        }
    }

    // r[verify vcf_header.string_map]
    #[test]
    fn string_map_pass_is_zero() {
        let mut builder = VcfHeader::builder().infos();
        builder
            .add_info("DP", InfoDef::new(Number::Count(1), ValueType::Integer, "Depth"))
            .unwrap();
        let header = builder.build().unwrap();

        assert_eq!(header.string_map().get("PASS"), Some(0));
        assert_eq!(header.string_map().get("DP"), Some(1));
    }

    // r[verify vcf_header.info_def]
    #[test]
    fn flag_requires_number_zero() {
        let mut builder = VcfHeader::builder().infos();
        let result =
            builder.add_info("DB", InfoDef::new(Number::Count(1), ValueType::Flag, "dbSNP"));
        assert!(result.is_err());

        // Number=0 + Flag should work
        let result =
            builder.add_info("DB", InfoDef::new(Number::Count(0), ValueType::Flag, "dbSNP"));
        assert!(result.is_ok());
    }

    // r[verify vcf_header.format_def]
    #[test]
    fn format_rejects_flag() {
        let mut builder = VcfHeader::builder().formats();
        let result =
            builder.add_format("X", FormatDef::new(Number::Count(0), ValueType::Flag, "bad"));
        assert!(result.is_err());
    }

    // r[verify vcf_header.serialization]
    #[test]
    fn serialization_format() {
        let mut builder = VcfHeader::builder();
        builder.add_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut builder = builder.infos();
        builder
            .add_info("DP", InfoDef::new(Number::Count(1), ValueType::Integer, "Depth"))
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();
        let header = builder.build().unwrap();

        let text = header.to_vcf_text();
        assert!(text.starts_with("##fileformat=VCFv4.5\n"));
        assert!(text.contains("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">"));
        assert!(text.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">"));
        assert!(text.contains("##contig=<ID=chr1,length=1000>"));
        assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"));
    }

    // r[verify vcf_header.serialization]
    /// A quote or backslash in a description is escaped, and attributes the
    /// model does not interpret follow the description in order.
    #[test]
    fn serialization_escapes_descriptions_and_keeps_attributes() {
        let mut builder = VcfHeader::builder();
        let mut contig = ContigDef::new(Some(10));
        contig.attributes.push(Attribute { key: "md5".into(), value: "abc".into(), quoted: false });
        builder.add_contig("chr1", contig).unwrap();
        let mut builder = builder.infos();
        let mut def = InfoDef::new(Number::Count(1), ValueType::Integer, r#"a "b" \c"#);
        def.attributes.push(Attribute { key: "Source".into(), value: "x,y".into(), quoted: true });
        builder.add_info("DP", def).unwrap();
        let text = builder.build().unwrap().to_vcf_text();
        assert!(
            text.contains(
                r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="a \"b\" \\c",Source="x,y">"#
            ),
            "{text}"
        );
        assert!(text.contains("##contig=<ID=chr1,length=10,md5=abc>\n"), "{text}");
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
        builder.add_contig("chr1", ContigDef::new(Some(1000))).unwrap();

        // Typestate enforces: Contigs → Filters → Infos → Formats → Samples
        let mut builder = builder.filters();
        builder.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let mut builder = builder.infos();
        builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
        let mut builder = builder.formats();
        builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();

        let header = builder.build();
        assert!(header.is_ok(), "canonical order should succeed");
    }
}
