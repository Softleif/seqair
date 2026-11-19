//! Format-agnostic record encoder with typed field keys.
//!
//! Typed field keys ([`InfoKey`], [`FormatKey`]) enforce value types at compile
//! time. The [`InfoEncoder`] and [`FormatEncoder`] traits provide object-safe
//! APIs for encoding INFO and FORMAT fields into either BCF binary or VCF text.
//!
//! # Field definitions
//!
//! [`InfoFieldDef`] and [`FormatFieldDef`] combine header metadata with
//! type-safe key resolution. The `Type` a definition declares comes from its
//! marker ([`FieldKind::VALUE_TYPE`]), so it cannot disagree with the values
//! the key encodes. Register them on [`VcfHeaderBuilder`](super::header::VcfHeaderBuilder)
//! to build the header and resolve keys in one step:
//!
//! ```rust
//! use seqair::vcf::{Number, VcfHeader, InfoInt, FormatGt};
//! use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut builder = VcfHeader::builder();
//! // advance past contigs/filters to the infos phase
//! let mut builder = builder.infos();
//!
//! // Returns InfoKey<Scalar<i32>> (aliased InfoInt) — can only encode i32 values
//! let dp_key: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")
//! )?;
//!
//! // advance to the formats phase
//! let mut builder = builder.formats();
//!
//! // Returns FormatKey<Gt> (aliased FormatGt) — can only encode Genotype slices
//! let gt_key: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", "Genotype")
//! )?;
//!
//! assert_eq!(dp_key.id().name(), "DP");
//! assert_eq!(gt_key.id().name(), "GT");
//! # Ok(())
//! # }
//! ```

use super::error::{KeyError, VcfError};
use super::header::{DictId, HeaderId, Number, ValueType};
use super::record::Genotype;
use seqair_types::SmolStr;
use std::marker::PhantomData;

// ── Value type markers (uninhabited) ───────────────────────────────────

// r[impl record_encoder.typed_keys]

/// Marker for a single scalar value: a field declared `Number=1`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scalar<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for a variable-length array. Reading one, an element missing
/// inside it (`1,.,3`) is an error ([`MissingElement`](super::MissingElement));
/// read through [`OptArr`] where elements may be missing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arr<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for an array with optional (missing) elements. It reads a field
/// of any `Number`, so it is also the lenient key: a `GQ` declared
/// `Number=.` where you expect one value reads with `.get(0)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptArr<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for a flag field (no value).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

/// Marker for a string field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Str {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

/// Marker for genotype (GT) format field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Gt {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

// ── Field kinds ────────────────────────────────────────────────────────

/// The element types a numeric marker ([`Scalar`], [`Arr`], [`OptArr`]) holds.
pub trait Numeric {
    /// `Integer` for `i32`, `Float` for `f32`.
    const VALUE_TYPE: ValueType;
}

impl Numeric for i32 {
    const VALUE_TYPE: ValueType = ValueType::Integer;
}

impl Numeric for f32 {
    const VALUE_TYPE: ValueType = ValueType::Float;
}

// r[impl vcf_values.field_kind]
/// What a marker declares and what declarations it can read: the single
/// source of a marker's `Type` and accepted `Number`s, for the writer's
/// [`InfoFieldDef`]/[`FormatFieldDef`] and for resolving keys from a parsed
/// header ([`VcfHeader::info_key`](super::VcfHeader::info_key)).
pub trait FieldKind {
    /// The `Type` a definition with this marker declares.
    const VALUE_TYPE: ValueType;

    /// True for [`Gt`]: the marker names only the FORMAT `GT` field.
    const GENOTYPE: bool = false;

    /// Whether a field declared with this `Type` can be read through the marker.
    fn reads_type(declared: ValueType) -> bool {
        declared == Self::VALUE_TYPE
    }

    /// Whether a field declared with this `Number` fits the marker's values.
    fn holds_number(_declared: Number) -> bool {
        true
    }
}

/// `Number=1`.
impl<T: Numeric> FieldKind for Scalar<T> {
    const VALUE_TYPE: ValueType = T::VALUE_TYPE;

    fn holds_number(declared: Number) -> bool {
        declared == Number::Count(1)
    }
}

impl<T: Numeric> FieldKind for Arr<T> {
    const VALUE_TYPE: ValueType = T::VALUE_TYPE;
}

impl<T: Numeric> FieldKind for OptArr<T> {
    const VALUE_TYPE: ValueType = T::VALUE_TYPE;
}

impl FieldKind for Flag {
    const VALUE_TYPE: ValueType = ValueType::Flag;
}

/// Writes `String`; reads `String` and `Character` (both are text on disk).
impl FieldKind for Str {
    const VALUE_TYPE: ValueType = ValueType::String;

    fn reads_type(declared: ValueType) -> bool {
        matches!(declared, ValueType::String | ValueType::Character)
    }
}

impl FieldKind for Gt {
    const VALUE_TYPE: ValueType = ValueType::String;
    const GENOTYPE: bool = true;
}

/// Which section a key names a field of.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Section {
    Info,
    Format,
}

/// Check that marker `V` may name the field `name` of `section` at all: [`Gt`]
/// names only the FORMAT `GT` field.
pub(crate) fn check_genotype<V: FieldKind>(section: Section, name: &str) -> Result<(), KeyError> {
    if V::GENOTYPE && (section != Section::Format || name != "GT") {
        return Err(KeyError::NotGenotype { name: SmolStr::from(name) });
    }
    Ok(())
}

// r[impl vcf_values.keys+2]
/// Check that marker `V` can name a field declared as `declared`/`number`.
pub(crate) fn check_kind<V: FieldKind>(
    section: Section,
    name: &str,
    declared: ValueType,
    number: Number,
) -> Result<(), KeyError> {
    check_genotype::<V>(section, name)?;
    if !V::reads_type(declared) {
        return Err(KeyError::TypeMismatch {
            name: SmolStr::from(name),
            declared,
            requested: V::VALUE_TYPE,
        });
    }
    if !V::holds_number(number) {
        return Err(KeyError::NumberMismatch { name: SmolStr::from(name), declared: number });
    }
    Ok(())
}

// ── Field identifiers ──────────────────────────────────────────────────

// r[impl record_encoder.field_id+2]
/// An INFO or FORMAT field as resolved from a header, once: everything a
/// writer encodes it by and a reader looks it up and decodes it by, so
/// neither consults the header per record. The marker `V` of the
/// [`InfoKey<V>`]/[`FormatKey<V>`] around it was checked against `typ` and
/// `number` when it was resolved.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FieldId {
    /// The BCF dictionary id: the BCF encoder's key, and the reader's slot
    /// in a record's index.
    pub(crate) id: DictId,
    /// The VCF text encoder's key, and what errors name.
    pub(crate) name: SmolStr,
    /// The declared `Type` (`String` for a field no header declares): what
    /// an [`Any`](super::reader::Any) read decodes the value as.
    pub(crate) typ: ValueType,
    /// The declared `Number`, which the writer checks a value's count against
    /// (`r[record_encoder.count_check+3]`).
    pub(crate) number: Number,
    /// The FORMAT `GT` field: an `Any` read decodes it as a genotype.
    pub(crate) genotype: bool,
    /// The header the id was resolved from (`r[vcf_values.header_brand]`).
    pub(crate) header: HeaderId,
}

impl FieldId {
    /// The field `name` of `section`, declared `typ`/`number`, with
    /// dictionary id `id` in header `header`.
    pub(crate) fn new(
        id: DictId,
        name: SmolStr,
        typ: ValueType,
        number: Number,
        section: Section,
        header: HeaderId,
    ) -> Self {
        let genotype = section == Section::Format && name == "GT";
        Self { id, name, typ, number, genotype, header }
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn dict_idx(&self) -> u32 {
        self.id.get()
    }
    pub fn number(&self) -> Number {
        self.number
    }

    // r[impl record_encoder.count_check+3]
    /// Check that `found` values fit the declared `Number` in a record of
    /// `n_allele` alleles. `ploidy` is the record's GT ploidy once GT is
    /// encoded; without it `G` is not checked.
    pub(crate) fn check_count(
        &self,
        n_allele: usize,
        ploidy: Option<usize>,
        found: usize,
    ) -> Result<(), VcfError> {
        let expected = match self.number {
            Number::Count(n) => usize::try_from(n).ok(),
            Number::AlternateBases => Some(n_allele.saturating_sub(1)),
            Number::ReferenceAlternateBases => Some(n_allele),
            Number::Genotypes => ploidy.and_then(|p| genotype_count(n_allele, p)),
            Number::Ploidy
            | Number::LocalAlternateBases
            | Number::LocalReferenceAlternateBases
            | Number::LocalGenotypes
            | Number::BaseModification
            | Number::Unknown => None,
        };
        match expected {
            Some(expected) if expected != found => Err(VcfError::CountMismatch {
                field: self.name.clone(),
                number: self.number,
                expected,
                found,
            }),
            _ => Ok(()),
        }
    }
}

/// The number of unordered genotypes of `ploidy` alleles drawn from
/// `n_allele`: `binomial(n_allele + ploidy - 1, ploidy)`. `None` if it
/// overflows, which no real record reaches.
fn genotype_count(n_allele: usize, ploidy: usize) -> Option<usize> {
    let mut count: usize = 1;
    for k in 1..=ploidy {
        // Exact at every step: the running product is binomial(n_allele + k - 1, k).
        count = count.checked_mul(n_allele.checked_add(k)?.checked_sub(1)?)?.checked_div(k)?;
    }
    Some(count)
}

// r[impl record_encoder.contig_id]
/// Resolved contig identifier carrying both integer tid and string name.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ContigId {
    pub(crate) tid: u32,
    pub(crate) name: SmolStr,
    /// The header the tid was resolved from (`r[vcf_values.header_brand]`).
    pub(crate) header: HeaderId,
}

impl ContigId {
    // r[impl vcf_header.contigs]
    pub fn name(&self) -> &SmolStr {
        &self.name
    }
    pub fn tid(&self) -> u32 {
        self.tid
    }
    /// The contig's place in its header's contig dictionary, from 0: the
    /// BCF `rid`.
    pub fn index(&self) -> usize {
        self.tid as usize
    }
}

// r[impl record_encoder.filter_id]
/// Resolved filter identifier carrying both dictionary index and name.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FilterId {
    pub(crate) dict_idx: u32,
    pub(crate) name: SmolStr,
    /// The header the index was resolved from; `None` for [`Self::PASS`],
    /// which is index 0 in every header (`r[vcf_header.filter_def]`).
    pub(crate) header: Option<HeaderId>,
}

impl FilterId {
    pub const PASS: FilterId =
        Self { dict_idx: 0, name: SmolStr::new_inline("PASS"), header: None };
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn dict_idx(&self) -> u32 {
        self.dict_idx
    }
}

// ── Typed field keys ───────────────────────────────────────────────────

// r[impl record_encoder.key_types+2]

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct InfoKey<V>(pub(crate) FieldId, pub(crate) PhantomData<V>);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FormatKey<V>(pub(crate) FieldId, pub(crate) PhantomData<V>);

impl<V> InfoKey<V> {
    pub fn id(&self) -> &FieldId {
        &self.0
    }
}

impl<V> FormatKey<V> {
    pub fn id(&self) -> &FieldId {
        &self.0
    }
}

pub type InfoInt = InfoKey<Scalar<i32>>;
pub type InfoFloat = InfoKey<Scalar<f32>>;
pub type InfoInts = InfoKey<Arr<i32>>;
pub type InfoFloats = InfoKey<Arr<f32>>;
pub type InfoFlag = InfoKey<Flag>;
pub type InfoString = InfoKey<Str>;
pub type InfoIntOpts = InfoKey<OptArr<i32>>;
pub type InfoFloatOpts = InfoKey<OptArr<f32>>;

pub type FormatGt = FormatKey<Gt>;
pub type FormatInt = FormatKey<Scalar<i32>>;
pub type FormatFloat = FormatKey<Scalar<f32>>;
pub type FormatInts = FormatKey<Arr<i32>>;
/// Fails on an element missing inside the array; see [`FormatFloatOpts`].
pub type FormatFloats = FormatKey<Arr<f32>>;
pub type FormatFloatOpts = FormatKey<OptArr<f32>>;
pub type FormatIntOpts = FormatKey<OptArr<i32>>;
pub type FormatString = FormatKey<Str>;

// ── InfoEncoder trait ──────────────────────────────────────────────────

// r[impl record_encoder.info_methods+3]
/// Object-safe trait for encoding INFO fields. Every method fails with
/// [`VcfError::ForeignKey`] for a key of another header
/// (`r[record_encoder.foreign_keys]`) and, for arrays, with
/// [`VcfError::CountMismatch`] for a count its `Number` does not allow
/// (`r[record_encoder.count_check+3]`), before writing anything.
pub trait InfoEncoder {
    // r[impl record_encoder.reserved_rejection+2]
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedIntValue`] if `value` is a reserved BCF pattern.
    fn info_int(&mut self, id: &FieldId, value: i32) -> Result<(), VcfError>;
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if `value` is a reserved BCF pattern.
    fn info_float(&mut self, id: &FieldId, value: f32) -> Result<(), VcfError>;
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedIntValue`] if any value is a reserved BCF pattern.
    fn info_ints(&mut self, id: &FieldId, values: &[i32]) -> Result<(), VcfError>;
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if any value is a reserved BCF pattern.
    fn info_floats(&mut self, id: &FieldId, values: &[f32]) -> Result<(), VcfError>;
    /// # Errors
    /// [`VcfError::ForeignKey`] if `id` is from another header.
    fn info_flag(&mut self, id: &FieldId) -> Result<(), VcfError>;
    /// # Errors
    /// [`VcfError::ForeignKey`] if `id` is from another header.
    fn info_string(&mut self, id: &FieldId, value: &str) -> Result<(), VcfError>;
    /// # Errors
    /// [`VcfError::ForeignKey`] if `id` is from another header,
    /// `CountMismatch` if the declared `Number` disagrees with the count.
    fn info_strings(&mut self, id: &FieldId, values: &[&str]) -> Result<(), VcfError>;
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedIntValue`] if any concrete value is reserved.
    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]) -> Result<(), VcfError>;
    // r[impl vcf_passthrough.opt_arrays]
    /// # Errors
    /// [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if any concrete value is reserved.
    fn info_float_opts(&mut self, id: &FieldId, values: &[Option<f32>]) -> Result<(), VcfError>;
    fn n_allele(&self) -> usize;
    fn n_alt(&self) -> usize;
}

// ── FormatEncoder trait ────────────────────────────────────────────────

// r[impl record_encoder.format_methods+3]
/// Object-safe trait for encoding FORMAT fields. Every method fails with
/// [`VcfError::ForeignKey`] for a key of another header and, for arrays, with
/// [`VcfError::CountMismatch`] for a sample whose count its `Number` does not
/// allow, before writing anything.
pub trait FormatEncoder {
    /// Encode a GT FORMAT field for all samples. Slice length MUST equal sample count.
    fn format_gt(&mut self, id: &FieldId, gts: &[Genotype]) -> Result<(), VcfError>;
    /// Encode a scalar integer FORMAT field for all samples.
    fn format_int(&mut self, id: &FieldId, values: &[i32]) -> Result<(), VcfError>;
    /// Encode a scalar float FORMAT field for all samples.
    fn format_float(&mut self, id: &FieldId, values: &[f32]) -> Result<(), VcfError>;
    /// Like [`format_int`](Self::format_int), but a sample may be missing
    /// (`None`): `.` in VCF text, the missing value of the column's integer
    /// type in BCF.
    fn format_int_opt(&mut self, id: &FieldId, values: &[Option<i32>]) -> Result<(), VcfError>;
    /// Like [`format_float`](Self::format_float), but a sample may be missing
    /// (`None`): `.` in VCF text, the float missing value in BCF.
    fn format_float_opt(&mut self, id: &FieldId, values: &[Option<f32>]) -> Result<(), VcfError>;
    /// Encode a multi-valued integer FORMAT field. `per_sample` holds one slice
    /// per sample (its length MUST equal the sample count); each inner slice is
    /// that sample's values. For BCF a single integer type is chosen across all
    /// samples' values (smallest that fits), the stored width is the longest
    /// inner slice, and shorter samples are padded with the integer
    /// end-of-vector sentinel for that type.
    fn format_ints(&mut self, id: &FieldId, per_sample: &[&[i32]]) -> Result<(), VcfError>;
    /// Encode a multi-valued float FORMAT field. `per_sample` holds one slice
    /// per sample (its length MUST equal the sample count); each inner slice is
    /// that sample's values. For BCF the stored width is the longest inner
    /// slice, shorter samples padded with the float end-of-vector sentinel.
    fn format_floats(&mut self, id: &FieldId, per_sample: &[&[f32]]) -> Result<(), VcfError>;
    /// Like [`format_floats`](Self::format_floats), but an element may be `None`:
    /// the BCF float missing value, `.` in VCF text.
    fn format_float_opts(
        &mut self,
        id: &FieldId,
        per_sample: &[&[Option<f32>]],
    ) -> Result<(), VcfError>;
    // r[impl vcf_passthrough.opt_arrays]
    /// Like [`format_ints`](Self::format_ints), but an element may be `None`:
    /// the missing value of the column's integer type in BCF, `.` in VCF text.
    fn format_int_opts(
        &mut self,
        id: &FieldId,
        per_sample: &[&[Option<i32>]],
    ) -> Result<(), VcfError>;
    /// Encode a string FORMAT field — one string per sample (its length MUST
    /// equal the sample count). For BCF the stored width is the longest string,
    /// shorter samples NUL-padded.
    fn format_string(&mut self, id: &FieldId, values: &[&str]) -> Result<(), VcfError>;
    fn n_allele(&self) -> usize;
    fn n_alt(&self) -> usize;
    /// Number of samples declared by the header.
    fn n_samples(&self) -> usize;
}

// ── Key encode methods ─────────────────────────────────────────────────

// r[impl record_encoder.key_encode+2]

impl InfoKey<Scalar<i32>> {
    /// Encode a scalar integer INFO value.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedIntValue`] if `value` is one
    /// of the BCF reserved integer patterns (see `r[bcf_encoder.reserved_bands]`).
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        value: i32,
    ) -> Result<(), VcfError> {
        enc.info_int(&self.0, value)
    }
}

impl InfoKey<Scalar<f32>> {
    /// Encode a scalar float INFO value.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if `value`'s bit
    /// pattern is one of the BCF reserved float patterns. `+inf` and quiet NaN are
    /// ordinary values and are accepted.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        value: f32,
    ) -> Result<(), VcfError> {
        enc.info_float(&self.0, value)
    }
}

impl InfoKey<Arr<i32>> {
    /// Encode an integer array INFO value.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedIntValue`] if any element is
    /// one of the BCF reserved integer patterns.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        values: &[i32],
    ) -> Result<(), VcfError> {
        enc.info_ints(&self.0, values)
    }
}

impl InfoKey<Arr<f32>> {
    /// Encode a float array INFO value.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if any element's
    /// bit pattern is one of the BCF reserved float patterns.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        values: &[f32],
    ) -> Result<(), VcfError> {
        enc.info_floats(&self.0, values)
    }
}

impl InfoKey<Flag> {
    /// Set a flag INFO field.
    ///
    /// # Errors
    /// Returns [`VcfError::ForeignKey`] if the key is from another header.
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized)) -> Result<(), VcfError> {
        enc.info_flag(&self.0)
    }
}

impl InfoKey<Str> {
    /// Encode a string INFO value.
    ///
    /// # Errors
    /// Returns [`VcfError::ForeignKey`] if the key is from another header.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        value: &str,
    ) -> Result<(), VcfError> {
        enc.info_string(&self.0, value)
    }

    // r[impl record_encoder.string_list]
    /// Encode a list of strings (`Number=.`, `A`, …): each element is
    /// percent-encoded, then they are joined with `,`, so an element may
    /// hold a comma. An empty list is `.`.
    ///
    /// # Errors
    /// [`VcfError::ForeignKey`] if the key is from another header,
    /// `CountMismatch` if the declared `Number` disagrees with the count.
    pub fn encode_list(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        values: &[&str],
    ) -> Result<(), VcfError> {
        enc.info_strings(&self.0, values)
    }
}

impl InfoKey<OptArr<i32>> {
    /// Encode an integer array INFO value with per-element missing markers.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedIntValue`] if any concrete
    /// element is one of the BCF reserved integer patterns. `None` elements are the
    /// missing marker by construction and are always accepted.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        values: &[Option<i32>],
    ) -> Result<(), VcfError> {
        enc.info_int_opts(&self.0, values)
    }
}

// r[impl vcf_passthrough.opt_arrays]
impl InfoKey<OptArr<f32>> {
    /// Encode a float array INFO value with per-element missing markers.
    ///
    /// # Errors
    /// Returns [`crate::vcf::VcfEncodeError::ReservedFloatValue`] if any concrete
    /// element is a reserved BCF float pattern; `None` elements are always accepted.
    pub fn encode(
        &self,
        enc: &mut (impl InfoEncoder + ?Sized),
        values: &[Option<f32>],
    ) -> Result<(), VcfError> {
        enc.info_float_opts(&self.0, values)
    }
}

// r[impl vcf_passthrough.opt_arrays]
impl FormatKey<OptArr<i32>> {
    /// Encode a multi-valued integer FORMAT field with per-element missing
    /// markers, one inner slice per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `per_sample.len() != n_samples`,
    /// and [`crate::vcf::VcfEncodeError::ReservedIntValue`] for a reserved concrete element.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        per_sample: &[&[Option<i32>]],
    ) -> Result<(), VcfError> {
        check_samples(enc.n_samples(), per_sample.len())?;
        enc.format_int_opts(&self.0, per_sample)
    }
}

// r[impl vcf_passthrough.opt_arrays]
impl FormatKey<OptArr<f32>> {
    /// Encode a multi-valued float FORMAT field with per-element missing
    /// markers, one inner slice per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `per_sample.len() != n_samples`,
    /// and [`crate::vcf::VcfEncodeError::ReservedFloatValue`] for a reserved concrete element.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        per_sample: &[&[Option<f32>]],
    ) -> Result<(), VcfError> {
        check_samples(enc.n_samples(), per_sample.len())?;
        enc.format_float_opts(&self.0, per_sample)
    }
}

fn check_samples(expected: usize, got: usize) -> Result<(), VcfError> {
    if expected == got { Ok(()) } else { Err(VcfError::SampleCountMismatch { expected, got }) }
}

impl FormatKey<Gt> {
    /// Encode a GT FORMAT field. `gts` must have exactly one entry per sample;
    /// samples may differ in ploidy (`r[record_encoder.mixed_ploidy]`).
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `gts.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        gts: &[Genotype],
    ) -> Result<(), VcfError> {
        if gts.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: gts.len(),
            });
        }
        enc.format_gt(&self.0, gts)
    }
}

impl FormatKey<Scalar<i32>> {
    /// Encode a scalar integer FORMAT field. `values` must have exactly one entry per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `values.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[i32],
    ) -> Result<(), VcfError> {
        if values.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: values.len(),
            });
        }
        enc.format_int(&self.0, values)
    }

    // r[impl record_encoder.format_scalar_opt]
    /// Encode a scalar integer FORMAT field where a sample may be missing:
    /// `None` is `.` in VCF text and the missing value of the column's
    /// integer width in BCF. `values` must have exactly one entry per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `values.len() != n_samples`,
    /// and [`crate::vcf::VcfEncodeError::ReservedIntValue`] for a concrete value
    /// BCF reserves.
    pub fn encode_opt(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[Option<i32>],
    ) -> Result<(), VcfError> {
        if values.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: values.len(),
            });
        }
        enc.format_int_opt(&self.0, values)
    }
}

impl FormatKey<Scalar<f32>> {
    /// Encode a scalar float FORMAT field. `values` must have exactly one entry per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `values.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[f32],
    ) -> Result<(), VcfError> {
        if values.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: values.len(),
            });
        }
        enc.format_float(&self.0, values)
    }

    // r[impl record_encoder.format_scalar_opt]
    /// Encode a scalar float FORMAT field where a sample may be missing:
    /// `None` is `.` in VCF text and the float missing value in BCF. `values`
    /// must have exactly one entry per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `values.len() != n_samples`,
    /// and [`crate::vcf::VcfEncodeError::ReservedFloatValue`] for a concrete
    /// value BCF reserves.
    pub fn encode_opt(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[Option<f32>],
    ) -> Result<(), VcfError> {
        if values.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: values.len(),
            });
        }
        enc.format_float_opt(&self.0, values)
    }
}

impl FormatKey<Arr<i32>> {
    /// Encode a multi-valued integer FORMAT field. `per_sample` must have exactly
    /// one slice per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `per_sample.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        per_sample: &[&[i32]],
    ) -> Result<(), VcfError> {
        if per_sample.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: per_sample.len(),
            });
        }
        enc.format_ints(&self.0, per_sample)
    }
}

impl FormatKey<Arr<f32>> {
    /// Encode a multi-valued float FORMAT field. `per_sample` must have exactly
    /// one slice per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `per_sample.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        per_sample: &[&[f32]],
    ) -> Result<(), VcfError> {
        if per_sample.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: per_sample.len(),
            });
        }
        enc.format_floats(&self.0, per_sample)
    }
}

// r[impl record_encoder.format_single_sample]

impl FormatKey<Arr<i32>> {
    /// Encode this multi-valued integer FORMAT field for a single-sample record,
    /// passing that one sample's values directly instead of wrapping them in a
    /// `&[&[i32]]`.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if the header declares anything
    /// other than exactly one sample.
    pub fn encode_single_sample(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[i32],
    ) -> Result<(), VcfError> {
        if enc.n_samples() != 1 {
            return Err(VcfError::SampleCountMismatch { expected: enc.n_samples(), got: 1 });
        }
        enc.format_ints(&self.0, &[values])
    }
}

impl FormatKey<Arr<f32>> {
    /// Encode this multi-valued float FORMAT field for a single-sample record,
    /// passing that one sample's values directly instead of wrapping them in a
    /// `&[&[f32]]`.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if the header declares anything
    /// other than exactly one sample.
    pub fn encode_single_sample(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[f32],
    ) -> Result<(), VcfError> {
        if enc.n_samples() != 1 {
            return Err(VcfError::SampleCountMismatch { expected: enc.n_samples(), got: 1 });
        }
        enc.format_floats(&self.0, &[values])
    }
}

impl FormatKey<OptArr<f32>> {
    /// Encode this multi-valued float FORMAT field for a single-sample record, with
    /// per-element missing markers.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if the header declares anything
    /// other than exactly one sample, and [`crate::vcf::VcfEncodeError::ReservedFloatValue`]
    /// if a concrete element is a reserved BCF float pattern. `None` elements are the
    /// missing marker by construction and are always accepted.
    pub fn encode_single_sample(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[Option<f32>],
    ) -> Result<(), VcfError> {
        if enc.n_samples() != 1 {
            return Err(VcfError::SampleCountMismatch { expected: enc.n_samples(), got: 1 });
        }
        enc.format_float_opts(&self.0, &[values])
    }
}

impl FormatKey<Str> {
    /// Encode a string FORMAT field. `values` must have exactly one string per sample.
    ///
    /// # Errors
    /// Returns [`VcfError::SampleCountMismatch`] if `values.len() != n_samples`.
    pub fn encode(
        &self,
        enc: &mut (impl FormatEncoder + ?Sized),
        values: &[&str],
    ) -> Result<(), VcfError> {
        if values.len() != enc.n_samples() {
            return Err(VcfError::SampleCountMismatch {
                expected: enc.n_samples(),
                got: values.len(),
            });
        }
        enc.format_string(&self.0, values)
    }
}

// ── Field definitions ──────────────────────────────────────────────────

// r[impl record_encoder.field_def+5]
// r[impl record_encoder.field_def_types]

/// An INFO field definition; its `Type` is the marker's ([`FieldKind::VALUE_TYPE`]),
/// and so is its `Number` where the marker fixes one: [`Scalar`] declares
/// `Number=1`, [`Flag`] `Number=0`, and their definitions are built with
/// [`new`](Self::new). [`Arr`], [`OptArr`] and [`Str`] definitions take a
/// `Number`, through [`with_number`](Self::with_number). The marker comes
/// from the binding's type:
///
/// ```
/// use seqair::vcf::Number;
/// use seqair::vcf::record_encoder::{Arr, Flag, InfoFieldDef, Scalar};
///
/// const DP: InfoFieldDef<Scalar<i32>> = InfoFieldDef::new("DP", "Depth");
/// const AD: InfoFieldDef<Arr<i32>> =
///     InfoFieldDef::with_number("AD", Number::ReferenceAlternateBases, "Allele depth");
/// const DB: InfoFieldDef<Flag> = InfoFieldDef::new("DB", "dbSNP");
/// // Without a typed binding, name the marker:
/// let _ = InfoFieldDef::<Scalar<f32>>::new("BQ", "Base quality");
/// ```
///
/// A scalar cannot be declared with another `Number`; one value per ALT is an
/// [`Arr`] with `Number=A`:
///
/// ```compile_fail,E0599
/// use seqair::vcf::Number;
/// use seqair::vcf::record_encoder::{InfoFieldDef, Scalar};
///
/// // ERROR: a Scalar's Number is 1; it has no `with_number`
/// let _ = InfoFieldDef::<Scalar<i32>>::with_number("DP", Number::Count(1), "Depth");
/// ```
///
/// Nor is there an INFO genotype:
///
/// ```compile_fail,E0599
/// use seqair::vcf::record_encoder::{Gt, InfoFieldDef};
///
/// let _ = InfoFieldDef::<Gt>::new("GT", "Genotype");
/// ```
pub struct InfoFieldDef<V> {
    pub(crate) name: &'static str,
    pub(crate) number: Number,
    pub(crate) description: &'static str,
    _marker: PhantomData<V>,
}

/// A FORMAT field definition; its `Type` is the marker's ([`FieldKind::VALUE_TYPE`]),
/// and so is its `Number` where the marker fixes one: [`Scalar`] and [`Gt`]
/// declare `Number=1`. [`Arr`], [`OptArr`] and [`Str`] definitions take a
/// `Number`.
pub struct FormatFieldDef<V> {
    pub(crate) name: &'static str,
    pub(crate) number: Number,
    pub(crate) description: &'static str,
    _marker: PhantomData<V>,
}

mod sealed {
    pub trait Sealed {}
    impl<T: super::Numeric> Sealed for super::Scalar<T> {}
    impl<T: super::Numeric> Sealed for super::Arr<T> {}
    impl<T: super::Numeric> Sealed for super::OptArr<T> {}
    impl Sealed for super::Str {}
    impl Sealed for super::Flag {}
    impl Sealed for super::Gt {}
}

/// An INFO marker whose `Number` the marker fixes: [`Scalar`] (`1`) and
/// [`Flag`] (`0`). Their definitions are built with [`InfoFieldDef::new`].
pub trait InfoFixedNumber: FieldKind + sealed::Sealed {
    /// The `Number` a definition with this marker declares.
    const NUMBER: Number;
}

/// A FORMAT marker whose `Number` the marker fixes: [`Scalar`] and [`Gt`]
/// (`1`). Their definitions are built with [`FormatFieldDef::new`].
pub trait FormatFixedNumber: FieldKind + sealed::Sealed {
    /// The `Number` a definition with this marker declares.
    const NUMBER: Number;
}

/// A marker whose definition takes a `Number`: [`Arr`], [`OptArr`] and
/// [`Str`]. Their definitions are built with `with_number`.
pub trait ChosenNumber: FieldKind + sealed::Sealed {}

impl<T: Numeric> InfoFixedNumber for Scalar<T> {
    const NUMBER: Number = Number::Count(1);
}
impl InfoFixedNumber for Flag {
    const NUMBER: Number = Number::Count(0);
}
impl<T: Numeric> FormatFixedNumber for Scalar<T> {
    const NUMBER: Number = Number::Count(1);
}
// No INFO `Gt`: the marker names only the FORMAT `GT` field.
impl FormatFixedNumber for Gt {
    const NUMBER: Number = Number::Count(1);
}
impl<T: Numeric> ChosenNumber for Arr<T> {}
impl<T: Numeric> ChosenNumber for OptArr<T> {}
// One string either way; `Number=.` is how real files declare a list,
// which `InfoKey<Str>::encode_list` writes.
impl ChosenNumber for Str {}

// r[impl record_encoder.field_def_const]
// r[impl record_encoder.field_def+5]
impl<V: InfoFixedNumber> InfoFieldDef<V> {
    /// A definition whose `Number` is the marker's: `Scalar` declares `1`,
    /// `Flag` `0`. The marker comes from the binding's type:
    /// `const DP: InfoFieldDef<Scalar<i32>> = InfoFieldDef::new("DP", "Depth");`.
    pub const fn new(name: &'static str, description: &'static str) -> Self {
        Self { name, number: V::NUMBER, description, _marker: PhantomData }
    }
}

impl<V: ChosenNumber> InfoFieldDef<V> {
    /// A definition of an array or string field, with its `Number`.
    pub const fn with_number(
        name: &'static str,
        number: Number,
        description: &'static str,
    ) -> Self {
        Self { name, number, description, _marker: PhantomData }
    }
}

impl<V: FormatFixedNumber> FormatFieldDef<V> {
    /// A definition whose `Number` is the marker's: `1` for `Scalar` and `Gt`.
    pub const fn new(name: &'static str, description: &'static str) -> Self {
        Self { name, number: V::NUMBER, description, _marker: PhantomData }
    }
}

impl<V: ChosenNumber> FormatFieldDef<V> {
    /// A definition of an array or string field, with its `Number`.
    pub const fn with_number(
        name: &'static str,
        number: Number,
        description: &'static str,
    ) -> Self {
        Self { name, number, description, _marker: PhantomData }
    }
}

pub struct FilterFieldDef {
    pub name: &'static str,
    pub description: &'static str,
}

impl FilterFieldDef {
    pub const fn new(name: &'static str, description: &'static str) -> Self {
        Self { name, description }
    }
}

// ── FieldDescription trait ─────────────────────────────────────────────

// r[impl record_encoder.field_description]
pub trait FieldDescription {
    fn name(&self) -> &str;
    fn number(&self) -> Number;
    fn value_type(&self) -> ValueType;
    fn description(&self) -> &str;
}

impl<V: FieldKind> FieldDescription for InfoFieldDef<V> {
    fn name(&self) -> &str {
        self.name
    }
    fn number(&self) -> Number {
        self.number
    }
    fn value_type(&self) -> ValueType {
        V::VALUE_TYPE
    }
    fn description(&self) -> &str {
        self.description
    }
}

impl<V: FieldKind> FieldDescription for FormatFieldDef<V> {
    fn name(&self) -> &str {
        self.name
    }
    fn number(&self) -> Number {
        self.number
    }
    fn value_type(&self) -> ValueType {
        V::VALUE_TYPE
    }
    fn description(&self) -> &str {
        self.description
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::OutputFormat;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::header::ContigDef;
    use crate::vcf::record::Genotype;
    use crate::vcf::unified::Writer;
    use seqair_types::{Base, Pos1, SmallVec};
    use std::sync::Arc;

    // r[verify record_encoder.field_def_const]
    #[test]
    fn field_defs_are_const() {
        const DP: InfoFieldDef<Scalar<i32>> =
            InfoFieldDef::<Scalar<i32>>::new("DP", "Combined depth");
        const GT: FormatFieldDef<Gt> = FormatFieldDef::<Gt>::new("GT", "Genotype");
        const LOW_DP: FilterFieldDef = FilterFieldDef::new("lowDp", "Low read depth");

        assert_eq!(DP.name, "DP");
        assert_eq!(GT.name, "GT");
        assert_eq!(LOW_DP.name, "lowDp");
    }

    // r[verify record_encoder.field_description]
    #[test]
    fn field_description_trait() {
        let dp = InfoFieldDef::<Scalar<i32>>::new("DP", "Combined depth");
        let desc: &dyn FieldDescription = &dp;
        assert_eq!(desc.name(), "DP");
        assert_eq!(desc.number(), Number::Count(1));
        assert_eq!(desc.value_type(), ValueType::Integer);
        assert_eq!(desc.description(), "Combined depth");
    }

    // r[verify record_encoder.field_def+5]
    // r[verify vcf_values.field_kind]
    /// The declared `Type` is the marker's, read through every path a
    /// definition reaches the header by: the def, the trait, the header line.
    #[test]
    fn a_definition_declares_its_markers_type() {
        assert_eq!(InfoFieldDef::<Scalar<i32>>::new("A", "").value_type(), ValueType::Integer);
        assert_eq!(InfoFieldDef::<Scalar<f32>>::new("A", "").value_type(), ValueType::Float);
        assert_eq!(
            InfoFieldDef::<Arr<i32>>::with_number("A", Number::Unknown, "").value_type(),
            ValueType::Integer
        );
        assert_eq!(
            InfoFieldDef::<OptArr<f32>>::with_number("A", Number::Unknown, "").value_type(),
            ValueType::Float
        );
        assert_eq!(InfoFieldDef::<Flag>::new("A", "").value_type(), ValueType::Flag);
        assert_eq!(
            InfoFieldDef::<Str>::with_number("A", Number::Count(1), "").value_type(),
            ValueType::String
        );
        assert_eq!(FormatFieldDef::<Gt>::new("GT", "").value_type(), ValueType::String);

        let mut builder = crate::vcf::header::VcfHeader::builder().infos();
        builder
            .register_info(&InfoFieldDef::<Arr<f32>>::with_number(
                "AF",
                Number::AlternateBases,
                "x",
            ))
            .unwrap();
        let header = builder.build().unwrap();
        assert!(header.to_vcf_text().contains("##INFO=<ID=AF,Number=A,Type=Float,"));
    }

    // r[verify record_encoder.field_def+5]
    /// A marker that fixes the `Number` declares it, and the declaration is
    /// one a key of the same marker resolves from a parsed header.
    #[test]
    fn a_definition_declares_its_markers_number() {
        fn check<V: FieldKind>(def: &dyn FieldDescription, number: Number) {
            assert_eq!(def.number(), number, "{}", def.name());
            assert!(V::holds_number(number), "{} must resolve back", def.name());
        }
        check::<Scalar<i32>>(&InfoFieldDef::<Scalar<i32>>::new("A", ""), Number::Count(1));
        check::<Scalar<f32>>(&FormatFieldDef::<Scalar<f32>>::new("A", ""), Number::Count(1));
        check::<Str>(&InfoFieldDef::<Str>::with_number("A", Number::Unknown, ""), Number::Unknown);
        check::<Str>(
            &FormatFieldDef::<Str>::with_number("A", Number::Count(1), ""),
            Number::Count(1),
        );
        check::<Flag>(&InfoFieldDef::<Flag>::new("A", ""), Number::Count(0));
        check::<Gt>(&FormatFieldDef::<Gt>::new("GT", ""), Number::Count(1));
        check::<Arr<i32>>(
            &InfoFieldDef::<Arr<i32>>::with_number("A", Number::Genotypes, ""),
            Number::Genotypes,
        );
        check::<OptArr<f32>>(
            &FormatFieldDef::<OptArr<f32>>::with_number("A", Number::Count(4), ""),
            Number::Count(4),
        );
    }

    // r[verify record_encoder.field_def+5]
    #[test]
    fn a_genotype_definition_names_only_format_gt() {
        // An INFO `Gt` definition has no constructor at all.
        let mut builder = crate::vcf::header::VcfHeader::builder().formats();
        let err = builder.register_format(&FormatFieldDef::<Gt>::new("GQ", "")).unwrap_err();
        assert!(
            matches!(err, crate::vcf::VcfHeaderError::Key(KeyError::NotGenotype { ref name }) if name == "GQ"),
            "{err:?}"
        );
        builder.register_format(&FormatFieldDef::<Gt>::new("GT", "")).unwrap();
    }

    // r[verify record_encoder.filter_id]
    #[test]
    fn filter_id_pass() {
        let pass = FilterId::PASS;
        assert_eq!(pass.dict_idx, 0);
        assert_eq!(pass.name(), "PASS");
    }

    // r[verify record_encoder.typed_keys]
    #[test]
    fn typed_keys_prevent_misuse() {
        let _int_key: InfoInt = InfoKey(
            FieldId::new(
                DictId::new(0).unwrap(),
                SmolStr::from("DP"),
                ValueType::Integer,
                Number::Count(1),
                Section::Info,
                HeaderId::next(),
            ),
            PhantomData,
        );
        let _float_key: InfoFloat = InfoKey(
            FieldId::new(
                DictId::new(1).unwrap(),
                SmolStr::from("BQ"),
                ValueType::Float,
                Number::Count(1),
                Section::Info,
                HeaderId::next(),
            ),
            PhantomData,
        );
        let _flag_key: InfoFlag = InfoKey(
            FieldId::new(
                DictId::new(2).unwrap(),
                SmolStr::from("CPG"),
                ValueType::Flag,
                Number::Count(0),
                Section::Info,
                HeaderId::next(),
            ),
            PhantomData,
        );
    }

    struct TestSetup {
        header: Arc<crate::vcf::header::VcfHeader>,
        contig: ContigId,
        dp_info: InfoInt,
        bq_info: InfoFloat,
        db_flag: InfoFlag,
        ad_info: InfoInts,
        gt_fmt: FormatGt,
        dp_fmt: FormatInt,
    }

    impl TestSetup {
        fn new() -> Self {
            let mut builder = crate::vcf::header::VcfHeader::builder();
            let contig = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
            let mut builder = builder.infos();
            let dp_info =
                builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
            let bq_info = builder
                .register_info(&InfoFieldDef::<Scalar<f32>>::new("BQ", "Base quality"))
                .unwrap();
            let db_flag = builder.register_info(&InfoFieldDef::<Flag>::new("DB", "dbSNP")).unwrap();
            let ad_info = builder
                .register_info(&InfoFieldDef::<Arr<i32>>::with_number(
                    "AD",
                    Number::ReferenceAlternateBases,
                    "Allele depth",
                ))
                .unwrap();
            let mut builder = builder.formats();
            let gt_fmt =
                builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
            let dp_fmt = builder
                .register_format(&FormatFieldDef::<Scalar<i32>>::new("DP", "Read depth"))
                .unwrap();
            let mut builder = builder.samples();
            builder.add_sample("S1").unwrap();
            let header = Arc::new(builder.build().unwrap());
            Self { header, contig, dp_info, bq_info, db_flag, ad_info, gt_fmt, dp_fmt }
        }
    }

    // r[verify record_encoder.register]
    // r[verify record_encoder.register_contig]
    // r[verify record_encoder.field_id+2]
    // r[verify record_encoder.contig_id]
    #[test]
    fn register_resolves_correct_dict_indices() {
        let setup = TestSetup::new();
        assert_eq!(setup.contig.tid(), 0);
        assert_eq!(setup.contig.name(), "chr1");

        let map = setup.header.string_map();
        assert_eq!(map.get("DP").unwrap(), setup.dp_info.id().dict_idx() as usize, "INFO DP");
        assert_eq!(map.get("BQ").unwrap(), setup.bq_info.id().dict_idx() as usize, "INFO BQ");
        assert_eq!(map.get("DB").unwrap(), setup.db_flag.id().dict_idx() as usize, "INFO DB");
        assert_eq!(map.get("GT").unwrap(), setup.gt_fmt.id().dict_idx() as usize, "FORMAT GT");

        // The declared Type and the GT flag, from registering and from resolving.
        let resolved = setup.header.format_key::<crate::vcf::reader::Any>("GT").unwrap();
        for (id, typ, genotype) in [
            (setup.dp_info.id(), ValueType::Integer, false),
            (setup.bq_info.id(), ValueType::Float, false),
            (setup.db_flag.id(), ValueType::Flag, false),
            (setup.gt_fmt.id(), ValueType::String, true),
            (&resolved.0, ValueType::String, true),
            (setup.dp_fmt.id(), ValueType::Integer, false),
        ] {
            assert_eq!((id.typ, id.genotype), (typ, genotype), "{}", id.name());
        }
        let dp = setup.header.info_key::<crate::vcf::reader::Any>("DP").unwrap();
        assert_eq!(dp.0, *setup.dp_info.id());
    }

    // r[verify record_encoder.field_def+5]
    // r[verify record_encoder.field_def_types]
    #[test]
    fn field_def_generic_keys_encode() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 42).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.filters+3]
    #[test]
    fn filter_fail_bcf_and_vcf_text() {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut builder = builder.filters();
        let low_dp = builder.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let header = Arc::new(builder.build().unwrap());

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let pos = Pos1::new(10).unwrap();

        let bcf_bytes = {
            let mut buf = Vec::new();
            let writer = Writer::new(&mut buf, OutputFormat::Bcf);
            let mut writer = writer.write_header(&header).unwrap();
            writer
                .begin_record(&contig, pos, &alleles, None)
                .unwrap()
                .filter_fail([&low_dp])
                .unwrap()
                .emit()
                .unwrap();
            writer.finish().unwrap();
            buf
        };
        assert!(!bcf_bytes.is_empty(), "BCF output must not be empty");

        let vcf_text = {
            let mut buf = Vec::new();
            let writer = Writer::new(&mut buf, OutputFormat::Vcf);
            let mut writer = writer.write_header(&header).unwrap();
            writer
                .begin_record(&contig, pos, &alleles, None)
                .unwrap()
                .filter_fail([&low_dp])
                .unwrap()
                .emit()
                .unwrap();
            writer.finish().unwrap();
            String::from_utf8(buf).unwrap()
        };
        let data_line = vcf_text.lines().find(|l| !l.starts_with('#')).unwrap();
        let filter_col = data_line.split('\t').nth(6).unwrap();
        assert_eq!(filter_col, "lowDp", "filter_fail must write the filter name");
    }

    // r[verify record_encoder.info_methods+3]
    // r[verify record_encoder.format_methods+3]
    // r[verify record_encoder.key_encode+2]
    // r[verify record_encoder.bcf_encoding]
    #[test]
    fn typestate_api_bcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        setup.bq_info.encode(&mut enc, 35.5).unwrap();
        setup.db_flag.encode(&mut enc).unwrap();
        setup.ad_info.encode(&mut enc, &[30, 20]).unwrap();
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[45]).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.vcf_encoding]
    #[test]
    fn typestate_emit_without_samples() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(cols.first(), Some(&"chr1"), "CHROM");
        assert_eq!(cols.get(1), Some(&"100"), "POS");
        assert_eq!(cols.get(6), Some(&"PASS"), "FILTER");
    }

    // r[verify record_encoder.filters+3]
    #[test]
    fn typestate_filter_fail() {
        let mut hb = crate::vcf::header::VcfHeader::builder();
        let contig = hb.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut hb = hb.filters();
        let ld = hb.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let header = Arc::new(hb.build().unwrap());
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        writer
            .begin_record(&contig, Pos1::new(10).unwrap(), &alleles, None)
            .unwrap()
            .filter_fail([&ld])
            .unwrap()
            .emit()
            .unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let filter_col = data_line.split('\t').nth(6).unwrap();
        assert_eq!(filter_col, "lowDp");

        for format in [OutputFormat::Vcf, OutputFormat::Bcf] {
            let mut buf = Vec::new();
            let mut writer = Writer::new(&mut buf, format).write_header(&header).unwrap();
            let begun = writer.begin_record(&contig, Pos1::new(10).unwrap(), &alleles, None);
            let err = begun.unwrap().filter_fail(std::iter::empty()).err().unwrap();
            assert!(matches!(err, VcfError::EmptyFilterList), "{format:?}: {err:?}");
        }
    }

    // r[verify record_encoder.vcf_encoding]
    #[test]
    fn typestate_no_filter() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        writer
            .begin_record(&setup.contig, Pos1::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass()
            .emit()
            .unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(data_line.split('\t').nth(6), Some("PASS"));
    }

    // r[verify record_encoder.bcf_encoding]
    // r[verify record_encoder.vcf_encoding]
    #[test]
    fn typestate_multiple_records() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        for pos in [1u32, 100, 200] {
            writer
                .begin_record(&setup.contig, Pos1::new(pos).unwrap(), &alleles, None)
                .unwrap()
                .filter_pass()
                .emit()
                .unwrap();
        }
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 3, "expected 3 data records");
    }

    // ── Multi-sample tests ───────────────────────────────────────────

    fn multi_sample_setup()
    -> (Arc<crate::vcf::header::VcfHeader>, ContigId, InfoInt, FormatGt, FormatInt) {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut builder = builder.infos();
        let dp_info =
            builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
        let mut builder = builder.formats();
        let gt_fmt = builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
        let dp_fmt = builder
            .register_format(&FormatFieldDef::<Scalar<i32>>::new("DP", "Read depth"))
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();
        builder.add_sample("S2").unwrap();
        builder.add_sample("S3").unwrap();
        let header = Arc::new(builder.build().unwrap());
        (header, contig, dp_info, gt_fmt, dp_fmt)
    }

    // r[verify record_encoder.format_methods+3]
    #[test]
    fn multi_sample_bcf() {
        let (header, contig, dp_info, gt_fmt, dp_fmt) = multi_sample_setup();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc =
            writer.begin_record(&contig, Pos1::new(100).unwrap(), &alleles, Some(30.0)).unwrap();
        let mut enc = enc.filter_pass();
        dp_info.encode(&mut enc, 150).unwrap();
        let mut enc = enc.begin_samples();
        gt_fmt
            .encode(
                &mut enc,
                &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
            )
            .unwrap();
        dp_fmt.encode(&mut enc, &[45, 52, 38]).unwrap();
        enc.emit().unwrap();

        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.format_methods+3]
    #[test]
    fn multi_sample_vcf_text() {
        let (header, contig, dp_info, gt_fmt, dp_fmt) = multi_sample_setup();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc =
            writer.begin_record(&contig, Pos1::new(100).unwrap(), &alleles, Some(30.0)).unwrap();
        let mut enc = enc.filter_pass();
        dp_info.encode(&mut enc, 150).unwrap();
        let mut enc = enc.begin_samples();
        gt_fmt
            .encode(
                &mut enc,
                &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
            )
            .unwrap();
        dp_fmt.encode(&mut enc, &[45, 52, 38]).unwrap();
        enc.emit().unwrap();

        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        // 8 fixed + FORMAT + 3 samples = 12 columns
        assert_eq!(cols.len(), 12, "expected 12 columns for 3-sample VCF");
        assert_eq!(cols[8], "GT:DP", "FORMAT column");
        assert_eq!(cols[9], "0/1:45", "sample 1");
        assert_eq!(cols[10], "0/0:52", "sample 2");
        assert_eq!(cols[11], "1/1:38", "sample 3");
    }

    // r[verify record_encoder.format_methods+3]
    #[test]
    fn single_sample_slice_api() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap();
        let mut enc = enc.filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        let mut enc = enc.begin_samples();
        // Single-sample callers pass 1-element slices
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[45]).unwrap();
        enc.emit().unwrap();

        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.format_single_sample]
    #[test]
    fn format_single_sample_helper_matches_nested_slice() {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut builder = builder.formats();
        let pl_fmt: FormatInts = builder
            .register_format(&FormatFieldDef::<Arr<i32>>::with_number(
                "PL",
                Number::Genotypes,
                "Phred likelihoods",
            ))
            .unwrap();
        let mut builder = builder.samples();
        builder.add_sample("S1").unwrap();
        let header = Arc::new(builder.build().unwrap());

        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&contig, Pos1::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // The flat helper must produce the same bytes as the explicit &[&[..]] form.
        pl_fmt.encode_single_sample(&mut enc, &[10, 20, 30]).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(cols[8], "PL");
        assert_eq!(cols[9], "10,20,30");
    }

    #[test]
    fn format_single_sample_helper_rejects_multi_sample() {
        let mut hb = crate::vcf::header::VcfHeader::builder();
        let contig = hb.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut hb = hb.formats();
        let af: FormatFloats = hb
            .register_format(&FormatFieldDef::<Arr<f32>>::with_number(
                "AF",
                Number::Unknown,
                "Allele freq",
            ))
            .unwrap();
        let mut hb = hb.samples();
        hb.add_sample("A").unwrap();
        hb.add_sample("B").unwrap();
        let header = Arc::new(hb.build().unwrap());

        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&contig, Pos1::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Header has 2 samples → the single-sample helper must refuse.
        let result = af.encode_single_sample(&mut enc, &[0.5]);
        assert!(matches!(result, Err(VcfError::SampleCountMismatch { .. })));
    }

    #[test]
    fn format_int_panics_on_sample_count_mismatch() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Pass 2 values for 1-sample header — must return Err
        let result = setup.dp_fmt.encode(&mut enc, &[45, 52]);
        assert!(result.is_err());
    }

    #[test]
    fn format_gt_panics_on_sample_count_mismatch() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Pass 3 genotypes for 1-sample header — must return Err
        let result = setup.gt_fmt.encode(
            &mut enc,
            &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
        );
        assert!(result.is_err());
    }

    // r[verify record_encoder.mixed_ploidy]
    // r[verify bcf_writer.gt_encoding]
    /// Haploid, diploid and triploid samples share one BCF GT column three
    /// wide, each shorter sample padded with int8 end-of-vector (`0x81`), and
    /// a genotype with no alleles is one missing call. The bytes are the BCF2
    /// table's: `(allele + 1) << 1 | phased`, a haploid call phased.
    #[test]
    fn mixed_ploidy_pads_with_end_of_vector() {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
        let mut builder = builder.formats();
        let gt_fmt = builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
        let mut builder = builder.samples();
        for name in ["S1", "S2", "S3", "S4"] {
            builder.add_sample(name).unwrap();
        }
        let header = builder.build().unwrap();
        let mut buf = Vec::new();
        let mut writer = Writer::new(&mut buf, OutputFormat::Bcf).write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&contig, Pos1::new(100).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass()
            .begin_samples();
        let triploid = Genotype {
            alleles: [Some(0), Some(1), Some(1)].into_iter().collect(),
            phased: [false, false].into_iter().collect(),
        };
        let none = Genotype { alleles: SmallVec::new(), phased: SmallVec::new() };
        gt_fmt
            .encode(&mut enc, &[Genotype::haploid(0), Genotype::unphased(0, 1), triploid, none])
            .unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();

        let mut data = Vec::new();
        crate::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(buf))
            .read_to_end(&mut data)
            .unwrap();
        let gt_idx = u8::try_from(header.string_map().get("GT").unwrap()).unwrap();
        let expected = [
            0x11, gt_idx, 0x31, // key, then 3 × int8
            0x03, 0x81, 0x81, // 0
            0x02, 0x04, 0x81, // 0/1
            0x02, 0x04, 0x04, // 0/1/1
            0x00, 0x81, 0x81, // no alleles: a missing call
        ];
        assert!(data.ends_with(&expected), "{:02x?}", &data[data.len() - expected.len()..]);
    }

    // ── INFO dedup tests ────────────────────────────────────────────

    // r[verify record_encoder.info_dedup+2]
    #[test]
    fn info_overwrite_last_value_wins_vcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        setup.dp_info.encode(&mut enc, 100).unwrap(); // overwrite
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let info_col = data_line.split('\t').nth(7).unwrap();
        assert_eq!(info_col, "DP=100", "duplicate INFO should keep only the last value");
    }

    // r[verify record_encoder.info_dedup+2]
    #[test]
    fn info_overwrite_first_of_multiple_vcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        setup.bq_info.encode(&mut enc, 35.5).unwrap();
        setup.dp_info.encode(&mut enc, 100).unwrap(); // overwrite first
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let info_col = data_line.split('\t').nth(7).unwrap();
        // In-place overwrite (htslib parity): DP keeps its original first position.
        assert_eq!(info_col, "DP=100;BQ=35.5", "overwritten field keeps its original position");
    }

    // r[verify record_encoder.info_dedup+2]
    #[test]
    fn info_overwrite_middle_vcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        setup.bq_info.encode(&mut enc, 35.5).unwrap();
        setup.db_flag.encode(&mut enc).unwrap();
        setup.bq_info.encode(&mut enc, 99.0).unwrap(); // overwrite middle
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let info_col = data_line.split('\t').nth(7).unwrap();
        // In-place overwrite (htslib parity): BQ keeps its middle position.
        assert_eq!(info_col, "DP=50;BQ=99;DB", "overwritten middle field keeps its position");
    }

    // r[verify record_encoder.info_dedup+2]
    #[test]
    fn info_overwrite_last_value_wins_bcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        setup.bq_info.encode(&mut enc, 35.5).unwrap();
        setup.dp_info.encode(&mut enc, 100).unwrap(); // overwrite
        enc.emit().unwrap();
        writer.finish().unwrap();

        // Parse with noodles to verify the BCF is valid and has the correct value
        use noodles::bcf;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoBufValue;
        let mut reader = bcf::io::Reader::new(std::io::Cursor::new(&buf));
        let header = reader.read_header().unwrap();
        let mut records = reader.record_bufs(&header);
        let rec = records.next().unwrap().unwrap();
        let dp =
            rec.info().get("DP").expect("DP should exist").expect("DP value should not be None");
        assert!(
            matches!(dp, &InfoBufValue::Integer(100)),
            "BCF DP should be overwritten to 100, got {dp:?}"
        );
    }

    // r[verify record_encoder.info_dedup+2]
    #[test]
    fn info_dedup_tracker_resets_between_records() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();

        // Record 1: write DP=50
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50).unwrap();
        enc.emit().unwrap();

        // Record 2: write DP=100 — should NOT be treated as a duplicate
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(200).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 100).unwrap();
        enc.emit().unwrap();

        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 2);
        assert_eq!(data_lines[0].split('\t').nth(7).unwrap(), "DP=50");
        assert_eq!(data_lines[1].split('\t').nth(7).unwrap(), "DP=100");
    }

    // ── FORMAT dedup tests ──────────────────────────────────────────

    // r[verify record_encoder.format_dedup+2]
    #[test]
    fn format_overwrite_last_value_wins_vcf() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[45]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[99]).unwrap(); // overwrite
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(cols[8], "GT:DP", "FORMAT column should list each key once");
        assert_eq!(cols[9], "0/1:99", "overwritten FORMAT DP should be 99");
    }

    // r[verify record_encoder.format_dedup+2]
    #[test]
    fn format_overwrite_first_of_multiple_vcf() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[45]).unwrap();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(1, 1)]).unwrap(); // overwrite first
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        // In-place overwrite (htslib parity): GT keeps its original first column.
        assert_eq!(cols[8], "GT:DP", "overwritten GT keeps its original position");
        assert_eq!(cols[9], "1/1:45", "sample data follows the unchanged key order");
    }

    // r[verify record_encoder.format_dedup+2]
    #[test]
    fn format_overwrite_bcf() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[45]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[99]).unwrap(); // overwrite
        enc.emit().unwrap();
        writer.finish().unwrap();

        // Verify BCF is parseable by noodles and structurally correct.
        // Value correctness is validated by format_overwrite_last_value_wins_vcf
        // and the bcftools round-trip test in tests/vcf_bcftools_roundtrip.rs.
        use noodles::bcf;
        use noodles::vcf::variant::record::Samples as _;
        let mut reader = bcf::io::Reader::new(std::io::Cursor::new(&buf));
        let header = reader.read_header().unwrap();
        let mut records = reader.record_bufs(&header);
        let rec = records.next().unwrap().unwrap();
        let keys: Vec<&str> = rec.samples().column_names(&header).map(|r| r.unwrap()).collect();
        assert_eq!(
            keys.iter().filter(|k| **k == "DP").count(),
            1,
            "BCF should have exactly one DP FORMAT key"
        );
    }
}
