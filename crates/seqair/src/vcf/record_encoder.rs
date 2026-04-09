//! Format-agnostic record encoder with typed field keys.
//!
//! The [`RecordEncoder`] trait provides a unified API for encoding VCF records
//! into both BCF binary and VCF text formats. Typed field keys ([`InfoKey`],
//! [`FormatKey`]) enforce value types at compile time.
//!
//! # Field definitions
//!
//! [`InfoFieldDef`] and [`FormatFieldDef`] combine header metadata with
//! type-safe key resolution. Register them on [`VcfHeaderBuilder`](super::header::VcfHeaderBuilder)
//! to build the header and resolve keys in one step:
//!
//! ```ignore
//! let dp = InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Depth");
//! let dp_key = builder.register_info(&dp)?;
//! // dp_key is InfoKey<Scalar<i32>> — can only encode i32 values
//! ```

use super::header::{Number, ValueType};
use super::record::Genotype;
use crate::vcf::alleles::Alleles;
use crate::vcf::error::VcfError;
use seqair_types::{One, Pos, SmolStr};
use std::marker::PhantomData;

// ── Value type markers (uninhabited) ───────────────────────────────────

// r[impl record_encoder.typed_keys]

/// Marker for a single scalar value.
pub enum Scalar<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for a variable-length array.
pub enum Arr<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for an array with optional (missing) elements.
pub enum OptArr<T> {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible, PhantomData<T>),
}

/// Marker for a flag field (no value).
pub enum Flag {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

/// Marker for a string field.
pub enum Str {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

/// Marker for genotype (GT) format field.
pub enum Gt {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

// ── Field identifiers ──────────────────────────────────────────────────

// r[impl record_encoder.field_id]
/// Base field identifier carrying both a BCF dictionary index and a string name.
///
/// The BCF encoder uses `dict_idx` for binary encoding; the VCF text encoder
/// uses `name` for text output. Both are resolved once at registration time.
#[derive(Debug, Clone)]
pub struct FieldId {
    pub(crate) dict_idx: u32,
    pub(crate) name: SmolStr,
}

impl FieldId {
    /// The field name (e.g., `"DP"`, `"GT"`).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The BCF dictionary index.
    pub fn dict_idx(&self) -> u32 {
        self.dict_idx
    }
}

// r[impl record_encoder.contig_id]
/// Resolved contig identifier carrying both integer tid and string name.
#[derive(Debug, Clone)]
pub struct ContigId {
    pub(crate) tid: u32,
    pub(crate) name: SmolStr,
}

impl ContigId {
    /// The contig name (e.g., `"chr1"`).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The integer tid (contig index in header order).
    pub fn tid(&self) -> u32 {
        self.tid
    }
}

// r[impl record_encoder.filter_id]
/// Resolved filter identifier carrying both dictionary index and name.
#[derive(Debug, Clone)]
pub struct FilterId {
    pub(crate) dict_idx: u32,
    pub(crate) name: SmolStr,
}

impl FilterId {
    /// The PASS filter (always BCF dictionary index 0).
    pub fn pass() -> Self {
        Self { dict_idx: 0, name: SmolStr::new_static("PASS") }
    }

    /// The filter name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The BCF dictionary index.
    pub fn dict_idx(&self) -> u32 {
        self.dict_idx
    }
}

// ── Typed field keys ───────────────────────────────────────────────────

// r[impl record_encoder.key_types]

/// Typed INFO field key. The value type marker `V` enforces that only the
/// correct value type can be encoded with this key.
#[derive(Debug, Clone)]
pub struct InfoKey<V>(pub(crate) FieldId, pub(crate) PhantomData<V>);

/// Typed FORMAT field key.
#[derive(Debug, Clone)]
pub struct FormatKey<V>(pub(crate) FieldId, pub(crate) PhantomData<V>);

impl<V> InfoKey<V> {
    /// The underlying field identifier.
    pub fn id(&self) -> &FieldId {
        &self.0
    }
}

impl<V> FormatKey<V> {
    /// The underlying field identifier.
    pub fn id(&self) -> &FieldId {
        &self.0
    }
}

/// Convenience type aliases for common key types.
pub type InfoInt = InfoKey<Scalar<i32>>;
pub type InfoFloat = InfoKey<Scalar<f32>>;
pub type InfoInts = InfoKey<Arr<i32>>;
pub type InfoFloats = InfoKey<Arr<f32>>;
pub type InfoFlag = InfoKey<Flag>;
pub type InfoString = InfoKey<Str>;
pub type InfoIntOpts = InfoKey<OptArr<i32>>;

pub type FormatGt = FormatKey<Gt>;
pub type FormatInt = FormatKey<Scalar<i32>>;
pub type FormatFloat = FormatKey<Scalar<f32>>;

// ── Key encode methods ─────────────────────────────────────────────────

// r[impl record_encoder.key_encode]

impl InfoKey<Scalar<i32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), value: i32) {
        enc.info_int(&self.0, value);
    }
}

impl InfoKey<Scalar<f32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), value: f32) {
        enc.info_float(&self.0, value);
    }
}

impl InfoKey<Arr<i32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), values: &[i32]) {
        enc.info_ints(&self.0, values);
    }
}

impl InfoKey<Arr<f32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), values: &[f32]) {
        enc.info_floats(&self.0, values);
    }
}

impl InfoKey<Flag> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized)) {
        enc.info_flag(&self.0);
    }
}

impl InfoKey<Str> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), value: &str) {
        enc.info_string(&self.0, value);
    }
}

impl InfoKey<OptArr<i32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), values: &[Option<i32>]) {
        enc.info_int_opts(&self.0, values);
    }
}

impl FormatKey<Gt> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), gt: &Genotype) {
        enc.format_gt(&self.0, gt);
    }
}

impl FormatKey<Scalar<i32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), value: i32) {
        enc.format_int(&self.0, value);
    }
}

impl FormatKey<Scalar<f32>> {
    pub fn encode(&self, enc: &mut (impl RecordEncoder + ?Sized), value: f32) {
        enc.format_float(&self.0, value);
    }
}

// ── Field definitions ──────────────────────────────────────────────────

// r[impl record_encoder.field_def]
// r[impl record_encoder.field_def_types]

/// INFO field definition carrying all metadata for header construction,
/// key resolution, and documentation. Parameterized by value type marker
/// for type-safe key resolution.
pub struct InfoFieldDef<V> {
    pub name: &'static str,
    pub number: Number,
    pub value_type: ValueType,
    pub description: &'static str,
    _marker: PhantomData<V>,
}

// r[impl record_encoder.field_def_const]
impl<V> InfoFieldDef<V> {
    pub const fn new(
        name: &'static str,
        number: Number,
        value_type: ValueType,
        description: &'static str,
    ) -> Self {
        Self { name, number, value_type, description, _marker: PhantomData }
    }
}

/// FORMAT field definition.
pub struct FormatFieldDef<V> {
    pub name: &'static str,
    pub number: Number,
    pub value_type: ValueType,
    pub description: &'static str,
    _marker: PhantomData<V>,
}

impl<V> FormatFieldDef<V> {
    pub const fn new(
        name: &'static str,
        number: Number,
        value_type: ValueType,
        description: &'static str,
    ) -> Self {
        Self { name, number, value_type, description, _marker: PhantomData }
    }
}

/// FILTER field definition.
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
/// Read access to field metadata for documentation generation.
pub trait FieldDescription {
    fn name(&self) -> &str;
    fn number(&self) -> Number;
    fn value_type(&self) -> ValueType;
    fn description(&self) -> &str;
}

impl<V> FieldDescription for InfoFieldDef<V> {
    fn name(&self) -> &str {
        self.name
    }
    fn number(&self) -> Number {
        self.number
    }
    fn value_type(&self) -> ValueType {
        self.value_type
    }
    fn description(&self) -> &str {
        self.description
    }
}

impl<V> FieldDescription for FormatFieldDef<V> {
    fn name(&self) -> &str {
        self.name
    }
    fn number(&self) -> Number {
        self.number
    }
    fn value_type(&self) -> ValueType {
        self.value_type
    }
    fn description(&self) -> &str {
        self.description
    }
}

// ── RecordEncoder trait ────────────────────────────────────────────────

// r[impl record_encoder.trait]
/// Format-agnostic record encoder. Implemented by both [`BcfRecordEncoder`](super::encoder::BcfRecordEncoder)
/// and `VcfRecordEncoder` to enable the same encoding code for VCF text and BCF binary output.
///
/// All INFO/FORMAT/filter methods are infallible (they write to in-memory buffers).
/// Only [`begin`](Self::begin) and [`emit`](Self::emit) return `Result`.
pub trait RecordEncoder {
    // r[impl record_encoder.begin]
    /// Begin a new record: clear per-record state, write/buffer fixed fields.
    fn begin(
        &mut self,
        contig: &ContigId,
        pos: Pos<One>,
        alleles: &Alleles,
        qual: Option<f32>,
    ) -> Result<(), VcfError>;

    // r[impl record_encoder.filters]
    /// Record passed all filters.
    fn filter_pass(&mut self);

    /// Record failed one or more filters.
    fn filter_fail(&mut self, filters: &[&FilterId]);

    // r[impl record_encoder.info_methods]
    /// Encode a scalar integer INFO field.
    fn info_int(&mut self, id: &FieldId, value: i32);
    /// Encode a scalar float INFO field.
    fn info_float(&mut self, id: &FieldId, value: f32);
    /// Encode an integer array INFO field.
    fn info_ints(&mut self, id: &FieldId, values: &[i32]);
    /// Encode a float array INFO field.
    fn info_floats(&mut self, id: &FieldId, values: &[f32]);
    /// Encode a flag INFO field (presence only, no value).
    fn info_flag(&mut self, id: &FieldId);
    /// Encode a string INFO field.
    fn info_string(&mut self, id: &FieldId, value: &str);
    /// Encode an integer array INFO field with optional (missing) elements.
    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]);

    // r[impl record_encoder.format_methods]
    /// Declare the number of samples. Must be called before any format method.
    fn begin_samples(&mut self, n: u32);
    /// Encode a genotype (GT) FORMAT field.
    fn format_gt(&mut self, id: &FieldId, gt: &Genotype);
    /// Encode a scalar integer FORMAT field.
    fn format_int(&mut self, id: &FieldId, value: i32);
    /// Encode a scalar float FORMAT field.
    fn format_float(&mut self, id: &FieldId, value: f32);

    // r[impl record_encoder.state_queries]
    /// Number of alleles (REF + ALTs) for the current record.
    fn n_allele(&self) -> usize;
    /// Number of alternate alleles for the current record.
    fn n_alt(&self) -> usize;

    // r[impl record_encoder.emit]
    /// Finalize and write the record to output. This is the only method that
    /// performs I/O and may return an error.
    fn emit(&mut self) -> Result<(), VcfError>;
}

// ── Custom type encoding traits ────────────────────────────────────────

// r[impl record_encoder.encode_info_trait]
// r[impl record_encoder.encode_dyn]
/// Trait for types that can encode themselves as one or more VCF INFO fields.
///
/// The associated `Key` type specifies which pre-resolved keys are needed.
/// For single-field types this is typically an [`InfoKey<V>`]; for types
/// that expand to multiple VCF fields (e.g., strand-specific OT/OB pairs)
/// this can be a tuple of keys.
pub trait EncodeInfo {
    type Key;
    fn encode_info(&self, enc: &mut dyn RecordEncoder, key: &Self::Key);
}

// r[impl record_encoder.encode_format_trait]
/// Trait for types that can encode themselves as a VCF FORMAT field.
///
/// Implementations that produce no output for certain values (e.g., unknown
/// methylation status) simply do not call any encoder methods — the FORMAT
/// key will not appear in the output.
pub trait EncodeFormat {
    type Key;
    fn encode_format(&self, enc: &mut dyn RecordEncoder, key: &Self::Key);
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify record_encoder.field_def_const]
    #[test]
    fn field_defs_are_const() {
        const DP: InfoFieldDef<Scalar<i32>> =
            InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Combined depth");
        const GT: FormatFieldDef<Gt> =
            FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "Genotype");
        const LOW_DP: FilterFieldDef = FilterFieldDef::new("lowDp", "Low read depth");

        assert_eq!(DP.name, "DP");
        assert_eq!(GT.name, "GT");
        assert_eq!(LOW_DP.name, "lowDp");
    }

    // r[verify record_encoder.field_description]
    #[test]
    fn field_description_trait() {
        let dp = InfoFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Combined depth",
        );
        let desc: &dyn FieldDescription = &dp;
        assert_eq!(desc.name(), "DP");
        assert_eq!(desc.number(), Number::Count(1));
        assert_eq!(desc.value_type(), ValueType::Integer);
        assert_eq!(desc.description(), "Combined depth");
    }

    // r[verify record_encoder.filter_id]
    #[test]
    fn filter_id_pass() {
        let pass = FilterId::pass();
        assert_eq!(pass.dict_idx, 0);
        assert_eq!(pass.name(), "PASS");
    }

    // r[verify record_encoder.typed_keys]
    #[test]
    fn typed_keys_prevent_misuse() {
        // This test verifies compile-time type safety by construction.
        // InfoKey<Scalar<i32>> can only encode i32 values.
        // InfoKey<Scalar<f32>> can only encode f32 values.
        // If you tried to call InfoInt::encode with a f32, it wouldn't compile.
        let _int_key: InfoInt =
            InfoKey(FieldId { dict_idx: 0, name: SmolStr::from("DP") }, PhantomData);
        let _float_key: InfoFloat =
            InfoKey(FieldId { dict_idx: 1, name: SmolStr::from("BQ") }, PhantomData);
        let _flag_key: InfoFlag =
            InfoKey(FieldId { dict_idx: 2, name: SmolStr::from("CPG") }, PhantomData);
    }

    // ── Registration + RecordEncoder integration tests ─────────────────

    use crate::vcf::bcf_writer::BcfWriter;
    use crate::vcf::header::ContigDef;
    use crate::vcf::record::{SampleValue, VcfRecordBuilder};
    use seqair_types::Base;
    use std::sync::Arc;

    /// Build a header using register_* methods for the RecordEncoder path.
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

            let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
            let dp_info = builder
                .register_info(&InfoFieldDef::new(
                    "DP",
                    Number::Count(1),
                    ValueType::Integer,
                    "Depth",
                ))
                .unwrap();
            let bq_info = builder
                .register_info(&InfoFieldDef::new(
                    "BQ",
                    Number::Count(1),
                    ValueType::Float,
                    "Base quality",
                ))
                .unwrap();
            let db_flag = builder
                .register_info(&InfoFieldDef::new("DB", Number::Count(0), ValueType::Flag, "dbSNP"))
                .unwrap();
            let ad_info = builder
                .register_info(&InfoFieldDef::new(
                    "AD",
                    Number::ReferenceAlternateBases,
                    ValueType::Integer,
                    "Allele depth",
                ))
                .unwrap();
            let gt_fmt = builder
                .register_format(&FormatFieldDef::new(
                    "GT",
                    Number::Count(1),
                    ValueType::String,
                    "Genotype",
                ))
                .unwrap();
            let dp_fmt = builder
                .register_format(&FormatFieldDef::new(
                    "DP",
                    Number::Count(1),
                    ValueType::Integer,
                    "Read depth",
                ))
                .unwrap();

            let header = Arc::new(builder.add_sample("S1").unwrap().build().unwrap());

            Self { header, contig, dp_info, bq_info, db_flag, ad_info, gt_fmt, dp_fmt }
        }
    }

    /// Encode a simple SNV record via RecordEncoder and return the raw BCF bytes.
    fn encode_via_record_encoder(setup: &TestSetup) -> Vec<u8> {
        let mut buf = Vec::new();
        let mut writer = BcfWriter::new(&mut buf, setup.header.clone(), false);
        writer.write_header().unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer.record_encoder();
        enc.begin(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0)).unwrap();
        enc.filter_pass();
        setup.dp_info.encode(&mut enc, 50);
        setup.bq_info.encode(&mut enc, 35.5);
        setup.db_flag.encode(&mut enc);
        setup.ad_info.encode(&mut enc, &[30, 20]);
        enc.begin_samples(1);
        setup.gt_fmt.encode(&mut enc, &Genotype::unphased(0, 1));
        setup.dp_fmt.encode(&mut enc, 45);
        enc.emit().unwrap();

        writer.finish().unwrap();
        buf
    }

    /// Encode the same record via VcfRecordBuilder and return the raw BCF bytes.
    fn encode_via_vcf_record(setup: &TestSetup) -> Vec<u8> {
        let mut buf = Vec::new();
        let mut writer = BcfWriter::new(&mut buf, setup.header.clone(), false);
        writer.write_header().unwrap();

        let record = VcfRecordBuilder::new(
            "chr1",
            Pos::<One>::new(100).unwrap(),
            Alleles::snv(Base::A, Base::T).unwrap(),
        )
        .qual(30.0)
        .filter_pass()
        .info_integer("DP", 50)
        .info_float("BQ", 35.5)
        .info_flag("DB")
        .info_integers("AD", &[30, 20])
        .format_keys(&["GT", "DP"])
        .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1)), SampleValue::Integer(45)])
        .build(&setup.header)
        .unwrap();

        writer.write_vcf_record(&record).unwrap();
        writer.finish().unwrap();
        buf
    }

    // r[verify record_encoder.bcf_impl]
    // r[verify record_encoder.record_path_equivalence]
    // r[verify record_encoder.begin]
    // r[verify record_encoder.emit]
    // r[verify record_encoder.filters]
    // r[verify record_encoder.info_methods]
    // r[verify record_encoder.format_methods]
    // r[verify record_encoder.state_queries]
    // r[verify record_encoder.key_encode]
    #[test]
    fn record_encoder_matches_vcf_record_path() {
        let setup = TestSetup::new();
        let encoder_bytes = encode_via_record_encoder(&setup);
        let record_bytes = encode_via_vcf_record(&setup);
        assert_eq!(
            encoder_bytes, record_bytes,
            "RecordEncoder and VcfRecord paths must produce identical BCF output"
        );
    }

    // r[verify record_encoder.register]
    // r[verify record_encoder.register_contig]
    // r[verify record_encoder.field_id]
    // r[verify record_encoder.contig_id]
    #[test]
    fn register_resolves_correct_dict_indices() {
        let setup = TestSetup::new();
        // PASS is always index 0
        // Filters come first, then INFO, then FORMAT in the string map
        assert_eq!(setup.contig.tid(), 0);
        assert_eq!(setup.contig.name(), "chr1");

        // Verify dict indices match what the header's string_map would assign
        let map = setup.header.string_map();
        assert_eq!(
            map.get("DP").unwrap() as u32,
            setup.dp_info.id().dict_idx(),
            "INFO DP dict_idx mismatch"
        );
        assert_eq!(
            map.get("BQ").unwrap() as u32,
            setup.bq_info.id().dict_idx(),
            "INFO BQ dict_idx mismatch"
        );
        assert_eq!(
            map.get("DB").unwrap() as u32,
            setup.db_flag.id().dict_idx(),
            "INFO DB dict_idx mismatch"
        );
        assert_eq!(
            map.get("GT").unwrap() as u32,
            setup.gt_fmt.id().dict_idx(),
            "FORMAT GT dict_idx mismatch"
        );
    }

    // r[verify record_encoder.trait]
    #[test]
    fn record_encoder_is_object_safe() {
        // Verify RecordEncoder can be used as a trait object
        fn accepts_dyn(_enc: &mut dyn RecordEncoder) {}
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let mut writer = BcfWriter::new(&mut buf, setup.header.clone(), false);
        writer.write_header().unwrap();
        let mut enc = writer.record_encoder();
        accepts_dyn(&mut enc);
    }

    // r[verify record_encoder.field_def]
    // r[verify record_encoder.field_def_types]
    // r[verify record_encoder.encode_info_trait]
    // r[verify record_encoder.encode_format_trait]
    // r[verify record_encoder.encode_dyn]
    #[test]
    fn custom_encode_traits_work() {
        // Simple custom type that encodes itself via EncodeInfo
        struct Depth(i32);
        impl EncodeInfo for Depth {
            type Key = InfoInt;
            fn encode_info(&self, enc: &mut dyn RecordEncoder, key: &Self::Key) {
                key.encode(enc, self.0);
            }
        }

        // Custom type for FORMAT
        struct Score(f32);
        impl EncodeFormat for Score {
            type Key = FormatFloat;
            fn encode_format(&self, enc: &mut dyn RecordEncoder, key: &Self::Key) {
                key.encode(enc, self.0);
            }
        }

        let setup = TestSetup::new();
        let depth = Depth(42);
        let mut buf = Vec::new();
        let mut writer = BcfWriter::new(&mut buf, setup.header.clone(), false);
        writer.write_header().unwrap();
        let mut enc = writer.record_encoder();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        enc.begin(&setup.contig, Pos::<One>::new(1).unwrap(), &alleles, None).unwrap();
        enc.filter_pass();
        depth.encode_info(&mut enc, &setup.dp_info);
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.key_types]
    // r[verify record_encoder.bcf_backwards_compat]
    #[test]
    fn old_handle_api_still_works() {
        use crate::vcf::encoder::{ContigHandle, FilterHandle, ScalarInfoHandle};

        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let mut writer = BcfWriter::new(&mut buf, setup.header.clone(), false);
        writer.write_header().unwrap();

        // Use old handle API directly
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer.record_encoder();
        alleles
            .begin_record(&mut enc, ContigHandle(0), Pos::<One>::new(100).unwrap(), Some(30.0))
            .unwrap();
        FilterHandle::PASS.encode(&mut enc);
        let dp_handle =
            ScalarInfoHandle::<i32> { dict_idx: setup.dp_info.id().dict_idx, _marker: PhantomData };
        dp_handle.encode(&mut enc, 50);
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::bcf_writer::BcfWriter;
    use crate::vcf::header::ContigDef;
    use crate::vcf::record::{SampleValue, VcfRecordBuilder};
    use proptest::prelude::*;
    use seqair_types::Base;
    use std::sync::Arc;

    fn test_header_and_keys() -> (
        Arc<crate::vcf::header::VcfHeader>,
        ContigId,
        InfoInt,
        InfoFloat,
        InfoFlag,
        FormatGt,
        FormatInt,
    ) {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let dp = builder
            .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
            .unwrap();
        let bq = builder
            .register_info(&InfoFieldDef::new("BQ", Number::Count(1), ValueType::Float, "BQ"))
            .unwrap();
        let db = builder
            .register_info(&InfoFieldDef::new("DB", Number::Count(0), ValueType::Flag, "dbSNP"))
            .unwrap();
        let gt = builder
            .register_format(&FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "GT"))
            .unwrap();
        let dp_fmt = builder
            .register_format(&FormatFieldDef::new("DP", Number::Count(1), ValueType::Integer, "DP"))
            .unwrap();
        let header = Arc::new(builder.add_sample("S1").unwrap().build().unwrap());
        (header, contig, dp, bq, db, gt, dp_fmt)
    }

    // r[verify record_encoder.record_path_equivalence]
    proptest! {
        #[test]
        fn bcf_record_encoder_matches_vcf_record_builder(
            pos in 1u32..1_000_000,
            qual in proptest::option::of(0.0f32..1000.0),
            dp_val in 0i32..10000,
            bq_val in 0.0f32..100.0,
            has_flag in proptest::bool::ANY,
            gt_a0 in 0u16..2,
            gt_a1 in 0u16..2,
            fmt_dp in 0i32..10000,
            alt_base_idx in 1u8..4, // 1=C, 2=G, 3=T (skip 0=A since ref is A)
        ) {
            let (header, contig, dp_key, bq_key, db_key, gt_key, dp_fmt_key) =
                test_header_and_keys();

            let alt_base = match alt_base_idx {
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            };
            let alleles = Alleles::snv(Base::A, alt_base).unwrap();
            let pos1 = Pos::<One>::new(pos).unwrap();
            let gt = Genotype::unphased(gt_a0, gt_a1);

            // Path 1: RecordEncoder
            let encoder_bytes = {
                let mut buf = Vec::new();
                let mut writer = BcfWriter::new(&mut buf, header.clone(), false);
                writer.write_header().unwrap();
                let mut enc = writer.record_encoder();
                enc.begin(&contig, pos1, &alleles, qual).unwrap();
                enc.filter_pass();
                dp_key.encode(&mut enc, dp_val);
                bq_key.encode(&mut enc, bq_val);
                if has_flag {
                    db_key.encode(&mut enc);
                }
                enc.begin_samples(1);
                gt_key.encode(&mut enc, &gt);
                dp_fmt_key.encode(&mut enc, fmt_dp);
                enc.emit().unwrap();
                writer.finish().unwrap();
                buf
            };

            // Path 2: VcfRecordBuilder
            let record_bytes = {
                let mut buf = Vec::new();
                let mut writer = BcfWriter::new(&mut buf, header.clone(), false);
                writer.write_header().unwrap();

                let mut builder = VcfRecordBuilder::new("chr1", pos1, alleles)
                    .filter_pass()
                    .info_integer("DP", dp_val)
                    .info_float("BQ", bq_val);
                if let Some(q) = qual {
                    builder = builder.qual(q);
                }
                if has_flag {
                    builder = builder.info_flag("DB");
                }
                let record = builder
                    .format_keys(&["GT", "DP"])
                    .add_sample(vec![
                        SampleValue::Genotype(gt),
                        SampleValue::Integer(fmt_dp),
                    ])
                    .build(&header)
                    .unwrap();

                writer.write_vcf_record(&record).unwrap();
                writer.finish().unwrap();
                buf
            };

            prop_assert_eq!(
                encoder_bytes,
                record_bytes,
                "RecordEncoder and VcfRecordBuilder BCF output must be identical"
            );
        }

        // r[verify record_encoder.vcf_text_encoder]
        // r[verify record_encoder.vcf_text_begin]
        // r[verify record_encoder.vcf_text_info]
        // r[verify record_encoder.vcf_text_format_accumulation]
        // r[verify record_encoder.vcf_text_emit]
        // r[verify record_encoder.vcf_text_buffer_reuse]
        // r[verify record_encoder.vcf_text_output]
        // r[verify record_encoder.record_path_equivalence]
        #[test]
        fn vcf_text_record_encoder_matches_vcf_record_builder(
            pos in 1u32..1_000_000,
            qual in proptest::option::of(0.0f32..1000.0),
            dp_val in 0i32..10000,
            bq_val in 0.0f32..100.0,
            has_flag in proptest::bool::ANY,
            gt_a0 in 0u16..2,
            gt_a1 in 0u16..2,
            fmt_dp in 0i32..10000,
            alt_base_idx in 1u8..4,
        ) {
            use crate::vcf::writer::VcfWriter;

            let (header, contig, dp_key, bq_key, db_key, gt_key, dp_fmt_key) =
                test_header_and_keys();

            let alt_base = match alt_base_idx {
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            };
            let alleles = Alleles::snv(Base::A, alt_base).unwrap();
            let pos1 = Pos::<One>::new(pos).unwrap();
            let gt = Genotype::unphased(gt_a0, gt_a1);

            // Path 1: VcfRecordEncoder
            let encoder_text = {
                let mut buf = Vec::new();
                let mut writer = VcfWriter::new(&mut buf, header.clone());
                writer.write_header().unwrap();
                let mut enc = writer.record_encoder();
                enc.begin(&contig, pos1, &alleles, qual).unwrap();
                enc.filter_pass();
                dp_key.encode(&mut enc, dp_val);
                bq_key.encode(&mut enc, bq_val);
                if has_flag {
                    db_key.encode(&mut enc);
                }
                enc.begin_samples(1);
                gt_key.encode(&mut enc, &gt);
                dp_fmt_key.encode(&mut enc, fmt_dp);
                enc.emit().unwrap();
                drop(enc);
                writer.finish().unwrap();
                String::from_utf8(buf).unwrap()
            };

            // Path 2: VcfRecordBuilder
            let record_text = {
                let mut buf = Vec::new();
                let mut writer = VcfWriter::new(&mut buf, header.clone());
                writer.write_header().unwrap();

                let mut builder = VcfRecordBuilder::new("chr1", pos1, alleles)
                    .filter_pass()
                    .info_integer("DP", dp_val)
                    .info_float("BQ", bq_val);
                if let Some(q) = qual {
                    builder = builder.qual(q);
                }
                if has_flag {
                    builder = builder.info_flag("DB");
                }
                let record = builder
                    .format_keys(&["GT", "DP"])
                    .add_sample(vec![
                        SampleValue::Genotype(gt),
                        SampleValue::Integer(fmt_dp),
                    ])
                    .build(&header)
                    .unwrap();

                writer.write_record(&record).unwrap();
                writer.finish().unwrap();
                String::from_utf8(buf).unwrap()
            };

            // Compare only data lines (skip header)
            let encoder_data: Vec<&str> = encoder_text.lines()
                .filter(|l| !l.starts_with('#'))
                .collect();
            let record_data: Vec<&str> = record_text.lines()
                .filter(|l| !l.starts_with('#'))
                .collect();
            prop_assert_eq!(
                encoder_data,
                record_data,
                "VcfRecordEncoder and VcfRecordBuilder VCF text output must be identical"
            );
        }
    }
}

#[cfg(test)]
mod cross_format_tests {
    use super::*;
    use crate::vcf::bcf_writer::BcfWriter;
    use crate::vcf::header::ContigDef;
    use crate::vcf::record::Genotype;
    use crate::vcf::writer::VcfWriter;
    use proptest::prelude::*;
    use seqair_types::Base;
    use std::sync::Arc;

    fn bcftools_available() -> bool {
        std::process::Command::new("bcftools")
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    fn setup() -> (
        Arc<crate::vcf::header::VcfHeader>,
        ContigId,
        InfoInt,
        InfoFloat,
        InfoFlag,
        InfoInts,
        FormatGt,
        FormatInt,
    ) {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let dp = builder
            .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
            .unwrap();
        let bq = builder
            .register_info(&InfoFieldDef::new("BQ", Number::Count(1), ValueType::Float, "BQ"))
            .unwrap();
        let db = builder
            .register_info(&InfoFieldDef::new("DB", Number::Count(0), ValueType::Flag, "dbSNP"))
            .unwrap();
        let ad = builder
            .register_info(&InfoFieldDef::new(
                "AD",
                Number::ReferenceAlternateBases,
                ValueType::Integer,
                "Allele depth",
            ))
            .unwrap();
        let gt = builder
            .register_format(&FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "GT"))
            .unwrap();
        let dp_fmt = builder
            .register_format(&FormatFieldDef::new(
                "DP",
                Number::Count(1),
                ValueType::Integer,
                "DP",
            ))
            .unwrap();
        let header = Arc::new(builder.add_sample("S1").unwrap().build().unwrap());
        (header, contig, dp, bq, db, ad, gt, dp_fmt)
    }

    /// Encode a record via RecordEncoder into BCF, convert to VCF text via bcftools,
    /// and compare to direct VCF text encoding.
    // r[verify record_encoder.vcf_bcf_equivalence]
    #[test]
    fn bcf_and_vcf_produce_equivalent_records_via_bcftools() {
        if !bcftools_available() {
            eprintln!("skipping: bcftools not found");
            return;
        }

        let (header, contig, dp_key, bq_key, db_key, ad_key, gt_key, dp_fmt_key) = setup();

        let pos = seqair_types::Pos::<seqair_types::One>::new(42).unwrap();
        let alleles = crate::vcf::alleles::Alleles::snv(Base::A, Base::T).unwrap();
        let gt = Genotype::unphased(0, 1);

        // Write BCF via RecordEncoder
        let bcf_path = std::env::temp_dir().join("record_encoder_test.bcf");
        {
            let file = std::fs::File::create(&bcf_path).unwrap();
            let mut writer = BcfWriter::new(file, header.clone(), false);
            writer.write_header().unwrap();
            let mut enc = writer.record_encoder();
            enc.begin(&contig, pos, &alleles, Some(30.0)).unwrap();
            enc.filter_pass();
            dp_key.encode(&mut enc, 50);
            bq_key.encode(&mut enc, 35.5);
            db_key.encode(&mut enc);
            ad_key.encode(&mut enc, &[30, 20]);
            enc.begin_samples(1);
            gt_key.encode(&mut enc, &gt);
            dp_fmt_key.encode(&mut enc, 45);
            enc.emit().unwrap();
            writer.finish().unwrap();
        }

        // Write VCF text via RecordEncoder
        let vcf_path = std::env::temp_dir().join("record_encoder_test.vcf");
        {
            let file = std::fs::File::create(&vcf_path).unwrap();
            let mut writer = VcfWriter::new(file, header.clone());
            writer.write_header().unwrap();
            let mut enc = writer.record_encoder();
            enc.begin(&contig, pos, &alleles, Some(30.0)).unwrap();
            enc.filter_pass();
            dp_key.encode(&mut enc, 50);
            bq_key.encode(&mut enc, 35.5);
            db_key.encode(&mut enc);
            ad_key.encode(&mut enc, &[30, 20]);
            enc.begin_samples(1);
            gt_key.encode(&mut enc, &gt);
            dp_fmt_key.encode(&mut enc, 45);
            enc.emit().unwrap();
            drop(enc);
            writer.finish().unwrap();
        }

        // Convert BCF → VCF via bcftools
        let bcftools_output = std::process::Command::new("bcftools")
            .args(["view", "-H", bcf_path.to_str().unwrap()])
            .output()
            .expect("failed to run bcftools");
        assert!(
            bcftools_output.status.success(),
            "bcftools view failed: {}",
            String::from_utf8_lossy(&bcftools_output.stderr)
        );
        let bcf_as_vcf = String::from_utf8(bcftools_output.stdout).unwrap();

        // Read VCF data lines
        let vcf_text = std::fs::read_to_string(&vcf_path).unwrap();
        let vcf_data: Vec<&str> = vcf_text.lines().filter(|l| !l.starts_with('#')).collect();
        let bcf_data: Vec<&str> = bcf_as_vcf.lines().filter(|l| !l.is_empty()).collect();

        assert_eq!(
            vcf_data, bcf_data,
            "VCF text and BCF (via bcftools) data lines must match"
        );

        // Cleanup
        let _ = std::fs::remove_file(&bcf_path);
        let _ = std::fs::remove_file(&vcf_path);
    }

    // r[verify record_encoder.vcf_bcf_equivalence]
    proptest! {
        #[test]
        fn bcf_vcf_equivalence_proptest(
            pos in 1u32..100_000,
            dp_val in 0i32..10000,
            bq_val in 0.0f32..100.0,
            has_flag in proptest::bool::ANY,
            gt_a0 in 0u16..2,
            gt_a1 in 0u16..2,
            fmt_dp in 0i32..10000,
            alt_base_idx in 1u8..4,
        ) {
            if !bcftools_available() {
                return Ok(());
            }

            let (header, contig, dp_key, bq_key, db_key, _ad_key, gt_key, dp_fmt_key) = setup();

            let alt_base = match alt_base_idx {
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            };
            let alleles = crate::vcf::alleles::Alleles::snv(Base::A, alt_base).unwrap();
            let pos1 = seqair_types::Pos::<seqair_types::One>::new(pos).unwrap();
            let gt = Genotype::unphased(gt_a0, gt_a1);

            // BCF via RecordEncoder → temp file → bcftools view
            let bcf_path = std::env::temp_dir().join(format!("proptest_bcf_{pos}.bcf"));
            {
                let file = std::fs::File::create(&bcf_path).unwrap();
                let mut writer = BcfWriter::new(file, header.clone(), false);
                writer.write_header().unwrap();
                let mut enc = writer.record_encoder();
                enc.begin(&contig, pos1, &alleles, Some(30.0)).unwrap();
                enc.filter_pass();
                dp_key.encode(&mut enc, dp_val);
                bq_key.encode(&mut enc, bq_val);
                if has_flag { db_key.encode(&mut enc); }
                enc.begin_samples(1);
                gt_key.encode(&mut enc, &gt);
                dp_fmt_key.encode(&mut enc, fmt_dp);
                enc.emit().unwrap();
                writer.finish().unwrap();
            }

            // VCF text via RecordEncoder → in-memory
            let vcf_text = {
                let mut buf = Vec::new();
                let mut writer = VcfWriter::new(&mut buf, header.clone());
                writer.write_header().unwrap();
                let mut enc = writer.record_encoder();
                enc.begin(&contig, pos1, &alleles, Some(30.0)).unwrap();
                enc.filter_pass();
                dp_key.encode(&mut enc, dp_val);
                bq_key.encode(&mut enc, bq_val);
                if has_flag { db_key.encode(&mut enc); }
                enc.begin_samples(1);
                gt_key.encode(&mut enc, &gt);
                dp_fmt_key.encode(&mut enc, fmt_dp);
                enc.emit().unwrap();
                drop(enc);
                writer.finish().unwrap();
                String::from_utf8(buf).unwrap()
            };

            let bcftools_output = std::process::Command::new("bcftools")
                .args(["view", "-H", bcf_path.to_str().unwrap()])
                .output()
                .expect("failed to run bcftools");
            let bcf_as_vcf = String::from_utf8(bcftools_output.stdout).unwrap();

            let vcf_data: Vec<&str> = vcf_text.lines().filter(|l| !l.starts_with('#')).collect();
            let bcf_data: Vec<&str> = bcf_as_vcf.lines().filter(|l| !l.is_empty()).collect();

            let _ = std::fs::remove_file(&bcf_path);

            prop_assert_eq!(vcf_data, bcf_data, "BCF↔VCF equivalence failed");
        }
    }
}
