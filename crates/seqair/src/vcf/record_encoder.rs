//! Format-agnostic record encoder with typed field keys.
//!
//! Typed field keys ([`InfoKey`], [`FormatKey`]) enforce value types at compile
//! time. The [`InfoEncoder`] and [`FormatEncoder`] traits provide object-safe
//! APIs for encoding INFO and FORMAT fields into either BCF binary or VCF text.
//!
//! # Field definitions
//!
//! [`InfoFieldDef`] and [`FormatFieldDef`] combine header metadata with
//! type-safe key resolution. Register them on [`VcfHeaderBuilder`](super::header::VcfHeaderBuilder)
//! to build the header and resolve keys in one step:
//!
//! ```
//! use seqair::vcf::{Number, ValueType, VcfHeader, InfoInt, FormatGt};
//! use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut builder = VcfHeader::builder();
//!
//! // Returns InfoKey<Scalar<i32>> (aliased InfoInt) — can only encode i32 values
//! let dp_key: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Depth")
//! )?;
//!
//! // Returns FormatKey<Gt> (aliased FormatGt) — can only encode Genotype slices
//! let gt_key: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//!
//! assert_eq!(dp_key.id().name(), "DP");
//! assert_eq!(gt_key.id().name(), "GT");
//! # Ok(())
//! # }
//! ```

use super::header::{Number, ValueType};
use super::record::Genotype;
use seqair_types::SmolStr;
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
#[derive(Debug, Clone)]
pub struct FieldId {
    pub(crate) dict_idx: u32,
    pub(crate) name: SmolStr,
}

impl FieldId {
    pub fn name(&self) -> &str {
        &self.name
    }
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
    pub fn name(&self) -> &str {
        &self.name
    }
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
    pub const PASS: FilterId = Self { dict_idx: 0, name: SmolStr::new_inline("PASS") };
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn dict_idx(&self) -> u32 {
        self.dict_idx
    }
}

// ── Typed field keys ───────────────────────────────────────────────────

// r[impl record_encoder.key_types]

#[derive(Debug, Clone)]
pub struct InfoKey<V>(pub(crate) FieldId, pub(crate) PhantomData<V>);

#[derive(Debug, Clone)]
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

pub type FormatGt = FormatKey<Gt>;
pub type FormatInt = FormatKey<Scalar<i32>>;
pub type FormatFloat = FormatKey<Scalar<f32>>;

// ── InfoEncoder trait ──────────────────────────────────────────────────

// r[impl record_encoder.info_methods]
/// Object-safe trait for encoding INFO fields.
pub trait InfoEncoder {
    fn info_int(&mut self, id: &FieldId, value: i32);
    fn info_float(&mut self, id: &FieldId, value: f32);
    fn info_ints(&mut self, id: &FieldId, values: &[i32]);
    fn info_floats(&mut self, id: &FieldId, values: &[f32]);
    fn info_flag(&mut self, id: &FieldId);
    fn info_string(&mut self, id: &FieldId, value: &str);
    fn info_int_opts(&mut self, id: &FieldId, values: &[Option<i32>]);
    fn n_allele(&self) -> usize;
    fn n_alt(&self) -> usize;
}

// ── FormatEncoder trait ────────────────────────────────────────────────

// r[impl record_encoder.format_methods]
/// Object-safe trait for encoding FORMAT fields.
pub trait FormatEncoder {
    /// Encode a GT FORMAT field for all samples. Slice length MUST equal sample count.
    fn format_gt(&mut self, id: &FieldId, gts: &[Genotype]);
    /// Encode a scalar integer FORMAT field for all samples.
    fn format_int(&mut self, id: &FieldId, values: &[i32]);
    /// Encode a scalar float FORMAT field for all samples.
    fn format_float(&mut self, id: &FieldId, values: &[f32]);
    fn n_allele(&self) -> usize;
    fn n_alt(&self) -> usize;
    /// Number of samples declared by the header.
    fn n_samples(&self) -> usize;
}

// ── Key encode methods ─────────────────────────────────────────────────

// r[impl record_encoder.key_encode]

impl InfoKey<Scalar<i32>> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), value: i32) {
        enc.info_int(&self.0, value);
    }
}

impl InfoKey<Scalar<f32>> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), value: f32) {
        enc.info_float(&self.0, value);
    }
}

impl InfoKey<Arr<i32>> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), values: &[i32]) {
        enc.info_ints(&self.0, values);
    }
}

impl InfoKey<Arr<f32>> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), values: &[f32]) {
        enc.info_floats(&self.0, values);
    }
}

impl InfoKey<Flag> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized)) {
        enc.info_flag(&self.0);
    }
}

impl InfoKey<Str> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), value: &str) {
        enc.info_string(&self.0, value);
    }
}

impl InfoKey<OptArr<i32>> {
    pub fn encode(&self, enc: &mut (impl InfoEncoder + ?Sized), values: &[Option<i32>]) {
        enc.info_int_opts(&self.0, values);
    }
}

impl FormatKey<Gt> {
    /// Encode a GT FORMAT field. `gts` must have exactly one entry per sample.
    ///
    /// # Panics
    /// Panics if `gts.len() != n_samples` (programming error — the caller controls the slice).
    pub fn encode(&self, enc: &mut (impl FormatEncoder + ?Sized), gts: &[Genotype]) {
        assert_eq!(
            gts.len(),
            enc.n_samples(),
            "FormatKey<Gt>::encode: expected {} genotypes (one per sample), got {}",
            enc.n_samples(),
            gts.len(),
        );
        enc.format_gt(&self.0, gts);
    }
}

impl FormatKey<Scalar<i32>> {
    /// Encode a scalar integer FORMAT field. `values` must have exactly one entry per sample.
    ///
    /// # Panics
    /// Panics if `values.len() != n_samples` (programming error — the caller controls the slice).
    pub fn encode(&self, enc: &mut (impl FormatEncoder + ?Sized), values: &[i32]) {
        assert_eq!(
            values.len(),
            enc.n_samples(),
            "FormatKey<Scalar<i32>>::encode: expected {} values (one per sample), got {}",
            enc.n_samples(),
            values.len(),
        );
        enc.format_int(&self.0, values);
    }
}

impl FormatKey<Scalar<f32>> {
    /// Encode a scalar float FORMAT field. `values` must have exactly one entry per sample.
    ///
    /// # Panics
    /// Panics if `values.len() != n_samples` (programming error — the caller controls the slice).
    pub fn encode(&self, enc: &mut (impl FormatEncoder + ?Sized), values: &[f32]) {
        assert_eq!(
            values.len(),
            enc.n_samples(),
            "FormatKey<Scalar<f32>>::encode: expected {} values (one per sample), got {}",
            enc.n_samples(),
            values.len(),
        );
        enc.format_float(&self.0, values);
    }
}

// ── Field definitions ──────────────────────────────────────────────────

// r[impl record_encoder.field_def]
// r[impl record_encoder.field_def_types]

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

// ── Custom type encoding traits ────────────────────────────────────────

// r[impl record_encoder.encode_info_trait]
// r[impl record_encoder.encode_dyn]
pub trait EncodeInfo {
    type Key;
    fn encode_info(&self, enc: &mut dyn InfoEncoder, key: &Self::Key);
}

// r[impl record_encoder.encode_format_trait]
pub trait EncodeFormat {
    type Key;
    fn encode_format(&self, enc: &mut dyn FormatEncoder, key: &Self::Key);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::OutputFormat;
    use crate::vcf::alleles::Alleles;
    use crate::vcf::header::ContigDef;
    use crate::vcf::record::Genotype;
    use crate::vcf::unified::Writer;
    use seqair_types::{Base, One, Pos};
    use std::sync::Arc;

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
        let pass = FilterId::PASS;
        assert_eq!(pass.dict_idx, 0);
        assert_eq!(pass.name(), "PASS");
    }

    // r[verify record_encoder.typed_keys]
    #[test]
    fn typed_keys_prevent_misuse() {
        let _int_key: InfoInt =
            InfoKey(FieldId { dict_idx: 0, name: SmolStr::from("DP") }, PhantomData);
        let _float_key: InfoFloat =
            InfoKey(FieldId { dict_idx: 1, name: SmolStr::from("BQ") }, PhantomData);
        let _flag_key: InfoFlag =
            InfoKey(FieldId { dict_idx: 2, name: SmolStr::from("CPG") }, PhantomData);
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

    // r[verify record_encoder.register]
    // r[verify record_encoder.register_contig]
    // r[verify record_encoder.field_id]
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
    }

    // r[verify record_encoder.encode_info_trait]
    // r[verify record_encoder.encode_dyn]
    // r[verify record_encoder.field_def]
    // r[verify record_encoder.field_def_types]
    #[test]
    fn typestate_custom_encode_traits() {
        struct Depth(i32);
        impl EncodeInfo for Depth {
            type Key = InfoInt;
            fn encode_info(&self, enc: &mut dyn InfoEncoder, key: &Self::Key) {
                key.encode(enc, self.0);
            }
        }

        #[allow(dead_code, reason = "only exists to verify EncodeFormat impl compiles")]
        struct Score(f32);
        impl EncodeFormat for Score {
            type Key = FormatFloat;
            fn encode_format(&self, enc: &mut dyn FormatEncoder, key: &Self::Key) {
                key.encode(enc, &[self.0]);
            }
        }

        let setup = TestSetup::new();
        let depth = Depth(42);
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        depth.encode_info(&mut enc, &setup.dp_info);
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.filters]
    #[test]
    fn filter_fail_bcf_and_vcf_text() {
        let mut builder = crate::vcf::header::VcfHeader::builder();
        let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let low_dp = builder.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let header = Arc::new(builder.build().unwrap());

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let pos = Pos::<One>::new(10).unwrap();

        let bcf_bytes = {
            let mut buf = Vec::new();
            let writer = Writer::new(&mut buf, OutputFormat::Bcf);
            let mut writer = writer.write_header(&header).unwrap();
            writer
                .begin_record(&contig, pos, &alleles, None)
                .unwrap()
                .filter_fail(&[&low_dp])
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
                .filter_fail(&[&low_dp])
                .emit()
                .unwrap();
            writer.finish().unwrap();
            String::from_utf8(buf).unwrap()
        };
        let data_line = vcf_text.lines().find(|l| !l.starts_with('#')).unwrap();
        let filter_col = data_line.split('\t').nth(6).unwrap();
        assert_eq!(filter_col, "lowDp", "filter_fail must write the filter name");
    }

    // r[verify record_encoder.info_methods]
    // r[verify record_encoder.format_methods]
    // r[verify record_encoder.key_encode]
    // r[verify record_encoder.bcf_impl]
    #[test]
    fn typestate_api_bcf() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50);
        setup.bq_info.encode(&mut enc, 35.5);
        setup.db_flag.encode(&mut enc);
        setup.ad_info.encode(&mut enc, &[30, 20]);
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
        setup.dp_fmt.encode(&mut enc, &[45]);
        enc.emit().unwrap();
        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.vcf_text_encoder]
    #[test]
    fn typestate_emit_without_samples() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50);
        enc.emit().unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let cols: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(cols.first(), Some(&"chr1"), "CHROM");
        assert_eq!(cols.get(1), Some(&"100"), "POS");
        assert_eq!(cols.get(6), Some(&"PASS"), "FILTER");
    }

    // r[verify record_encoder.filters]
    #[test]
    fn typestate_filter_fail() {
        let mut hb = crate::vcf::header::VcfHeader::builder();
        let contig = hb.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let ld = hb.register_filter(&FilterFieldDef::new("lowDp", "Low depth")).unwrap();
        let header = Arc::new(hb.build().unwrap());
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        writer
            .begin_record(&contig, Pos::<One>::new(10).unwrap(), &alleles, None)
            .unwrap()
            .filter_fail(&[&ld])
            .emit()
            .unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        let filter_col = data_line.split('\t').nth(6).unwrap();
        assert_eq!(filter_col, "lowDp");
    }

    // r[verify record_encoder.vcf_text_encoder]
    #[test]
    fn typestate_no_filter() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        writer
            .begin_record(&setup.contig, Pos::<One>::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass()
            .emit()
            .unwrap();
        writer.finish().unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(data_line.split('\t').nth(6), Some("PASS"));
    }

    // r[verify record_encoder.encode_info_trait]
    // r[verify record_encoder.encode_format_trait]
    #[test]
    fn typestate_info_and_format_encoder_are_object_safe() {
        fn accepts_info_enc(_enc: &mut dyn InfoEncoder) {}
        fn accepts_format_enc(_enc: &mut dyn FormatEncoder) {}

        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut filtered = writer
            .begin_record(&setup.contig, Pos::<One>::new(1).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        accepts_info_enc(&mut filtered);
        let mut with_samples = filtered.begin_samples();
        accepts_format_enc(&mut with_samples);
        with_samples.emit().unwrap();
        writer.finish().unwrap();
    }

    // r[verify record_encoder.bcf_impl]
    // r[verify record_encoder.vcf_text_encoder]
    #[test]
    fn typestate_multiple_records() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        for pos in [1u32, 100, 200] {
            writer
                .begin_record(&setup.contig, Pos::<One>::new(pos).unwrap(), &alleles, None)
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
        let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
        let dp_info = builder
            .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
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
        let header = Arc::new(
            builder
                .add_sample("S1")
                .unwrap()
                .add_sample("S2")
                .unwrap()
                .add_sample("S3")
                .unwrap()
                .build()
                .unwrap(),
        );
        (header, contig, dp_info, gt_fmt, dp_fmt)
    }

    // r[verify record_encoder.format_methods]
    #[test]
    fn multi_sample_bcf() {
        let (header, contig, dp_info, gt_fmt, dp_fmt) = multi_sample_setup();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap();
        let mut enc = enc.filter_pass();
        dp_info.encode(&mut enc, 150);
        let mut enc = enc.begin_samples();
        gt_fmt.encode(
            &mut enc,
            &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
        );
        dp_fmt.encode(&mut enc, &[45, 52, 38]);
        enc.emit().unwrap();

        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    // r[verify record_encoder.format_methods]
    #[test]
    fn multi_sample_vcf_text() {
        let (header, contig, dp_info, gt_fmt, dp_fmt) = multi_sample_setup();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap();
        let mut enc = enc.filter_pass();
        dp_info.encode(&mut enc, 150);
        let mut enc = enc.begin_samples();
        gt_fmt.encode(
            &mut enc,
            &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
        );
        dp_fmt.encode(&mut enc, &[45, 52, 38]);
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

    // r[verify record_encoder.format_methods]
    #[test]
    fn single_sample_slice_api() {
        let setup = TestSetup::new();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();

        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap();
        let mut enc = enc.filter_pass();
        setup.dp_info.encode(&mut enc, 50);
        let mut enc = enc.begin_samples();
        // Single-sample callers pass 1-element slices
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
        setup.dp_fmt.encode(&mut enc, &[45]);
        enc.emit().unwrap();

        writer.finish().unwrap();
        assert!(!buf.is_empty());
    }

    #[test]
    #[should_panic(expected = "expected 1 values (one per sample), got 2")]
    fn format_int_panics_on_sample_count_mismatch() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Pass 2 values for 1-sample header — must panic
        setup.dp_fmt.encode(&mut enc, &[45, 52]);
    }

    #[test]
    #[should_panic(expected = "expected 1 genotypes (one per sample), got 3")]
    fn format_gt_panics_on_sample_count_mismatch() {
        let setup = TestSetup::new(); // 1 sample
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&setup.contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Pass 3 genotypes for 1-sample header — must panic
        setup.gt_fmt.encode(
            &mut enc,
            &[Genotype::unphased(0, 1), Genotype::unphased(0, 0), Genotype::unphased(1, 1)],
        );
    }

    #[test]
    #[should_panic(expected = "mixed ploidy not supported")]
    fn format_gt_panics_on_mixed_ploidy() {
        let (header, contig, _dp_info, gt_fmt, _dp_fmt) = multi_sample_setup();
        let mut buf = Vec::new();
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&contig, Pos::<One>::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        // Haploid + diploid + diploid — must panic
        gt_fmt.encode(
            &mut enc,
            &[Genotype::haploid(0), Genotype::unphased(0, 1), Genotype::unphased(1, 1)],
        );
    }
}
