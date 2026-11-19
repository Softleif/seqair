//! Typed access (`r[vcf_values.access]`): what each key marker reads a stored
//! value as, and the untyped [`DynValue`] for generic tools.

use super::bcf_value::{TypedValue, gt_allele};
use super::error::FieldError;
use super::file::FileShared;
use super::text_value::{TextGt, TextNumber, TextVec, text_str};
use super::values::{GenotypeRef, GtAlleles, OptValues, Values};
use crate::vcf::header::ValueType;
use crate::vcf::record_encoder::{
    Arr, FieldKind, Flag, FormatKey, Gt, InfoKey, OptArr, Scalar, Str,
};

/// Most alleles a genotype may have (`r[vcf_reader.limits+2]`).
pub(crate) const MAX_PLOIDY: usize = 256;

/// One field's value as a record stores it, not yet decoded.
#[derive(Debug, Clone, Copy)]
pub struct Stored<'a>(pub(crate) StoredRepr<'a>);

/// Which format the value is in, so each marker picks its decoder.
#[derive(Debug, Clone, Copy)]
pub(crate) enum StoredRepr<'a> {
    Bcf(TypedValue<'a>),
    /// A VCF text value: `None` for an INFO key without `=`, which holds
    /// nothing (BCF's empty vector); an empty value is `.` (missing).
    Text(Option<&'a [u8]>),
}

/// The text of a value as it parses: an empty value is `.`.
fn text_value(value: &[u8]) -> &[u8] {
    if value.is_empty() { b"." } else { value }
}

fn text_vec<'a, T: TextNumber>(
    value: Option<&'a [u8]>,
    cx: DecodeCx<'_>,
) -> Result<TextVec<'a, T>, FieldError> {
    let Some(v) = value else { return Ok(TextVec::empty()) };
    let v = TextVec::new(text_value(v))?;
    // r[impl vcf_values.text_numbers+2]
    if v.has_reserved()
        && let Some(file) = cx.file
    {
        file.warn_reserved_int(cx.name);
    }
    Ok(v)
}

/// How to decode: what a `Gt` needs from the file, and what an [`Any`]
/// needs from the field's definition.
#[derive(Debug, Clone, Copy)]
pub struct DecodeCx<'f> {
    /// Whether a genotype's first allele carries its own phase bit (VCF
    /// 4.4+, `r[vcf_record.gt_first_phase+3]`), from the file's version.
    pub(crate) first_phase_is_read: bool,
    /// The field's declared `Type`, from its key ([`FieldId`](crate::vcf::record_encoder::FieldId))
    /// or the field table: what [`Any`] decodes as.
    pub(crate) declared: ValueType,
    /// Whether the field is FORMAT `GT`, which [`Any`] decodes as a genotype.
    pub(crate) genotype: bool,
    /// The file to warn of a reserved integer read as missing, and the field's name.
    pub(crate) file: Option<&'f FileShared>,
    pub(crate) name: &'f str,
}

mod sealed {
    pub trait Sealed {}
}

// r[impl vcf_values.access+2]
/// A key marker a reader can decode, and what it decodes to:
///
/// | marker | `Out` | absent or `.` |
/// |---|---|---|
/// | `Scalar<i32>`, `Scalar<f32>` | `Option<T>` | `None` |
/// | `Arr<T>` | `Option<Values<'_, T>>` | `None`; a missing element inside is an error |
/// | `OptArr<T>` | `Option<OptValues<'_, T>>` | `None` |
/// | `Flag` | `bool` | `false` |
/// | `Str` | `Option<&str>` | `None` |
/// | `Gt` | `Option<GenotypeRef<'_>>` | `None` |
///
/// A field absent from the record and a field present as `.` both read as
/// absent: VCF lets a sample drop trailing FORMAT fields.
pub trait ReadKind: FieldKind + sealed::Sealed {
    type Out<'a>;

    #[doc(hidden)]
    fn absent<'a>() -> Self::Out<'a>;

    #[doc(hidden)]
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError>;

    /// `true` if every sample's value decodes, as the type of the whole
    /// column (`column`, every sample's elements, `per_sample` each) alone
    /// shows; `false` when each sample must be decoded to tell.
    #[doc(hidden)]
    fn column_decodes(_column: Stored<'_>, _per_sample: usize) -> bool {
        false
    }
}

fn ints<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<OptValues<'a, i32>, FieldError> {
    Ok(match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(v.ints()?),
        StoredRepr::Text(t) => OptValues::text(text_vec(t, cx)?),
    })
}

fn floats<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<OptValues<'a, f32>, FieldError> {
    Ok(match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(v.floats()?),
        StoredRepr::Text(t) => OptValues::text(text_vec(t, cx)?),
    })
}

macro_rules! numeric_kinds {
    ($($t:ty => $read:ident),*) => {$(
        impl sealed::Sealed for Scalar<$t> {}
        impl ReadKind for Scalar<$t> {
            type Out<'a> = Option<$t>;

            fn absent<'a>() -> Self::Out<'a> {
                None
            }

            fn column_decodes(column: Stored<'_>, per_sample: usize) -> bool {
                matches!(column.0, StoredRepr::Bcf(v) if per_sample <= 1 && v.reads_as::<$t>())
            }

            fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value, cx)?;
                match v.len() {
                    0 => Ok(None),
                    1 => Ok(v.get(0).flatten()),
                    len => Err(FieldError::NotScalar { len }),
                }
            }
        }

        impl sealed::Sealed for Arr<$t> {}
        impl ReadKind for Arr<$t> {
            type Out<'a> = Option<Values<'a, $t>>;

            fn absent<'a>() -> Self::Out<'a> {
                None
            }

            fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value, cx)?;
                if v.is_missing() {
                    return Ok(None);
                }
                Ok(Some(v.all_present()?))
            }
        }

        impl sealed::Sealed for OptArr<$t> {}
        impl ReadKind for OptArr<$t> {
            type Out<'a> = Option<OptValues<'a, $t>>;

            fn absent<'a>() -> Self::Out<'a> {
                None
            }

            fn column_decodes(column: Stored<'_>, _: usize) -> bool {
                matches!(column.0, StoredRepr::Bcf(v) if v.reads_as::<$t>())
            }

            fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value, cx)?;
                Ok((!v.is_missing()).then_some(v))
            }
        }
    )*};
}

numeric_kinds!(i32 => ints, f32 => floats);

impl sealed::Sealed for Flag {}
impl ReadKind for Flag {
    type Out<'a> = bool;

    fn absent<'a>() -> Self::Out<'a> {
        false
    }

    fn column_decodes(_: Stored<'_>, _: usize) -> bool {
        true
    }

    /// Present is set, whatever value a writer gave it.
    fn decode<'a>(_: Stored<'a>, _: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        Ok(true)
    }
}

impl sealed::Sealed for Str {}
impl ReadKind for Str {
    type Out<'a> = Option<&'a str>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn decode<'a>(value: Stored<'a>, _: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        let s = match value.0 {
            StoredRepr::Bcf(v) => v.str()?,
            StoredRepr::Text(t) => text_str(t.unwrap_or_default())?,
        };
        Ok((!s.is_empty() && s != ".").then_some(s))
    }
}

impl sealed::Sealed for Gt {}
impl ReadKind for Gt {
    type Out<'a> = Option<GenotypeRef<'a>>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn column_decodes(column: Stored<'_>, per_sample: usize) -> bool {
        matches!(column.0, StoredRepr::Bcf(v)
            if per_sample <= MAX_PLOIDY && v.every_int(|x| gt_allele(x).is_ok()))
    }

    // r[impl vcf_values.gt+2]
    /// No alleles, or a lone missing one, is `.`; a missing diploid `./.`
    /// keeps its ploidy.
    #[inline]
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        let gt = match value.0 {
            StoredRepr::Bcf(v) => {
                check_ploidy(v.count())?;
                let gt = v.gt(cx.first_phase_is_read)?;
                return Ok((!gt.is_missing()).then_some(GenotypeRef(GtAlleles::Bcf(gt))));
            }
            StoredRepr::Text(t) => {
                let t = text_value(t.unwrap_or_default());
                // Every allele but the first follows a separator; so may the first.
                let sep = |b: &u8| *b == b'/' || *b == b'|';
                let seps = t.iter().filter(|b| sep(b)).count();
                check_ploidy(seps.saturating_add(usize::from(!t.first().is_some_and(sep))))?;
                GtAlleles::Text(TextGt::new(t, cx.first_phase_is_read)?)
            }
        };
        let mut alleles = gt.clone();
        Ok(match (alleles.next(), alleles.next()) {
            (None, _) | (Some(super::values::GtAllele { allele: None, .. }), None) => None,
            _ => Some(GenotypeRef(gt)),
        })
    }
}

// r[impl vcf_values.access+2]
/// Marker for a field of any declared type, read as its definition says:
/// `Out` is `Option<DynValue>`, `None` when the record does not hold the
/// field and `Some` (maybe [`DynValue::is_missing`]) when it does. For
/// "present whatever the type" (`info.get(&trid)?.is_some_and(|v| !v.is_missing())`)
/// and for tools that handle fields alike.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Any {
    #[doc(hidden)]
    _Uninhabited(std::convert::Infallible),
}

/// Reads every `Type` and `Number`; declares nothing a writer could encode
/// (an `Any` key has no `encode`).
impl FieldKind for Any {
    const VALUE_TYPE: ValueType = ValueType::String;

    fn reads_type(_: ValueType) -> bool {
        true
    }
}

/// An INFO key of any type.
pub type InfoAny = InfoKey<Any>;
/// A FORMAT key of any type.
pub type FormatAny = FormatKey<Any>;

impl sealed::Sealed for Any {}
impl ReadKind for Any {
    type Out<'a> = Option<DynValue<'a>>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn decode<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        DynValue::decode(value, cx.declared, cx.genotype, cx).map(Some)
    }
}

// r[impl vcf_reader.limits+2]
/// Refuse a genotype of more than [`MAX_PLOIDY`] alleles before reading it.
fn check_ploidy(n: usize) -> Result<(), FieldError> {
    if n > MAX_PLOIDY {
        return Err(FieldError::Ploidy { ploidy: n, max: MAX_PLOIDY });
    }
    Ok(())
}

// r[impl vcf_values.dynamic.missing]
/// `v`, or no elements if it is missing as a whole: VCF text spells an
/// empty vector `.`, so the formats agree only on this.
fn whole<T: super::values::Number>(v: OptValues<'_, T>) -> OptValues<'_, T> {
    if v.is_missing() { OptValues::text(TextVec::empty()) } else { v }
}

// r[impl vcf_values.dynamic+2]
/// A field's value read by its header declaration, for tools that handle
/// every field alike (`convert`, `view`).
#[derive(Debug, Clone)]
pub enum DynValue<'a> {
    /// `Type=Integer`; missing elements are `None`, and a vector missing as
    /// a whole (`.`, BCF's empty vector or lone missing value) is empty.
    Int(OptValues<'a, i32>),
    /// `Type=Float`, as `Int`.
    Float(OptValues<'a, f32>),
    /// `Type=Flag` (INFO only), present.
    Flag,
    /// `Type=String` or `Character`; `None` for `.`.
    Str(Option<&'a str>),
    /// The FORMAT `GT` field; `None` for `.`.
    Gt(Option<GenotypeRef<'a>>),
}

impl<'a> DynValue<'a> {
    // r[impl vcf_values.dynamic+2]
    /// Whether the value is `.`: an empty number vector (VCF text's `.`,
    /// BCF's empty or lone missing vector), a missing string or genotype.
    /// A flag is never missing.
    pub fn is_missing(&self) -> bool {
        match self {
            Self::Int(v) => v.len() == 0,
            Self::Float(v) => v.len() == 0,
            Self::Str(s) => s.is_none(),
            Self::Gt(gt) => gt.is_none(),
            Self::Flag => false,
        }
    }

    /// Decode `value` as `declared` says; `genotype` for the FORMAT GT field.
    pub(crate) fn decode(
        value: Stored<'a>,
        declared: ValueType,
        genotype: bool,
        cx: DecodeCx<'_>,
    ) -> Result<Self, FieldError> {
        if genotype {
            return Ok(Self::Gt(Gt::decode(value, cx)?));
        }
        Ok(match declared {
            ValueType::Integer => Self::Int(whole(ints(value, cx)?)),
            ValueType::Float => Self::Float(whole(floats(value, cx)?)),
            ValueType::Flag => Self::Flag,
            ValueType::String | ValueType::Character => Self::Str(Str::decode(value, cx)?),
        })
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::arithmetic_side_effects, reason = "test sizes below 100")]

    use super::*;
    use crate::vcf::bcf_encoding::encode_type_byte;
    use crate::vcf::reader::bcf_value::FormatColumn;
    use hegel::prelude::*;

    /// Whether `column_decodes` vouches for `column`, checked against
    /// decoding every sample.
    fn vouches<V: ReadKind>(column: &FormatColumn<'_>, n: usize) -> bool {
        let (all, per_sample) = column.all();
        let cx = DecodeCx {
            first_phase_is_read: true,
            declared: ValueType::String,
            genotype: false,
            file: None,
            name: "",
        };
        let vouched = V::column_decodes(Stored(StoredRepr::Bcf(all)), per_sample);
        if vouched {
            for i in 0..n {
                let value = Stored(StoredRepr::Bcf(column.sample(i).unwrap()));
                assert!(V::decode(value, cx).is_ok(), "sample {i}");
            }
        }
        vouched
    }

    // r[verify vcf_fields.samples]
    /// A column the type check vouches for decodes for every sample, on
    /// columns of every type with arbitrary elements (sentinels included).
    #[hegel::test]
    fn a_vouched_column_decodes_for_every_sample(tc: TestCase) {
        let code = tc.draw(gs::sampled_from(vec![0u8, 1, 2, 3, 5, 7]));
        let per_sample = tc.draw(gs::integers::<usize>().max_value(3));
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
        let width = [0, 1, 2, 4, 0, 4, 0, 1][usize::from(code)];
        let mut bytes = Vec::new();
        encode_type_byte(&mut bytes, per_sample, code);
        let payload = tc.draw(
            gs::vecs(gs::sampled_from(vec![0u8, 1, 2, 3, 4, 0x7F, 0x80, 0x81, 0xFF]))
                .min_size(per_sample * width * n)
                .max_size(per_sample * width * n),
        );
        bytes.extend(payload);
        let (column, _) = FormatColumn::parse(&bytes, n).unwrap();
        let all = [
            vouches::<Gt>(&column, n),
            vouches::<Scalar<i32>>(&column, n),
            vouches::<OptArr<i32>>(&column, n),
            vouches::<Scalar<f32>>(&column, n),
            vouches::<OptArr<f32>>(&column, n),
            vouches::<Flag>(&column, n),
        ];
        if code == 1 && per_sample <= 1 {
            assert!(all[1] && all[2]);
        }
    }

    /// Diploid int8 genotypes, haploid ones padded, and missing ones are vouched for.
    #[test]
    fn common_gt_columns_are_vouched_for() {
        let bytes = [0x21, 2, 4, 3, 0x81, 0, 0, 0x80, 0x81];
        let (column, _) = FormatColumn::parse(&bytes, 4).unwrap();
        assert!(vouches::<Gt>(&column, 4));
        let bad = [0x21, 2, 0xFB];
        let (column, _) = FormatColumn::parse(&bad, 1).unwrap();
        assert!(!vouches::<Gt>(&column, 1));
    }
}
