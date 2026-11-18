//! Typed access (`r[vcf_values.access]`): what each key marker reads a stored
//! value as, and the untyped [`DynValue`] for generic tools.

use super::bcf_value::TypedValue;
use super::error::FieldError;
use super::text_value::{TextGt, TextNumber, TextVec, text_str};
use super::values::{GenotypeRef, GtAlleles, OptValues, Values};
use crate::vcf::header::ValueType;
use crate::vcf::record_encoder::{
    Arr, FieldKind, Flag, FormatKey, Gt, InfoKey, OptArr, Scalar, Str,
};

/// Most alleles a genotype may have (`r[vcf_reader.limits]`).
pub(crate) const MAX_PLOIDY: usize = 256;

/// One field's value as a record stores it, not yet decoded.
#[derive(Debug, Clone, Copy)]
pub struct Stored<'a>(pub(crate) StoredRepr<'a>);

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

fn text_vec<T: TextNumber>(value: Option<&[u8]>) -> Result<TextVec<'_, T>, FieldError> {
    Ok(match value {
        None => TextVec::empty(),
        Some(v) => TextVec::new(text_value(v))?,
    })
}

/// How to decode: what a `Gt` needs from the file, and what an [`Any`]
/// needs from the field's definition.
#[derive(Debug, Clone, Copy)]
pub struct DecodeCx {
    pub(crate) first_phase_is_read: bool,
    /// The field's `Type`; looked up only for a marker that reads it.
    pub(crate) declared: ValueType,
    /// Whether the field is FORMAT `GT`.
    pub(crate) genotype: bool,
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

    /// Whether decoding needs the field's definition ([`Any`]).
    #[doc(hidden)]
    const DYNAMIC: bool = false;

    #[doc(hidden)]
    fn absent<'a>() -> Self::Out<'a>;

    #[doc(hidden)]
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx) -> Result<Self::Out<'a>, FieldError>;
}

fn ints(value: Stored<'_>) -> Result<OptValues<'_, i32>, FieldError> {
    Ok(match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(v.ints()?),
        StoredRepr::Text(t) => OptValues::text(text_vec(t)?),
    })
}

fn floats(value: Stored<'_>) -> Result<OptValues<'_, f32>, FieldError> {
    Ok(match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(v.floats()?),
        StoredRepr::Text(t) => OptValues::text(text_vec(t)?),
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

            fn decode<'a>(value: Stored<'a>, _: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value)?;
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

            fn decode<'a>(value: Stored<'a>, _: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value)?;
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

            fn decode<'a>(value: Stored<'a>, _: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
                let v = $read(value)?;
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

    /// Present is set, whatever value a writer gave it.
    fn decode<'a>(_: Stored<'a>, _: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
        Ok(true)
    }
}

impl sealed::Sealed for Str {}
impl ReadKind for Str {
    type Out<'a> = Option<&'a str>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn decode<'a>(value: Stored<'a>, _: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
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

    // r[impl vcf_values.gt+2]
    /// No alleles, or a lone missing one, is `.`; a missing diploid `./.`
    /// keeps its ploidy.
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
        let gt = match value.0 {
            StoredRepr::Bcf(v) => {
                check_ploidy(v.count())?;
                GtAlleles::Bcf(v.gt(cx.first_phase_is_read)?)
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
    const DYNAMIC: bool = true;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn decode<'a>(value: Stored<'a>, cx: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
        DynValue::decode(value, cx.declared, cx.genotype, cx).map(Some)
    }
}

// r[impl vcf_reader.limits]
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
        cx: DecodeCx,
    ) -> Result<Self, FieldError> {
        if genotype {
            return Ok(Self::Gt(Gt::decode(value, cx)?));
        }
        Ok(match declared {
            ValueType::Integer => Self::Int(whole(ints(value)?)),
            ValueType::Float => Self::Float(whole(floats(value)?)),
            ValueType::Flag => Self::Flag,
            ValueType::String | ValueType::Character => Self::Str(Str::decode(value, cx)?),
        })
    }
}
