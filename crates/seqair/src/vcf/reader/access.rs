//! Typed access (`r[vcf_values.access]`): what each key marker reads a stored
//! value as, and the untyped [`DynValue`] for generic tools.

use super::bcf_value::TypedValue;
use super::error::FieldError;
use super::values::{GenotypeRef, GtAlleles, OptValues, Values};
use crate::vcf::header::ValueType;
use crate::vcf::record_encoder::{Arr, FieldKind, Flag, Gt, OptArr, Scalar, Str};

/// Most alleles a genotype may have (`r[vcf_reader.limits]`).
pub(crate) const MAX_PLOIDY: usize = 256;

/// One field's value as a record stores it, not yet decoded.
#[derive(Debug, Clone, Copy)]
pub struct Stored<'a>(pub(crate) StoredRepr<'a>);

#[derive(Debug, Clone, Copy)]
pub(crate) enum StoredRepr<'a> {
    Bcf(TypedValue<'a>),
}

/// How to decode: what a `Gt` needs from the file.
#[derive(Debug, Clone, Copy)]
pub struct DecodeCx {
    pub(crate) first_phase_is_read: bool,
}

mod sealed {
    pub trait Sealed {}
}

// r[impl vcf_values.access]
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
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx) -> Result<Self::Out<'a>, FieldError>;
}

fn typed(value: Stored<'_>) -> TypedValue<'_> {
    match value.0 {
        StoredRepr::Bcf(v) => v,
    }
}

fn ints(value: Stored<'_>) -> Result<OptValues<'_, i32>, FieldError> {
    Ok(OptValues::bcf(typed(value).ints()?))
}

fn floats(value: Stored<'_>) -> Result<OptValues<'_, f32>, FieldError> {
    Ok(OptValues::bcf(typed(value).floats()?))
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
        let s = typed(value).str()?;
        Ok((!s.is_empty() && s != ".").then_some(s))
    }
}

impl sealed::Sealed for Gt {}
impl ReadKind for Gt {
    type Out<'a> = Option<GenotypeRef<'a>>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    // r[impl vcf_values.gt]
    /// No alleles, or a lone missing one, is `.`; a missing diploid `./.`
    /// keeps its ploidy.
    fn decode<'a>(value: Stored<'a>, cx: DecodeCx) -> Result<Self::Out<'a>, FieldError> {
        let v = typed(value);
        if v.count() > MAX_PLOIDY {
            return Err(FieldError::Ploidy { ploidy: v.count(), max: MAX_PLOIDY });
        }
        let gt = GtAlleles::Bcf(v.gt(cx.first_phase_is_read)?);
        let mut alleles = gt.clone();
        Ok(match (alleles.next(), alleles.next()) {
            (None, _) | (Some(super::values::GtAllele { allele: None, .. }), None) => None,
            _ => Some(GenotypeRef(gt)),
        })
    }
}

// r[impl vcf_values.dynamic]
/// A field's value read by its header declaration, for tools that handle
/// every field alike (`convert`, `view`).
#[derive(Debug, Clone)]
pub enum DynValue<'a> {
    /// `Type=Integer`; missing elements are `None`.
    Int(OptValues<'a, i32>),
    /// `Type=Float`; missing elements are `None`.
    Float(OptValues<'a, f32>),
    /// `Type=Flag` (INFO only), present.
    Flag,
    /// `Type=String` or `Character`; `None` for `.`.
    Str(Option<&'a str>),
    /// The FORMAT `GT` field; `None` for `.`.
    Gt(Option<GenotypeRef<'a>>),
}

impl<'a> DynValue<'a> {
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
            ValueType::Integer => Self::Int(ints(value)?),
            ValueType::Float => Self::Float(floats(value)?),
            ValueType::Flag => Self::Flag,
            ValueType::String | ValueType::Character => Self::Str(Str::decode(value, cx)?),
        })
    }
}
