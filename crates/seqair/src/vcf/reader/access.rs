//! Typed access (`r[vcf_values.access]`): what each key marker reads a stored
//! value as, and the untyped [`DynValue`] for generic tools.
//!
//! Reading is two steps. [`Valid::check`] checks stored values once — one
//! value, or every sample of a BCF FORMAT column at once — and only it
//! fails. [`ReadKind::read`] then builds the marker's output from a
//! [`Valid`] value without failing: the readers underneath (`BcfVec`,
//! `BcfGt`, `TextVec`, `TextGt`) are total, right on checked bytes, and read
//! an element a check would reject as missing.

use super::bcf_value::{BcfGt, BcfVec, TypedValue};
use super::error::FieldError;
use super::file::FileShared;
use super::text_value::{TextGt, TextVec, text_str};
use super::values::{GenotypeRef, GtAlleles, Number, OptValues, Values};
use crate::vcf::header::ValueType;
use crate::vcf::record_encoder::{
    Arr, FieldKind, Flag, FormatKey, Gt, InfoKey, OptArr, Scalar, Str,
};
use std::marker::PhantomData;

/// Most alleles a genotype may have (`r[vcf_reader.limits+2]`).
pub(crate) const MAX_PLOIDY: usize = 256;

/// One field's values as a record stores them, not yet checked.
#[derive(Debug, Clone, Copy)]
pub struct Stored<'a>(pub(crate) StoredRepr<'a>);

/// Which format the values are in, so each marker picks its check and reader.
#[derive(Debug, Clone, Copy)]
pub(crate) enum StoredRepr<'a> {
    /// One BCF value, or every sample's of a FORMAT column.
    Bcf(TypedValue<'a>),
    /// A VCF text value: `None` for an INFO key without `=`, which holds
    /// nothing (BCF's empty vector); an empty value is `.` (missing).
    Text(Option<&'a [u8]>),
}

/// The text of a value as it parses: an empty value is `.`.
fn text_value(value: &[u8]) -> &[u8] {
    if value.is_empty() { b"." } else { value }
}

/// What reading a checked value needs, from the field's key or its
/// definition in the file's field table.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadCx {
    /// Whether a genotype's first allele carries its own phase bit (VCF
    /// 4.4+, `r[vcf_record.gt_first_phase+3]`), from the file's version.
    pub(crate) first_phase_is_read: bool,
    /// The field's declared `Type`: what [`Any`] reads it as.
    pub(crate) declared: ValueType,
    /// Whether the field is FORMAT `GT`, which [`Any`] reads as a genotype.
    pub(crate) genotype: bool,
}

/// How to check and read: what the readers need, and where to warn.
#[derive(Debug, Clone, Copy)]
pub struct DecodeCx<'f> {
    pub(crate) read: ReadCx,
    /// The file to warn of a reserved integer read as missing, and the field's name.
    pub(crate) file: Option<&'f FileShared>,
    pub(crate) name: &'f str,
}

/// One value that marker `V`'s check passed: what [`ReadKind::read`] reads
/// without failing. Made by [`Valid::check`], or by a column for a sample
/// whose whole column that check passed.
pub struct Valid<'a, V: ?Sized> {
    pub(crate) value: Stored<'a>,
    pub(crate) cx: ReadCx,
    _v: PhantomData<fn() -> V>,
}

impl<V: ?Sized> Clone for Valid<'_, V> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<V: ?Sized> Copy for Valid<'_, V> {}

impl<'a, V: ReadKind + ?Sized> Valid<'a, V> {
    /// Check `values`: one value, or every sample's of a BCF column, failing
    /// with the first sample's error as reading them one by one would.
    #[inline]
    pub(crate) fn check(values: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self, FieldError> {
        // A column of no samples holds no value to fail, whatever its type.
        if !matches!(values.0, StoredRepr::Bcf(v) if v.n() == 0) {
            V::check(values, cx)?;
        }
        Ok(Self::checked(values, cx.read))
    }

    /// `value`, which the check of it or of the column holding it passed.
    #[inline]
    pub(crate) fn checked(value: Stored<'a>, cx: ReadCx) -> Self {
        Self { value, cx, _v: PhantomData }
    }

    /// The value read as `V`.
    #[inline]
    pub(crate) fn read(self) -> V::Out<'a> {
        V::read(self)
    }
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

    /// Check `values` (one, or every sample's of a BCF column, in sample
    /// order); `Ok` means [`Self::read`] reads each as stored.
    #[doc(hidden)]
    fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError>;

    #[doc(hidden)]
    fn read(value: Valid<'_, Self>) -> Self::Out<'_>;

    /// Check one value and read it. Markers whose check builds the value
    /// anyway (numbers, strings) read it in that same pass, and derive
    /// `check` from this.
    #[doc(hidden)]
    #[inline]
    fn get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        Ok(Valid::<Self>::check(value, cx)?.read())
    }
}

/// Run `check` on each value of `values`: each sample of a BCF column, in
/// order, or the one value.
fn each<'a>(
    values: Stored<'a>,
    mut check: impl FnMut(Stored<'a>) -> Result<(), FieldError>,
) -> Result<(), FieldError> {
    match values.0 {
        StoredRepr::Bcf(v) if v.n() != 1 => {
            v.samples().try_for_each(|s| check(Stored(StoredRepr::Bcf(s))))
        }
        _ => check(values),
    }
}

/// One value's numbers, checked: the type of BCF values, every element of
/// VCF text (warning of a reserved integer read as missing).
fn numbers<'a, T: Number>(
    value: Stored<'a>,
    cx: DecodeCx<'_>,
) -> Result<OptValues<'a, T>, FieldError> {
    Ok(match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(v.numbers()?),
        StoredRepr::Text(None) => OptValues::text(TextVec::empty()),
        // r[impl vcf_values.text_numbers+2]
        StoredRepr::Text(Some(t)) => {
            let (v, reserved) = TextVec::parse(text_value(t))?;
            if reserved && let Some(file) = cx.file {
                file.warn_reserved_int(cx.name);
            }
            OptValues::text(v)
        }
    })
}

/// One value's numbers, read as checked.
fn read_numbers<T: Number>(value: Stored<'_>) -> OptValues<'_, T> {
    match value.0 {
        StoredRepr::Bcf(v) => OptValues::bcf(BcfVec::new(v)),
        StoredRepr::Text(None) => OptValues::text(TextVec::empty()),
        StoredRepr::Text(Some(t)) => OptValues::text(TextVec::new(text_value(t))),
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

macro_rules! numeric_kinds {
    ($($t:ty),*) => {$(
        impl sealed::Sealed for Scalar<$t> {}
        impl ReadKind for Scalar<$t> {
            type Out<'a> = Option<$t>;

            fn absent<'a>() -> Self::Out<'a> {
                None
            }

            fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError> {
                each(values, |v| Self::get(v, cx).map(drop))
            }

            fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
                read_numbers::<$t>(value.value).next().flatten()
            }

            fn get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = numbers::<$t>(value, cx)?;
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

            fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError> {
                each(values, |v| Self::get(v, cx).map(drop))
            }

            fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
                let v = read_numbers(value.value);
                (!v.is_missing()).then(|| Values::checked(v))
            }

            fn get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = numbers::<$t>(value, cx)?;
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

            fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError> {
                each(values, |v| Self::get(v, cx).map(drop))
            }

            fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
                let v = read_numbers(value.value);
                (!v.is_missing()).then_some(v)
            }

            fn get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
                let v = numbers::<$t>(value, cx)?;
                Ok((!v.is_missing()).then_some(v))
            }
        }
    )*};
}

numeric_kinds!(i32, f32);

impl sealed::Sealed for Flag {}
impl ReadKind for Flag {
    type Out<'a> = bool;

    fn absent<'a>() -> Self::Out<'a> {
        false
    }

    fn check(_: Stored<'_>, _: DecodeCx<'_>) -> Result<(), FieldError> {
        Ok(())
    }

    /// Present is set, whatever value a writer gave it.
    fn read(_: Valid<'_, Self>) -> Self::Out<'_> {
        true
    }
}

/// One value as a string, `None` for `.` or nothing.
fn string(value: Stored<'_>) -> Result<Option<&str>, FieldError> {
    let s = match value.0 {
        StoredRepr::Bcf(v) => v.str()?,
        StoredRepr::Text(t) => text_str(t.unwrap_or_default())?,
    };
    Ok((!s.is_empty() && s != ".").then_some(s))
}

impl sealed::Sealed for Str {}
impl ReadKind for Str {
    type Out<'a> = Option<&'a str>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError> {
        each(values, |v| Self::get(v, cx).map(drop))
    }

    /// UTF-8 has no unchecked reading without `unsafe`: checked again, and
    /// bytes the check would reject read as `None`.
    fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
        string(value.value).ok().flatten()
    }

    fn get<'a>(value: Stored<'a>, _: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        string(value)
    }
}

impl sealed::Sealed for Gt {}
impl ReadKind for Gt {
    type Out<'a> = Option<GenotypeRef<'a>>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    // r[impl vcf_values.gt+2]
    /// The ploidy before the alleles: a BCF column's declared count, a text
    /// genotype's separators.
    fn check(values: Stored<'_>, _: DecodeCx<'_>) -> Result<(), FieldError> {
        match values.0 {
            StoredRepr::Bcf(v) => {
                check_ploidy(v.count())?;
                Ok(v.check_gt()?)
            }
            StoredRepr::Text(t) => {
                let t = text_value(t.unwrap_or_default());
                // Every allele but the first follows a separator; so may the first.
                let sep = |b: &u8| *b == b'/' || *b == b'|';
                let seps = t.iter().filter(|b| sep(b)).count();
                check_ploidy(seps.saturating_add(usize::from(!t.first().is_some_and(sep))))?;
                Ok(TextGt::check(t)?)
            }
        }
    }

    // r[impl vcf_values.gt+2]
    /// No alleles, or a lone missing one, is `.`; a missing diploid `./.`
    /// keeps its ploidy.
    #[inline(always)]
    fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
        let first_phase_is_read = value.cx.first_phase_is_read;
        let gt = match value.value.0 {
            StoredRepr::Bcf(v) => GtAlleles::Bcf(BcfGt::new(v, first_phase_is_read)),
            StoredRepr::Text(t) => {
                GtAlleles::Text(TextGt::new(text_value(t.unwrap_or_default()), first_phase_is_read))
            }
        };
        (!gt.is_missing()).then_some(GenotypeRef(gt))
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

/// Re-mark a checked value: `Any` checked it as the marker its declaration names.
fn as_kind<V>(value: Valid<'_, Any>) -> Valid<'_, V> {
    Valid { value: value.value, cx: value.cx, _v: PhantomData }
}

impl sealed::Sealed for Any {}
impl ReadKind for Any {
    type Out<'a> = Option<DynValue<'a>>;

    fn absent<'a>() -> Self::Out<'a> {
        None
    }

    fn check(values: Stored<'_>, cx: DecodeCx<'_>) -> Result<(), FieldError> {
        if cx.read.genotype {
            return Gt::check(values, cx);
        }
        each(values, |v| dyn_get(v, cx).map(drop))
    }

    fn read(value: Valid<'_, Self>) -> Self::Out<'_> {
        Some(value.read_dyn())
    }

    fn get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<Self::Out<'a>, FieldError> {
        dyn_get(value, cx).map(Some)
    }
}

/// One value checked and read as its declaration says (see
/// [`Valid::read_dyn`]).
pub(crate) fn dyn_get<'a>(value: Stored<'a>, cx: DecodeCx<'_>) -> Result<DynValue<'a>, FieldError> {
    if cx.read.genotype {
        return Gt::get(value, cx).map(DynValue::Gt);
    }
    Ok(match cx.read.declared {
        ValueType::Integer => DynValue::Int(whole(numbers(value, cx)?)),
        ValueType::Float => DynValue::Float(whole(numbers(value, cx)?)),
        ValueType::Flag => DynValue::Flag,
        ValueType::String | ValueType::Character => DynValue::Str(Str::get(value, cx)?),
    })
}

impl<'a> Valid<'a, Any> {
    // r[impl vcf_values.dynamic.missing]
    /// The value as its declaration says. A number vector missing as a
    /// whole is empty: VCF text spells an empty vector `.`, so the formats
    /// agree only on this.
    pub(crate) fn read_dyn(self) -> DynValue<'a> {
        if self.cx.genotype {
            return DynValue::Gt(Gt::read(as_kind(self)));
        }
        match self.cx.declared {
            ValueType::Integer => DynValue::Int(whole(read_numbers(self.value))),
            ValueType::Float => DynValue::Float(whole(read_numbers(self.value))),
            ValueType::Flag => DynValue::Flag,
            ValueType::String | ValueType::Character => DynValue::Str(Str::read(as_kind(self))),
        }
    }
}

/// `v`, or no elements if it is missing as a whole.
fn whole<T: Number>(v: OptValues<'_, T>) -> OptValues<'_, T> {
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

impl DynValue<'_> {
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
}
