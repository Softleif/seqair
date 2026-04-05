//! VCF/BCF writing. Construct a [`VcfHeader`] via its builder, create
//! [`VcfRecord`]s with type-safe [`Alleles`], and write to VCF text or BCF binary.

pub mod alleles;
pub mod error;
pub mod header;

pub use alleles::Alleles;
pub use error::{AllelesError, VcfEncodeError, VcfError, VcfHeaderError};
pub use header::{
    ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader, VcfHeaderBuilder,
};
