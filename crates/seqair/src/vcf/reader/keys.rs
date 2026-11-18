//! Resolving keys through a file rather than its header
//! (`r[vcf_values.file_keys]`, `r[vcf_values.optional_keys]`).

use super::fields::optional;
use super::file::{FileShared, VcfFile};
use crate::vcf::error::KeyError;
use crate::vcf::record_encoder::{FieldKind, FormatKey, InfoKey};

impl FileShared {
    pub(crate) fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.header.info_key(name)
    }

    pub(crate) fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.header.format_key(name)
    }
}

impl VcfFile {
    // r[impl vcf_values.file_keys]
    /// Resolve an INFO key for reading this file, checking marker `V`
    /// against the declaration: `let dp: InfoInt = file.info_key("DP")?;`.
    ///
    /// # Errors
    /// [`KeyError::UnknownField`] for a field the file does not have,
    /// `TypeMismatch`/`NumberMismatch`/`NotGenotype` for a marker that
    /// cannot read it.
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.0.info_key(name)
    }

    // r[impl vcf_values.file_keys]
    /// Resolve a FORMAT key for reading this file (see [`Self::info_key`]).
    ///
    /// # Errors
    /// As [`Self::info_key`].
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.0.format_key(name)
    }

    // r[impl vcf_values.optional_keys]
    /// An INFO key, or `None` if the file does not have the field: for an
    /// optional annotation, one `?` and one `Option`.
    ///
    /// # Errors
    /// As [`Self::info_key`], but for `UnknownField`.
    pub fn info_key_opt<V: FieldKind>(&self, name: &str) -> Result<Option<InfoKey<V>>, KeyError> {
        optional(self.info_key(name))
    }

    // r[impl vcf_values.optional_keys]
    /// A FORMAT key, or `None` if the file does not have the field.
    ///
    /// # Errors
    /// As [`Self::format_key`], but for `UnknownField`.
    pub fn format_key_opt<V: FieldKind>(
        &self,
        name: &str,
    ) -> Result<Option<FormatKey<V>>, KeyError> {
        optional(self.format_key(name))
    }
}
