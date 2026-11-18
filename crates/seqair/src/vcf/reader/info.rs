//! [`Info`]: the INFO column, indexed on first access (`r[vcf_fields.info]`).

use super::access::{DecodeCx, DynValue, ReadKind, Stored, StoredRepr};
use super::bcf_value::TypedValue;
use super::entries::{BcfEntries, LazyIndex, index_bcf, index_text_info};
use super::error::{KeyKind, RecordAt, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::file::FileShared;
use super::index::{DictId, FieldIndex};
use super::raw::RawRecord;
use crate::vcf::header::{InfoDef, ValueType, VcfHeader};
use crate::vcf::record_encoder::{FieldId, InfoKey};
use seqair_types::SmolStr;
use std::sync::Arc;

// r[impl vcf_values.header_brand]
/// Fail with `ForeignKey` unless `key` was resolved from `file`'s header.
pub(crate) fn check_brand(key: &FieldId, file: &FileShared) -> Result<(), VcfReadError> {
    if key.header == file.header.id() {
        Ok(())
    } else {
        Err(VcfReadError::ForeignKey { name: key.name.clone() })
    }
}

/// One field's value in a section's bytes: a BCF typed value, or the VCF
/// text that follows an INFO key (empty, or `=` and the value) or one
/// sample's value of a FORMAT field.
pub(crate) fn stored<'a>(
    bytes: &'a [u8],
    file: &FileShared,
    part: RecordPart,
    at: RecordAt,
) -> Result<Stored<'a>, VcfReadError> {
    if !file.is_text() {
        let (value, _) = TypedValue::parse(bytes)
            .map_err(|source| at.error(file, RecordError::Value { part, source }))?;
        return Ok(Stored(StoredRepr::Bcf(value)));
    }
    Ok(Stored(StoredRepr::Text(match part {
        RecordPart::Info => bytes.split_first().map(|(_, value)| value),
        _ => Some(bytes),
    })))
}

/// Decode `value`, the field `name`'s, as `V`.
pub(crate) fn decode<'a, V: ReadKind>(
    value: Stored<'a>,
    name: &SmolStr,
    file: &FileShared,
    at: RecordAt,
) -> Result<V::Out<'a>, VcfReadError> {
    let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
    V::decode(value, cx)
        .map_err(|source| at.error(file, RecordError::Field { name: name.clone(), source }))
}

/// Decode `value`, the field `name`'s, as `declared` says.
pub(crate) fn dyn_decode<'a>(
    value: Stored<'a>,
    name: &str,
    declared: ValueType,
    genotype: bool,
    file: &FileShared,
    at: RecordAt,
) -> Result<DynValue<'a>, VcfReadError> {
    let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
    DynValue::decode(value, declared, genotype, cx)
        .map_err(|source| at.error(file, RecordError::Field { name: name.into(), source }))
}

/// Index an INFO column as its file's format lays it out.
fn build(
    index: &mut FieldIndex,
    bytes: &[u8],
    n: usize,
    file: &FileShared,
) -> Result<(), RecordError> {
    if file.is_text() {
        index_text_info(index, bytes, file)
    } else {
        index_bcf(index, BcfEntries::new(bytes, n, None), KeyKind::Info, file)
    }
}

/// A record's INFO column, borrowed; indexed on first access.
#[derive(Debug, Clone, Copy)]
pub struct InfoRef<'a> {
    bytes: &'a [u8],
    /// BCF's `n_info`.
    n: usize,
    index: &'a LazyIndex,
    file: &'a FileShared,
    at: RecordAt,
}

impl<'a> InfoRef<'a> {
    pub(crate) fn new(
        bytes: &'a [u8],
        n: usize,
        index: &'a LazyIndex,
        file: &'a FileShared,
        at: RecordAt,
    ) -> Self {
        Self { bytes, n, index, file, at }
    }

    fn index(&self) -> Result<&'a FieldIndex, VcfReadError> {
        let Self { bytes, n, file, at, .. } = *self;
        self.index.get(|index| build(index, bytes, n, file)).map_err(|e| at.error(file, e.clone()))
    }

    fn value(&self, span: std::ops::Range<usize>) -> Result<Stored<'a>, VcfReadError> {
        let bytes = self.bytes.get(span).unwrap_or_default();
        stored(bytes, self.file, RecordPart::Info, self.at)
    }

    // r[impl vcf_values.access]
    /// The field `key` names, decoded as its marker says; absent or `.` is
    /// the marker's absent value (`None`, or `false` for a flag).
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key of another header, and
    /// [`VcfReadError::Record`] for malformed data.
    pub fn get<V: ReadKind>(&self, key: &InfoKey<V>) -> Result<V::Out<'a>, VcfReadError> {
        check_brand(&key.0, self.file)?;
        let Some(id) = DictId::new(key.0.dict_idx) else { return Ok(V::absent()) };
        match self.index()?.get(id) {
            Some(entry) => decode::<V>(self.value(entry.span())?, &key.0.name, self.file, self.at),
            None => Ok(V::absent()),
        }
    }

    // r[impl vcf_values.dynamic]
    /// Every declared field present, in record order (a repeated key once),
    /// with its declaration and its value read as declared.
    ///
    /// # Errors
    /// [`VcfReadError::Record`] if the column is malformed; each item, if
    /// its value is.
    pub fn iter(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<(&'a str, &'a InfoDef, DynValue<'a>), VcfReadError>> + use<'a>,
        VcfReadError,
    > {
        let this = *self;
        let header = &self.file.header;
        Ok(self.index()?.entries().iter().filter_map(move |e| {
            let name = header.string_map().name(e.id.get() as usize)?;
            let (name, def) = header.infos().get_key_value(name)?;
            let value = this
                .value(e.span())
                .and_then(|v| dyn_decode(v, name, def.typ, false, this.file, this.at));
            Some(value.map(|v| (name.as_str(), def, v)))
        }))
    }
}

// r[impl vcf_fields.info]
// r[impl vcf_fields.owned]
/// A record's INFO column, owned; its index is built on first access, so a
/// record rejected on its position or FILTER never pays for it.
#[derive(Debug, Clone)]
pub struct Info {
    bytes: Vec<u8>,
    n: usize,
    file: Arc<FileShared>,
    at: RecordAt,
    index: LazyIndex,
}

impl Info {
    fn view(&self) -> InfoRef<'_> {
        InfoRef::new(&self.bytes, self.n, &self.index, &self.file, self.at)
    }

    /// The field `key` names (see [`InfoRef::get`]).
    ///
    /// # Errors
    /// As [`InfoRef::get`].
    pub fn get<V: ReadKind>(&self, key: &InfoKey<V>) -> Result<V::Out<'_>, VcfReadError> {
        self.view().get(key)
    }

    /// Every declared field present (see [`InfoRef::iter`]).
    ///
    /// # Errors
    /// As [`InfoRef::iter`].
    pub fn iter(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<(&str, &InfoDef, DynValue<'_>), VcfReadError>>,
        VcfReadError,
    > {
        self.view().iter()
    }
}

impl Fields for Info {
    type Plan = ();

    fn plan(_: &VcfHeader, _: &ReadOptions) -> Result<(), PlanError> {
        Ok(())
    }

    fn needs(_: &()) -> Sections {
        Sections::INFO
    }

    fn new(plan: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut info = Self {
            bytes: Vec::new(),
            n: 0,
            file: Arc::clone(raw.file),
            at: raw.at(),
            index: LazyIndex::default(),
        };
        info.fill(plan, raw)?;
        Ok(info)
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        let (bytes, n) = raw.info_bytes()?;
        self.bytes.clear();
        self.bytes.extend_from_slice(bytes);
        self.n = n;
        self.at = raw.at();
        self.index.reset();
        raw.keep_file(&mut self.file);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::reader::{Locus, VcfFile};
    use crate::vcf::record_encoder::InfoInt;

    // r[verify vcf_fields.lazy]
    /// Filling a record does not index its INFO; the first access does, and
    /// the next fill drops that index (keeping its buffers).
    #[test]
    fn info_is_indexed_on_first_access_only() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
        let file = VcfFile::open(path).unwrap();
        let dp: InfoInt = file.header().info_key("DP").unwrap();
        let mut reader = file.reader::<(Locus, Info)>().open().unwrap();
        let (_, info) = reader.next_record().unwrap().unwrap();
        assert!(!info.index.is_built());
        assert_eq!(info.get(&dp).unwrap(), Some(140));
        assert!(info.index.is_built());
        let (_, info) = reader.next_record().unwrap().unwrap();
        assert!(!info.index.is_built());
        assert_eq!(info.get(&dp).unwrap(), Some(60));
    }
}
