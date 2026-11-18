//! [`Info`]: the INFO column, indexed on first access (`r[vcf_fields.info]`).

use super::access::{DecodeCx, DynValue, ReadKind, Stored, StoredRepr};
use super::bcf_value::TypedValue;
use super::entries::{BcfEntries, LazyIndex, index_bcf, scan};
use super::error::{KeyKind, RecordAt, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::file::FileShared;
use super::index::{DictId, FieldIndex};
use super::raw::RawRecord;
use crate::vcf::header::{InfoDef, VcfHeader};
use crate::vcf::record_encoder::{FieldId, InfoKey};
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

/// Decode the typed value in `bytes` as `V`.
pub(crate) fn decode<'a, V: ReadKind>(
    bytes: &'a [u8],
    key: &FieldId,
    part: RecordPart,
    file: &FileShared,
    at: RecordAt,
) -> Result<V::Out<'a>, VcfReadError> {
    let (value, _) =
        TypedValue::parse(bytes).map_err(|source| at.error(RecordError::Value { part, source }))?;
    let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
    V::decode(Stored(StoredRepr::Bcf(value)), cx)
        .map_err(|source| at.error(RecordError::Field { name: key.name.clone(), source }))
}

/// Where an INFO view finds a field.
#[derive(Debug, Clone, Copy)]
enum Source<'a> {
    Indexed(&'a Info),
    Scan { bytes: &'a [u8], n: usize, file: &'a FileShared, at: RecordAt },
}

/// A record's INFO column, borrowed.
#[derive(Debug, Clone, Copy)]
pub struct InfoRef<'a>(Source<'a>);

impl<'a> InfoRef<'a> {
    /// A view that walks the BCF entries on each access.
    pub(crate) fn scan(bytes: &'a [u8], n: usize, file: &'a FileShared, at: RecordAt) -> Self {
        Self(Source::Scan { bytes, n, file, at })
    }

    // r[impl vcf_values.access]
    /// The field `key` names, decoded as its marker says; absent or `.` is
    /// the marker's absent value (`None`, or `false` for a flag).
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key of another header, and
    /// [`VcfReadError::Record`] for malformed data.
    pub fn get<V: ReadKind>(&self, key: &InfoKey<V>) -> Result<V::Out<'a>, VcfReadError> {
        let (bytes, file, at, span) = match self.0 {
            Source::Indexed(info) => {
                let Some(file) = info.file.as_deref() else { return Ok(V::absent()) };
                check_brand(&key.0, file)?;
                let span = DictId::new(key.0.dict_idx)
                    .map(|id| Ok::<_, VcfReadError>(info.index(file)?.get(id).map(|e| e.span())))
                    .transpose()?
                    .flatten();
                (info.bytes.as_slice(), file, info.at, span)
            }
            Source::Scan { bytes, n, file, at } => {
                check_brand(&key.0, file)?;
                let span = DictId::new(key.0.dict_idx)
                    .map(|id| scan(BcfEntries::new(bytes, n, None), id))
                    .transpose()
                    .map_err(|source| {
                        at.error(RecordError::Value { part: RecordPart::Info, source })
                    })?
                    .flatten();
                (bytes, file, at, span)
            }
        };
        match span {
            Some(span) => {
                decode::<V>(bytes.get(span).unwrap_or_default(), &key.0, RecordPart::Info, file, at)
            }
            None => Ok(V::absent()),
        }
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
    file: Option<Arc<FileShared>>,
    at: RecordAt,
    index: LazyIndex,
}

impl Info {
    fn index(&self, file: &FileShared) -> Result<&FieldIndex, VcfReadError> {
        let entries = BcfEntries::new(&self.bytes, self.n, None);
        self.index
            .get(|index| index_bcf(index, entries, KeyKind::Info, file))
            .map_err(|e| self.at.error(e.clone()))
    }

    /// The field `key` names (see [`InfoRef::get`]).
    ///
    /// # Errors
    /// As [`InfoRef::get`].
    pub fn get<V: ReadKind>(&self, key: &InfoKey<V>) -> Result<V::Out<'_>, VcfReadError> {
        InfoRef(Source::Indexed(self)).get(key)
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
        impl Iterator<Item = Result<(&str, &InfoDef, DynValue<'_>), VcfReadError>>,
        VcfReadError,
    > {
        let entries = match self.file.as_deref() {
            Some(file) => self.index(file)?.entries(),
            None => &[],
        };
        Ok(entries.iter().filter_map(move |e| {
            let file = self.file.as_deref()?;
            let name = file.header.string_map().name(e.id.get() as usize)?;
            let (name, def) = file.header.infos().get_key_value(name)?;
            Some(
                dyn_value(
                    self.bytes.get(e.span()).unwrap_or_default(),
                    name,
                    def.typ,
                    false,
                    file,
                    self.at,
                )
                .map(|v| (name.as_str(), def, v)),
            )
        }))
    }
}

/// A value read as its declaration says.
pub(crate) fn dyn_value<'a>(
    bytes: &'a [u8],
    name: &str,
    typ: crate::vcf::ValueType,
    genotype: bool,
    file: &FileShared,
    at: RecordAt,
) -> Result<DynValue<'a>, VcfReadError> {
    let part = if genotype { RecordPart::Format } else { RecordPart::Info };
    let (value, _) =
        TypedValue::parse(bytes).map_err(|source| at.error(RecordError::Value { part, source }))?;
    let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
    DynValue::decode(Stored(StoredRepr::Bcf(value)), typ, genotype, cx)
        .map_err(|source| at.error(RecordError::Field { name: name.into(), source }))
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
        let mut info =
            Self { bytes: Vec::new(), n: 0, file: None, at: raw.at(), index: LazyIndex::default() };
        info.fill(plan, raw)?;
        Ok(info)
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        let (bytes, n) = raw.bcf_info()?;
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
