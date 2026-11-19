//! [`Info`]: the INFO column, indexed on first access (`r[vcf_fields.info]`),
//! and what the INFO and FORMAT views share: the brand check and decoding a
//! stored value through a key ([`get`], [`check`]) or a field table
//! definition ([`get_dyn`]).

use super::access::{DecodeCx, DynValue, ReadCx, ReadKind, Stored, StoredRepr, Valid, dyn_get};
use super::bcf_value::TypedValue;
use super::error::FieldError;
use super::error::{RecordAt, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanContext, PlanError, Sections};
use super::file::FileShared;
use super::index::{DictId, Entry, SectionIndex};
use super::raw::RawRecord;
use crate::vcf::header::{InfoDef, ValueType};
use crate::vcf::record_encoder::{FieldId, InfoKey};
use std::cell::RefCell;
use std::sync::Arc;

// r[impl vcf_values.header_brand+2]
/// Fail with `ForeignKey` unless `key` was resolved from `file`'s header.
pub(crate) fn check_brand(key: &FieldId, file: &FileShared) -> Result<(), VcfReadError> {
    if key.header == file.header.id() || key.header == file.fields.id {
        Ok(())
    } else {
        Err(VcfReadError::ForeignKey { name: key.name.clone() })
    }
}

/// A record's INFO column as its file lays it out; what a view scans and
/// reads values from.
#[derive(Debug, Clone, Copy)]
pub(crate) enum InfoSource<'a> {
    /// BCF: `n` entries, each a typed key and a typed value.
    Bcf { bytes: &'a [u8], n: usize },
    /// VCF text: `;`-separated `key` or `key=value` entries.
    Text(&'a [u8]),
}

impl<'a> InfoSource<'a> {
    /// `bytes` as `file` lays INFO out; `n` is BCF's `n_info`.
    pub(crate) fn new(bytes: &'a [u8], n: usize, file: &FileShared) -> Self {
        if file.is_text() { Self::Text(bytes) } else { Self::Bcf { bytes, n } }
    }

    /// The value of `entry`: a BCF typed value, or the VCF text that follows
    /// its key, which the scan keeps with its `=` so that a key without one
    /// (holding nothing) and `key=` (holding `.`) differ.
    #[inline]
    fn value(self, entry: Entry) -> Result<Stored<'a>, RecordError> {
        match self {
            Self::Bcf { bytes, .. } => {
                let bytes = bytes.get(entry.span()).unwrap_or_default();
                let (value, _) = TypedValue::parse(bytes)
                    .map_err(|source| RecordError::Value { part: RecordPart::Info, source })?;
                Ok(Stored(StoredRepr::Bcf(value)))
            }
            Self::Text(bytes) => {
                let after_key = bytes.get(entry.span()).unwrap_or_default();
                Ok(Stored(StoredRepr::Text(after_key.split_first().map(|(_, value)| value))))
            }
        }
    }
}

/// How to check and read the field `name`: its key's, or its definition's
/// in the file's field table.
pub(crate) fn decode_cx<'f>(
    file: &'f FileShared,
    name: &'f str,
    declared: ValueType,
    genotype: bool,
) -> DecodeCx<'f> {
    let first_phase_is_read = file.first_phase_is_read;
    DecodeCx { read: ReadCx { first_phase_is_read, declared, genotype }, file: Some(file), name }
}

/// How to check and read the field `key` names.
pub(crate) fn key_cx<'f>(key: &'f FieldId, file: &'f FileShared) -> DecodeCx<'f> {
    decode_cx(file, &key.name, key.typ, key.genotype)
}

/// Check `values` as `V`, an error saying which field and record.
#[inline]
pub(crate) fn check<'a, V: ReadKind>(
    values: Stored<'a>,
    cx: DecodeCx<'_>,
    file: &FileShared,
    at: RecordAt,
) -> Result<Valid<'a, V>, VcfReadError> {
    Valid::check(values, cx).map_err(|source| field_error(cx, file, at, source))
}

/// Check and read one value as `V`, an error saying which field and record.
#[inline]
pub(crate) fn get<'a, V: ReadKind>(
    value: Stored<'a>,
    cx: DecodeCx<'_>,
    file: &FileShared,
    at: RecordAt,
) -> Result<V::Out<'a>, VcfReadError> {
    V::get(value, cx).map_err(|source| field_error(cx, file, at, source))
}

/// Check and read one value as its declaration says.
pub(crate) fn get_dyn<'a>(
    value: Stored<'a>,
    cx: DecodeCx<'_>,
    file: &FileShared,
    at: RecordAt,
) -> Result<DynValue<'a>, VcfReadError> {
    dyn_get(value, cx).map_err(|source| field_error(cx, file, at, source))
}

#[cold]
fn field_error(
    cx: DecodeCx<'_>,
    file: &FileShared,
    at: RecordAt,
    source: FieldError,
) -> VcfReadError {
    at.error(file, RecordError::Field { name: cx.name.into(), source })
}

// r[impl vcf_reader.debug]
/// `fields` as a map of name to value; a malformed section or value is
/// shown where it is found, and ends the map.
pub(crate) fn debug_fields<'a, D: 'a>(
    f: &mut std::fmt::Formatter<'_>,
    fields: impl Iterator<Item = Result<(&'a str, D, DynValue<'a>), VcfReadError>>,
) -> std::fmt::Result {
    let mut map = f.debug_map();
    for field in fields {
        match field {
            Ok((name, _, value)) => map.entry(&name, &value),
            Err(e) => return map.entry(&"error", &e.to_string()).finish(),
        };
    }
    map.finish()
}

/// A record's INFO column, borrowed; indexed on first access.
#[derive(Clone, Copy)]
pub struct InfoRef<'a> {
    info: InfoSource<'a>,
    /// The record's INFO index: the component's own, or the reader's for a
    /// [`RawRecord`] view. Borrowed for one method at a time.
    index: &'a RefCell<SectionIndex>,
    /// The field table that maps keys to ids, and the warn-once state.
    file: &'a FileShared,
    /// Where the record is, for errors.
    at: RecordAt,
}

/// The fields present, by name.
impl std::fmt::Debug for InfoRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Info ")?;
        debug_fields(f, self.iter())
    }
}

impl<'a> InfoRef<'a> {
    pub(crate) fn new(
        bytes: &'a [u8],
        n: usize,
        index: &'a RefCell<SectionIndex>,
        file: &'a FileShared,
        at: RecordAt,
    ) -> Self {
        Self { info: InfoSource::new(bytes, n, file), index, file, at }
    }

    /// The first entry of `id`, indexing no further than that.
    #[inline]
    fn find(&self, id: DictId) -> Result<Option<Entry>, VcfReadError> {
        let mut index = self.index.borrow_mut();
        index.find(id, self.info, self.file).map_err(|e| self.at.error(self.file, e.clone()))
    }

    #[inline]
    fn value(&self, entry: Entry) -> Result<Stored<'a>, VcfReadError> {
        self.info.value(entry).map_err(|e| self.at.error(self.file, e))
    }

    // r[impl vcf_values.access+2]
    /// The field `key` names, decoded as its marker says; absent or `.` is
    /// the marker's absent value (`None`, or `false` for a flag).
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key of another header, and
    /// [`VcfReadError::Record`] for malformed data.
    pub fn get<V: ReadKind>(&self, key: &InfoKey<V>) -> Result<V::Out<'a>, VcfReadError> {
        check_brand(&key.0, self.file)?;
        match self.find(key.0.id)? {
            Some(entry) => {
                get::<V>(self.value(entry)?, key_cx(&key.0, self.file), self.file, self.at)
            }
            None => Ok(V::absent()),
        }
    }

    // r[impl vcf_values.access+2]
    /// Whether the record holds the field `key` names, whatever its value:
    /// absent is `!has`, and `.` is `has` with an absent [`get`](Self::get).
    ///
    /// # Errors
    /// As [`Self::get`].
    pub fn has<V>(&self, key: &InfoKey<V>) -> Result<bool, VcfReadError> {
        check_brand(&key.0, self.file)?;
        Ok(self.find(key.0.id)?.is_some())
    }

    // r[impl vcf_values.dynamic+2]
    /// Every field present, in record order (a repeated key once), with its
    /// name, its definition in the file's field table and its value read as
    /// that says. A column that does not index is the first item's error.
    pub fn iter(
        &self,
    ) -> impl Iterator<Item = Result<(&'a str, &'a InfoDef, DynValue<'a>), VcfReadError>> + use<'a>
    {
        self.iter_ids().map(|field| field.map(|(_, name, def, value)| (name, def, value)))
    }

    /// [`Self::iter`], each field with its id in the file's field table.
    pub(crate) fn iter_ids(
        &self,
    ) -> impl Iterator<Item = Result<(DictId, &'a str, &'a InfoDef, DynValue<'a>), VcfReadError>> + use<'a>
    {
        let this = *self;
        let Self { info, index, file, at } = this;
        let finished = index.borrow_mut().finish(info, file).map_err(RecordError::clone);
        let (error, len) = match finished {
            Ok(()) => (None, index.borrow().entries().len()),
            Err(e) => (Some(Err(at.error(file, e))), 0),
        };
        error.into_iter().chain((0..len).filter_map(move |i| {
            let e = index.borrow().entries().get(i).copied()?;
            let (name, def) = file.fields.info_def(&file.header, e.id)?;
            let cx = decode_cx(file, name, def.typ, false);
            let value = this.value(e).and_then(|v| get_dyn(v, cx, file, at));
            Some(value.map(|v| (e.id, name, def, v)))
        }))
    }
}

// r[impl vcf_fields.info+2]
// r[impl vcf_fields.owned+2]
/// A record's INFO column, owned; its index is built on first access, so a
/// record rejected on its position or FILTER never pays for it.
#[derive(Clone)]
pub struct Info {
    /// The INFO column as stored (BCF entries or VCF text).
    bytes: Vec<u8>,
    /// BCF's `n_info`; 0 for VCF text.
    n: usize,
    file: Arc<FileShared>,
    /// Where the record is, for errors.
    at: RecordAt,
    /// The index of `bytes`, reset in [`Fields::fill`] and built on access.
    index: RefCell<SectionIndex>,
}

impl std::fmt::Debug for Info {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.view().fmt(f)
    }
}

impl Info {
    /// The file the record was read from.
    pub(crate) fn file(&self) -> &Arc<FileShared> {
        &self.file
    }

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

    /// Whether the record holds the field `key` names (see [`InfoRef::has`]).
    ///
    /// # Errors
    /// As [`InfoRef::has`].
    pub fn has<V>(&self, key: &InfoKey<V>) -> Result<bool, VcfReadError> {
        self.view().has(key)
    }

    /// Every field present (see [`InfoRef::iter`]).
    pub fn iter(
        &self,
    ) -> impl Iterator<Item = Result<(&str, &InfoDef, DynValue<'_>), VcfReadError>> {
        self.view().iter()
    }

    /// Every field present, with its id (see [`InfoRef::iter`]).
    pub(crate) fn iter_ids(
        &self,
    ) -> impl Iterator<Item = Result<(DictId, &str, &InfoDef, DynValue<'_>), VcfReadError>> {
        self.view().iter_ids()
    }
}

impl Fields for Info {
    type Plan = ();

    fn plan(_: &PlanContext<'_>) -> Result<(), PlanError> {
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
            index: RefCell::default(),
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
        self.index.get_mut().reset();
        raw.keep_file(&mut self.file);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::reader::{Locus, VcfFile};
    use crate::vcf::record_encoder::InfoInt;

    // r[verify vcf_fields.lazy+4]
    /// Filling a record does not index its INFO; the first access does, and
    /// the next fill drops that index (keeping its buffers).
    #[test]
    fn info_is_indexed_on_first_access_only() {
        let is_built = |info: &Info| !info.index.borrow().entries().is_empty();
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
        let file = VcfFile::open(path).unwrap();
        let dp: InfoInt = file.header().info_key("DP").unwrap();
        let mut reader = file.reader::<(Locus, Info)>().open().unwrap();
        let (_, info) = reader.next_record().unwrap().unwrap();
        assert!(!is_built(info));
        assert_eq!(info.get(&dp).unwrap(), Some(140));
        assert!(is_built(info));
        let (_, info) = reader.next_record().unwrap().unwrap();
        assert!(!is_built(info));
        assert_eq!(info.get(&dp).unwrap(), Some(60));
    }
}
