//! FORMAT values: one sample ([`OnlySample`], [`Sample`]) or all
//! ([`Samples`]), by sample or by column (`r[vcf_fields.only_sample]`,
//! `r[vcf_fields.sample]`, `r[vcf_fields.samples]`).

use super::access::{DynValue, ReadKind, Stored, StoredRepr};
use super::bcf_value::{BcfValueError, FormatColumn, read_typed_int};
use super::entries::{
    BcfEntries, Cursor, FormatIndex, LazyIndex, Scan, index_bcf, index_text_format,
};
use super::error::{KeyKind, RecordAt, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanContext, PlanError, Sections};
use super::file::FileShared;
use super::index::{DictId, Entry};
use super::info::{check_brand, debug_fields, decode, dyn_decode};
use super::raw::RawRecord;
use crate::vcf::header::{FormatDef, SampleIdx};
use crate::vcf::record_encoder::{FieldId, FormatKey};
use seqair_types::{SmallVec, SmolStr};
use std::marker::PhantomData;
use std::ops::Range;
use std::sync::Arc;

/// Where a record's FORMAT data comes from, as the reader found it.
#[derive(Debug, Clone, Copy)]
pub(crate) enum FormatSource<'a> {
    /// The individual section: per field, its key, its type descriptor,
    /// then each sample's elements.
    Bcf { indiv: &'a [u8], n_fmt: usize, n_samples: usize },
    /// The FORMAT column and the TAB-separated columns of the samples read.
    Text { keys: &'a [u8], samples: &'a [u8], n_samples: usize },
}

/// FORMAT data, owned: BCF's individual section (all samples, or one
/// sample's elements of each field), or VCF text's FORMAT column, a TAB,
/// and the sample columns read.
#[derive(Debug, Clone)]
struct FormatData {
    bytes: Vec<u8>,
    /// BCF: the number of fields.
    n_fmt: usize,
    /// VCF text: the FORMAT column's length.
    keys_len: usize,
    n_samples: usize,
    file: Arc<FileShared>,
    at: RecordAt,
    index: LazyIndex<FormatIndex>,
    /// The error copying one sample's values out of a multi-sample BCF
    /// record met, returned by every access (`r[vcf_fields.lazy+2]`).
    failed: Option<RecordError>,
}

impl FormatData {
    fn new(raw: &RawRecord<'_>) -> Self {
        Self {
            bytes: Vec::new(),
            n_fmt: 0,
            keys_len: 0,
            n_samples: 0,
            file: Arc::clone(raw.file),
            at: raw.at(),
            index: LazyIndex::default(),
            failed: None,
        }
    }

    fn begin(&mut self, raw: &RawRecord<'_>, n_samples: usize) {
        self.bytes.clear();
        self.n_fmt = 0;
        self.keys_len = 0;
        self.n_samples = n_samples;
        self.at = raw.at();
        self.index.reset();
        self.failed = None;
        raw.keep_file(&mut self.file);
    }

    /// Copy the FORMAT data of `sample`, or of every sample.
    fn fill(&mut self, raw: &RawRecord<'_>, sample: Option<SampleIdx>) -> Result<(), VcfReadError> {
        match raw.format_source(sample)? {
            FormatSource::Bcf { indiv, n_fmt, n_samples } => match sample {
                // One sample's values are the whole section: copy it as it
                // is, leaving malformed data to the index, which fails on access.
                None | Some(_) if n_samples == 1 => {
                    self.begin(raw, 1);
                    self.bytes.extend_from_slice(indiv);
                    self.n_fmt = n_fmt;
                }
                None => {
                    self.begin(raw, n_samples);
                    self.bytes.extend_from_slice(indiv);
                    self.n_fmt = n_fmt;
                }
                // r[impl vcf_fields.lazy+2]
                // Malformed data fails the accesses, not the read, as it
                // does for every other shape of FORMAT data.
                Some(sample) => {
                    self.begin(raw, 1);
                    self.n_fmt = n_fmt;
                    if let Err(e) =
                        copy_one_sample(&mut self.bytes, indiv, n_fmt, n_samples, sample.index())
                    {
                        self.bytes.clear();
                        self.failed = Some(e);
                    }
                }
            },
            FormatSource::Text { keys, samples, n_samples } => {
                self.begin(raw, n_samples);
                self.keys_len = keys.len();
                self.bytes.extend_from_slice(keys);
                if !keys.is_empty() {
                    self.bytes.push(b'\t');
                    self.bytes.extend_from_slice(samples);
                }
            }
        }
        Ok(())
    }

    fn view(&self) -> FormatView<'_> {
        let (keys, samples) = if self.file.is_text() {
            let samples = self.bytes.get(self.keys_len.saturating_add(1)..).unwrap_or_default();
            (self.bytes.get(..self.keys_len).unwrap_or_default(), samples)
        } else {
            (&[][..], &self.bytes[..])
        };
        FormatView {
            keys,
            samples,
            n_fmt: self.n_fmt,
            n_samples: self.n_samples,
            index: &self.index,
            failed: self.failed.as_ref(),
            file: &self.file,
            at: self.at,
        }
    }
}

/// Copy only sample `sample`'s values from a BCF individual section: each
/// field's key and descriptor, then that sample's elements.
fn copy_one_sample(
    out: &mut Vec<u8>,
    indiv: &[u8],
    n_fmt: usize,
    n_samples: usize,
    sample: usize,
) -> Result<(), RecordError> {
    let malformed = |source| RecordError::Value { part: RecordPart::Format, source };
    let mut at = 0usize;
    for _ in 0..n_fmt {
        let rest = indiv.get(at..).unwrap_or_default();
        let (_, key_len) = read_typed_int(rest).map_err(malformed)?;
        let (column, len) = FormatColumn::parse(rest.get(key_len..).unwrap_or_default(), n_samples)
            .map_err(malformed)?;
        let head = key_len.saturating_add(column.header_len());
        out.extend_from_slice(rest.get(..head).unwrap_or_default());
        out.extend_from_slice(column.sample_bytes(sample).unwrap_or_default());
        at = at.saturating_add(key_len).saturating_add(len);
    }
    let len = indiv.len().saturating_sub(at);
    if len > 0 {
        return Err(RecordError::TrailingBytes { part: RecordPart::Format, len });
    }
    Ok(())
}

/// FORMAT data as a view reads it, indexed on first access.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FormatView<'a> {
    /// VCF text: the FORMAT column.
    keys: &'a [u8],
    /// BCF: the individual section. VCF text: the sample columns.
    samples: &'a [u8],
    n_fmt: usize,
    n_samples: usize,
    index: &'a LazyIndex<FormatIndex>,
    /// Copying the data failed: every access returns this.
    failed: Option<&'a RecordError>,
    file: &'a FileShared,
    at: RecordAt,
}

impl<'a> FormatView<'a> {
    /// The view of a record's FORMAT data straight from the reader.
    pub(crate) fn new(
        source: FormatSource<'a>,
        index: &'a LazyIndex<FormatIndex>,
        file: &'a FileShared,
        at: RecordAt,
    ) -> Self {
        let (keys, samples, n_fmt, n_samples) = match source {
            FormatSource::Bcf { indiv, n_fmt, n_samples } => (&[][..], indiv, n_fmt, n_samples),
            FormatSource::Text { keys, samples, n_samples } => (keys, samples, 0, n_samples),
        };
        Self { keys, samples, n_fmt, n_samples, index, failed: None, file, at }
    }

    fn index(&self) -> Result<&'a FormatIndex, VcfReadError> {
        let Self { keys, samples, n_fmt, n_samples, file, at, .. } = *self;
        if let Some(e) = self.failed {
            return Err(at.error(file, e.clone()));
        }
        self.index
            .get(|index: &mut FormatIndex, cursor: &mut Cursor, _| {
                if file.is_text() {
                    index_text_format(index, keys, samples, n_samples, file).map(|()| Scan::End)
                } else {
                    let entries = BcfEntries::new(samples, n_fmt, Some(n_samples), *cursor);
                    index_bcf(&mut index.keys, entries, cursor, None, KeyKind::Format, file)
                }
            })
            .map_err(|e| at.error(file, e.clone()))
    }

    fn malformed(&self, source: BcfValueError) -> VcfReadError {
        self.at.error(self.file, RecordError::Value { part: RecordPart::Format, source })
    }

    /// The values of the field an index entry names.
    fn field_at(&self, index: &'a FormatIndex, entry: Entry) -> Result<FieldAt<'a>, VcfReadError> {
        if self.file.is_text() {
            return Ok(FieldAt::Text {
                values: &index.values,
                samples: self.samples,
                n_keys: index.n_keys,
                k: entry.start as usize,
                n_samples: self.n_samples,
            });
        }
        let bytes = self.samples.get(entry.span()).unwrap_or_default();
        let (column, _) =
            FormatColumn::parse(bytes, self.n_samples).map_err(|e| self.malformed(e))?;
        Ok(FieldAt::Bcf(column))
    }

    /// Whether the record has field `key`.
    fn has(&self, key: &FieldId) -> Result<bool, VcfReadError> {
        check_brand(key, self.file)?;
        let Some(id) = DictId::new(key.dict_idx) else { return Ok(false) };
        Ok(self.index()?.keys.get(id).is_some())
    }

    /// Field `key`'s values, if the record has it.
    fn field(&self, key: &FieldId) -> Result<Option<FieldAt<'a>>, VcfReadError> {
        check_brand(key, self.file)?;
        let Some(id) = DictId::new(key.dict_idx) else { return Ok(None) };
        let index = self.index()?;
        index.keys.get(id).map(|e| self.field_at(index, e)).transpose()
    }
}

/// One FORMAT field's values, sample by sample.
#[derive(Debug, Clone, Copy)]
enum FieldAt<'a> {
    Bcf(FormatColumn<'a>),
    Text { values: &'a [Range<u32>], samples: &'a [u8], n_keys: usize, k: usize, n_samples: usize },
}

impl<'a> FieldAt<'a> {
    /// Sample `i`'s value, `None` past the last sample.
    fn sample(&self, i: usize) -> Option<Stored<'a>> {
        match *self {
            Self::Bcf(column) => column.sample(i).map(|v| Stored(StoredRepr::Bcf(v))),
            Self::Text { values, samples, n_keys, k, n_samples } => {
                if i >= n_samples {
                    return None;
                }
                let r = values.get(i.checked_mul(n_keys)?.checked_add(k)?)?;
                let value = samples.get(r.start as usize..r.end as usize).unwrap_or_default();
                Some(Stored(StoredRepr::Text(Some(value))))
            }
        }
    }
}

// r[impl vcf_fields.sample_idx]
/// Fail with `ForeignKey` unless `sample` is of `file`'s header.
fn check_sample(sample: Option<SampleIdx>, file: &FileShared) -> Result<(), VcfReadError> {
    match sample {
        Some(s) if s.header != file.header.id() => {
            Err(VcfReadError::ForeignKey { name: SmolStr::from(format!("sample {}", s.index())) })
        }
        _ => Ok(()),
    }
}

/// One sample's FORMAT values.
#[derive(Clone, Copy)]
pub struct SampleRef<'a> {
    view: FormatView<'a>,
    column: usize,
    /// The handle it was asked for by, checked against the record's header.
    sample: Option<SampleIdx>,
}

// r[impl vcf_reader.debug]
/// The fields present, by name.
impl std::fmt::Debug for SampleRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Sample ")?;
        debug_fields(f, self.iter())
    }
}

impl<'a> SampleRef<'a> {
    pub(crate) fn new(view: FormatView<'a>, column: usize, sample: Option<SampleIdx>) -> Self {
        Self { view, column, sample }
    }

    // r[impl vcf_values.access+2]
    /// The sample's value of the field `key` names, decoded as its marker
    /// says; absent or `.` is the marker's absent value.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key or sample of another header,
    /// and [`VcfReadError::Record`] for malformed data.
    pub fn get<V: ReadKind>(&self, key: &FormatKey<V>) -> Result<V::Out<'a>, VcfReadError> {
        let file = self.view.file;
        check_sample(self.sample, file)?;
        let Some(field) = self.view.field(&key.0)? else { return Ok(V::absent()) };
        match field.sample(self.column) {
            Some(value) => decode::<V>(value, &key.0, KeyKind::Format, file, self.view.at),
            None => Ok(V::absent()),
        }
    }

    // r[impl vcf_values.access+2]
    /// Whether the record's FORMAT column names the field `key` names,
    /// whatever this sample's value: absent is `!has`, and `.` (or a
    /// trailing field the sample dropped) is `has` with an absent
    /// [`get`](Self::get). Needs only this sample's section.
    ///
    /// # Errors
    /// As [`Self::get`].
    pub fn has<V>(&self, key: &FormatKey<V>) -> Result<bool, VcfReadError> {
        check_sample(self.sample, self.view.file)?;
        self.view.has(&key.0)
    }

    // r[impl vcf_values.dynamic+2]
    /// Every FORMAT field present, in record order (a repeated key once),
    /// with its name, its definition in the file's field table and this
    /// sample's value read as that says. FORMAT data that does not index is
    /// the first item's error.
    pub fn iter(
        &self,
    ) -> impl Iterator<Item = Result<(&'a str, &'a FormatDef, DynValue<'a>), VcfReadError>> + use<'a>
    {
        let this = *self;
        let index = check_sample(self.sample, self.view.file).and_then(|()| self.view.index());
        let (error, entries) = match index {
            Ok(index) => (None, Some(index)),
            Err(e) => (Some(Err(e)), None),
        };
        let fields = entries.into_iter().flat_map(move |index| {
            index.keys.entries().iter().map(move |&e| this.dyn_entry(index, e))
        });
        error.into_iter().chain(fields)
    }

    fn dyn_entry(
        &self,
        index: &'a FormatIndex,
        entry: Entry,
    ) -> Result<(&'a str, &'a FormatDef, DynValue<'a>), VcfReadError> {
        let (file, at) = (self.view.file, self.view.at);
        let id = i32::try_from(entry.id.get()).unwrap_or(i32::MAX);
        let unknown =
            || at.error(file, RecordError::UnknownKey { index: id, kind: KeyKind::Format });
        let (name, def) = file.fields.format_def(&file.header, entry.id).ok_or_else(unknown)?;
        let value = self.view.field_at(index, entry)?.sample(self.column).ok_or_else(unknown)?;
        let value = dyn_decode(value, name, def.typ, name == "GT", file, at)?;
        Ok((name, def, value))
    }
}

/// One FORMAT field of a record with its id in the file's field table, read
/// sample by sample as its definition says: what a passthrough copy writes.
pub(crate) struct DynColumn<'a> {
    pub(crate) id: DictId,
    pub(crate) name: &'a str,
    pub(crate) def: &'a FormatDef,
    field: FieldAt<'a>,
    file: &'a FileShared,
    at: RecordAt,
}

impl<'a> DynColumn<'a> {
    /// Sample `i`'s value; `None` past the last sample.
    pub(crate) fn sample(&self, i: usize) -> Result<Option<DynValue<'a>>, VcfReadError> {
        let Some(value) = self.field.sample(i) else { return Ok(None) };
        let genotype = self.name == "GT";
        dyn_decode(value, self.name, self.def.typ, genotype, self.file, self.at).map(Some)
    }
}

/// Every sample's FORMAT values, borrowed.
#[derive(Clone, Copy)]
pub struct SamplesRef<'a> {
    view: FormatView<'a>,
}

/// How many samples, and the fields the record has; not every value.
impl std::fmt::Debug for SamplesRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut d = f.debug_struct("Samples");
        d.field("len", &self.len());
        match self.view.index() {
            Ok(index) => {
                let names: SmallVec<&str, 8> = index
                    .keys
                    .entries()
                    .iter()
                    .filter_map(|e| self.view.file.fields.name(&self.view.file.header, e.id))
                    .collect();
                d.field("fields", &names)
            }
            Err(e) => d.field("error", &e.to_string()),
        };
        d.finish()
    }
}

impl<'a> SamplesRef<'a> {
    pub(crate) fn new(view: FormatView<'a>) -> Self {
        Self { view }
    }

    /// One sample's view. A sample of another header fails on access with
    /// [`VcfReadError::ForeignKey`].
    pub fn sample(&self, sample: SampleIdx) -> SampleRef<'a> {
        SampleRef { view: self.view, column: sample.index(), sample: Some(sample) }
    }

    // r[impl vcf_fields.samples]
    /// One field across all samples, in sample order; `None` if the record
    /// does not have the field.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key of another header;
    /// [`VcfReadError::Record`] if any sample's value is malformed (every
    /// sample is checked here, so iterating cannot fail).
    pub fn column<V: ReadKind>(
        &self,
        key: &FormatKey<V>,
    ) -> Result<Option<Column<'a, V>>, VcfReadError> {
        let Some(field) = self.view.field(&key.0)? else { return Ok(None) };
        let (file, at, n) = (self.view.file, self.view.at, self.view.n_samples);
        let checked = match field {
            FieldAt::Bcf(column) => {
                let (all, per_sample) = column.all();
                V::column_decodes(Stored(StoredRepr::Bcf(all)), per_sample)
            }
            FieldAt::Text { .. } => false,
        };
        if !checked {
            for i in 0..n {
                if let Some(value) = field.sample(i) {
                    decode::<V>(value, &key.0, KeyKind::Format, file, at)?;
                }
            }
        }
        let key = key.0.clone();
        Ok(Some(Column { field, next: 0, end: n, key, file, at, _v: PhantomData }))
    }

    /// Whether the record's FORMAT column names `key` (see [`SampleRef::has`]).
    ///
    /// # Errors
    /// As [`SampleRef::has`].
    pub fn has<V>(&self, key: &FormatKey<V>) -> Result<bool, VcfReadError> {
        self.view.has(&key.0)
    }

    /// How many samples there are.
    pub fn len(&self) -> usize {
        self.view.n_samples
    }

    /// Every FORMAT field present, in record order (a repeated key once).
    pub(crate) fn columns(
        &self,
    ) -> Result<impl Iterator<Item = Result<DynColumn<'a>, VcfReadError>> + use<'a>, VcfReadError>
    {
        let view = self.view;
        let index = view.index()?;
        Ok(index.keys.entries().iter().map(move |&entry| {
            let (file, at) = (view.file, view.at);
            let index_id = i32::try_from(entry.id.get()).unwrap_or(i32::MAX);
            let (name, def) = file.fields.format_def(&file.header, entry.id).ok_or_else(|| {
                at.error(file, RecordError::UnknownKey { index: index_id, kind: KeyKind::Format })
            })?;
            let field = view.field_at(index, entry)?;
            Ok(DynColumn { id: entry.id, name, def, field, file, at })
        }))
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// One FORMAT field across all samples, in sample order.
pub struct Column<'a, V> {
    field: FieldAt<'a>,
    next: usize,
    end: usize,
    key: FieldId,
    file: &'a FileShared,
    at: RecordAt,
    _v: PhantomData<fn() -> V>,
}

impl<V> std::fmt::Debug for Column<'_, V> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Column")
            .field("name", &self.key.name())
            .field("left", &self.end.saturating_sub(self.next))
            .finish()
    }
}

impl<V> Clone for Column<'_, V> {
    fn clone(&self) -> Self {
        Self { key: self.key.clone(), _v: PhantomData, ..*self }
    }
}

impl<'a, V: ReadKind> Column<'a, V> {
    /// Sample `sample`'s value, wherever iteration is.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a sample of another header.
    pub fn get(&self, sample: SampleIdx) -> Result<V::Out<'a>, VcfReadError> {
        check_sample(Some(sample), self.file)?;
        match self.field.sample(sample.index()) {
            Some(value) => decode::<V>(value, &self.key, KeyKind::Format, self.file, self.at),
            None => Ok(V::absent()),
        }
    }
}

impl<'a, V: ReadKind> Iterator for Column<'a, V> {
    type Item = V::Out<'a>;

    fn next(&mut self) -> Option<V::Out<'a>> {
        if self.next >= self.end {
            return None;
        }
        let value = self.field.sample(self.next)?;
        self.next = self.next.saturating_add(1);
        let decoded = decode::<V>(value, &self.key, KeyKind::Format, self.file, self.at);
        debug_assert!(decoded.is_ok(), "column() checked every sample");
        Some(decoded.unwrap_or_else(|_| V::absent()))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.end.saturating_sub(self.next);
        (n, Some(n))
    }
}

impl<V: ReadKind> ExactSizeIterator for Column<'_, V> {}

// r[impl vcf_fields.only_sample]
/// The FORMAT values of a single-sample file, read without a sample handle.
/// Fails to plan unless the header has exactly one sample: there is no
/// "first sample" reader that would silently pick one of many.
#[derive(Clone)]
pub struct OnlySample(FormatData);

impl std::fmt::Debug for OnlySample {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.view().fmt(f)
    }
}

impl OnlySample {
    fn view(&self) -> SampleRef<'_> {
        SampleRef { view: self.0.view(), column: 0, sample: None }
    }

    /// The sample's value of the field `key` names (see [`SampleRef::get`]).
    ///
    /// # Errors
    /// As [`SampleRef::get`].
    pub fn get<V: ReadKind>(&self, key: &FormatKey<V>) -> Result<V::Out<'_>, VcfReadError> {
        self.view().get(key)
    }

    /// Whether the record's FORMAT column names `key` (see [`SampleRef::has`]).
    ///
    /// # Errors
    /// As [`SampleRef::has`].
    pub fn has<V>(&self, key: &FormatKey<V>) -> Result<bool, VcfReadError> {
        self.view().has(key)
    }

    /// Every FORMAT field present (see [`SampleRef::iter`]).
    pub fn iter(
        &self,
    ) -> impl Iterator<Item = Result<(&str, &FormatDef, DynValue<'_>), VcfReadError>> {
        self.view().iter()
    }
}

impl Fields for OnlySample {
    type Plan = SampleIdx;

    fn plan(cx: &PlanContext<'_>) -> Result<SampleIdx, PlanError> {
        let mut samples = cx.header().sample_indices();
        match (samples.len(), samples.next()) {
            (1, Some(sample)) => Ok(sample),
            (n, _) => Err(PlanError::NotSingleSample { n }),
        }
    }

    fn needs(sample: &SampleIdx) -> Sections {
        Sections::sample(*sample)
    }

    fn new(sample: &SampleIdx, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::new(&raw);
        data.fill(&raw, Some(*sample))?;
        Ok(Self(data))
    }

    fn fill(&mut self, sample: &SampleIdx, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.0.fill(&raw, Some(*sample))
    }
}

// r[impl vcf_fields.sample+2]
/// The FORMAT values of the sample the reader's builder named with
/// `.sample(name)` (or of the only sample, if there is one).
#[derive(Clone)]
pub struct Sample {
    data: FormatData,
    name: SmolStr,
}

impl std::fmt::Debug for Sample {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ", self.name)?;
        self.view().fmt(f)
    }
}

impl Sample {
    fn view(&self) -> SampleRef<'_> {
        SampleRef { view: self.data.view(), column: 0, sample: None }
    }

    /// The sample's name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The sample's value of the field `key` names (see [`SampleRef::get`]).
    ///
    /// # Errors
    /// As [`SampleRef::get`].
    pub fn get<V: ReadKind>(&self, key: &FormatKey<V>) -> Result<V::Out<'_>, VcfReadError> {
        self.view().get(key)
    }

    /// Whether the record's FORMAT column names `key` (see [`SampleRef::has`]).
    ///
    /// # Errors
    /// As [`SampleRef::has`].
    pub fn has<V>(&self, key: &FormatKey<V>) -> Result<bool, VcfReadError> {
        self.view().has(key)
    }

    /// Every FORMAT field present (see [`SampleRef::iter`]).
    pub fn iter(
        &self,
    ) -> impl Iterator<Item = Result<(&str, &FormatDef, DynValue<'_>), VcfReadError>> {
        self.view().iter()
    }
}

impl Fields for Sample {
    type Plan = (SampleIdx, SmolStr);

    fn plan(cx: &PlanContext<'_>) -> Result<Self::Plan, PlanError> {
        let sample = cx.selected_sample()?;
        let name = cx.header().samples().get(sample.index()).cloned().unwrap_or_default();
        Ok((sample, name))
    }

    fn needs((sample, _): &Self::Plan) -> Sections {
        Sections::sample(*sample)
    }

    fn new((sample, name): &Self::Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::new(&raw);
        data.fill(&raw, Some(*sample))?;
        Ok(Self { data, name: name.clone() })
    }

    fn fill(&mut self, (sample, _): &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.data.fill(&raw, Some(*sample))
    }
}

/// Every sample's FORMAT values, by sample ([`sample`](Self::sample)) or by
/// field ([`column`](Self::column)).
#[derive(Clone)]
pub struct Samples(FormatData);

impl std::fmt::Debug for Samples {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.view().fmt(f)
    }
}

impl Samples {
    fn view(&self) -> SamplesRef<'_> {
        SamplesRef { view: self.0.view() }
    }

    /// One sample's view (see [`SamplesRef::sample`]).
    pub fn sample(&self, sample: SampleIdx) -> SampleRef<'_> {
        self.view().sample(sample)
    }

    /// One field across all samples (see [`SamplesRef::column`]).
    ///
    /// # Errors
    /// As [`SamplesRef::column`].
    pub fn column<'s, V: ReadKind>(
        &'s self,
        key: &FormatKey<V>,
    ) -> Result<Option<Column<'s, V>>, VcfReadError> {
        self.view().column(key)
    }

    /// Whether the record's FORMAT column names `key` (see [`SampleRef::has`]).
    ///
    /// # Errors
    /// As [`SampleRef::has`].
    pub fn has<V>(&self, key: &FormatKey<V>) -> Result<bool, VcfReadError> {
        self.view().has(key)
    }

    /// How many samples there are.
    pub fn len(&self) -> usize {
        self.0.n_samples
    }

    /// Every FORMAT field present (see [`SamplesRef::columns`]).
    pub(crate) fn columns(
        &self,
    ) -> Result<impl Iterator<Item = Result<DynColumn<'_>, VcfReadError>>, VcfReadError> {
        SamplesRef::new(self.0.view()).columns()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Fields for Samples {
    type Plan = ();

    fn plan(_: &PlanContext<'_>) -> Result<(), PlanError> {
        Ok(())
    }

    fn needs(_: &()) -> Sections {
        Sections::SAMPLES
    }

    fn new(_: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::new(&raw);
        data.fill(&raw, None)?;
        Ok(Self(data))
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.0.fill(&raw, None)
    }
}
