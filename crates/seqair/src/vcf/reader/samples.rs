//! FORMAT values: one sample ([`OnlySample`], [`Sample`]) or all
//! ([`Samples`]), by sample or by column (`r[vcf_fields.only_sample]`,
//! `r[vcf_fields.sample]`, `r[vcf_fields.samples]`).

use super::access::{DecodeCx, DynValue, ReadKind, Stored, StoredRepr};
use super::bcf_value::{FormatColumn, TypedValue, read_typed_int};
use super::entries::{BcfEntries, LazyIndex, index_bcf, scan};
use super::error::{KeyKind, RecordAt, RecordError, RecordPart, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, Sections};
use super::file::FileShared;
use super::index::{DictId, Entry, FieldIndex};
use super::info::check_brand;
use super::raw::RawRecord;
use crate::vcf::header::{FormatDef, SampleIdx, VcfHeader};
use crate::vcf::record_encoder::{FieldId, FormatKey};
use seqair_types::SmolStr;
use std::marker::PhantomData;
use std::ops::Range;
use std::sync::Arc;

/// FORMAT values of some samples in BCF's layout: per field, its key, its
/// type descriptor, then each sample's elements.
#[derive(Debug, Clone, Default)]
struct FormatData {
    bytes: Vec<u8>,
    n_fmt: usize,
    n_samples: usize,
    file: Option<Arc<FileShared>>,
    at: RecordAt,
    index: LazyIndex,
}

impl FormatData {
    fn begin(&mut self, raw: &RawRecord<'_>, n_fmt: usize, n_samples: usize) {
        self.bytes.clear();
        self.n_fmt = n_fmt;
        self.n_samples = n_samples;
        self.at = raw.at();
        self.index.reset();
        raw.keep_file(&mut self.file);
    }

    fn fill_all(&mut self, raw: &RawRecord<'_>) -> Result<(), VcfReadError> {
        let (indiv, n_fmt, n_samples) = raw.bcf_samples(None)?;
        self.begin(raw, n_fmt, n_samples);
        self.bytes.extend_from_slice(indiv);
        Ok(())
    }

    /// Copy only sample `sample`'s values: each field's key and descriptor,
    /// then that sample's elements.
    fn fill_one(&mut self, raw: &RawRecord<'_>, sample: SampleIdx) -> Result<(), VcfReadError> {
        let (indiv, n_fmt, n_samples) = raw.bcf_samples(Some(sample))?;
        self.begin(raw, n_fmt, 1);
        let malformed =
            |source| raw.at().error(RecordError::Value { part: RecordPart::Format, source });
        let mut at = 0usize;
        for _ in 0..n_fmt {
            let rest = indiv.get(at..).unwrap_or_default();
            let (_, key_len) = read_typed_int(rest).map_err(malformed)?;
            let (column, len) =
                FormatColumn::parse(rest.get(key_len..).unwrap_or_default(), n_samples)
                    .map_err(malformed)?;
            let head = key_len.saturating_add(column.header_len());
            self.bytes.extend_from_slice(rest.get(..head).unwrap_or_default());
            self.bytes.extend_from_slice(column.sample_bytes(sample.index()).unwrap_or_default());
            at = at.saturating_add(key_len).saturating_add(len);
        }
        let len = indiv.len().saturating_sub(at);
        if len > 0 {
            return Err(raw
                .at()
                .error(RecordError::TrailingBytes { part: RecordPart::Format, len }));
        }
        Ok(())
    }

    fn index(&self, file: &FileShared) -> Result<&FieldIndex, VcfReadError> {
        let entries = BcfEntries::new(&self.bytes, self.n_fmt, Some(self.n_samples));
        self.index
            .get(|index| index_bcf(index, entries, KeyKind::Format, file))
            .map_err(|e| self.at.error(e.clone()))
    }

    fn source(&self) -> Source<'_> {
        Source::Indexed(self)
    }
}

/// Where a FORMAT view finds a field.
#[derive(Debug, Clone, Copy)]
enum Source<'a> {
    Indexed(&'a FormatData),
    Scan { bytes: &'a [u8], n_fmt: usize, n_samples: usize, file: &'a FileShared, at: RecordAt },
}

impl<'a> Source<'a> {
    fn bytes(&self) -> &'a [u8] {
        match *self {
            Self::Indexed(data) => &data.bytes,
            Self::Scan { bytes, .. } => bytes,
        }
    }

    fn n_samples(&self) -> usize {
        match *self {
            Self::Indexed(data) => data.n_samples,
            Self::Scan { n_samples, .. } => n_samples,
        }
    }

    fn file(&self) -> Option<&'a FileShared> {
        match *self {
            Self::Indexed(data) => data.file.as_deref(),
            Self::Scan { file, .. } => Some(file),
        }
    }

    fn at(&self) -> RecordAt {
        match *self {
            Self::Indexed(data) => data.at,
            Self::Scan { at, .. } => at,
        }
    }

    fn malformed(&self, source: super::bcf_value::BcfValueError) -> VcfReadError {
        self.at().error(RecordError::Value { part: RecordPart::Format, source })
    }

    /// The span of field `id`'s column, if the record has it.
    fn span(&self, file: &FileShared, id: &FieldId) -> Result<Option<Range<usize>>, VcfReadError> {
        let Some(id) = DictId::new(id.dict_idx) else { return Ok(None) };
        match *self {
            Self::Indexed(data) => Ok(data.index(file)?.get(id).map(Entry::span)),
            Self::Scan { bytes, n_fmt, n_samples, .. } => {
                scan(BcfEntries::new(bytes, n_fmt, Some(n_samples)), id)
                    .map_err(|e| self.malformed(e))
            }
        }
    }

    /// Field `key`'s column, if the record has it.
    fn column(
        &self,
        file: &FileShared,
        key: &FieldId,
    ) -> Result<Option<FormatColumn<'a>>, VcfReadError> {
        check_brand(key, file)?;
        let Some(span) = self.span(file, key)? else { return Ok(None) };
        let bytes = self.bytes().get(span).unwrap_or_default();
        let (column, _) =
            FormatColumn::parse(bytes, self.n_samples()).map_err(|e| self.malformed(e))?;
        Ok(Some(column))
    }

    /// The fields present, each as its id and column span.
    fn entries(&self, file: &'a FileShared) -> Result<Entries<'a>, VcfReadError> {
        Ok(match *self {
            Self::Indexed(data) => Entries::Indexed(data.index(file)?.entries().iter()),
            Self::Scan { bytes, n_fmt, n_samples, .. } => {
                Entries::Scan(BcfEntries::new(bytes, n_fmt, Some(n_samples)))
            }
        })
    }
}

enum Entries<'a> {
    Indexed(std::slice::Iter<'a, Entry>),
    Scan(BcfEntries<'a>),
}

/// Decode one sample's elements as `V`.
fn decode<'a, V: ReadKind>(
    value: TypedValue<'a>,
    name: &SmolStr,
    file: &FileShared,
    at: RecordAt,
) -> Result<V::Out<'a>, VcfReadError> {
    let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
    V::decode(Stored(StoredRepr::Bcf(value)), cx)
        .map_err(|source| at.error(RecordError::Field { name: name.clone(), source }))
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
#[derive(Debug, Clone, Copy)]
pub struct SampleRef<'a> {
    source: Source<'a>,
    column: usize,
    /// The handle it was asked for by, checked against the record's header.
    sample: Option<SampleIdx>,
}

impl<'a> SampleRef<'a> {
    // r[impl vcf_values.access]
    /// The sample's value of the field `key` names, decoded as its marker
    /// says; absent or `.` is the marker's absent value.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a key or sample of another header,
    /// and [`VcfReadError::Record`] for malformed data.
    pub fn get<V: ReadKind>(&self, key: &FormatKey<V>) -> Result<V::Out<'a>, VcfReadError> {
        let Some(file) = self.source.file() else { return Ok(V::absent()) };
        check_sample(self.sample, file)?;
        let Some(column) = self.source.column(file, &key.0)? else { return Ok(V::absent()) };
        match column.sample(self.column) {
            Some(value) => decode::<V>(value, &key.0.name, file, self.source.at()),
            None => Ok(V::absent()),
        }
    }

    // r[impl vcf_values.dynamic]
    /// Every declared FORMAT field present, in record order (from an owned
    /// component, a repeated key once), with its declaration and this
    /// sample's value read as declared.
    ///
    /// # Errors
    /// [`VcfReadError::Record`] if the FORMAT data is malformed; each item,
    /// if its value is.
    pub fn iter(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<(&'a str, &'a FormatDef, DynValue<'a>), VcfReadError>> + use<'a>,
        VcfReadError,
    > {
        let file = self.source.file();
        if let Some(file) = file {
            check_sample(self.sample, file)?;
        }
        let entries = match file {
            Some(file) => Some(self.source.entries(file)?),
            None => None,
        };
        let this = *self;
        let mut entries = entries;
        Ok(std::iter::from_fn(move || {
            let file = file?;
            let (index, span) = match entries.as_mut()? {
                Entries::Indexed(it) => {
                    let e = it.next()?;
                    (i32::try_from(e.id.get()).unwrap_or(i32::MAX), e.span())
                }
                Entries::Scan(it) => match it.next()? {
                    Ok(entry) => entry,
                    Err(e) => return Some(Err(this.source.malformed(e))),
                },
            };
            Some(this.dyn_entry(file, index, span))
        }))
    }

    fn dyn_entry(
        &self,
        file: &'a FileShared,
        index: i32,
        span: Range<usize>,
    ) -> Result<(&'a str, &'a FormatDef, DynValue<'a>), VcfReadError> {
        let at = self.source.at();
        let unknown = || at.error(RecordError::UnknownKey { index, kind: KeyKind::Format });
        let name = usize::try_from(index)
            .ok()
            .and_then(|i| file.header.string_map().name(i))
            .ok_or_else(unknown)?;
        let (name, def) = file.header.formats().get_key_value(name).ok_or_else(unknown)?;
        let bytes = self.source.bytes().get(span).unwrap_or_default();
        let (column, _) = FormatColumn::parse(bytes, self.source.n_samples())
            .map_err(|e| self.source.malformed(e))?;
        let value = column.sample(self.column).ok_or_else(unknown)?;
        let cx = DecodeCx { first_phase_is_read: file.first_phase_is_read };
        let value = DynValue::decode(Stored(StoredRepr::Bcf(value)), def.typ, name == "GT", cx)
            .map_err(|source| at.error(RecordError::Field { name: name.clone(), source }))?;
        Ok((name.as_str(), def, value))
    }
}

/// Every sample's FORMAT values, borrowed.
#[derive(Debug, Clone, Copy)]
pub struct SamplesRef<'a> {
    source: Source<'a>,
}

impl<'a> SamplesRef<'a> {
    /// The samples of a BCF record, walking its FORMAT entries on each access.
    pub(crate) fn scan(
        bytes: &'a [u8],
        n_fmt: usize,
        n_samples: usize,
        file: &'a FileShared,
        at: RecordAt,
    ) -> Self {
        Self { source: Source::Scan { bytes, n_fmt, n_samples, file, at } }
    }

    /// One sample's view. A sample of another header fails on access with
    /// [`VcfReadError::ForeignKey`].
    pub fn sample(&self, sample: SampleIdx) -> SampleRef<'a> {
        SampleRef { source: self.source, column: sample.index(), sample: Some(sample) }
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
        let Some(file) = self.source.file() else { return Ok(None) };
        let Some(column) = self.source.column(file, &key.0)? else { return Ok(None) };
        let at = self.source.at();
        let n = self.source.n_samples();
        for i in 0..n {
            if let Some(value) = column.sample(i) {
                decode::<V>(value, &key.0.name, file, at)?;
            }
        }
        let name = key.0.name.clone();
        Ok(Some(Column { column, next: 0, end: n, name, file, at, _v: PhantomData }))
    }

    /// How many samples there are.
    pub fn len(&self) -> usize {
        self.source.n_samples()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// One FORMAT field across all samples, in sample order.
#[derive(Debug)]
pub struct Column<'a, V> {
    column: FormatColumn<'a>,
    next: usize,
    end: usize,
    name: SmolStr,
    file: &'a FileShared,
    at: RecordAt,
    _v: PhantomData<fn() -> V>,
}

impl<V> Clone for Column<'_, V> {
    fn clone(&self) -> Self {
        Self { name: self.name.clone(), _v: PhantomData, ..*self }
    }
}

impl<'a, V: ReadKind> Column<'a, V> {
    /// Sample `sample`'s value, wherever iteration is.
    ///
    /// # Errors
    /// [`VcfReadError::ForeignKey`] for a sample of another header.
    pub fn get(&self, sample: SampleIdx) -> Result<V::Out<'a>, VcfReadError> {
        check_sample(Some(sample), self.file)?;
        match self.column.sample(sample.index()) {
            Some(value) => decode::<V>(value, &self.name, self.file, self.at),
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
        let value = self.column.sample(self.next)?;
        self.next = self.next.saturating_add(1);
        let decoded = decode::<V>(value, &self.name, self.file, self.at);
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
#[derive(Debug, Clone)]
pub struct OnlySample(FormatData);

impl OnlySample {
    fn view(&self) -> SampleRef<'_> {
        SampleRef { source: self.0.source(), column: 0, sample: None }
    }

    /// The sample's value of the field `key` names (see [`SampleRef::get`]).
    ///
    /// # Errors
    /// As [`SampleRef::get`].
    pub fn get<V: ReadKind>(&self, key: &FormatKey<V>) -> Result<V::Out<'_>, VcfReadError> {
        self.view().get(key)
    }

    /// Every FORMAT field present (see [`SampleRef::iter`]).
    ///
    /// # Errors
    /// As [`SampleRef::iter`].
    pub fn iter(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<(&str, &FormatDef, DynValue<'_>), VcfReadError>>,
        VcfReadError,
    > {
        self.view().iter()
    }
}

impl Fields for OnlySample {
    type Plan = SampleIdx;

    fn plan(header: &VcfHeader, _: &ReadOptions) -> Result<SampleIdx, PlanError> {
        let mut samples = header.sample_indices();
        match (samples.len(), samples.next()) {
            (1, Some(sample)) => Ok(sample),
            (n, _) => Err(PlanError::NotSingleSample { n }),
        }
    }

    fn needs(sample: &SampleIdx) -> Sections {
        Sections::sample(*sample)
    }

    fn new(sample: &SampleIdx, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::default();
        data.fill_one(&raw, *sample)?;
        Ok(Self(data))
    }

    fn fill(&mut self, sample: &SampleIdx, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.0.fill_one(&raw, *sample)
    }
}

// r[impl vcf_fields.sample]
/// The FORMAT values of the sample the reader's builder named with
/// `.sample(name)` (or of the only sample, if there is one).
#[derive(Debug, Clone)]
pub struct Sample {
    data: FormatData,
    name: SmolStr,
}

impl Sample {
    fn view(&self) -> SampleRef<'_> {
        SampleRef { source: self.data.source(), column: 0, sample: None }
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

    /// Every FORMAT field present (see [`SampleRef::iter`]).
    ///
    /// # Errors
    /// As [`SampleRef::iter`].
    pub fn iter(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<(&str, &FormatDef, DynValue<'_>), VcfReadError>>,
        VcfReadError,
    > {
        self.view().iter()
    }
}

impl Fields for Sample {
    type Plan = (SampleIdx, SmolStr);

    fn plan(header: &VcfHeader, opts: &ReadOptions) -> Result<Self::Plan, PlanError> {
        let sample = match opts.sample() {
            Some(name) => header
                .sample(name)
                .map_err(|_| PlanError::UnknownSample { name: SmolStr::from(name) })?,
            None => {
                let mut samples = header.sample_indices();
                match (samples.len(), samples.next()) {
                    (1, Some(sample)) => sample,
                    (n, _) => return Err(PlanError::SampleNotSelected { n }),
                }
            }
        };
        let name = header.samples().get(sample.index()).cloned().unwrap_or_default();
        Ok((sample, name))
    }

    fn needs((sample, _): &Self::Plan) -> Sections {
        Sections::sample(*sample)
    }

    fn new((sample, name): &Self::Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::default();
        data.fill_one(&raw, *sample)?;
        Ok(Self { data, name: name.clone() })
    }

    fn fill(&mut self, (sample, _): &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.data.fill_one(&raw, *sample)
    }
}

/// Every sample's FORMAT values, by sample ([`sample`](Self::sample)) or by
/// field ([`column`](Self::column)).
#[derive(Debug, Clone)]
pub struct Samples(FormatData);

impl Samples {
    fn view(&self) -> SamplesRef<'_> {
        SamplesRef { source: self.0.source() }
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

    /// How many samples there are.
    pub fn len(&self) -> usize {
        self.0.n_samples
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Fields for Samples {
    type Plan = ();

    fn plan(_: &VcfHeader, _: &ReadOptions) -> Result<(), PlanError> {
        Ok(())
    }

    fn needs(_: &()) -> Sections {
        Sections::SAMPLES
    }

    fn new(_: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut data = FormatData::default();
        data.fill_all(&raw)?;
        Ok(Self(data))
    }

    fn fill(&mut self, _: &(), raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        self.0.fill_all(&raw)
    }
}
