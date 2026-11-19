//! Passthrough writing (`docs/spec/5-vcf-passthrough.md`): the records a
//! reader yields, copied into a [`Writer`] through a [`FieldMap`] and open
//! for edits as a [`RecordEdit`]. The reading guide's last section
//! ([`crate::vcf::reader`], "Writing what you read") shows a conversion and
//! a rewrite of GT and PS.
//!
//! A copy is re-encoded from the reader's values, not copied as bytes, so
//! every output follows the writer's rules (and is what bcftools writes for
//! the same records), and editing a copied field is encoding it again.

use super::error::{FieldMapError, KeyKind, VcfError};
use super::header::{Number, ValueType, VcfHeader};
use super::reader::{DictId, DynValue, FilterStatus, FullRecord, VcfFile};
use super::record::Genotype;
use super::record_encoder::{ContigId, FieldId, FilterId, FormatEncoder, InfoEncoder};
use super::unified::{Escape, Ready, RecordEdit, Writer};
use seqair_types::{SmallVec, SmolStr};
use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use tracing::warn;

/// Where an input id's field goes in the output.
#[derive(Debug, Clone)]
enum Target {
    /// To this output field. Its `Number` is `.`: copied values are not
    /// counted (`r[vcf_passthrough.copy.values]`).
    Copy(FieldId),
    /// Nowhere: the input header declares it, the output removed it.
    Removed,
    /// Not defined in this section by the input: looked up by name.
    ByName,
}

// r[impl vcf_passthrough.field_map]
/// The translation from an input file's fields, filters and contigs to an
/// output header's, built and checked once ([`FieldMap::new`]) and used by
/// [`Writer::copy_record`] for every record.
///
/// Fields are matched by name in their section: the output header may be the
/// input's ([`VcfFile::header`]), or one derived from it with
/// [`VcfHeader::to_builder`] that removes, adds or replaces definitions.
/// A definition the output removed is left out of every copy.
#[derive(Debug)]
pub struct FieldMap {
    input: VcfFile,
    output: Arc<VcfHeader>,
    /// By input dictionary id (the header's and, for VCF text, the reserved).
    info: Vec<Target>,
    format: Vec<Target>,
    /// By input dictionary id: the output FILTER, `None` where removed.
    filters: Vec<Option<FilterId>>,
    /// By input contig tid (the header's).
    contigs: Vec<Option<ContigId>>,
    drop_undeclared: bool,
    warned: AtomicBool,
}

impl FieldMap {
    // r[impl vcf_passthrough.field_map]
    // r[impl vcf_passthrough.field_map.compatible]
    // r[impl vcf_passthrough.field_map.contigs]
    // r[impl vcf_passthrough.field_map.samples]
    /// The map from `input`'s fields to `output`'s.
    ///
    /// # Errors
    /// [`FieldMapError::Incompatible`] for a field both define with another
    /// `Type` or `Number`, [`FieldMapError::MissingContig`] for an input
    /// contig the output lacks, and [`FieldMapError::SampleMismatch`] unless
    /// both have the same samples in the same order.
    pub fn new(input: &VcfFile, output: &VcfHeader) -> Result<Self, FieldMapError> {
        Self::build(input, Arc::new(output.clone()))
    }

    fn build(input: &VcfFile, output: Arc<VcfHeader>) -> Result<Self, FieldMapError> {
        let in_header = input.header();
        // r[impl vcf_passthrough.field_map.removed+2]
        // A name removed from a header derived from the input's is a new
        // field if declared again, never the input's.
        let replaced = |kind, name: &str| output.was_removed_since(in_header.id(), kind, name);
        let (ins, outs) = (in_header.samples(), output.samples());
        if ins != outs {
            let index = ins.iter().zip(outs).take_while(|(a, b)| a == b).count();
            return Err(FieldMapError::SampleMismatch {
                index,
                input: ins.get(index).cloned(),
                output: outs.get(index).cloned(),
            });
        }
        // By tid, holes included: a contig dictionary may skip indices (`IDX=`).
        let mut contigs = vec![None; in_header.contig_map().len()];
        for (contig, _) in in_header.contigs() {
            let out = output
                .contig(contig.name())
                .map_err(|_| FieldMapError::MissingContig { name: contig.name().clone() })?;
            if let Some(slot) = contigs.get_mut(contig.tid() as usize) {
                *slot = Some(out);
            }
        }
        let table = &input.0.fields;
        let ids = (0..table.base()).filter_map(|i| DictId::new(u32::try_from(i).ok()?));
        let (mut info, mut format) = (Vec::new(), Vec::new());
        for id in ids {
            let declared = |kind| table.declared(kind, id).is_some();
            info.push(match table.info_def(in_header, id) {
                Some((name, def)) if declared(KeyKind::Info) => {
                    let declared = in_header.infos().contains_key(name);
                    if replaced(KeyKind::Info, name) {
                        Target::Removed
                    } else {
                        target(&output, KeyKind::Info, name, def.typ, def.number, declared)?
                    }
                }
                _ => Target::ByName,
            });
            format.push(match table.format_def(in_header, id) {
                Some((name, def)) if declared(KeyKind::Format) => {
                    let declared = in_header.formats().contains_key(name);
                    if replaced(KeyKind::Format, name) {
                        Target::Removed
                    } else {
                        target(&output, KeyKind::Format, name, def.typ, def.number, declared)?
                    }
                }
                _ => Target::ByName,
            });
        }
        let mut filters = vec![None; in_header.string_map().len()];
        for name in in_header.filters().keys() {
            if let (Ok(input_id), Ok(out)) = (in_header.filter(name), output.filter(name))
                && !replaced(KeyKind::Filter, name)
                && let Some(slot) = filters.get_mut(input_id.dict_idx() as usize)
            {
                *slot = Some(out);
            }
        }
        Ok(Self {
            input: input.clone(),
            output,
            info,
            format,
            filters,
            contigs,
            drop_undeclared: false,
            warned: AtomicBool::new(false),
        })
    }

    // r[impl vcf_passthrough.undeclared+2]
    /// Leave the keys and filters neither the input header nor the output
    /// header declares out of every copy, with one `warn!` naming the first,
    /// instead of failing the copy with [`VcfError::UndeclaredField`]. A
    /// record that would then fail no filter still fails the copy so: it
    /// must not pass for filters nobody declared.
    #[must_use]
    pub fn drop_undeclared(mut self) -> Self {
        self.drop_undeclared = true;
        self
    }

    /// The input file's header.
    pub fn input(&self) -> &VcfHeader {
        self.input.header()
    }

    /// The output header the map writes for.
    pub fn output(&self) -> &VcfHeader {
        &self.output
    }

    // r[impl vcf_passthrough.undeclared+2]
    /// What to do with a field or filter no header declares: `Ok(None)` drops it.
    fn undeclared<T>(&self, kind: KeyKind, name: &str) -> Result<Option<T>, VcfError> {
        if !self.drop_undeclared {
            return Err(VcfError::UndeclaredField { kind, name: SmolStr::from(name) });
        }
        if !self.warned.swap(true, Ordering::Relaxed) {
            warn!(
                %kind,
                name,
                "dropping a field neither header declares from the copies (further ones are not \
                 reported)"
            );
        }
        Ok(None)
    }

    /// The output field input field `id` (named `name`) of `kind` goes to,
    /// `None` if it is left out.
    fn field(
        &self,
        kind: KeyKind,
        id: DictId,
        name: &str,
    ) -> Result<Option<Cow<'_, FieldId>>, VcfError> {
        let table = if kind == KeyKind::Info { &self.info } else { &self.format };
        match table.get(id.get() as usize) {
            Some(Target::Copy(out)) => return Ok(Some(Cow::Borrowed(out))),
            Some(Target::Removed) => return Ok(None),
            Some(Target::ByName) | None => {}
        }
        // r[impl vcf_passthrough.undeclared+2]
        // Neither header defines it here: the reader read it as a String.
        let def = match kind {
            KeyKind::Info => self.output.infos().get(name).map(|d| (d.typ, d.number)),
            _ => self.output.formats().get(name).map(|d| (d.typ, d.number)),
        };
        let Some((typ, number)) = def else { return self.undeclared(kind, name) };
        if !matches!(typ, ValueType::String | ValueType::Character) {
            return Err(FieldMapError::Incompatible {
                kind,
                name: SmolStr::from(name),
                input_type: ValueType::String,
                input_number: Number::Unknown,
                output_type: typ,
                output_number: number,
            }
            .into());
        }
        Ok(output_id(&self.output, kind, name).map(Cow::Owned))
    }

    // r[impl vcf_passthrough.field_map.contigs]
    fn contig(&self, contig: &ContigId) -> Result<ContigId, VcfError> {
        if let Some(Some(out)) = self.contigs.get(contig.tid() as usize) {
            return Ok(out.clone());
        }
        self.output
            .contig(contig.name())
            .map_err(|_| VcfError::UndeclaredContig { name: contig.name().clone() })
    }

    // r[impl vcf_passthrough.field_map.removed+2]
    // r[impl vcf_passthrough.undeclared+2]
    /// The output filters a record failing `failed` fails.
    fn failed_filters(
        &self,
        record: &FullRecord,
        failed: super::reader::FailedFilters<'_>,
    ) -> Result<SmallVec<FilterId, 4>, VcfError> {
        let mut out = SmallVec::new();
        // `FailedFilters` leaves out the names the header does not declare:
        // then read the column as written.
        let filters = record.site().filters();
        let Some(names) = filters.names_with_undeclared() else {
            out.extend(failed.filter_map(|f| self.filters.get(f.dict_idx() as usize)?.clone()));
            return Ok(out);
        };
        let mut dropped = None;
        for name in names {
            let utf8 = std::str::from_utf8(name).ok();
            let mapped = match utf8.map(|n| (n, self.input().filter(n))) {
                Some((_, Ok(id))) => self.filters.get(id.dict_idx() as usize).cloned().flatten(),
                Some((n, Err(_))) if let Ok(id) = self.output.filter(n) => Some(id),
                _ => {
                    let name = utf8.map_or_else(
                        || SmolStr::from(name.escape_ascii().to_string()),
                        SmolStr::from,
                    );
                    let mapped = self.undeclared(KeyKind::Filter, &name)?;
                    dropped.get_or_insert(name);
                    mapped
                }
            };
            out.extend(mapped);
        }
        // r[impl vcf_passthrough.undeclared+2]
        // Dropping every filter it failed would make it pass for filters
        // nobody declared.
        match dropped {
            Some(name) if out.is_empty() => {
                Err(VcfError::UndeclaredField { kind: KeyKind::Filter, name })
            }
            _ => Ok(out),
        }
    }
}

/// The output field `name` of `kind`, with `Number=.` so a copy is not counted.
fn output_id(output: &VcfHeader, kind: KeyKind, name: &str) -> Option<FieldId> {
    let id = match kind {
        KeyKind::Info => output.info_key::<super::reader::Any>(name).ok()?.0,
        _ => output.format_key::<super::reader::Any>(name).ok()?.0,
    };
    Some(FieldId { number: Number::Unknown, ..id })
}

/// Where the input field `name` of `kind`, defined with `typ` and `number`,
/// goes; `declared`: the input header declares it (else VCF reserves it).
fn target(
    output: &VcfHeader,
    kind: KeyKind,
    name: &str,
    typ: ValueType,
    number: Number,
    declared: bool,
) -> Result<Target, FieldMapError> {
    let def = match kind {
        KeyKind::Info => output.infos().get(name).map(|d| (d.typ, d.number)),
        _ => output.formats().get(name).map(|d| (d.typ, d.number)),
    };
    let Some((output_type, output_number)) = def else {
        return Ok(if declared { Target::Removed } else { Target::ByName });
    };
    if (output_type, output_number) != (typ, number) {
        return Err(FieldMapError::Incompatible {
            kind,
            name: SmolStr::from(name),
            input_type: typ,
            input_number: number,
            output_type,
            output_number,
        });
    }
    Ok(output_id(output, kind, name).map_or(Target::Removed, Target::Copy))
}

/// Buffers a copy gathers a field's values in, kept by the writer and
/// reused across records.
#[derive(Debug, Default)]
pub(crate) struct CopyScratch {
    ints: Vec<Option<i32>>,
    floats: Vec<Option<f32>>,
    /// Where each sample's values end in `ints` or `floats`.
    ends: Vec<usize>,
    gts: Vec<Genotype>,
}

/// The header a writer wrote, and the map its last implicit copy used.
#[derive(Debug, Default)]
pub(crate) struct ImplicitMap {
    output: Option<Arc<VcfHeader>>,
    map: Option<Arc<FieldMap>>,
}

impl ImplicitMap {
    pub(crate) fn new(output: &VcfHeader) -> Self {
        Self { output: Some(Arc::new(output.clone())), map: None }
    }
}

impl<W: Write> Writer<W, Ready> {
    // r[impl vcf_passthrough.implicit_map+2]
    /// Begin a copy of `record` into this writer, whose header is the
    /// record's (a conversion) or was derived from it with
    /// [`VcfHeader::to_builder`]: the [`FieldMap`] is built from the two on
    /// the first record and kept. See [`copy_record_with`](Self::copy_record_with)
    /// for what a copy writes and how to edit it.
    ///
    /// # Errors
    /// [`VcfError::ForeignRecord`] for a record of an unrelated header (build
    /// a [`FieldMap`] for those), [`VcfError::FieldMap`] if the map cannot be
    /// built, and what [`copy_record_with`](Self::copy_record_with) fails with.
    pub fn copy_record(&mut self, record: &FullRecord) -> Result<RecordEdit<'_>, VcfError> {
        let implicit = &mut self.record.implicit;
        let map = match &implicit.map {
            Some(map) if map.input().id() == record.header_id() => Arc::clone(map),
            _ => {
                let output = implicit.output.clone().ok_or(VcfError::HeaderNotWritten)?;
                let id = record.header_id();
                if output.id() != id && !output.derives_from(id) {
                    return Err(VcfError::ForeignRecord);
                }
                let map = Arc::new(FieldMap::build(&record.file(), output)?);
                implicit.map = Some(Arc::clone(&map));
                map
            }
        };
        self.copy_record_with(&map, record)
    }

    // r[impl vcf_passthrough.copy]
    // r[impl vcf_passthrough.copy.equivalence]
    /// Begin a copy of `record`, read from `map`'s input, translated to this
    /// writer's header, `map`'s output: every column, INFO and FORMAT field,
    /// re-encoded as the writer encodes any record. The returned
    /// [`RecordEdit`] takes edits — `key.encode(&mut edit, value)` replaces a
    /// field where it stands or appends it, `remove_info`/`remove_format`
    /// drop one — and writes the record on `emit()`.
    ///
    /// # Errors
    /// [`VcfError::ForeignFieldMap`] unless this writer wrote `map`'s output
    /// header, [`VcfError::ForeignRecord`] unless `record` was read with
    /// `map`'s input header, [`VcfError::TelomereRecord`] for a record at
    /// POS 0, [`VcfError::UndeclaredField`] and [`VcfError::UndeclaredContig`]
    /// for what no header declares, [`VcfError::Read`] for a malformed value,
    /// and what encoding the values can fail with.
    pub fn copy_record_with(
        &mut self,
        map: &FieldMap,
        record: &FullRecord,
    ) -> Result<RecordEdit<'_>, VcfError> {
        if self.header != Some(map.output.id()) {
            return Err(VcfError::ForeignFieldMap);
        }
        if record.header_id() != map.input().id() {
            return Err(VcfError::ForeignRecord);
        }
        let locus = record.locus();
        let contig = map.contig(locus.contig())?;
        let pos = locus
            .start()
            .pos1()
            .ok_or_else(|| VcfError::TelomereRecord { contig: contig.name().clone() })?;
        let site = record.site();
        let alleles = site.to_alleles()?;
        let enc = self.begin_record(&contig, pos, &alleles, site.qual())?.ids(site.ids())?;
        let enc = match site.filters().status() {
            FilterStatus::Pass => enc.filter_pass(),
            FilterStatus::Missing => enc.no_filter(),
            FilterStatus::Failed(failed) => {
                let failed = map.failed_filters(record, failed)?;
                // r[impl vcf_passthrough.field_map.removed+2]
                if failed.is_empty() { enc.filter_pass() } else { enc.filter_fail(&failed)? }
            }
        };
        let mut enc = enc.begin_samples().into_edit();
        for field in record.info().iter_ids() {
            let (id, name, _, value) = field?;
            if let Some(out) = map.field(KeyKind::Info, id, name)? {
                copy_info(&mut enc, &out, value)?;
            }
        }
        if !map.output.samples().is_empty() {
            for column in record.samples().columns()? {
                let column = column?;
                if let Some(out) = map.field(KeyKind::Format, column.id, column.name)? {
                    copy_format(&mut enc, &out, &column, record.samples().len())?;
                }
            }
        }
        Ok(enc)
    }
}

// r[impl vcf_passthrough.copy.encoder]
// r[impl vcf_passthrough.copy.values]
fn copy_info(enc: &mut RecordEdit<'_>, id: &FieldId, value: DynValue<'_>) -> Result<(), VcfError> {
    let mut scratch = std::mem::take(&mut enc.record.copy);
    let result = match value {
        DynValue::Flag => enc.info_flag(id),
        DynValue::Str(s) => enc.info_str(id, s.unwrap_or(".").as_bytes(), Escape::Structure(b';')),
        DynValue::Int(values) => {
            scratch.ints.clear();
            scratch.ints.extend(values);
            if scratch.ints.is_empty() {
                scratch.ints.push(None);
            }
            enc.info_int_opts(id, &scratch.ints)
        }
        DynValue::Float(values) => {
            scratch.floats.clear();
            scratch.floats.extend(values);
            if scratch.floats.is_empty() {
                scratch.floats.push(None);
            }
            enc.info_float_opts(id, &scratch.floats)
        }
        // INFO is never decoded as a genotype.
        DynValue::Gt(_) => Ok(()),
    };
    enc.record.copy = scratch;
    result
}

// r[impl vcf_passthrough.copy.encoder]
// r[impl vcf_passthrough.copy.values]
fn copy_format(
    enc: &mut RecordEdit<'_>,
    id: &FieldId,
    column: &super::reader::DynColumn<'_>,
    n_samples: usize,
) -> Result<(), VcfError> {
    let mut s = std::mem::take(&mut enc.record.copy);
    let result = copy_format_with(enc, id, column, n_samples, &mut s);
    enc.record.copy = s;
    result
}

fn copy_format_with(
    enc: &mut RecordEdit<'_>,
    id: &FieldId,
    column: &super::reader::DynColumn<'_>,
    n: usize,
    s: &mut CopyScratch,
) -> Result<(), VcfError> {
    if column.name == "GT" {
        s.gts.clear();
        for i in 0..n {
            s.gts.push(match column.sample(i)? {
                Some(DynValue::Gt(Some(gt))) => gt.to_genotype(),
                _ => Genotype { alleles: SmallVec::new(), phased: SmallVec::new() },
            });
        }
        return enc.format_gt(id, &s.gts);
    }
    match column.def.typ {
        ValueType::Integer => {
            gather(column, n, &mut s.ints, &mut s.ends, |v, out| {
                if let DynValue::Int(v) = v {
                    out.extend(v);
                }
            })?;
            enc.format_int_opts(id, &slices(&s.ints, &s.ends))
        }
        ValueType::Float => {
            gather(column, n, &mut s.floats, &mut s.ends, |v, out| {
                if let DynValue::Float(v) = v {
                    out.extend(v);
                }
            })?;
            enc.format_float_opts(id, &slices(&s.floats, &s.ends))
        }
        ValueType::String | ValueType::Character => {
            let mut strings: SmallVec<&str, 8> = SmallVec::with_capacity(n);
            for i in 0..n {
                strings.push(match column.sample(i)? {
                    Some(DynValue::Str(Some(v))) => v,
                    _ => "",
                });
            }
            enc.format_str(id, &strings, Escape::Structure(b':'))
        }
        // A header declares no FORMAT Flag.
        ValueType::Flag => Ok(()),
    }
}

/// Every sample's elements into `values`, each sample's end into `ends`; a
/// sample without any is one missing element.
fn gather<T: Copy>(
    column: &super::reader::DynColumn<'_>,
    n: usize,
    values: &mut Vec<Option<T>>,
    ends: &mut Vec<usize>,
    elements: impl Fn(DynValue<'_>, &mut Vec<Option<T>>),
) -> Result<(), VcfError> {
    values.clear();
    ends.clear();
    for i in 0..n {
        let start = values.len();
        if let Some(v) = column.sample(i)? {
            elements(v, values);
        }
        if values.len() == start {
            values.push(None);
        }
        ends.push(values.len());
    }
    Ok(())
}

/// The per-sample slices of `values` that `ends` delimits.
fn slices<'v, T>(values: &'v [T], ends: &[usize]) -> SmallVec<&'v [T], 8> {
    let starts = std::iter::once(0).chain(ends.iter().copied());
    starts.zip(ends).map(|(start, &end)| values.get(start..end).unwrap_or_default()).collect()
}
