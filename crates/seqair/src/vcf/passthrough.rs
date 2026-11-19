//! Passthrough writing: copying the records a reader yields into a
//! [`Writer`], translated through a [`FieldMap`] and open for edits as a
//! [`RecordEdit`] (`docs/spec/5-vcf-passthrough.md`).

use super::error::{FieldMapError, KeyKind, VcfError};
use super::header::{Number, ValueType, VcfHeader};
use super::reader::{DictId, DynValue, FilterStatus, FullRecord, VcfFile};
use super::record::Genotype;
use super::record_encoder::{ContigId, FieldId, FilterId, FormatEncoder, InfoEncoder};
use super::unified::{Edit, Escape, Ready, RecordEdit, Writer};
use seqair_types::{SmallVec, SmolStr};
use std::borrow::Cow;
use std::io::Write;
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
    output: VcfHeader,
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
        let in_header = input.header();
        let (ins, outs) = (in_header.samples(), output.samples());
        if ins != outs {
            let index = ins.iter().zip(outs).take_while(|(a, b)| a == b).count();
            return Err(FieldMapError::SampleMismatch {
                index,
                input: ins.get(index).cloned(),
                output: outs.get(index).cloned(),
            });
        }
        let contigs = in_header
            .contigs()
            .map(|(contig, _)| {
                output.contig(contig.name()).map(Some).map_err(|_| FieldMapError::MissingContig {
                    name: contig.name().clone(),
                })
            })
            .collect::<Result<_, _>>()?;
        let table = &input.0.fields;
        let ids = (0..table.base()).filter_map(|i| DictId::new(u32::try_from(i).ok()?));
        let (mut info, mut format) = (Vec::new(), Vec::new());
        for id in ids {
            let declared = |kind| table.declared(kind, id).is_some();
            info.push(match table.info_def(in_header, id) {
                Some((name, def)) if declared(KeyKind::Info) => {
                    let declared = in_header.infos().contains_key(name);
                    target(output, KeyKind::Info, name, def.typ, def.number, declared)?
                }
                _ => Target::ByName,
            });
            format.push(match table.format_def(in_header, id) {
                Some((name, def)) if declared(KeyKind::Format) => {
                    let declared = in_header.formats().contains_key(name);
                    target(output, KeyKind::Format, name, def.typ, def.number, declared)?
                }
                _ => Target::ByName,
            });
        }
        let mut filters = vec![None; in_header.string_map().len()];
        for name in in_header.filters().keys() {
            if let (Ok(input_id), Ok(out)) = (in_header.filter(name), output.filter(name))
                && let Some(slot) = filters.get_mut(input_id.dict_idx() as usize)
            {
                *slot = Some(out);
            }
        }
        Ok(Self {
            input: input.clone(),
            output: output.clone(),
            info,
            format,
            filters,
            contigs,
            drop_undeclared: false,
            warned: AtomicBool::new(false),
        })
    }

    // r[impl vcf_passthrough.undeclared]
    /// Leave the keys and filters neither the input header nor the output
    /// header declares out of every copy, with one `warn!` naming the first,
    /// instead of failing the copy with [`VcfError::UndeclaredField`].
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

    // r[impl vcf_passthrough.undeclared]
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
        // r[impl vcf_passthrough.undeclared]
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

    // r[impl vcf_passthrough.field_map.removed]
    // r[impl vcf_passthrough.undeclared]
    /// The output filters a record failing `failed` fails.
    fn failed_filters(
        &self,
        record: &FullRecord,
        failed: super::reader::FailedFilters<'_>,
    ) -> Result<SmallVec<FilterId, 4>, VcfError> {
        let mut out = SmallVec::new();
        if !self.input.0.is_text() {
            out.extend(failed.filter_map(|f| self.filters.get(f.dict_idx() as usize)?.clone()));
            return Ok(out);
        }
        // `FailedFilters` leaves out the names the header does not declare:
        // read the column as written.
        let column = record.raw_bytes().split(|&b| b == b'\t').nth(6).unwrap_or_default();
        for name in column.split(|&b| b == b';') {
            let name = std::str::from_utf8(name).map_err(|_| VcfError::UndeclaredField {
                kind: KeyKind::Filter,
                name: SmolStr::from(name.escape_ascii().to_string()),
            })?;
            let mapped = match self.input().filter(name) {
                Ok(id) => self.filters.get(id.dict_idx() as usize).cloned().flatten(),
                Err(_) => match self.output.filter(name) {
                    Ok(id) => Some(id),
                    Err(_) => self.undeclared(KeyKind::Filter, name)?,
                },
            };
            out.extend(mapped);
        }
        Ok(out)
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

impl<W: Write> Writer<W, Ready> {
    // r[impl vcf_passthrough.copy]
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
    pub fn copy_record(
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
                // r[impl vcf_passthrough.field_map.removed]
                if failed.is_empty() {
                    enc.filter_pass()
                } else {
                    enc.filter_fail(&failed)?
                }
            }
        };
        let mut enc = enc.begin_samples().into_state::<Edit>();
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
        DynValue::Str(s) => enc.info_str(id, s.unwrap_or("."), Escape::Structure(b';')),
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
