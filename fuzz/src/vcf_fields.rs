//! Reading every field of every record of a VCF/BCF file, shared by the
//! VCF fuzz targets (`r[vcf_reader.fuzz]`).

use seqair::vcf::{
    DynValue, Fields, Flag, FullRecord, Info, Locus, OptArr, PlanContext, PlanError, RawRecord,
    Sample, Sections, Str, ValueType, VcfFile, VcfHeader, VcfReadError,
};

/// At most this many records per input, so a run stays fast.
const MAX_RECORDS: usize = 2000;

fn touch(value: &DynValue<'_>) -> usize {
    match value {
        DynValue::Int(v) => v.clone().flatten().count(),
        DynValue::Float(v) => v.clone().flatten().count(),
        DynValue::Flag => 1,
        DynValue::Str(s) => s.map_or(0, str::len),
        DynValue::Gt(g) => {
            g.as_ref().map_or(0, |g| g.alleles().count() + g.to_genotype().alleles.len())
        }
    }
}

/// Every INFO and FORMAT field of the header, read through a typed key too.
fn typed(header: &VcfHeader, rec: &FullRecord) -> usize {
    let mut n = 0;
    for (name, def) in header.infos() {
        let info = rec.info();
        n += match def.typ {
            ValueType::Integer => header
                .info_key::<OptArr<i32>>(name)
                .ok()
                .and_then(|k| info.get(&k).ok())
                .flatten()
                .map_or(0, Iterator::count),
            ValueType::Float => header
                .info_key::<OptArr<f32>>(name)
                .ok()
                .and_then(|k| info.get(&k).ok())
                .flatten()
                .map_or(0, Iterator::count),
            ValueType::Flag => header
                .info_key::<Flag>(name)
                .ok()
                .and_then(|k| info.get(&k).ok())
                .map_or(0, usize::from),
            ValueType::String | ValueType::Character => header
                .info_key::<Str>(name)
                .ok()
                .and_then(|k| info.get(&k).ok())
                .flatten()
                .map_or(0, str::len),
        };
    }
    for (name, def) in header.formats() {
        let samples = rec.samples();
        n += match (name.as_str(), def.typ) {
            ("GT", _) => header.format_key(name).ok().and_then(|k: seqair::vcf::FormatGt| {
                samples.column(&k).ok().flatten().map(|c| c.flatten().count())
            }),
            (_, ValueType::Integer) => header
                .format_key::<OptArr<i32>>(name)
                .ok()
                .and_then(|k| samples.column(&k).ok().flatten().map(|c| c.flatten().count())),
            (_, ValueType::Float) => header
                .format_key::<OptArr<f32>>(name)
                .ok()
                .and_then(|k| samples.column(&k).ok().flatten().map(|c| c.flatten().count())),
            _ => header
                .format_key::<Str>(name)
                .ok()
                .and_then(|k| samples.column(&k).ok().flatten().map(|c| c.flatten().count())),
        }
        .unwrap_or(0);
    }
    n
}

/// Read every field of every record (up to [`MAX_RECORDS`]) through
/// `FullRecord`, untyped and through typed keys, and a `Locus` reader.
pub fn read_all(file: &VcfFile) -> usize {
    let header = file.header();
    let Ok(mut reader) = file.reader::<FullRecord>().open() else { return 0 };
    let mut n = 0;
    for _ in 0..MAX_RECORDS {
        let rec = match reader.next_record() {
            Ok(Some(rec)) => rec,
            Ok(None) | Err(_) => break,
        };
        let locus = rec.locus();
        n += locus.contig().name().len() + usize::from(locus.start() <= locus.end());
        let site = rec.site();
        n += site.ids().count() + site.alleles().reference().len() + site.alleles().alts().count();
        n += site.to_alleles().map_or(0, |a| a.ref_text().len());
        n += site.qual().map_or(0, |q| usize::from(q.is_nan()));
        n += match site.filters().status() {
            seqair::vcf::FilterStatus::Failed(f) => f.count(),
            _ => 1,
        };
        n += rec.info().iter().map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
        for sample in header.sample_indices() {
            let fields = rec.samples().sample(sample).iter();
            n += fields.map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
        }
        n += typed(header, rec) + rec.raw_bytes().len();
    }
    n
}

/// Everything [`RawRecord`]'s own views reach, the way a user type reads.
#[derive(Clone, Debug)]
struct Raw(usize);

impl Fields for Raw {
    type Plan = ();
    fn plan(_: &PlanContext<'_>) -> Result<(), PlanError> {
        Ok(())
    }
    fn needs(_: &()) -> Sections {
        Sections::SITE | Sections::INFO | Sections::SAMPLES
    }
    fn new(_: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let mut n = raw.locus().map_or(0, |l| usize::from(l.start() <= l.end()));
        n += raw.site()?.alleles().alts().count();
        n += raw.info()?.iter().map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
        for sample in raw.header().sample_indices().take(4) {
            n += raw
                .sample(sample)?
                .iter()
                .map(|f| f.map_or(0, |(_, _, v)| touch(&v)))
                .sum::<usize>();
        }
        Ok(Self(n))
    }
}

/// Read the records as narrower types do: `Locus` alone, `Info` with one
/// named `Sample` (which copies only its column), and through
/// [`RawRecord`]'s views.
pub fn read_narrow(file: &VcfFile) -> usize {
    let mut n = 0;
    if let Ok(reader) = file.reader::<Locus>().open() {
        n += reader.into_iter().take(MAX_RECORDS).map_while(Result::ok).count();
    }
    if let Some(name) = file.header().samples().last()
        && let Ok(reader) = file.reader::<(Info, Sample)>().sample(name.clone()).open()
    {
        for rec in reader.into_iter().take(MAX_RECORDS).map_while(Result::ok) {
            n += rec.1.iter().map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
        }
    }
    if let Ok(reader) = file.reader::<Raw>().open() {
        n +=
            reader.into_iter().take(MAX_RECORDS).filter_map(Result::ok).map(|r| r.0).sum::<usize>();
    }
    n
}
