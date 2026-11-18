//! A VCF/BCF file opened from bytes and every field of every record read
//! through `FullRecord` (`r[vcf_reader.fuzz]`): nothing panics, however the
//! file is malformed.
//!
//! The first byte picks the framing: even, the rest is BGZF-compressed here
//! first (so mutations reach the records, not the deflate stream); odd, the
//! rest is the file as is.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::io::BgzfWriter;
use seqair::vcf::{DynValue, Flag, FullRecord, OptArr, Str, ValueType, VcfFile, VcfHeader};

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

fn read_all(file: &VcfFile) -> usize {
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
        if let Ok(fields) = rec.info().iter() {
            n += fields.map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
        }
        for sample in header.sample_indices() {
            if let Ok(fields) = rec.samples().sample(sample).iter() {
                n += fields.map(|f| f.map_or(0, |(_, _, v)| touch(&v))).sum::<usize>();
            }
        }
        n += typed(header, rec) + rec.raw_bytes().len();
    }
    n
}

fuzz_target!(|data: &[u8]| {
    let Some((&mode, rest)) = data.split_first() else { return };
    let bytes = if mode % 2 == 0 {
        let mut w = BgzfWriter::new(Vec::new());
        if w.write_all(rest).is_err() {
            return;
        }
        match w.finish() {
            Ok(b) => b,
            Err(_) => return,
        }
    } else {
        rest.to_vec()
    };
    if let Ok(file) = VcfFile::from_bytes(bytes, None) {
        std::hint::black_box(read_all(&file));
    }
});
