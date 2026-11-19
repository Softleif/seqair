//! A VCF/BCF file and its index from bytes, queried and partitioned
//! (`r[vcf_reader.fuzz+2]`, `docs/spec/2-vcf-1-reader.md` "Region queries"):
//! nothing panics, however file or index is malformed, and on a file that
//! reads whole a query yields only records the overlap rule selects.
//!
//! Layout: mode byte (even: the file is BGZF-compressed here first),
//! `u32` LE index length, the index (TBI or CSI), a region (contig byte,
//! `u32` LE start, `u32` LE length, `u32::MAX` for the contig's end; flags
//! byte: bit 0 `StartsIn`, bits 2.. the partition count), then the file.
#![no_main]

// r[verify vcf_reader.fuzz+2]

use libfuzzer_sys::fuzz_target;
use seqair::io::BgzfWriter;
use seqair::vcf::{Locus, Match, Region, VcfFile, VcfPos};
use seqair_types::Pos0;
use std::num::NonZeroUsize;

/// At most this many records per read, so a run stays fast.
const MAX_RECORDS: usize = 2000;

fn take<'a>(data: &mut &'a [u8], n: usize) -> Option<&'a [u8]> {
    let (head, rest) = data.split_at_checked(n)?;
    *data = rest;
    Some(head)
}

fn u32_le(data: &mut &[u8]) -> Option<u32> {
    take(data, 4)?.first_chunk::<4>().map(|b| u32::from_le_bytes(*b))
}

fn read(reader: &mut seqair::vcf::Reader<Locus>) -> Vec<Locus> {
    let mut out = Vec::new();
    while out.len() < MAX_RECORDS {
        match reader.next_record() {
            Ok(Some(locus)) => out.push(locus.clone()),
            _ => break,
        }
    }
    out
}

fuzz_target!(|data: &[u8]| {
    let mut data = data;
    let Some(&[mode]) = take(&mut data, 1) else { return };
    let Some(len) = u32_le(&mut data) else { return };
    let Some(index) = take(&mut data, len as usize) else { return };
    let Some(&[contig]) = take(&mut data, 1) else { return };
    let (Some(start), Some(span)) = (u32_le(&mut data), u32_le(&mut data)) else { return };
    let Some(&[flags]) = take(&mut data, 1) else { return };
    let bytes = if mode % 2 == 0 {
        let mut w = BgzfWriter::new(Vec::new());
        if w.write_all(data).is_err() {
            return;
        }
        match w.finish() {
            Ok(b) => b,
            Err(_) => return,
        }
    } else {
        data.to_vec()
    };
    let Ok(file) = VcfFile::from_bytes(bytes, Some(index)) else { return };
    let names: Vec<String> = file.header().contigs().map(|(id, _)| id.name().to_string()).collect();
    let Some(name) = names.get(usize::from(contig) % names.len().max(1)) else { return };
    let Some(start) = Pos0::new(start % (1 << 30)) else { return };
    let region = match span {
        u32::MAX => Region::new(name.as_str(), start..),
        span => match Pos0::new(start.as_u32().saturating_add(span % (1 << 20))) {
            Some(end) => Region::new(name.as_str(), start..end),
            None => return,
        },
    };
    let Ok(region) = region else { return };
    let matching = if flags & 1 == 0 { Match::Overlapping } else { Match::StartsIn };
    if let Ok(mut q) = file.reader::<Locus>().matching(matching).query(&region) {
        let got = read(&mut q);
        // Every record a query yields is on its contig; one it takes by
        // its start starts before the region ends.
        for locus in &got {
            assert_eq!(locus.contig().name(), name.as_str());
            if let (VcfPos::At(s), seqair::vcf::RegionEnd::Before(end)) =
                (locus.start(), region.end())
            {
                assert!(s.to_zero_based() < end, "{locus:?} starts past {region}");
            }
        }
        if q.requery(&region).is_ok() {
            std::hint::black_box(read(&mut q));
        }
    }
    let n = NonZeroUsize::new(usize::from(flags >> 2)).unwrap_or(NonZeroUsize::MIN);
    if let Ok(partitions) = file.partitions(n) {
        for p in partitions.iter().take(8) {
            if let Ok(mut r) = file.reader::<Locus>().partition(p) {
                std::hint::black_box(read(&mut r));
            }
        }
    }
});
