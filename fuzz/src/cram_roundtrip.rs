//! Write records as CRAM with seqair's writer, decode them with seqair's
//! stream reader, and check each against what was written, under the
//! normalizations `r[cram_writer.fidelity]` lists. Shared by the CRAM
//! write fuzz targets.

use seqair::bam::RawBamRecord;
use seqair::bam::RecordStore;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::cram::stream::CramStreamReader;
use seqair::cram::writer::{CramVersion, CramWriteError, CramWriterBuilder};
use seqair::fasta::IndexedFastaReader;
use seqair::sam::format::format_sam_record;
use seqair_types::IupacBase;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// A FASTA of `contigs` (name, bases) with its `.fai`, in `dir`.
pub fn write_fasta(dir: &Path, contigs: &[(String, Vec<u8>)]) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut fasta = Vec::new();
    let mut fai = String::new();
    for (name, bases) in contigs {
        fasta.push(b'>');
        fasta.extend_from_slice(name.as_bytes());
        fasta.push(b'\n');
        let offset = fasta.len();
        fasta.extend_from_slice(bases);
        fasta.push(b'\n');
        let n = bases.len();
        writeln!(fai, "{name}\t{n}\t{offset}\t{n}\t{}", n + 1).unwrap();
    }
    std::fs::write(&path, fasta).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), fai).unwrap();
    path
}

/// How a record should read back.
#[derive(Debug, PartialEq)]
struct Expected {
    qname: Vec<u8>,
    flags: u16,
    tid: Option<u32>,
    pos: Option<u32>,
    mapq: u8,
    cigar: Vec<(CigarOpType, u32)>,
    seq: Vec<IupacBase>,
    qual: Option<Vec<u8>>,
    mate: (Option<u32>, Option<u32>, i32),
    /// Tags other than `MD`, `NM` and `RG`, in order.
    tags: Vec<Vec<u8>>,
    rg: Option<Vec<u8>>,
}

/// One aux tag's head and value bytes, and the rest.
fn split_tags(aux: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut out = Vec::new();
    for tag in seqair::bam::aux::Aux::new(aux).iter_checked() {
        let Ok((name, value)) = tag else { break };
        // NM may be regenerated as another integer type: compare its value.
        let text = match (&name, value.as_i64()) {
            (b"NM", Some(n)) => format!("{n}"),
            _ => format!("{value:?}"),
        };
        out.push((name.to_vec(), text.into_bytes()));
    }
    out
}

/// `=`/`X` as `M`, zero-length operations dropped, neighbours of a kind
/// merged.
fn normalized_cigar(cigar: impl Iterator<Item = CigarOp>) -> Vec<(CigarOpType, u32)> {
    let mut out: Vec<(CigarOpType, u32)> = Vec::new();
    for op in cigar.filter(|op| op.len() > 0) {
        let ty = match op.op_type() {
            CigarOpType::SeqMatch | CigarOpType::SeqMismatch => CigarOpType::Match,
            ty => ty,
        };
        match out.last_mut() {
            Some((last, len)) if *last == ty => *len += op.len(),
            _ => out.push((ty, op.len())),
        }
    }
    out
}

fn expected(r: &RawBamRecord, normalize: bool) -> (Expected, Option<Vec<u8>>, Option<Vec<u8>>) {
    let unmapped = r.flags().is_unmapped();
    let mut tags = Vec::new();
    let (mut md, mut nm, mut rg) = (None, None, None);
    for (name, value) in split_tags(r.aux().as_bytes()) {
        match &name[..] {
            b"MD" => md = md.or(Some(value)),
            b"NM" => nm = nm.or(Some(value)),
            b"RG" if rg.is_none() => rg = Some(value),
            _ => tags.push([name, value].concat()),
        }
    }
    let mut seq = Vec::new();
    r.iupac_bases_into(&mut seq);
    let mate = r.mate_fields();
    let paired = r.flags().is_paired();
    let e = Expected {
        qname: r.qname().to_vec(),
        flags: r.flags().raw(),
        tid: r.tid().map(|t| t.as_u32()),
        pos: r.pos().map(|p| p.as_u32()),
        mapq: if unmapped && normalize { 0 } else { r.mapq().as_byte() },
        cigar: if unmapped && normalize {
            Vec::new()
        } else if normalize {
            normalized_cigar(r.cigar().iter())
        } else {
            r.cigar().iter().map(|op| (op.op_type(), op.len())).collect()
        },
        seq,
        qual: r.qual().map(|q| seqair_types::BaseQuality::slice_to_bytes(q).to_vec()),
        mate: (
            mate.tid.filter(|_| paired || !normalize).map(|t| t.as_u32()),
            mate.pos.map(|p| p.as_u32()),
            mate.template_len.map_or(0, i32::from),
        ),
        tags,
        rg,
    };
    (e, md, nm)
}

/// Write `records` (as `write` or `write_raw` would take them; here as raw
/// records) to a CRAM against `fasta`. Returns the file and the records the
/// writer accepted, in order.
pub fn write(
    header: &BamHeader,
    fasta: &Path,
    records: &[RawBamRecord],
    threads: usize,
    records_per_slice: u32,
    version: CramVersion,
) -> (Vec<u8>, Vec<RawBamRecord>) {
    let reference = IndexedFastaReader::open(fasta).unwrap();
    let mut writer = match CramWriterBuilder::to_writer(Vec::new(), header, reference)
        .compression_threads(threads)
        .records_per_slice(records_per_slice)
        .version(version)
        .write_index(true)
        .build()
    {
        Ok(w) => w,
        // A header the FASTA disagrees with.
        Err(_) => return (Vec::new(), Vec::new()),
    };
    let mut accepted = Vec::new();
    for r in records {
        match writer.write_raw(r) {
            Ok(()) => accepted.push(r.clone()),
            Err(CramWriteError::Io { .. } | CramWriteError::Poisoned) => {
                panic!("a Vec sink cannot fail")
            }
            // Refused: the writer stays usable.
            Err(_) => {}
        }
    }
    let (cram, index) = writer.finish().unwrap();
    // The index is well-formed and covers every slice it lists.
    let mut crai = Vec::new();
    index.unwrap().write_crai(&mut crai).unwrap();
    (cram, accepted)
}

/// The records of `cram`, decoded by seqair against `fasta`.
pub fn read(cram: &[u8], fasta: &Path) -> (BamHeader, Vec<RawBamRecord>) {
    let mut reader = CramStreamReader::new(cram).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap();
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    (reader.into_header(), records)
}

/// Write, read back, compare each record with what was written; then write
/// what was read and read it again, which must give the same records (the
/// normalizations are a fixed point). Every record read formats as SAM.
pub fn check(
    header: &BamHeader,
    fasta: &Path,
    records: &[RawBamRecord],
    threads: usize,
    rps: u32,
    version: CramVersion,
) {
    let (cram, accepted) = write(header, fasta, records, threads, rps, version);
    if cram.is_empty() {
        return;
    }
    let (cram_header, decoded) = read(&cram, fasta);
    assert_eq!(decoded.len(), accepted.len(), "every accepted record reads back");
    for (i, (input, output)) in accepted.iter().zip(&decoded).enumerate() {
        let (want, want_md, want_nm) = expected(input, true);
        let (got, got_md, got_nm) = expected(output, false);
        assert_eq!(got, want, "record {i}");
        // A stored MD/NM reads back as it was; a missing one may be
        // generated.
        if want_md.is_some() {
            assert_eq!(got_md, want_md, "record {i} MD");
        }
        if want_nm.is_some() {
            assert_eq!(got_nm, want_nm, "record {i} NM");
        }
        let mut sam = Vec::new();
        format_sam_record(output, &cram_header, &mut sam).unwrap();
    }

    store_writes_as_records(&cram_header, fasta, &decoded, threads, rps, version);

    let (again, accepted_again) = write(&cram_header, fasta, &decoded, threads, rps, version);
    assert_eq!(accepted_again.len(), decoded.len(), "decoded records are all writable");
    let (_, decoded_again) = read(&again, fasta);
    for (i, (a, b)) in decoded.iter().zip(&decoded_again).enumerate() {
        assert_eq!(a.as_payload(), b.as_payload(), "record {i} is not a fixed point");
    }
}

/// `write_store` of the placed records writes what `write_store_record`
/// writes for each in turn, in two stores so the first leaves a slice open.
fn store_writes_as_records(
    header: &BamHeader,
    fasta: &Path,
    records: &[RawBamRecord],
    threads: usize,
    rps: u32,
    version: CramVersion,
) {
    let mut stores = [RecordStore::new(), RecordStore::new()];
    let half = records.len() / 2;
    for (i, r) in records.iter().enumerate() {
        let Some(store) = stores.get_mut(usize::from(i >= half)) else { continue };
        if store.push_raw(r.as_payload(), header, &mut ()).is_err() {
            return;
        }
    }
    let write = |by_store: bool, threads: usize| {
        let reference = IndexedFastaReader::open(fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), header, reference)
            .compression_threads(threads)
            .records_per_slice(rps)
            .version(version)
            .write_index(true)
            .build()
            .unwrap();
        for store in &stores {
            if by_store {
                writer.write_store(store).unwrap();
            } else {
                for idx in store.indices() {
                    writer.write_store_record(store, idx).unwrap();
                }
            }
        }
        let (cram, index) = writer.finish().unwrap();
        let mut crai = Vec::new();
        index.unwrap().write_crai(&mut crai).unwrap();
        (cram, crai)
    };
    assert!(write(true, threads) == write(false, 0), "write_store differs");
}
