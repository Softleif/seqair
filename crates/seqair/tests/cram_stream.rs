//! `CramStreamReader` and `CramContainerDecoder` against samtools and the
//! indexed CRAM reader.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use seqair::bam::stream::BamStreamReader;
use seqair::bam::{RawBamRecord, RecordStore};
use seqair::cram::container::{ContainerHeader, RefSeqId};
use seqair::cram::reader::IndexedCramReader;
use seqair::cram::stream::{CramStreamError, CramStreamReader};
use seqair::reader::{AlignmentStream, AlignmentStreamError, ResolveTid, StreamFormat};
use seqair::sam::format::format_sam_record;
use seqair_types::Pos0;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data")).join(name)
}

/// Every CRAM fixture with its reference.
fn fixtures() -> Vec<(PathBuf, PathBuf)> {
    let test_fa = data("test.fasta.gz");
    let fidelity_fa = data("cram-stream/fidelity.fa");
    let ce = data("hts-specs/cram/ce.fa");
    let mut out = vec![
        (data("test.cram"), test_fa.clone()),
        (data("test_v30.cram"), test_fa.clone()),
        (data("test_gzip.cram"), test_fa),
        (data("unplaced_multiref.cram"), data("unplaced_multiref.fa")),
    ];
    for f in [
        "fidelity_v31",
        "fidelity_v30",
        "fidelity_slices",
        "fidelity_multiref",
        "fidelity_lossy_names",
        "fidelity_embed_ref2",
        "fidelity_no_ref",
        "onlyu",
    ] {
        out.push((data(&format!("cram-stream/{f}.cram")), fidelity_fa.clone()));
    }
    for level in 1..=4 {
        out.push((data(&format!("hts-specs/cram/3.1/passed/level-{level}.cram")), ce.clone()));
    }
    out
}

/// A reader that hands out at most `chunk` bytes per call, like a pipe, and
/// has no `Seek`.
struct Trickle<R> {
    inner: R,
    chunk: usize,
}

impl<R: Read> Read for Trickle<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = buf.len().min(self.chunk);
        self.inner.read(&mut buf[..n])
    }
}

fn stream_records(cram: &Path, fasta: &Path, chunk: usize) -> (Vec<RawBamRecord>, usize) {
    stream_records_with(cram, fasta, chunk, true)
}

fn stream_records_with(
    cram: &Path,
    fasta: &Path,
    chunk: usize,
    decode_md: bool,
) -> (Vec<RawBamRecord>, usize) {
    let file = std::fs::File::open(cram).unwrap();
    let mut reader = CramStreamReader::new(Trickle { inner: file, chunk }).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap().decode_md(decode_md);
    let mut records = Vec::new();
    let mut containers = 0;
    while let Some(container) = reader.next_container().unwrap() {
        containers += 1;
        let before = records.len();
        decoder.decode(&container, &mut records).unwrap();
        assert_eq!((records.len() - before) as u64, container.num_records());
    }
    assert!(reader.next_container().unwrap().is_none(), "the end stays the end");
    (records, containers)
}

fn sam_lines(records: &[RawBamRecord], cram: &Path) -> Vec<String> {
    let file = std::fs::File::open(cram).unwrap();
    let header = CramStreamReader::new(file).unwrap().into_header();
    records
        .iter()
        .map(|r| {
            let mut out = Vec::new();
            format_sam_record(r, &header, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        })
        .collect()
}

fn samtools(args: &[&str], cram: &Path, fasta: &Path) -> String {
    let out = Command::new("samtools")
        .args(args)
        .arg("--reference")
        .arg(fasta)
        .arg(cram)
        .output()
        .expect("samtools on PATH");
    assert!(out.status.success(), "samtools {args:?} {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

// r[verify cram.stream]
// r[verify cram.stream.container_header]
// r[verify cram.stream.framing]
// r[verify cram.stream.decoder]
// r[verify cram.stream.all_records]
// r[verify cram.stream.reference_window]
// r[verify cram.stream.raw_records]
/// Records that differ from samtools, as (fixture, the read name samtools
/// prints, why); an empty fixture is every fixture, an empty name every
/// record (`r[cram.stream.fidelity+2]`).
/// Each entry must differ somewhere, so fixing a gap means removing it.
const KNOWN_MISMATCHES: &[(&str, &str, &str)] = &[(
    "fidelity_lossy_names",
    "",
    "RN=false: detached records lose their stored names, attached ones get no generated \
         name, and the attached halfUnmapped pair gets TLEN ±20 where htslib zeroes it",
)];

fn known(fixture: &str, qname: &str) -> Option<usize> {
    KNOWN_MISMATCHES
        .iter()
        .position(|(f, q, _)| (f.is_empty() || *f == fixture) && (q.is_empty() || *q == qname))
}

/// samtools' records for `cram` as BAM payloads, via `samtools view -b`.
fn samtools_payloads(cram: &Path, fasta: &Path, args: &[&str]) -> Vec<Vec<u8>> {
    let out = Command::new("samtools")
        .args(["view", "--no-PG", "-b", "--reference"])
        .arg(fasta)
        .args(args)
        .arg(cram)
        .output()
        .expect("samtools on PATH");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mut reader = BamStreamReader::new(&out.stdout[..]).unwrap();
    let mut payloads = Vec::new();
    loop {
        let mut payload = Vec::new();
        if !reader.read_payload_into(&mut payload).unwrap() {
            break;
        }
        payloads.push(payload);
    }
    payloads
}

// r[verify cram.stream.fidelity+2]
// r[verify cram.stream.raw_records]
// r[verify cram.stream.md_nm]
// r[verify cram.stream.md_nm.option]
// r[verify cram.stream.iupac]
/// Each fixture, streamed through a reader that returns 1–7 bytes per call,
/// gives the records `samtools view -b` writes, byte for byte (`bin`
/// included), and formats as `samtools view` prints them; with
/// `decode_md(false)` the same holds against
/// `--input-fmt-option decode_md=0`. The exceptions are [`KNOWN_MISMATCHES`].
#[test]
fn streamed_records_match_samtools() {
    let mut seen_known = std::collections::BTreeSet::new();
    let mut unknown = Vec::new();
    for (cram, fasta) in fixtures() {
        let name = cram.file_stem().unwrap().display().to_string();
        let chunk = 1 + cram.as_os_str().len() % 7;
        let count = samtools(&["view", "-c"], &cram, &fasta);
        for (decode_md, args) in
            [(true, &[][..]), (false, &["--input-fmt-option", "decode_md=0"][..])]
        {
            let (records, _) = stream_records_with(&cram, &fasta, chunk, decode_md);
            let ours = sam_lines(&records, &cram);
            assert_eq!(ours.len(), count.trim().parse::<usize>().unwrap(), "{name}");
            let theirs = samtools(&[&["view", "--no-PG"][..], args].concat(), &cram, &fasta);
            let theirs: Vec<&str> = theirs.lines().collect();
            let payloads = samtools_payloads(&cram, &fasta, args);
            assert_eq!(theirs.len(), ours.len(), "{name}");
            assert_eq!(payloads.len(), ours.len(), "{name}");

            for (i, record) in records.iter().enumerate() {
                if ours[i] == theirs[i] && record.as_payload() == payloads[i].as_slice() {
                    continue;
                }
                let qname = theirs[i].split('\t').next().unwrap_or_default();
                match known(&name, qname) {
                    Some(entry) => {
                        seen_known.insert(entry);
                    }
                    None => unknown.push(format!(
                        "{name} record {i} (decode_md {decode_md}):\n ours     {}\n samtools {}",
                        ours[i], theirs[i]
                    )),
                }
            }
        }
    }
    assert!(unknown.is_empty(), "records differ from samtools:\n{}", unknown.join("\n"));
    for (entry, (fixture, qname, why)) in KNOWN_MISMATCHES.iter().enumerate() {
        assert!(seen_known.contains(&entry), "{fixture} {qname} ({why}) matches samtools now");
    }
}

/// The multi-reference fixture really has multi-reference containers, and
/// the slices one has several.
#[test]
fn fixtures_have_the_layouts_they_are_named_for() {
    let containers = |name: &str| -> Vec<ContainerHeader> {
        let bytes = std::fs::read(data(name)).unwrap();
        let mut at = 26;
        let mut out = Vec::new();
        while at < bytes.len() {
            let h = ContainerHeader::parse(&bytes[at..]).unwrap();
            at += h.header_size + h.length as usize;
            out.push(h);
        }
        out
    };
    let multi = containers("cram-stream/fidelity_multiref.cram");
    assert!(multi.iter().any(|c| c.ref_seq_id == RefSeqId::Multi));
    assert!(
        multi.iter().any(|c| c.ref_seq_id == RefSeqId::Unmapped || c.ref_seq_id == RefSeqId::Multi)
    );
    let sliced = containers("cram-stream/fidelity_slices.cram");
    assert!(sliced.iter().filter(|c| c.num_records > 0).count() > 2);
    assert!(sliced.iter().any(|c| c.landmarks.len() > 1));
}

/// The mapped records of every reference equal what the indexed reader
/// returns for a whole-reference query, in the same order. The indexed
/// reader generates no MD/NM (`r[cram.stream.md_nm.option]`).
#[test]
fn mapped_records_match_the_indexed_reader() {
    for (cram, fasta) in fixtures() {
        if !cram.with_extension("cram.crai").exists()
            && !PathBuf::from(format!("{}.crai", cram.display())).exists()
        {
            continue;
        }
        let (records, _) = stream_records_with(&cram, &fasta, 4096, false);
        let mut indexed = IndexedCramReader::open(&cram, &fasta).unwrap();
        let header = indexed.header().clone();
        let mut store = RecordStore::new();
        for t in 0..u32::try_from(header.target_count()).unwrap() {
            let tid = t.resolve_tid(&header).unwrap();
            indexed.fetch_into(tid, (Pos0::ZERO..=Pos0::MAX).into(), &mut store).unwrap();
            let streamed: Vec<&RawBamRecord> =
                records.iter().filter(|r| r.tid() == Some(tid) && r.pos().is_some()).collect();
            let fetched: Vec<_> = store.indices().map(|i| store.record(i).unwrap()).collect();
            assert_eq!(streamed.len(), fetched.len(), "{} tid {t}", cram.display());
            for (s, f) in streamed.iter().zip(&fetched) {
                assert_eq!(s.qname(), f.qname(), "{}", cram.display());
                assert_eq!(s.pos(), Some(f.pos));
                assert_eq!(s.flags(), f.flags);
                assert_eq!(s.mapq(), f.mapq);
                assert_eq!(s.mate_fields(), f.mate_fields);
                assert_eq!(s.cigar().collect::<Vec<_>>(), f.cigar());
                assert_eq!(s.bases(), f.seq());
                assert_eq!(s.aux().as_bytes(), f.aux());
            }
        }
    }
}

/// `AlignmentStream` detects CRAM, decodes it with the reference it was
/// given, and names the problem without one.
#[test]
fn alignment_stream_reads_cram() {
    let (cram, fasta) = (data("cram-stream/fidelity_slices.cram"), data("cram-stream/fidelity.fa"));
    let (want, _) = stream_records(&cram, &fasta, 4096);

    let file = std::fs::File::open(&cram).unwrap();
    let mut stream = AlignmentStream::builder(Trickle { inner: file, chunk: 3 })
        .reference(&fasta)
        .build()
        .unwrap();
    assert_eq!(stream.format(), StreamFormat::Cram);
    let mut got = Vec::new();
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record).unwrap() {
        got.push(record.as_payload().to_vec());
    }
    let want: Vec<_> = want.iter().map(|r| r.as_payload().to_vec()).collect();
    assert_eq!(got, want);

    let file = std::fs::File::open(&cram).unwrap();
    assert!(matches!(
        AlignmentStream::from_reader(file),
        Err(AlignmentStreamError::CramNeedsReference)
    ));
}

// r[verify cram.stream.framing]
/// A CRAM cut inside a container is a typed error that poisons the stream;
/// one cut cleanly before the EOF container ends with a warning.
#[test]
fn truncation_is_typed() {
    let cram = data("cram-stream/fidelity_slices.cram");
    let bytes = std::fs::read(&cram).unwrap();
    // Container boundaries, so a cut can land on one or inside one.
    let mut bounds = vec![26usize];
    while *bounds.last().unwrap() < bytes.len() {
        let at = *bounds.last().unwrap();
        let h = ContainerHeader::parse(&bytes[at..]).unwrap();
        bounds.push(at + h.header_size + h.length as usize);
    }
    let eof_start = bounds[bounds.len() - 2];

    let read_all = |cut: usize| -> Result<usize, CramStreamError> {
        let mut reader = CramStreamReader::new(&bytes[..cut])?;
        let mut n = 0;
        loop {
            match reader.next_container() {
                Ok(Some(_)) => n += 1,
                Ok(None) => return Ok(n),
                Err(e) => {
                    assert!(matches!(reader.next_container(), Err(CramStreamError::Poisoned)));
                    return Err(e);
                }
            }
        }
    };
    let whole = read_all(bytes.len()).unwrap();
    assert_eq!(read_all(eof_start).unwrap(), whole, "a missing EOF container is not an error");
    for &b in &bounds[2..bounds.len() - 1] {
        assert!(read_all(b).is_ok());
        for cut in [b - 1, b - 20, b + 3] {
            if cut > bounds[1] && cut < eof_start {
                assert!(
                    matches!(read_all(cut), Err(CramStreamError::TruncatedContainer { .. })),
                    "cut at {cut}"
                );
            }
        }
    }
    assert!(matches!(
        CramStreamReader::new(&bytes[..20]),
        Err(CramStreamError::TruncatedFileDefinition)
    ));
    assert!(matches!(
        CramStreamReader::new(&bytes[..bounds[1] - 1]),
        Err(CramStreamError::TruncatedContainer { container: 0 })
    ));
}

// r[verify cram.stream.container]
/// Containers decode on other threads, in any order, to the same records.
#[test]
fn containers_decode_on_worker_threads() {
    let (cram, fasta) = (data("test.cram"), data("test.fasta.gz"));
    let (want, containers) = stream_records(&cram, &fasta, 1 << 16);
    assert!(containers > 1);

    let mut reader = CramStreamReader::new(std::fs::File::open(&cram).unwrap()).unwrap();
    let decoder = reader.decoder(&fasta).unwrap();
    let mut all = Vec::new();
    while let Some(c) = reader.next_container().unwrap() {
        all.push(c);
    }
    let decoded: Vec<Vec<RawBamRecord>> = std::thread::scope(|s| {
        let handles: Vec<_> = all
            .iter()
            .rev()
            .map(|c| {
                let mut d = decoder.fork().unwrap();
                s.spawn(move || {
                    let mut out = Vec::new();
                    d.decode(c, &mut out).unwrap();
                    out
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).rev().collect()
    });
    let got: Vec<_> = decoded.iter().flatten().map(|r| r.as_payload()).collect();
    let want: Vec<_> = want.iter().map(|r| r.as_payload()).collect();
    assert_eq!(got, want);
}
