//! `CramWriter` against seqair's CRAM reader and samtools.
#![cfg(feature = "cram")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_lossless,
    clippy::format_collect,
    clippy::range_plus_one,
    clippy::manual_is_multiple_of,
    reason = "test code"
)]

use seqair::bam::RawBamRecord;
use seqair::bam::stream::BamStreamReader;
use seqair::cram::stream::CramStreamReader;
use seqair::cram::writer::{CramProfile, CramVersion, CramWriterBuilder};
use seqair::fasta::IndexedFastaReader;
use seqair::io::CompressionLevel;
use seqair::sam::format::format_sam_record;
use seqair_types::BaseQuality;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data").join(name)
}

fn samtools(args: &[&str]) -> String {
    let out = Command::new("samtools").args(args).output().expect("samtools on PATH");
    assert!(out.status.success(), "samtools {args:?}: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// The records `samtools view` prints for `path`.
fn samtools_lines(path: &Path, fasta: &Path) -> Vec<String> {
    samtools(&["view", "--no-PG", "-T", fasta.to_str().unwrap(), path.to_str().unwrap()])
        .lines()
        .map(String::from)
        .collect()
}

/// The records seqair's stream reader decodes from `cram`, as SAM lines.
fn seqair_lines(cram: &Path, fasta: &Path) -> Vec<String> {
    let mut reader = CramStreamReader::new(std::fs::File::open(cram).unwrap()).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap();
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    let header = reader.into_header();
    records
        .iter()
        .map(|r| {
            let mut out = Vec::new();
            format_sam_record(r, &header, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        })
        .collect()
}

/// A record's core fields, as two readers can be compared on them.
#[derive(Debug, PartialEq)]
struct Core {
    name: Vec<u8>,
    flags: u16,
    tid: Option<usize>,
    pos: Option<usize>,
    mapq: Option<u8>,
    cigar: Vec<(char, usize)>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    mate_tid: Option<usize>,
    mate_pos: Option<usize>,
    tlen: i32,
}

/// The records of `cram` as noodles reads them, a reader with nothing of
/// htslib's in it.
fn noodles_cores(cram: &Path, fasta: &Path) -> Vec<Core> {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta)
        .expect("open indexed FASTA");
    let repository = noodles::fasta::Repository::new(
        noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
    );
    let mut reader = noodles::cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram)
        .expect("noodles opens the CRAM");
    let header = reader.read_header().expect("noodles reads the header");
    reader
        .records(&header)
        .map(|r| {
            let r = r.expect("noodles reads every record");
            let unmapped = r.flags().is_unmapped();
            Core {
                name: r.name().map(|n| n.to_vec()).unwrap_or_default(),
                flags: u16::from(r.flags()),
                tid: r.reference_sequence_id(),
                pos: r.alignment_start().map(|p| usize::from(p) - 1),
                mapq: (!unmapped).then(|| r.mapping_quality().map_or(255, u8::from)),
                cigar: r
                    .cigar()
                    .as_ref()
                    .iter()
                    .map(|op| {
                        let c = match op.kind() {
                            Kind::Match => 'M',
                            Kind::Insertion => 'I',
                            Kind::Deletion => 'D',
                            Kind::Skip => 'N',
                            Kind::SoftClip => 'S',
                            Kind::HardClip => 'H',
                            Kind::Pad => 'P',
                            Kind::SequenceMatch => '=',
                            Kind::SequenceMismatch => 'X',
                        };
                        (c, op.len())
                    })
                    .collect(),
                // noodles copies a soft-masked reference's lower case into
                // SEQ; seqair and htslib upper-case it.
                seq: r.sequence().as_ref().to_ascii_uppercase(),
                qual: r.quality_scores().as_ref().to_vec(),
                // htslib and seqair drop the mate reference of an unpaired
                // record (`r[cram_writer.fidelity]`); noodles keeps it.
                mate_tid: r.mate_reference_sequence_id().filter(|_| r.flags().is_segmented()),
                mate_pos: r.mate_alignment_start().map(|p| usize::from(p) - 1),
                tlen: r.template_length(),
            }
        })
        .collect()
}

/// The same fields of seqair's records, and whether one runs past its
/// reference's end.
fn seqair_cores(cram: &Path, fasta: &Path) -> (Vec<Core>, bool) {
    let mut reader = CramStreamReader::new(std::fs::File::open(cram).unwrap()).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap();
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    let header = reader.into_header();
    let past_end = records.iter().any(|r| {
        let (Some(tid), Some(end)) = (r.tid(), r.end_pos()) else { return false };
        !r.flags().is_unmapped() && end.as_u64() >= header.target_len(tid).unwrap_or(u64::MAX)
    });
    let cores = records
        .iter()
        .map(|r| {
            let mut seq = Vec::new();
            r.decode_seq_into(&mut seq);
            let mate = r.mate_fields();
            let unmapped = r.flags().is_unmapped();
            Core {
                name: r.qname().to_vec(),
                flags: r.flags().raw(),
                tid: r.tid().map(|t| t.as_usize()),
                pos: r.pos().map(|p| p.as_usize()),
                mapq: (!unmapped).then(|| r.mapq().as_byte()),
                cigar: r
                    .cigar()
                    .iter()
                    .map(|op| {
                        (char::from(b"MIDNSHP=X"[usize::from(op.op_code())]), op.len() as usize)
                    })
                    .collect(),
                seq,
                qual: r.qual().map(|q| BaseQuality::slice_to_bytes(q).to_vec()).unwrap_or_default(),
                mate_tid: mate.tid.map(|t| t.as_usize()),
                mate_pos: mate.pos.map(|p| p.as_usize()),
                tlen: mate.template_len.map_or(0, i32::from),
            }
        })
        .collect();
    (cores, past_end)
}

/// noodles reads `cram` as seqair does, record for record; but noodles
/// (0.102) panics on a read that runs past its reference's end, so a file
/// with one is not given to it.
fn assert_noodles_reads_it_as_seqair(cram: &Path, fasta: &Path) {
    let (ours, past_end) = seqair_cores(cram, fasta);
    if past_end {
        return;
    }
    let theirs = noodles_cores(cram, fasta);
    assert_eq!(ours.len(), theirs.len(), "noodles reads as many records");
    for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
        assert_eq!(a, b, "record {i}, read by noodles");
    }
}

fn stream_bam(path: &Path) -> (seqair::bam::BamHeader, Vec<RawBamRecord>) {
    let mut reader = BamStreamReader::new(std::fs::File::open(path).unwrap()).unwrap();
    let records = reader.records().collect::<Result<_, _>>().unwrap();
    (reader.into_header(), records)
}

// r[verify cram_writer.test_roundtrip]
// r[verify cram_writer.fidelity]
// r[verify cram_writer.test_size]
// r[verify cram_writer.file_definition+2]
// r[verify cram_writer.header_container]
// r[verify cram_writer.eof]
// r[verify cram_writer.compression+3]
// r[verify cram_writer.builder]
// r[verify cram_writer.finish]
// r[verify cram_writer.scope+2]
/// The fixture BAM written with `write_raw` reads back, through samtools and
/// through seqair, as samtools' own CRAM of it reads back, in CRAM 3.0 and
/// 3.1, and is no larger than samtools' of the same version.
// r[verify cram_writer.version+2]
// r[verify cram_writer.profile]
// r[verify cram_writer.rans_nx16]
// r[verify cram_writer.tok3]
#[test]
fn fixture_reads_back_as_samtools_cram_does() {
    let bam = data("test.bam");
    let fasta = data("test.fasta.gz");
    let (header, records) = stream_bam(&bam);
    for (version, name) in [
        (CramVersion::V3_0, "3.0"),
        (CramVersion::V3_1(CramProfile::Normal), "3.1"),
        (CramVersion::V3_1(CramProfile::Archive), "3.1"),
    ] {
        let archive = version == CramVersion::V3_1(CramProfile::Archive);
        let dir = tempfile::tempdir().unwrap();
        let ours = dir.path().join("ours.cram");
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let mut writer =
            CramWriterBuilder::to_path(&ours, &header, reference).version(version).build().unwrap();
        for record in &records {
            writer.write_raw(record).unwrap();
        }
        writer.finish().unwrap();
        // The file definition says which minor version.
        assert_eq!(std::fs::read(&ours).unwrap()[5], u8::from(version != CramVersion::V3_0));

        let theirs = dir.path().join("theirs.cram");
        let version_option = format!("version={name}");
        let profile_option = if archive { "archive" } else { "normal" };
        samtools(&[
            "view",
            "-C",
            "--output-fmt-option",
            &version_option,
            "--output-fmt-option",
            profile_option,
            "-T",
            fasta.to_str().unwrap(),
            "-o",
            theirs.to_str().unwrap(),
            bam.to_str().unwrap(),
        ]);
        let expected = samtools_lines(&theirs, &fasta);
        let through_samtools = samtools_lines(&ours, &fasta);
        assert_eq!(through_samtools.len(), expected.len());
        for (i, (ours, theirs)) in through_samtools.iter().zip(&expected).enumerate() {
            assert_eq!(ours, theirs, "CRAM {name} record {i}");
        }
        assert!(seqair_lines(&ours, &fasta) == expected, "seqair reads it as samtools does");
        assert_noodles_reads_it_as_seqair(&ours, &fasta);

        samtools(&["quickcheck", ours.to_str().unwrap()]);
        let size = |p: &Path| std::fs::metadata(p).unwrap().len();
        eprintln!("CRAM {name}: ours {} bytes, samtools {} bytes", size(&ours), size(&theirs));
        assert!(size(&ours) <= size(&theirs), "larger than samtools' CRAM {name}");
        if archive {
            // The archive codecs are what the file holds: fqzcomp (`f`) for
            // the qualities, the arithmetic coder (`a`/`A`) somewhere, as
            // samtools reports them.
            let sizes = samtools(&["cram-size", ours.to_str().unwrap()]);
            let method = |series: &str| {
                sizes.lines().find_map(|l| {
                    let f: Vec<&str> = l.split_whitespace().collect();
                    (f.get(6) == Some(&series)).then(|| f.get(5).copied().unwrap_or("").to_owned())
                })
            };
            assert!(method("QS").is_some_and(|m| m.contains('f')), "{sizes}");
            assert!(
                sizes
                    .lines()
                    .any(|l| l.split_whitespace().nth(5).is_some_and(|m| m.contains(['a', 'A']))),
                "{sizes}"
            );
        }

        // Several containers, so the pool reorders some.
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
            .version(version)
            .compression_threads(4)
            .build()
            .unwrap();
        for record in &records {
            writer.write_raw(record).unwrap();
        }
        let threaded = writer.finish().unwrap().0;
        assert!(threaded[26..] == std::fs::read(&ours).unwrap()[26..], "threads change the file");
    }
}

// r[verify cram_writer.embedded_reference]
// r[verify cram_writer.test_size]
/// The fixture BAM against a FASTA without its lambda contig reads back as
/// samtools' `embed_ref=2` CRAM of it does, through samtools, seqair and
/// noodles, and is no larger.
#[test]
fn a_fixture_contig_the_fasta_lacks_embeds_its_reference() {
    let bam = data("test.bam");
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("no_lambda.fa");
    let all = data("test.fasta.gz");
    let text = samtools(&["faidx", all.to_str().unwrap(), "chr19", "2kb_3_Unmodified"]);
    std::fs::write(&fasta, text).unwrap();
    samtools(&["faidx", fasta.to_str().unwrap()]);
    let (header, records) = stream_bam(&bam);
    let ours = dir.path().join("ours.cram");
    let reference = IndexedFastaReader::open(&fasta).unwrap();
    let mut writer = CramWriterBuilder::to_path(&ours, &header, reference).build().unwrap();
    for record in &records {
        writer.write_raw(record).unwrap();
    }
    writer.finish().unwrap();
    assert!(embeds_a_reference(&ours));

    let theirs = dir.path().join("theirs.cram");
    let fa = fasta.to_str().unwrap();
    let out = theirs.to_str().unwrap();
    let options = ["--output-fmt-option", "embed_ref=2"];
    samtools(
        &[&["view", "-C", "-T", fa, "-o", out][..], &options, &[bam.to_str().unwrap()]].concat(),
    );
    let expected = samtools_lines(&theirs, &fasta);
    assert!(samtools_lines(&ours, &fasta) == expected, "samtools reads it as its own");
    assert!(seqair_lines(&ours, &fasta) == expected, "seqair reads it as samtools does");
    assert_noodles_reads_it_as_seqair(&ours, &fasta);
    let size = |p: &Path| std::fs::metadata(p).unwrap().len();
    eprintln!("ours {} bytes, samtools embed_ref=2 {} bytes", size(&ours), size(&theirs));
    assert!(size(&ours) <= size(&theirs), "larger than samtools' CRAM");
}

// r[verify cram_writer.compression.learning]
// r[verify cram_writer.threads+2]
/// With many slices, most try only the two codecs the trial slices chose:
/// the file reads back as the exhaustive search's does, is within 1% of its
/// size, and threads do not change a byte.
#[test]
fn learned_codecs_keep_the_records_and_the_size() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = data("test.fasta.gz");
    let (header, records) = stream_bam(&data("test.bam"));
    // Every 150th record: 196 records in slices of 1, 197 slices, so trials
    // at 0, 70 and 140 and each block small enough for a debug build.
    let records: Vec<_> = records.into_iter().step_by(150).collect();
    for version in [CramVersion::V3_0, CramVersion::V3_1(CramProfile::Normal)] {
        let write = |exhaustive: bool, threads: usize| {
            let reference = IndexedFastaReader::open(&fasta).unwrap();
            let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
                .version(version)
                .records_per_slice(1)
                .exhaustive_codec_search(exhaustive)
                .compression_threads(threads)
                .build()
                .unwrap();
            for record in &records {
                writer.write_raw(record).unwrap();
            }
            writer.finish().unwrap().0
        };
        let learned = write(false, 0);
        let exhaustive = write(true, 0);
        assert!(write(false, 3) == learned, "threads change the learned file");
        let (ours, theirs) = (dir.path().join("learned.cram"), dir.path().join("all.cram"));
        std::fs::write(&ours, &learned).unwrap();
        std::fs::write(&theirs, &exhaustive).unwrap();
        assert_eq!(samtools_lines(&ours, &fasta), samtools_lines(&theirs, &fasta));
        assert!(seqair_lines(&ours, &fasta) == samtools_lines(&theirs, &fasta));
        eprintln!("{version:?}: learned {} bytes, exhaustive {}", learned.len(), exhaustive.len());
        assert!(learned.len() * 100 <= exhaustive.len() * 101, "{version:?} grew more than 1%");
        if version == CramVersion::V3_1(CramProfile::Normal) {
            // Some learned slice made another choice than the full search:
            // the learning ran.
            assert!(learned != exhaustive, "no slice learned");
        }
    }
}

/// A `.crai`'s lines, decompressed.
fn crai_lines(path: &Path) -> Vec<String> {
    use std::io::Read as _;
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines().map(String::from).collect()
}

// r[verify cram_writer.slicing+2]
/// The archive profile cuts slices of up to 100,000 records: the fixture's
/// 29,358 records take one slice per reference, where the normal profile
/// cuts chr19 every 10,000.
#[test]
fn archive_slices_are_larger() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = data("test.fasta.gz");
    let (header, records) = stream_bam(&data("test.bam"));
    let slices = |version| {
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
            .version(version)
            .write_index(true)
            .build()
            .unwrap();
        for record in &records {
            writer.write_raw(record).unwrap();
        }
        let (_, index) = writer.finish().unwrap();
        let path = dir.path().join("index.crai");
        index.unwrap().write_crai(std::fs::File::create(&path).unwrap()).unwrap();
        crai_lines(&path).len()
    };
    let normal = slices(CramVersion::V3_1(CramProfile::Normal));
    let archive = slices(CramVersion::V3_1(CramProfile::Archive));
    let references =
        records.iter().filter_map(|r| r.tid()).collect::<std::collections::BTreeSet<_>>().len();
    let unplaced = usize::from(records.iter().any(|r| r.tid().is_none()));
    assert_eq!(archive, references + unplaced);
    assert!(normal > archive, "{normal} vs {archive}");
}

// r[verify cram_writer.crai+2]
/// The index the writer builds is the one `samtools index` builds for its
/// output, and seqair's indexed reader queries the file with it.
#[test]
fn the_index_is_samtools_index() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = data("test.fasta.gz");
    let (header, records) = stream_bam(&data("test.bam"));
    let cram = dir.path().join("ours.cram");
    let reference = IndexedFastaReader::open(&fasta).unwrap();
    let mut writer =
        CramWriterBuilder::to_path(&cram, &header, reference).write_index(true).build().unwrap();
    for record in &records {
        writer.write_raw(record).unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let ours = dir.path().join("ours.cram.crai");
    index.unwrap().write_crai(std::fs::File::create(&ours).unwrap()).unwrap();

    let theirs = dir.path().join("theirs.crai");
    samtools(&["index", "-o", theirs.to_str().unwrap(), cram.to_str().unwrap()]);
    assert_eq!(crai_lines(&ours), crai_lines(&theirs));

    let mut reader = seqair::cram::reader::IndexedCramReader::open(&cram, &fasta).unwrap();
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid("chr19").unwrap();
    let span =
        seqair_types::Pos0::new(6_103_000).unwrap()..=seqair_types::Pos0::new(6_104_000).unwrap();
    let n = reader.fetch_into(tid, span.into(), &mut store).unwrap();
    let expected = samtools(&[
        "view",
        "-c",
        "-T",
        fasta.to_str().unwrap(),
        cram.to_str().unwrap(),
        "chr19:6103001-6104001",
    ]);
    assert_eq!(n.to_string(), expected.trim());
}

// ── Property: drawn reads against samtools ──────────────────────────────

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::{AuxData, RecordStore};
use seqair_types::IupacBase;
use std::fmt::Write as _;

/// Several, so that slices of few records switch the writer to
/// multi-reference slices (`r[cram_writer.multi_ref]`).
const CONTIGS: [(&str, usize); 4] = [("chr1", 3000), ("chr2", 1500), ("chr3", 700), ("chr4", 400)];

/// xorshift, for bases: their values shrink to nothing useful.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn pick(&mut self, from: &[u8]) -> u8 {
        from[(self.next() % from.len() as u64) as usize]
    }
}

/// A reference with `N` runs, IUPAC codes and soft-masked (lower case)
/// stretches, as real ones have.
fn reference_bases(rng: &mut Rng) -> Vec<Vec<u8>> {
    CONTIGS
        .iter()
        .map(|&(_, len)| {
            let mut lower = false;
            (0..len)
                .map(|_| {
                    if rng.next() % 200 == 0 {
                        lower = !lower;
                    }
                    let b = match rng.next() % 100 {
                        0 => b'N',
                        1 => rng.pick(b"RYKM"),
                        // A letter that is no IUPAC code: the decoders write
                        // it into MD differently. (Not `*` or `-`: htslib
                        // uppercases those into other bytes and its MD5s
                        // then disagree with the spec's.)
                        2 if rng.next() % 10 == 0 => b'X',
                        _ => rng.pick(b"ACGT"),
                    };
                    if lower { b.to_ascii_lowercase() } else { b }
                })
                .collect()
        })
        .collect()
}

fn write_fasta(dir: &Path, refs: &[Vec<u8>]) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut text = String::new();
    // An empty one is left out.
    for ((name, _), bases) in CONTIGS.iter().zip(refs).filter(|(_, b)| !b.is_empty()) {
        writeln!(text, ">{name}").unwrap();
        for line in bases.chunks(60) {
            writeln!(text, "{}", std::str::from_utf8(line).unwrap()).unwrap();
        }
    }
    std::fs::write(&path, text).unwrap();
    samtools(&["faidx", path.to_str().unwrap()]);
    path
}

const READ_GROUPS: [&str; 2] = ["rg1", "rg2"];

fn header_text() -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").unwrap();
    }
    for rg in READ_GROUPS {
        writeln!(text, "@RG\tID:{rg}\tSM:s").unwrap();
    }
    text
}

/// One SAM line and its sort key.
#[derive(Debug)]
struct Line {
    key: (usize, usize),
    text: String,
}

/// A read: mapped with substitutions, IUPAC and `=` bases, indels, clips,
/// skips and padding, sometimes past its reference's end; or unmapped,
/// placed or not. SEQ `*`, missing qualities, and tags of every type.
#[hegel::composite]
fn read(tc: &TestCase, refs: Vec<Vec<u8>>, name: String, paired_primary: bool) -> Line {
    let mut rng = Rng(tc.draw(gs::integers::<u64>()) | 1);
    let unmapped = tc.draw(gs::integers::<u8>().max_value(9)) == 0;
    // Any FLAG but the unmapped bit, which `unmapped` sets. htslib's decode
    // adds 0x1 to a mate-reverse (0x20) read, which seqair keeps as it is.
    let mut flag = tc.draw(gs::integers::<u16>().max_value(0xfff)) & !0x4;
    if paired_primary {
        flag = (flag | 0x1) & !0x900;
    }
    if flag & 0x1 == 0 {
        flag &= !0x20;
    }
    let contig = tc.draw(gs::integers::<usize>().max_value(CONTIGS.len() - 1));
    let len = CONTIGS[contig].1;
    // Near the end often, so some reads run past it.
    let pos = if tc.draw(gs::booleans()) {
        len - 1 - tc.draw(gs::integers::<usize>().max_value(40))
    } else {
        tc.draw(gs::integers::<usize>().max_value(len - 1))
    };

    let (rname, pos1, mapq, cigar, mut seq, key) = if unmapped {
        flag |= 0x4;
        let n = tc.draw(gs::integers::<usize>().max_value(40));
        let seq: Vec<u8> = (0..n).map(|_| rng.pick(b"ACGTNACGTRYKM")).collect();
        if tc.draw(gs::booleans()) {
            // Some aligners give unmapped reads a MAPQ, which CRAM drops.
            let mapq = tc.draw(gs::sampled_from(vec![0, 0, 37]));
            (CONTIGS[contig].0.to_string(), pos + 1, mapq, "*".to_string(), seq, (contig, pos))
        } else {
            let pos1 = tc.draw(gs::sampled_from(vec![0, 1, 1, 17]));
            ("*".to_string(), pos1, 0, "*".to_string(), seq, (usize::MAX, 0))
        }
    } else {
        let mapq = tc.draw(gs::integers::<u8>());
        let reference = &refs[contig];
        let mut ops: Vec<(char, usize)> = Vec::new();
        if tc.draw(gs::booleans()) {
            ops.push(('H', 1 + tc.draw(gs::integers::<usize>().max_value(4))));
        }
        if tc.draw(gs::booleans()) {
            ops.push(('S', 1 + tc.draw(gs::integers::<usize>().max_value(8))));
        }
        if tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
            ops.push((
                tc.draw(gs::sampled_from(vec!['I', 'D'])),
                1 + tc.draw(gs::integers::<usize>().max_value(3)),
            ));
        }
        for i in 0..1 + tc.draw(gs::integers::<usize>().max_value(3)) {
            if i > 0 {
                let op = tc.draw(gs::sampled_from(vec!['I', 'I', 'D', 'N', 'P']));
                let n =
                    1 + tc.draw(gs::integers::<usize>().max_value(if op == 'N' { 200 } else { 4 }));
                ops.push((op, n));
            }
            let m = tc.draw(gs::sampled_from(vec!['M', 'M', '=', 'X']));
            ops.push((m, 1 + tc.draw(gs::integers::<usize>().max_value(40))));
        }
        if tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
            ops.push((
                tc.draw(gs::sampled_from(vec!['I', 'D'])),
                1 + tc.draw(gs::integers::<usize>().max_value(3)),
            ));
        }
        if tc.draw(gs::booleans()) {
            ops.push(('S', 1 + tc.draw(gs::integers::<usize>().max_value(8))));
        }
        if tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
            ops.push(('H', 1 + tc.draw(gs::integers::<usize>().max_value(4))));
        }
        let substitutions = tc.draw(gs::integers::<u64>().max_value(30));
        let mut seq = Vec::new();
        let mut g = pos;
        for &(op, n) in &ops {
            match op {
                'M' | '=' | 'X' => {
                    for _ in 0..n {
                        let r = reference.get(g).map(u8::to_ascii_uppercase);
                        let b = match r {
                            Some(r) if rng.next() % 100 >= substitutions => r,
                            _ => rng.pick(b"ACGTNACGT=RYKMN"),
                        };
                        seq.push(b);
                        g += 1;
                    }
                }
                'I' | 'S' => seq.extend((0..n).map(|_| rng.pick(b"ACGTN"))),
                'D' | 'N' => g += n,
                _ => {}
            }
        }
        let cigar: String = ops.iter().map(|(op, n)| format!("{n}{op}")).collect();
        (CONTIGS[contig].0.to_string(), pos + 1, mapq, cigar, seq, (contig, pos))
    };
    if !unmapped && tc.draw(gs::integers::<u8>().max_value(19)) == 0 {
        seq.clear();
    }
    let qual = if seq.is_empty() || tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
        "*".to_string()
    } else {
        seq.iter().map(|_| char::from(b'!' + (rng.next() % 41) as u8)).collect()
    };
    let seq = if seq.is_empty() { "*".to_string() } else { String::from_utf8(seq).unwrap() };

    let (rnext, pnext, tlen) = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => ("*", 0, 0),
        1 => (
            "=",
            1 + tc.draw(gs::integers::<usize>().max_value(1000)),
            tc.draw(gs::integers::<i32>().min_value(-i32::MAX)),
        ),
        _ => (
            "chr2",
            1 + tc.draw(gs::integers::<usize>().max_value(1000)),
            tc.draw(gs::integers::<i32>().min_value(-i32::MAX)),
        ),
    };
    let rnext = if rname == "*" && rnext == "=" { "*" } else { rnext };

    let mut tags = Vec::new();
    if tc.draw(gs::booleans()) {
        tags.push(format!("RG:Z:{}", tc.draw(gs::sampled_from(vec!["rg1", "rg2", "elsewhere"]))));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("NM:i:{}", tc.draw(gs::integers::<u32>().max_value(300))));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("MD:Z:{}", tc.draw(gs::integers::<u32>().max_value(300))));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("XS:i:{}", tc.draw(gs::integers::<i32>())));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("XA:A:{}", tc.draw(gs::sampled_from(vec!['a', 'Z', '!']))));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("XF:f:{}", tc.draw(gs::integers::<i16>()) as f32 / 8.0));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("ZZ:Z:text {}", tc.draw(gs::integers::<u16>())));
    }
    if tc.draw(gs::booleans()) {
        tags.push("XH:H:1AE301".to_string());
    }
    if tc.draw(gs::booleans()) {
        let n = tc.draw(gs::integers::<usize>().max_value(4));
        let values: String = (0..n).map(|i| format!(",{}", i as i32 * 1000 - 1500)).collect();
        tags.push(format!("XB:B:s{values}"));
    }
    if tc.draw(gs::booleans()) {
        tags.push(format!("XC:B:c,{},-3", tc.draw(gs::integers::<i8>())));
    }
    if tc.draw(gs::booleans()) {
        tags.push("XG:B:f,1.5,-0.25".to_string());
    }
    // Same tags in another order make another tag line.
    for i in (1..tags.len()).rev() {
        tags.swap(i, (rng.next() % (i as u64 + 1)) as usize);
    }
    let mut text = format!(
        "{name}\t{flag}\t{rname}\t{pos1}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
    );
    for tag in tags {
        text.push('\t');
        text.push_str(&tag);
    }
    text.push('\n');
    Line { key, text }
}

/// A SAM line with `MD` and `NM` moved to the end, sorted: the decoders put
/// a tag they generate after the stored ones, and htslib keeps `MD`/`NM`
/// stored where seqair drops them. Every other tag keeps its place.
fn normalized(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    let (fixed, tags) = fields.split_at(11.min(fields.len()));
    let moved = |t: &&str| ["MD:", "NM:"].iter().any(|p| t.starts_with(p));
    let (mut tail, kept): (Vec<&str>, Vec<&str>) = tags.iter().copied().partition(moved);
    tail.sort_unstable();
    [fixed, &kept, &tail].concat().join("\t")
}

/// How the writer is configured in a test.
#[derive(Debug, Clone, Copy)]
struct Options {
    records_per_slice: u32,
    version: CramVersion,
    /// The FASTA the CRAM is written against lacks `chr4`, the last
    /// reference, which a multi-reference slice of the others ends before.
    lacks_chr4: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self { records_per_slice: 10_000, version: CramVersion::default(), lacks_chr4: false }
    }
}

/// The records of `bam` written through each of the writer's inputs, and
/// with and without threads, which must all give the same bytes. Returns
/// the file and its `.crai`.
fn write_each_way(dir: &Path, bam: &Path, fasta: &Path, options: Options) -> (PathBuf, PathBuf) {
    type Writer = seqair::cram::writer::CramWriter<Vec<u8>>;
    let (header, records) = stream_bam(bam);
    // The comparisons between inputs store blocks raw: the records' bytes
    // are what they compare, and the codec search dominates a debug build.
    let write_at = |level: CompressionLevel, threads: usize, f: &dyn Fn(&mut Writer)| {
        let reference = IndexedFastaReader::open(fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
            .write_index(true)
            .compression_level(level)
            .compression_threads(threads)
            .records_per_slice(options.records_per_slice)
            .version(options.version)
            .build()
            .unwrap();
        f(&mut writer);
        let (bytes, index) = writer.finish().unwrap();
        let mut crai = Vec::new();
        index.unwrap().write_crai(&mut crai).unwrap();
        (bytes, crai)
    };
    let records = &records;
    let raw_of = |placed_only: bool| {
        move |w: &mut Writer| {
            for r in records.iter().filter(|r| !placed_only || r.tid().is_some()) {
                w.write_raw(r).unwrap();
            }
        }
    };
    let write =
        |threads: usize, f: &dyn Fn(&mut Writer)| write_at(CompressionLevel::Level0, threads, f);
    let raw = write(0, &raw_of(false));
    let owned = write(0, &|w| {
        for r in records {
            let mut seq = Vec::new();
            r.iupac_bases_into(&mut seq);
            let owned = OwnedBamRecord {
                tid: r.tid(),
                pos: r.pos(),
                mapq: r.mapq(),
                flags: r.flags(),
                mate_fields: r.mate_fields(),
                qname: r.qname().to_vec(),
                cigar: r.cigar().iter().collect(),
                seq,
                qual: r.qual().map(<[_]>::to_vec).unwrap_or_default(),
                aux: AuxData::from_bytes(r.aux().as_bytes().to_vec()),
            };
            w.write(&owned).unwrap();
        }
    });
    assert!(raw == owned, "write and write_raw differ");
    // r[verify cram_writer.threads+2]
    assert!(write(3, &raw_of(false)) == raw, "the threaded writer's file or index differs");
    // The store keeps placed records only.
    let mut store = RecordStore::new();
    let placed: Vec<_> = records
        .iter()
        .filter(|r| r.tid().is_some())
        .map(|r| store.push_raw(r.as_payload(), &header, &mut ()).unwrap().unwrap())
        .collect();
    let stored = write(0, &|w| {
        for &idx in &placed {
            w.write_store_record(&store, idx).unwrap();
        }
    });
    assert!(write(0, &raw_of(true)) == stored, "write_store_record and write_raw differ");
    // r[verify cram_writer.store]
    let one_by_one = |w: &mut Writer| {
        for &idx in &placed {
            w.write_store_record(&store, idx).unwrap();
        }
    };
    let whole = |w: &mut Writer| w.write_store(&store).unwrap();
    assert!(write(3, &whole) == stored, "write_store and write_store_record differ");
    let level = CompressionLevel::default();
    assert!(
        write_at(level, 2, &whole) == write_at(level, 0, &one_by_one),
        "write_store and write_store_record differ, compressed"
    );
    let (cram, index) = write_at(CompressionLevel::default(), 0, &raw_of(false));
    let path = dir.join("ours.cram");
    std::fs::write(&path, cram).unwrap();
    let crai = dir.join("ours.cram.crai");
    std::fs::write(&crai, index).unwrap();
    (path, crai)
}

// r[verify cram_writer.test_roundtrip]
// r[verify cram_writer.inputs]
// r[verify cram_writer.features]
// r[verify cram_writer.substitution+2]
// r[verify cram_writer.past_reference_end]
// r[verify cram_writer.md_nm+2]
// r[verify cram_writer.embedded_reference]
// r[verify cram_writer.record.read_group]
// r[verify cram_writer.record.unmapped]
// r[verify cram_writer.record.mate]
// r[verify cram_writer.mates+2]
// r[verify cram_writer.record.order]
// r[verify cram_writer.record.flags]
// r[verify cram_writer.record.position]
// r[verify cram_writer.record.quality]
// r[verify cram_writer.record.mapq]
// r[verify cram_writer.slicing+2]
// r[verify cram_writer.container_header+2]
// r[verify cram_writer.slice_header+2]
// r[verify cram_writer.slice_extent+2]
// r[verify cram_writer.preservation_map+2]
// r[verify cram_writer.data_series+2]
/// Drawn reads written by seqair read back as samtools' CRAM 3.0 of them
/// does, through samtools and through seqair; and every input writes the
/// same bytes.
#[hegel::test(test_cases = pinned::cases(100))]
fn drawn_reads_read_back_as_samtools_cram_does(tc: TestCase) {
    let refs = reference_bases(&mut Rng(tc.draw(gs::integers::<u64>()) | 1));
    let n = tc.draw(gs::integers::<usize>().max_value(60));
    // Pairs of reads share a name, and samtools fixmate makes their mate
    // fields agree, so the writer attaches them.
    let pairs = tc.draw(gs::booleans());
    let name = |i: usize| if pairs { format!("p{}", i / 2) } else { format!("r{i}") };
    // fixmate only fixes pairs of paired primary reads.
    let primary_pairs = pairs && tc.draw(gs::booleans());
    let mut lines: Vec<Line> =
        (0..n).map(|i| tc.draw(read(refs.clone(), name(i), primary_pairs))).collect();
    lines.sort_by_key(|l| l.key);
    let mut text: String = lines.iter().map(|l| l.text.as_str()).collect();
    if pairs && let Some(fixed) = fixmate(&text) {
        text = fixed;
    }
    // Small slices sometimes, so records cross slice boundaries.
    let options = Options {
        records_per_slice: tc.draw(gs::sampled_from(vec![10_000, 10_000, 1, 2, 5])),
        version: tc.draw(gs::sampled_from(vec![
            CramVersion::V3_0,
            CramVersion::V3_1(CramProfile::Normal),
            CramVersion::V3_1(CramProfile::Archive),
        ])),
        lacks_chr4: tc.draw(gs::booleans()),
    };
    let calmd = tc.draw(gs::booleans());
    let compared = assert_reads_back_as_samtools_cram_does(&refs, &text, calmd, options);
    tc.assume(compared);
}

/// `records` (SAM, no header) with their mate fields made consistent by
/// `samtools fixmate`, coordinate-sorted again; `None` where samtools
/// fails (fixmate breaks some reads that run past the reference's end).
fn fixmate(records: &str) -> Option<String> {
    let dir = tempfile::tempdir().unwrap();
    let path = |name: &str| dir.path().join(name).to_str().unwrap().to_owned();
    std::fs::write(path("in.sam"), header_text() + records).unwrap();
    let ok = |args: &[&str]| Command::new("samtools").args(args).output().unwrap().status.success();
    let (sam, names, fixed, sorted) =
        (path("in.sam"), path("names.bam"), path("fixed.bam"), path("sorted.bam"));
    let steps: [&[&str]; 3] = [
        &["sort", "-n", "-o", &names, &sam],
        &["fixmate", &names, &fixed],
        &["sort", "-o", &sorted, &fixed],
    ];
    if !steps.iter().all(|args| ok(args)) {
        return None;
    }
    Some(
        samtools(&["view", "--no-PG", &path("sorted.bam")])
            .lines()
            .map(|l| format!("{l}\n"))
            .collect(),
    )
}

/// `sam` (records, no header) against `refs`, through seqair and through
/// samtools: both CRAMs read back alike. With `calmd`, the records get
/// correct `MD`/`NM` first.
///
/// Returns `false`, having checked only that samtools reads seqair's file,
/// when samtools cannot read its own: a deletion past the reference's end
/// without `MD`/`NM` (`r[cram_writer.md_nm.past_reference_end]`).
fn assert_reads_back_as_samtools_cram_does(
    refs: &[Vec<u8>],
    records: &str,
    calmd: bool,
    options: Options,
) -> bool {
    let dir = tempfile::tempdir().unwrap();
    let full = tempfile::tempdir().unwrap();
    let full_fasta = write_fasta(full.path(), refs);
    // Without chr4, its slices embed a reference built from their reads.
    let fasta = if options.lacks_chr4 {
        let mut lacking = refs.to_vec();
        lacking[3].clear();
        write_fasta(dir.path(), &lacking)
    } else {
        full_fasta.clone()
    };
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, header_text() + records).unwrap();
    let mut bam = dir.path().join("in.bam");
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    if calmd {
        // Correct MD/NM, which the writer drops and the decode regenerates.
        let out = Command::new("samtools")
            .args(["calmd", "-b", bam.to_str().unwrap(), full_fasta.to_str().unwrap()])
            .output()
            .unwrap();
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        bam = dir.path().join("calmd.bam");
        std::fs::write(&bam, out.stdout).unwrap();
    }

    let (ours, crai) = write_each_way(dir.path(), &bam, &fasta, options);
    // r[verify cram_writer.crai+2]
    let samtools_crai = dir.path().join("samtools.crai");
    samtools(&["index", "-o", samtools_crai.to_str().unwrap(), ours.to_str().unwrap()]);
    assert_eq!(crai_lines(&crai), crai_lines(&samtools_crai), "the index is samtools index's");
    // samtools' own CRAM, the oracle of the normalizations, against every
    // reference: given a FASTA that lacks one, htslib may write the whole
    // file without a reference.
    let theirs = dir.path().join("theirs.cram");
    samtools(&[
        "view",
        "-C",
        "--output-fmt-option",
        "version=3.0",
        "-T",
        full_fasta.to_str().unwrap(),
        "-o",
        theirs.to_str().unwrap(),
        bam.to_str().unwrap(),
    ]);
    let through_samtools = samtools_lines(&ours, &fasta);
    assert_noodles_reads_it_as_seqair(&ours, &fasta);
    // Where the reference has a letter that is no IUPAC code, seqair's
    // decode writes it into a generated MD as N, htslib's as it is
    // (`r[cram.stream.iupac+1]`).
    if refs.iter().flatten().all(|b| IupacBase::ASCII.contains(&b.to_ascii_uppercase())) {
        assert_eq!(
            seqair_lines(&ours, &fasta),
            through_samtools,
            "seqair and samtools read it alike"
        );
    }
    let theirs_read = Command::new("samtools")
        .args(["view", "--no-PG", "-T", full_fasta.to_str().unwrap(), theirs.to_str().unwrap()])
        .output()
        .unwrap();
    if !theirs_read.status.success() {
        let stderr = String::from_utf8_lossy(&theirs_read.stderr);
        assert!(stderr.contains("beyond slice reference extents"), "{stderr}");
        return false;
    }
    // A chr4 record without MD/NM reads back without them when chr4's slice
    // embeds its reference (`r[cram_writer.embedded_reference.md_nm]`),
    // where htslib generates them from the FASTA's chr4.
    let input = samtools_lines(&bam, &full_fasta);
    let expected: Vec<String> = String::from_utf8(theirs_read.stdout)
        .unwrap()
        .lines()
        .zip(&input)
        .map(|(line, input)| {
            let on_chr4 = line.split('\t').nth(2) == Some("chr4");
            if !(options.lacks_chr4 && on_chr4) {
                return normalized(line);
            }
            let generated = |tag: &str| {
                line.contains(&format!("\t{tag}:")) && !input.contains(&format!("\t{tag}:"))
            };
            let (md, nm) = (generated("MD"), generated("NM"));
            let kept: Vec<&str> = line
                .split('\t')
                .filter(|t| !(md && t.starts_with("MD:") || nm && t.starts_with("NM:")))
                .collect();
            normalized(&kept.join("\t"))
        })
        .collect();
    let normalized_ours: Vec<String> = through_samtools.iter().map(|l| normalized(l)).collect();
    assert_eq!(normalized_ours, expected);
    true
}

/// Reads the property test found.
#[test]
fn found_cases_read_back_as_samtools_cram_does() {
    let refs = reference_bases(&mut Rng(1));
    for records in [
        "r0\t4\tchr1\t1\t0\t*\t*\t0\t0\t*\t*\n",
        // Past the end without a deletion: MD/NM generated, as htslib does.
        "r1\t0\tchr1\t3000\t0\t2M\t*\t0\t0\tAA\t*\n",
    ] {
        assert!(assert_reads_back_as_samtools_cram_does(&refs, records, false, Options::default()));
    }
}

// r[verify cram_writer.md_nm.past_reference_end]
/// A deletion past the reference's end: samtools cannot read the CRAM it
/// writes of it, and reads seqair's, with the record as it was.
#[test]
fn a_deletion_past_the_reference_end_reads_back_through_samtools() {
    let refs = reference_bases(&mut Rng(1));
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path(), &refs);
    let last = refs[0][2999].to_ascii_uppercase() as char;
    let records = [
        format!("bare\t0\tchr1\t3000\t9\t1M3D\t*\t0\t0\t{last}\tI\n"),
        format!("tagged\t0\tchr1\t3000\t9\t1M3D\t*\t0\t0\t{last}\tI\tMD:Z:1^XYZ\tNM:i:3\n"),
        format!("md_only\t0\tchr1\t3000\t9\t1M3D\t*\t0\t0\t{last}\tI\tMD:Z:1^XYZ\n"),
    ];
    let text: String = records.concat();
    assert!(!assert_reads_back_as_samtools_cram_does(&refs, &text, false, Options::default()));
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, header_text() + &text).unwrap();
    let bam = dir.path().join("in.bam");
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    let (ours, _) = write_each_way(dir.path(), &bam, &fasta, Options::default());
    let read: Vec<String> = samtools_lines(&ours, &fasta);
    let expected: Vec<String> = text.lines().map(String::from).collect();
    assert_eq!(read, expected);
}

/// The header of a BAM whose text is `text` and whose binary reference list
/// is [`CONTIGS`].
fn bam_header_with_text(text: &str) -> seqair::bam::BamHeader {
    let mut bytes = b"BAM\x01".to_vec();
    bytes.extend_from_slice(&(text.len() as i32).to_le_bytes());
    bytes.extend_from_slice(text.as_bytes());
    bytes.extend_from_slice(&(CONTIGS.len() as i32).to_le_bytes());
    for (name, len) in CONTIGS {
        bytes.extend_from_slice(&(name.len() as i32 + 1).to_le_bytes());
        bytes.extend_from_slice(name.as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(&(len as i32).to_le_bytes());
    }
    let mut bgzf = seqair::io::BgzfWriter::new(Vec::new());
    bgzf.write_all(&bytes).unwrap();
    let bam = bgzf.finish().unwrap();
    BamStreamReader::new(std::io::Cursor::new(bam)).unwrap().into_header()
}

// r[verify cram_writer.header_text]
// r[verify cram_writer.header_m5]
/// A BAM header whose text lists no references gets one `@SQ` line per
/// target, with `M5`, after a leading `@HD` or else first.
#[test]
fn a_header_without_sq_lines_gets_them() {
    let refs = reference_bases(&mut Rng(3));
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path(), &refs);
    for (text, before) in [
        ("", ""),
        ("@HD\tVN:1.6\n@CO\tc\n", "@HD\tVN:1.6\n"),
        ("@CO\tc\n@HD\tVN:1.6\n", ""),
        ("@RG\tID:a\n", ""),
    ] {
        let header = bam_header_with_text(text);
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference).build().unwrap();
        let (cram, _) = writer.finish().unwrap();
        let written = CramStreamReader::new(&cram[..]).unwrap().into_header();
        let written = written.header_text();
        let sq: String = CONTIGS
            .iter()
            .zip(&refs)
            .map(|((name, len), bases)| {
                let m5 = format!("{:x}", md5::compute(bases.to_ascii_uppercase()));
                format!("@SQ\tSN:{name}\tLN:{len}\tM5:{m5}\n")
            })
            .collect();
        let rest = text.strip_prefix(before).unwrap();
        assert_eq!(written, format!("{before}{sq}{rest}"), "{text:?}");
    }
}

/// Whether `samtools cram-size` finds a slice of `cram` that embeds its
/// reference.
fn embeds_a_reference(cram: &Path) -> bool {
    samtools(&["cram-size", "-v", cram.to_str().unwrap()]).contains("embedded_ref")
}

// r[verify cram_writer.embedded_reference]
// r[verify cram_writer.embedded_reference.md_nm]
/// Records on a reference the FASTA lacks are written against one built
/// from them, and read back through samtools, seqair and noodles as they
/// were: an `MD` that fits its CIGAR gives the bases (and is dropped), one
/// that does not stays, and a record without `MD`/`NM` reads back without.
#[test]
fn a_reference_the_fasta_lacks_is_embedded() {
    let refs = reference_bases(&mut Rng(5));
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path(), &refs);
    let text = header_text().replacen("@RG", "@SQ\tSN:chrX\tLN:500\n@RG", 1);
    let records = [
        "a\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII\tXS:i:1",
        "b\t0\tchrX\t3\t60\t2S3M1I2M2D1M1N2M\t*\t0\t0\tTTACGNACGAA\tIIIIIIIIIII\tNM:i:3\tMD:Z:0A",
        "b2\t0\tchrX\t4\t60\t4M2D2M\t*\t0\t0\tCGTAAC\tIIIIII\tMD:Z:1T2^GG2\tNM:i:3",
        "c\t81\tchrX\t9\t7\t1H3M1P1M\t=\t3\t-7\tACGT\t*",
        "d\t133\tchrX\t9\t0\t*\t=\t9\t0\tAC\t##",
    ];
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, text + &records.join("\n") + "\n").unwrap();
    let bam = dir.path().join("in.bam");
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    let (ours, _) = write_each_way(dir.path(), &bam, &fasta, Options::default());
    let through_samtools = samtools_lines(&ours, &fasta);
    // chr1 has a reference: its record gains MD/NM on read.
    assert!(through_samtools[0].starts_with(records[0]));
    assert_eq!(through_samtools[1..], records[1..]);
    assert_eq!(seqair_lines(&ours, &fasta), through_samtools);
    assert_noodles_reads_it_as_seqair(&ours, &fasta);
    assert!(embeds_a_reference(&ours));
}

// r[verify cram_writer.reference_free+2]
/// Records on a reference the FASTA lacks with no aligned base to build one
/// from, or too few (under 1% of the extent), are written without a
/// reference.
#[test]
fn a_reference_the_fasta_lacks_without_aligned_bases_is_written_without_it() {
    let refs = reference_bases(&mut Rng(5));
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path(), &refs);
    let text = header_text().replacen("@RG", "@SQ\tSN:chrX\tLN:500\n@RG", 1);
    let records = [
        "c\t0\tchrX\t9\t7\t2S1I\t*\t0\t0\tACG\t*",
        "d\t0\tchrX\t9\t7\t3M\t*\t0\t0\t*\t*",
        "e\t4\tchrX\t9\t0\t*\t*\t0\t0\tAC\t##",
    ];
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, text.clone() + &records.join("\n") + "\n").unwrap();
    let bam = dir.path().join("in.bam");
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    let (ours, _) = write_each_way(dir.path(), &bam, &fasta, Options::default());
    let through_samtools = samtools_lines(&ours, &fasta);
    assert_eq!(through_samtools, records);
    assert_eq!(seqair_lines(&ours, &fasta), through_samtools);
    assert_noodles_reads_it_as_seqair(&ours, &fasta);
    assert!(!embeds_a_reference(&ours));

    // Aligned bases under 1% of the extent: 2 × 2 over 450 positions.
    let sparse = ["f\t0\tchrX\t9\t7\t2M\t*\t0\t0\tAC\t*", "g\t0\tchrX\t457\t7\t2M\t*\t0\t0\tGT\t*"];
    std::fs::write(&sam, text + &sparse.join("\n") + "\n").unwrap();
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    let (ours, _) = write_each_way(dir.path(), &bam, &fasta, Options::default());
    let through_samtools = samtools_lines(&ours, &fasta);
    assert_eq!(through_samtools, sparse);
    assert!(!embeds_a_reference(&ours));
}

// r[verify cram_writer.multi_ref]
// r[verify cram_writer.multi_ref.records]
// r[verify cram_writer.crai+2]
// r[verify cram_writer.slicing+2]
/// One read each on chr1 and chr2 makes two small slices; the next slices
/// hold several references, chr3's reads and chr4's, until a full slice of
/// chr4's alone switches multi-reference slices off again. Read back through
/// samtools, seqair and noodles as samtools' own CRAM, indexed as
/// `samtools index` indexes it.
#[test]
fn few_reads_on_many_references_share_slices() {
    let refs = reference_bases(&mut Rng(9));
    let read = |name: String, contig: usize, pos: usize| {
        let (rname, len) = CONTIGS[contig];
        let pos = pos % (len - 10);
        let bases = std::str::from_utf8(&refs[contig][pos..pos + 10]).unwrap().to_ascii_uppercase();
        format!("{name}\t0\t{rname}\t{}\t60\t10M\t*\t0\t0\t{bases}\t*\n", pos + 1)
    };
    let mut text = read("a".into(), 0, 5) + &read("b".into(), 1, 5);
    text += &read("c1".into(), 2, 5);
    text += &read("c2".into(), 2, 9);
    for i in 0..80 {
        text += &read(format!("d{i}"), 3, 3 * i);
    }
    text += "u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\nu2\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t*\n";
    // Slices of at most 40 records are small below 4.
    let options = Options { records_per_slice: 40, ..Options::default() };
    assert!(assert_reads_back_as_samtools_cram_does(&refs, &text, true, options));

    let dir = tempfile::tempdir().unwrap();
    let fasta = write_fasta(dir.path(), &refs);
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, header_text() + &text).unwrap();
    let bam = dir.path().join("in.bam");
    samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
    let (_, crai) = write_each_way(dir.path(), &bam, &fasta, options);
    // Reference id, and the container's offset.
    let entries: Vec<(String, String)> = crai_lines(&crai)
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_owned(), f[3].to_owned())
        })
        .collect();
    let refs_by_container: Vec<Vec<&str>> = entries
        .chunk_by(|a, b| a.1 == b.1)
        .map(|c| c.iter().map(|(r, _)| r.as_str()).collect())
        .collect();
    assert_eq!(
        refs_by_container,
        [
            vec!["0"],
            vec!["1"],
            // Multi-reference from here, until 40 records of chr4 alone.
            vec!["2", "3"],
            vec!["3"],
            vec!["3"],
            // Unplaced records are never in one.
            vec!["-1"],
        ]
    );
}

/// The block sizes `samtools cram-size -v` reports for `series` in `cram`.
fn series_bytes(cram: &Path, series: &str) -> usize {
    samtools(&["cram-size", "-v", cram.to_str().unwrap()])
        .lines()
        .filter(|l| l.starts_with("BLOCK") && l.split_whitespace().last() == Some(series))
        .filter_map(|l| l.split_whitespace().nth(2)?.parse::<usize>().ok())
        .sum()
}

// r[verify cram_writer.mates+2]
/// Pairs the decode rebuilds are attached and read back as they were
/// through samtools, seqair and noodles; others detached. Attachment shows
/// as a `CF` block: the first of a pair has `CF` 0x5, the second 0x1, where
/// every detached record here has 0x3, a constant, which takes no block.
#[test]
fn mate_pairs_attach_where_every_decoder_rebuilds_them() {
    let refs = reference_bases(&mut Rng(7));
    let base =
        |at: usize, n: usize| String::from_utf8(refs[0][at..at + n].to_ascii_uppercase()).unwrap();
    let cases: [(&str, Vec<String>, bool); 7] = [
        (
            "proper pair",
            vec![
                format!("a\t99\tchr1\t101\t60\t10M\t=\t201\t110\t{}\t*", base(100, 10)),
                format!("a\t147\tchr1\t201\t60\t10M\t=\t101\t-110\t{}\t*", base(200, 10)),
            ],
            true,
        ),
        (
            "placed unmapped mate",
            vec![
                format!("b\t73\tchr1\t101\t60\t10M\t=\t101\t0\t{}\t*", base(100, 10)),
                "b\t133\tchr1\t101\t0\t*\t=\t101\t0\tACGT\t*".to_string(),
            ],
            false,
        ),
        (
            "same start and end",
            vec![
                format!("c\t99\tchr1\t101\t60\t10M\t=\t101\t10\t{}\t*", base(100, 10)),
                format!("c\t147\tchr1\t101\t60\t10M\t=\t101\t-10\t{}\t*", base(100, 10)),
            ],
            true,
        ),
        (
            "SEQ *",
            vec![
                "d\t99\tchr1\t101\t60\t10M\t=\t201\t110\t*\t*".to_string(),
                "d\t147\tchr1\t201\t60\t10M\t=\t101\t-110\t*\t*".to_string(),
            ],
            true,
        ),
        (
            "unplaced pair",
            vec![
                "e\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\t*".to_string(),
                "e\t141\t*\t0\t0\t*\t*\t0\t0\tACGT\t*".to_string(),
            ],
            false,
        ),
        (
            "a TLEN the decode computes otherwise",
            vec![
                format!("f\t99\tchr1\t101\t60\t10M\t=\t201\t500\t{}\t*", base(100, 10)),
                format!("f\t147\tchr1\t201\t60\t10M\t=\t101\t-500\t{}\t*", base(200, 10)),
            ],
            false,
        ),
        (
            "three records of a name",
            vec![
                format!("g\t99\tchr1\t101\t60\t10M\t=\t201\t110\t{}\t*", base(100, 10)),
                format!("g\t2147\tchr1\t150\t60\t10M\t=\t101\t0\t{}\t*", base(149, 10)),
                format!("g\t147\tchr1\t201\t60\t10M\t=\t101\t-110\t{}\t*", base(200, 10)),
            ],
            false,
        ),
    ];
    for (what, records, attaches) in cases {
        let text: String = records.iter().map(|r| format!("{r}\n")).collect();
        let dir = tempfile::tempdir().unwrap();
        let fasta = write_fasta(dir.path(), &refs);
        let sam = dir.path().join("in.sam");
        std::fs::write(&sam, header_text() + &text).unwrap();
        let bam = dir.path().join("in.bam");
        samtools(&["view", "-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]);
        let (ours, _) = write_each_way(dir.path(), &bam, &fasta, Options::default());
        assert_eq!(series_bytes(&ours, "CF") > 0, attaches, "{what}");
        let read: Vec<String> = samtools_lines(&ours, &fasta)
            .iter()
            .map(|l| l.split('\t').take(11).collect::<Vec<_>>().join("\t"))
            .collect();
        assert_eq!(read, records, "{what}");
        assert_eq!(seqair_lines(&ours, &fasta).len(), records.len());
        assert_noodles_reads_it_as_seqair(&ours, &fasta);
    }
}

// r[verify cram_writer.store]
/// `write_store` writes the bytes and index `write_store_record` writes one
/// record at a time: whatever the threads, slice size, version and profile,
/// whether the FASTA lacks a contig (embedded references, held records),
/// and when the records come in two stores, the first leaving a slice open
/// that the second completes.
// Each case writes the file twice: slow by hegel's measure, on purpose.
#[hegel::test(
    test_cases = pinned::cases(8),
    suppress_health_check = [hegel::HealthCheck::TooSlow]
)]
fn write_store_writes_what_write_store_record_writes(tc: TestCase) {
    use seqair::bam::RecordStore;
    let dir = tempfile::tempdir().unwrap();
    let fasta = if tc.draw(gs::booleans()) {
        data("test.fasta.gz")
    } else {
        let fasta = dir.path().join("no_lambda.fa");
        let all = data("test.fasta.gz");
        let text = samtools(&["faidx", all.to_str().unwrap(), "chr19", "2kb_3_Unmodified"]);
        std::fs::write(&fasta, text).unwrap();
        samtools(&["faidx", fasta.to_str().unwrap()]);
        fasta
    };
    let (header, records) = stream_bam(&data("test.bam"));
    // Every 13th record, so a debug build writes quickly; small slices make
    // trial slices, many slices and multi-reference ones.
    let step = 13;
    let records: Vec<_> = records.into_iter().step_by(step).collect();
    let mut stores = [RecordStore::new(), RecordStore::new()];
    let split = tc.draw(gs::integers::<usize>().max_value(records.len()));
    for (i, r) in records.iter().enumerate() {
        let store = &mut stores[usize::from(i >= split)];
        store.push_raw(r.as_payload(), &header, &mut ()).unwrap();
    }
    let options = (
        tc.draw(gs::sampled_from(vec![1_u32, 7, 15, 10_000])),
        tc.draw(gs::sampled_from(vec![1_usize, 3])),
        tc.draw(gs::sampled_from(vec![
            CramVersion::V3_0,
            CramVersion::V3_1(CramProfile::Normal),
            CramVersion::V3_1(CramProfile::Archive),
        ])),
    );
    let (records_per_slice, threads, version) = options;
    let write = |by_store: bool, threads: usize| {
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), &header, reference)
            .write_index(true)
            .records_per_slice(records_per_slice)
            .compression_threads(threads)
            .version(version)
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
        let (bytes, index) = writer.finish().unwrap();
        let mut crai = Vec::new();
        index.unwrap().write_crai(&mut crai).unwrap();
        (bytes, crai)
    };
    let expected = write(false, 0);
    assert!(write(true, threads) == expected, "{options:?}, split {split}");
}

// r[verify cram_writer.store]
/// A record `write_store` refuses ends the call: the records before it are
/// written as `write_store_record` writes them, it and the later ones not.
#[test]
fn write_store_stops_at_a_refused_record() {
    use seqair::bam::RecordStore;
    use seqair::cram::writer::CramWriteError;
    let (header, records) = stream_bam(&data("test.bam"));
    let fasta = data("test.fasta.gz");
    let placed: Vec<_> = records.iter().filter(|r| r.tid().is_some()).collect();
    let store_of = |records: &mut dyn Iterator<Item = &&RawBamRecord>| {
        let mut store = RecordStore::new();
        for r in records {
            store.push_raw(r.as_payload(), &header, &mut ()).unwrap();
        }
        store
    };
    let write = |header: &seqair::bam::BamHeader,
                 threads: usize,
                 by_store: bool,
                 stores: &[RecordStore],
                 errors: &mut Vec<CramWriteError>| {
        let reference = IndexedFastaReader::open(&fasta).unwrap();
        let mut writer = CramWriterBuilder::to_writer(Vec::new(), header, reference)
            .records_per_slice(7)
            .compression_threads(threads)
            .write_index(true)
            .build()
            .unwrap();
        for store in stores {
            let result = if by_store {
                writer.write_store(store)
            } else {
                // One by one up to the first refusal, as `write_store` stops.
                store.indices().try_for_each(|idx| writer.write_store_record(store, idx))
            };
            errors.extend(result.err());
        }
        let (bytes, index) = writer.finish().unwrap();
        let mut crai = Vec::new();
        index.unwrap().write_crai(&mut crai).unwrap();
        (bytes, crai)
    };
    // Records 1..301 in order, then record 0: unsorted; then a store after.
    let unsorted = [
        store_of(&mut placed.iter().skip(1).take(300).chain([&placed[0]])),
        store_of(&mut placed.iter().skip(400).take(50)),
    ];
    // A wrong `M5` for the second reference: its first record is refused,
    // after the slice before it was handed out; and so are the later ones.
    let text = header.header_text().replacen(
        "SN:2kb_3_Unmodified\tLN:2018",
        "SN:2kb_3_Unmodified\tLN:2018\tM5:00000000000000000000000000000000",
        1,
    );
    let wrong_m5 = seqair::bam::BamHeader::from_sam_text(&text).unwrap();
    let all = [store_of(&mut placed.iter())];
    for threads in [1, 3] {
        for (header, stores, refused) in
            [(&header, &unsorted[..], "Unsorted"), (&wrong_m5, &all[..], "ReferenceMd5Mismatch")]
        {
            let (mut ours, mut theirs) = (Vec::new(), Vec::new());
            let by_store = write(header, threads, true, stores, &mut ours);
            assert_eq!(ours.len(), 1, "{ours:?}");
            assert!(format!("{:?}", ours[0]).starts_with(refused), "{ours:?}");
            let one_by_one = write(header, 0, false, stores, &mut theirs);
            assert_eq!(format!("{ours:?}"), format!("{theirs:?}"));
            assert!(by_store == one_by_one, "{refused}, {threads} threads");
        }
    }
}
