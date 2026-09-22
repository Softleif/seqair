//! The CRAM shapes samtools can be *told* to write, on generated data.
//!
//! `cram_bam_equivalence_properties.rs` generates reads and compares CRAM
//! against BAM, but only ever for one contig and one set of writer options —
//! the default single-ref, external-reference, one-slice-per-container layout.
//! The layouts that make CRAM awkward to decode are the other ones: a slice
//! whose records come from several references (`ref_seq_id == -2`), a slice
//! carrying its own reference bases, a file cut into many small slices so the
//! CRAI has to actually route a query. Those were fixture-only.
//!
//! Every property here writes the *same* generated reads as BAM and as CRAM
//! and reads both with seqair: BAM is the oracle for CRAM, and the generated
//! reads are the oracle for BAM, so a misreading shared by neither can hide.
//! What varies is the writer option matrix — version, `embed_ref`,
//! `seqs_per_slice`, `slices_per_container`, `multi_seq_per_slice` — and, for
//! the query properties, the region asked for.
//!
//! The awkward layouts are not hoped for: each property that exists to cover
//! one asserts, per case, that the file it just wrote really has it, by
//! parsing the container/slice headers back out of the bytes.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(clippy::print_stdout, reason = "a skipped test says why on stdout")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test code with known small values"
)]

use hegel::prelude::*;
use seqair::bam::{Pos0, RecordStore, RejectUnmapped};
use seqair::cram::block::{self, ContentType};
use seqair::cram::container::ContainerHeader;
use seqair::cram::index::{CraiEntry, CramIndex};
use seqair::cram::slice::SliceHeader;
use seqair::reader::Readers;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Every contig is this long. Short contigs are the point: samtools decides to
/// merge references into one slice based on how few records it has for each.
const CONTIG_LEN: u32 = 300;
/// The Sanger quality range a generated QUAL string is drawn from.
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
/// The CRAM file definition: `CRAM` + major + minor + 20 bytes of file ID.
const FILE_DEFINITION_LEN: usize = 26;

// ── The generated sample ─────────────────────────────────────────────

/// One generated read, in the form the SAM line and the assertions both need.
#[derive(Debug, Clone)]
struct GenRead {
    contig: usize,
    /// 1-based, as SAM writes it.
    pos: u32,
    cigar: Vec<(u32, char)>,
    seq: String,
    qual: String,
    flags: u32,
    mapq: u32,
}

impl GenRead {
    fn cigar_text(&self) -> String {
        let mut out = String::new();
        for (len, op) in &self.cigar {
            write!(out, "{len}{op}").expect("writing to a String cannot fail");
        }
        out
    }

    /// Reference bases consumed — `M` and `D`, not the clips or insertions.
    fn ref_span(&self) -> u32 {
        self.cigar.iter().filter(|(_, op)| matches!(op, 'M' | 'D')).map(|(l, _)| l).sum()
    }

    fn start0(&self) -> u32 {
        self.pos - 1
    }

    /// Last reference position covered, 0-based inclusive.
    fn last0(&self) -> u32 {
        self.start0() + self.ref_span() - 1
    }

    fn overlaps(&self, start0: u32, end0: u32) -> bool {
        self.start0() <= end0 && self.last0() >= start0
    }
}

/// A reference, the reads aligned to it, and some reads aligned to nothing.
#[derive(Debug, Clone)]
struct Sample {
    /// Bases per contig; contig `i` is named `c{i}`.
    contigs: Vec<String>,
    reads: Vec<GenRead>,
    /// `(seq, qual)` for reads with flag 0x4 and no reference.
    unmapped: Vec<(String, String)>,
}

impl Sample {
    fn contig_name(i: usize) -> String {
        format!("c{i}")
    }

    fn reads_on(&self, contig: usize) -> impl Iterator<Item = &GenRead> {
        self.reads.iter().filter(move |r| r.contig == contig)
    }
}

fn bases(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet("ACGT").min_size(len).max_size(len)
}

/// A CIGAR valid by construction: matches at both ends of the aligned body,
/// at most one indel between them, soft clips outside. Kept short so several
/// reads fit on a 300-base contig without crowding out the end positions.
#[hegel::composite]
fn arb_cigar(tc: &TestCase) -> Vec<(u32, char)> {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let mut body = vec![(tc.draw_silent(n(8, 30)), 'M')];
    if tc.draw_silent(gs::booleans()) {
        body.push((tc.draw_silent(n(1, 5)), tc.draw_silent(gs::sampled_from(&['I', 'D']))));
        body.push((tc.draw_silent(n(8, 30)), 'M'));
    }

    let mut ops = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 6)), 'S'));
    }
    ops.extend(body);
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 6)), 'S'));
    }
    ops
}

/// A read whose aligned bases start out equal to the reference, with a handful
/// of positions overwritten: CRAM stores an aligned base as a match or as a
/// substitution depending on that comparison, so a read that differs
/// everywhere would only ever exercise the substitution path.
fn arb_read(tc: &TestCase, contig: usize, reference: &str) -> GenRead {
    let ops = tc.draw_silent(arb_cigar());
    let span: u32 = ops.iter().filter(|(_, op)| matches!(op, 'M' | 'D')).map(|(l, _)| l).sum();
    let last_start = CONTIG_LEN - span + 1;
    let pos = tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(last_start));

    let ref_bytes = reference.as_bytes();
    let mut seq: Vec<u8> = Vec::new();
    let mut ref_cursor = (pos - 1) as usize;
    for &(len, op) in &ops {
        match op {
            'M' => {
                seq.extend_from_slice(&ref_bytes[ref_cursor..ref_cursor + len as usize]);
                ref_cursor += len as usize;
            }
            'D' => ref_cursor += len as usize,
            'I' | 'S' => seq.extend_from_slice(tc.draw_silent(bases(len as usize)).as_bytes()),
            _ => unreachable!("arb_cigar emits only M/I/D/S"),
        }
    }
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(3)) {
        let at = tc.draw_silent(gs::integers::<usize>().max_value(seq.len() - 1));
        seq[at] = tc.draw_silent(gs::sampled_from(b"ACGT"));
    }

    let qual = tc.draw_silent(qual_string(seq.len()));
    GenRead {
        contig,
        pos,
        cigar: ops,
        seq: String::from_utf8(seq).expect("ACGT is ASCII"),
        qual,
        flags: tc.draw_silent(gs::sampled_from(&[0u32, 16])),
        mapq: tc.draw_silent(gs::integers::<u32>().max_value(60)),
    }
}

/// More contigs, fewer reads on each: the shape that gives samtools a reason
/// to merge references into one slice.
#[hegel::composite]
fn arb_many_contig_sample(tc: &TestCase) -> Sample {
    arb_sample_sized(tc, 3, 6, 1, 2, 0, 1)
}

/// Guarantees reads that align nowhere, so the file gets an unmapped slice.
#[hegel::composite]
fn arb_sample_with_unmapped(tc: &TestCase) -> Sample {
    arb_sample_sized(tc, 2, 3, 1, 3, 1, 3)
}

fn qual_string(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet(QUAL_CHARS).min_size(len).max_size(len)
}

/// Several short contigs, each with at least one read. "At least one" is what
/// makes the multi-reference properties reachable by construction: a slice can
/// only span two references if two references have records.
#[hegel::composite]
fn arb_sample(tc: &TestCase) -> Sample {
    arb_sample_sized(tc, 2, 4, 1, 3, 0, 2)
}

/// `arb_sample` with the shape knobs exposed, for the properties that need a
/// particular one (many contigs, or unmapped reads, or neither).
fn arb_sample_sized(
    tc: &TestCase,
    min_contigs: usize,
    max_contigs: usize,
    min_reads: usize,
    max_reads: usize,
    min_unmapped: usize,
    max_unmapped: usize,
) -> Sample {
    let n_contigs =
        tc.draw_silent(gs::integers::<usize>().min_value(min_contigs).max_value(max_contigs));
    let contigs: Vec<String> =
        (0..n_contigs).map(|_| tc.draw_silent(bases(CONTIG_LEN as usize))).collect();

    let mut reads = Vec::new();
    for (i, reference) in contigs.iter().enumerate() {
        let n = tc.draw_silent(gs::integers::<usize>().min_value(min_reads).max_value(max_reads));
        for _ in 0..n {
            reads.push(arb_read(tc, i, reference));
        }
    }

    let n_unmapped =
        tc.draw_silent(gs::integers::<usize>().min_value(min_unmapped).max_value(max_unmapped));
    let unmapped = (0..n_unmapped)
        .map(|_| {
            let len = tc.draw_silent(gs::integers::<usize>().min_value(20).max_value(40));
            (tc.draw_silent(bases(len)), tc.draw_silent(qual_string(len)))
        })
        .collect();

    Sample { contigs, reads, unmapped }
}

fn sam_text(sample: &Sample) -> String {
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for i in 0..sample.contigs.len() {
        writeln!(sam, "@SQ\tSN:{}\tLN:{CONTIG_LEN}", Sample::contig_name(i))
            .expect("writing to a String cannot fail");
    }
    for (i, r) in sample.reads.iter().enumerate() {
        writeln!(
            sam,
            "r{i}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            r.flags,
            Sample::contig_name(r.contig),
            r.pos,
            r.mapq,
            r.cigar_text(),
            r.seq,
            r.qual
        )
        .expect("writing to a String cannot fail");
    }
    for (i, (seq, qual)) in sample.unmapped.iter().enumerate() {
        writeln!(sam, "u{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
            .expect("writing to a String cannot fail");
    }
    sam
}

// ── The samtools writer option matrix ────────────────────────────────

/// The CRAM writer options this file varies. Each maps to one
/// `--output-fmt-option` argument to `samtools view -C`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CramOpts {
    /// `3.0` or `3.1`. See `cram_v2_is_out_of_scope` for why not `2.1`.
    version: &'static str,
    /// 0 = external reference only, 1 = embed the reference bases,
    /// 2 = embed a consensus computed from the reads.
    embed_ref: u8,
    seqs_per_slice: u32,
    /// More than one slice per container gives the CRAI several entries at the
    /// same container offset, distinguished only by their slice landmark.
    slices_per_container: u32,
    /// Whether one slice may hold records from more than one reference.
    /// `None` leaves the decision to samtools, which turns it on by itself
    /// when a container would otherwise hold very few records per reference.
    multi_seq: Option<bool>,
}

impl CramOpts {
    fn args(self) -> Vec<String> {
        let mut args = vec![format!("version={}", self.version)];
        if self.embed_ref > 0 {
            args.push(format!("embed_ref={}", self.embed_ref));
        }
        args.push(format!("seqs_per_slice={}", self.seqs_per_slice));
        args.push(format!("slices_per_container={}", self.slices_per_container));
        if let Some(multi_seq) = self.multi_seq {
            args.push(format!("multi_seq_per_slice={}", u8::from(multi_seq)));
        }
        args
    }

    fn label(self) -> String {
        format!(
            "v{} embed_ref={} seqs_per_slice={} slices_per_container={} multi_seq={:?}",
            self.version,
            self.embed_ref,
            self.seqs_per_slice,
            self.slices_per_container,
            self.multi_seq
        )
    }
}

#[hegel::composite]
fn arb_opts(tc: &TestCase) -> CramOpts {
    let embed_ref = tc.draw_silent(gs::sampled_from(&[0u8, 1]));
    let multi_seq = arb_multi_seq(tc, embed_ref);
    CramOpts {
        version: tc.draw_silent(gs::sampled_from(&["3.0", "3.1"])),
        embed_ref,
        seqs_per_slice: tc.draw_silent(gs::sampled_from(&[1u32, 2, 5, 10_000])),
        slices_per_container: arb_slices_per_container(tc, multi_seq),
        multi_seq,
    }
}

/// Several slices per container are only drawn when multi-reference slices are
/// ruled out. A multi-ref container with two slices for the same reference
/// decodes the second one's tail as `N` — see
/// `multi_ref_containers_use_every_slices_reference_range` — and leaving that
/// in the matrix would fail every property here for the same one reason.
fn arb_slices_per_container(tc: &TestCase, multi_seq: Option<bool>) -> u32 {
    if multi_seq == Some(false) { tc.draw_silent(gs::sampled_from(&[1u32, 2])) } else { 1 }
}

/// `embed_ref` and `multi_seq_per_slice` together make htslib fall back to
/// *no-reference* mode — it warns about changing from `embed_ref` to `no_ref`
/// — and seqair cannot currently decode a no-reference read that has an
/// insertion; see `no_ref_cram_decodes_reads_with_insertions`. Until that is
/// fixed, pin `multi_seq_per_slice=0` whenever the reference is embedded, so
/// this matrix stays about the layouts it means to cover.
fn arb_multi_seq(tc: &TestCase, embed_ref: u8) -> Option<bool> {
    if embed_ref > 0 {
        return Some(false);
    }
    tc.draw_silent(gs::sampled_from(&[None, Some(false), Some(true)]))
}

// ── Running samtools ─────────────────────────────────────────────────

fn run(cmd: &mut Command, what: &str) {
    let output = cmd.stdout(Stdio::null()).output().unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(output.status.success(), "{what} failed: {}", String::from_utf8_lossy(&output.stderr));
}

/// Write the reference and index it, so samtools can compute the M5 tags the
/// CRAM header carries and seqair can verify them on read.
fn write_reference(dir: &Path, sample: &Sample) -> PathBuf {
    let fasta = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fasta).expect("create FASTA");
    for (i, contig) in sample.contigs.iter().enumerate() {
        writeln!(f, ">{}", Sample::contig_name(i)).expect("write");
        for chunk in contig.as_bytes().chunks(60) {
            f.write_all(chunk).expect("write");
            f.write_all(b"\n").expect("write");
        }
    }
    drop(f);
    run(Command::new("samtools").arg("faidx").arg(&fasta), "samtools faidx");
    fasta
}

/// Sort the generated SAM into an indexed BAM. This is the oracle file.
fn write_bam(dir: &Path, sample: &Sample) -> PathBuf {
    let sam_path = dir.join("generated.sam");
    std::fs::File::create(&sam_path)
        .expect("create SAM")
        .write_all(sam_text(sample).as_bytes())
        .expect("write SAM");

    let bam = dir.join("generated.bam");
    run(Command::new("samtools").args(["sort", "-o"]).arg(&bam).arg(&sam_path), "samtools sort");
    run(Command::new("samtools").arg("index").arg(&bam), "samtools index (bam)");
    bam
}

/// Transcode the oracle BAM to CRAM under `opts`, and index it.
fn write_cram(dir: &Path, bam: &Path, fasta: &Path, opts: CramOpts) -> PathBuf {
    let cram = dir.join(format!("generated_{}.cram", opts.version));
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "-C"]);
    for opt in opts.args() {
        cmd.arg("--output-fmt-option").arg(opt);
    }
    cmd.arg("-T").arg(fasta).arg("-o").arg(&cram).arg(bam);
    run(&mut cmd, &format!("samtools view -C ({})", opts.label()));
    run(Command::new("samtools").arg("index").arg(&cram), "samtools index (cram)");
    cram
}

// ── Reading back with seqair ─────────────────────────────────────────

/// A decoded record, in the fields both formats must agree on.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Decoded {
    qname: Vec<u8>,
    pos: i64,
    flags: u16,
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    /// Phred values, not the ASCII SAM writes.
    qual: Vec<u8>,
}

/// Everything seqair returns for `[start0, end0]` on one contig.
///
/// `RejectUnmapped` is the customizer the other CRAM comparison tests use, and
/// it is what makes BAM a fair oracle here: a BAM index cannot reach an
/// unplaced read from a contig query at all, while a CRAM multi-ref slice
/// decodes its unplaced records alongside the mapped ones and hands them to
/// `filter_raw` to decide on (`r[cram.index.unmapped]`). Rejecting them puts
/// both formats on the same footing.
fn fetch(path: &Path, fasta: &Path, contig: usize, start0: u32, end0: u32) -> Vec<Decoded> {
    let mut readers = Readers::open_customized(path, fasta, RejectUnmapped).expect("open");
    let name = Sample::contig_name(contig);
    let tid = readers.header().tid(&name).expect("contig in header");
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, Pos0::new(start0).unwrap(), Pos0::new(end0).unwrap(), &mut store)
        .expect("fetch");
    decode_store(&store)
}

fn decode_store<E: Default>(store: &RecordStore<E>) -> Vec<Decoded> {
    store
        .indices()
        .map(|idx| {
            let rec = store.record(idx).unwrap();
            Decoded {
                qname: rec.qname().to_vec(),
                pos: rec.pos.as_i64(),
                flags: rec.flags.raw(),
                mapq: rec.mapq,
                cigar: rec.cigar().iter().map(ToString::to_string).collect::<String>(),
                seq: rec.seq().iter().map(|b| **b).collect(),
                qual: seqair_types::BaseQuality::slice_to_bytes(rec.qual()).to_vec(),
            }
        })
        .collect()
}

/// Every record on every contig, contig by contig, from one open handle.
fn fetch_all(path: &Path, fasta: &Path, sample: &Sample) -> Vec<Decoded> {
    let mut readers = Readers::open_customized(path, fasta, RejectUnmapped).expect("open");
    let mut store = RecordStore::new();
    let mut out = Vec::new();
    for contig in 0..sample.contigs.len() {
        let name = Sample::contig_name(contig);
        let tid = readers.header().tid(&name).expect("contig in header");
        store.clear();
        readers
            .fetch_into(tid, Pos0::ZERO, Pos0::new(CONTIG_LEN - 1).unwrap(), &mut store)
            .expect("fetch");
        out.extend(decode_store(&store));
    }
    out
}

/// The qnames the generator gave the reads a region should return, sorted.
fn expected_qnames(sample: &Sample, contig: usize, start0: u32, end0: u32) -> Vec<String> {
    let mut names: Vec<String> = sample
        .reads
        .iter()
        .enumerate()
        .filter(|(_, r)| r.contig == contig && r.overlaps(start0, end0))
        .map(|(i, _)| format!("r{i}"))
        .collect();
    names.sort();
    names
}

fn qnames_of(decoded: &[Decoded]) -> Vec<String> {
    let mut names: Vec<String> =
        decoded.iter().map(|d| String::from_utf8(d.qname.clone()).expect("ASCII")).collect();
    names.sort();
    names
}

/// Check every decoded record against the read the generator wrote, matched on
/// the `r<index>` qname `sam_text` assigned (`samtools sort` reorders them).
fn assert_matches_generator(decoded: &[Decoded], sample: &Sample, label: &str) {
    for got in decoded {
        let name = std::str::from_utf8(&got.qname).expect("ASCII");
        let i: usize = name.trim_start_matches('r').parse().expect("qnames are r<index>");
        let want = &sample.reads[i];
        assert_eq!(got.pos, i64::from(want.pos) - 1, "{label}/{name}: pos");
        assert_eq!(u32::from(got.flags), want.flags, "{label}/{name}: flags");
        assert_eq!(u32::from(got.mapq), want.mapq, "{label}/{name}: mapq");
        assert_eq!(got.cigar, want.cigar_text(), "{label}/{name}: cigar");
        assert_eq!(got.seq, want.seq.as_bytes(), "{label}/{name}: seq");
        // SAM writes QUAL as `phred + 33`; seqair stores the Phred value.
        let want_qual: Vec<u8> = want.qual.bytes().map(|b| b - 33).collect();
        assert_eq!(got.qual, want_qual, "{label}/{name}: qual");
    }
}

// ── Looking at what samtools actually wrote ──────────────────────────

/// Every slice header in the file, read straight out of the bytes.
///
/// Walks containers from the file definition, and inside each container jumps
/// to every landmark — which is exactly the offset `r[cram.index.crai_per_slice]`
/// says a CRAI `slice_offset` repeats — and parses the slice header block there.
fn slice_headers(cram: &Path) -> Vec<SliceHeader> {
    let data = std::fs::read(cram).expect("read CRAM");
    let mut headers = Vec::new();
    let mut offset = FILE_DEFINITION_LEN;
    // The first container holds the SAM header, not slices.
    let mut is_header_container = true;

    while offset < data.len() {
        let container = match ContainerHeader::parse(&data[offset..]) {
            Ok(c) => c,
            Err(e) => panic!("container header at {offset}: {e:?}"),
        };
        if container.is_eof() {
            break;
        }
        let data_start = offset + container.header_size;
        if is_header_container {
            is_header_container = false;
            offset = data_start + container.length as usize;
            continue;
        }
        for &landmark in &container.landmarks {
            let at = data_start + landmark as usize;
            let (blk, _) = block::parse_block(&data[at..]).expect("slice header block");
            assert_eq!(blk.content_type, ContentType::SliceHeader, "landmark is a slice header");
            headers.push(SliceHeader::parse(&blk.data).expect("parse slice header"));
        }
        offset = data_start + container.length as usize;
    }

    headers
}

fn crai_entries(cram: &Path) -> Vec<CraiEntry> {
    let crai = cram.with_extension("cram.crai");
    CramIndex::from_path(&crai).expect("read CRAI").entries().to_vec()
}

// ── Properties ───────────────────────────────────────────────────────

// r[verify unified.fetch_equivalence]
// r[verify cram.record.sequence]
// r[verify cram.record.cigar_reconstruction]
// r[verify cram.record.quality]
// r[verify cram.scope.versions]
// r[verify cram.slice.embedded_ref]
/// Across the writer option matrix, CRAM decodes to exactly what BAM decodes,
/// and to exactly what the generator wrote.
///
/// The options change how the records are laid out on disk — how many slices,
/// whether the reference travels with them, which container each lands in —
/// and none of that is allowed to change a single decoded field.
#[hegel::test(test_cases = 48)]
fn cram_matches_bam_across_the_writer_option_matrix(tc: TestCase) {
    let sample = tc.draw(arb_sample().print_as_debug());
    let opts = tc.draw(arb_opts().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let from_bam = fetch_all(&bam, &fasta, &sample);
    let from_cram = fetch_all(&cram, &fasta, &sample);

    assert_eq!(from_cram.len(), sample.reads.len(), "{}: record count", opts.label());
    assert_eq!(from_cram, from_bam, "{}: CRAM and BAM disagree", opts.label());
    assert_matches_generator(&from_cram, &sample, &opts.label());

    tc.event(format!("version {}", opts.version));
    tc.event(format!("embed_ref={}", opts.embed_ref));
    tc.event_value("slices", slice_headers(&cram).len() as f64);
    tc.event_value("reads", sample.reads.len() as f64);
}

// r[verify cram.slice.multi_ref]
// r[verify cram.slice.multi_ref_skip]
// r[verify cram.index.multi_ref_slices]
// r[verify unified.fetch_equivalence]
/// A slice holding records from several references decodes per reference.
///
/// With `multi_seq_per_slice=1` and every contig carrying at least one read,
/// samtools puts records from different references in one slice and marks it
/// `ref_seq_id == -2`; the `RI` series then decides which reference each
/// record belongs to. The assertion that such a slice exists is part of the
/// property: without it this is just the matrix test again.
#[hegel::test(test_cases = 24)]
fn multi_reference_slices_decode_per_reference(tc: TestCase) {
    let sample = tc.draw(arb_many_contig_sample().print_as_debug());
    let opts = CramOpts {
        version: tc.draw(gs::sampled_from(&["3.0", "3.1"])),
        embed_ref: 0,
        seqs_per_slice: 10_000,
        slices_per_container: 1,
        multi_seq: Some(true),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let headers = slice_headers(&cram);
    let multi_ref = headers.iter().filter(|h| h.ref_seq_id == -2).count();
    assert!(
        multi_ref > 0,
        "{}: no multi-ref slice was written; slice ref ids were {:?}",
        opts.label(),
        headers.iter().map(|h| h.ref_seq_id).collect::<Vec<_>>()
    );

    // Per contig, not just in bulk: a multi-ref slice is decoded once per
    // query and must hand back only the records for the queried reference.
    for contig in 0..sample.contigs.len() {
        let from_cram = fetch(&cram, &fasta, contig, 0, CONTIG_LEN - 1);
        let from_bam = fetch(&bam, &fasta, contig, 0, CONTIG_LEN - 1);
        let want = sample.reads_on(contig).count();
        assert_eq!(from_cram.len(), want, "contig {contig}: record count");
        assert_eq!(from_cram, from_bam, "contig {contig}: CRAM and BAM disagree");
        assert_matches_generator(&from_cram, &sample, &format!("c{contig}"));
    }

    tc.event_value("multi-ref slices", multi_ref as f64);
    tc.event_value("contigs", sample.contigs.len() as f64);
}

// r[verify cram.slice.embedded_ref]
// r[verify unified.fetch_equivalence]
/// A slice that carries its own reference bases decodes to the same records.
///
/// `embed_ref=1` stores the reference segment the slice spans; `embed_ref=2`
/// stores a consensus computed from the reads, which is *not* the FASTA — so
/// the decoder must reconstruct against the embedded block and not against the
/// reference it was handed.
#[hegel::test(test_cases = 24)]
fn embedded_reference_slices_decode_to_the_same_records(tc: TestCase) {
    let sample = tc.draw(arb_sample().print_as_debug());
    let opts = CramOpts {
        version: tc.draw(gs::sampled_from(&["3.0", "3.1"])),
        embed_ref: 1,
        seqs_per_slice: tc.draw(gs::sampled_from(&[2u32, 10_000])),
        slices_per_container: tc.draw(gs::sampled_from(&[1u32, 2])),
        multi_seq: Some(false),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let headers = slice_headers(&cram);
    let embedded = headers.iter().filter(|h| h.embedded_reference >= 0).count();
    assert!(
        embedded > 0,
        "{}: no slice embedded a reference; embedded_reference ids were {:?}",
        opts.label(),
        headers.iter().map(|h| h.embedded_reference).collect::<Vec<_>>()
    );

    let from_bam = fetch_all(&bam, &fasta, &sample);
    let from_cram = fetch_all(&cram, &fasta, &sample);
    assert_eq!(from_cram, from_bam, "{}: CRAM and BAM disagree", opts.label());
    assert_matches_generator(&from_cram, &sample, &opts.label());

    tc.event(format!("embed_ref={}", opts.embed_ref));
    tc.event_value("slices with an embedded reference", embedded as f64);
}

/// A slice whose embedded reference is a *consensus* rather than the FASTA.
///
/// `embed_ref=2` makes samtools compute a consensus from the reads and store
/// that as the slice's reference; htslib then MD5s the consensus into the
/// slice header. seqair compares that MD5 against the FASTA regardless of
/// whether the slice brought its own reference, so every such file whose
/// consensus differs from the reference — which is any file with low coverage
/// or mismatching reads — fails to open with `ReferenceMd5Mismatch`.
///
/// The reads below are one per contig with substitutions, so the consensus is
/// guaranteed to differ. Un-ignore this once the MD5 check learns to skip
/// slices with `embedded_reference >= 0`.
// r[verify cram.slice.embedded_ref]
// r[verify cram.edge.reference_mismatch]
#[test]
#[ignore = "seqair MD5-checks the FASTA even when the slice embeds its own reference; \
            embed_ref=2 files whose consensus differs from the reference cannot be opened"]
fn embed_ref_consensus_slices_decode_to_the_same_records() {
    let reference = "ACGT".repeat(CONTIG_LEN as usize / 4);
    let mut seq: Vec<u8> = reference.as_bytes()[10..60].to_vec();
    for at in [3usize, 17, 29, 41] {
        seq[at] = if seq[at] == b'A' { b'C' } else { b'A' };
    }
    let sample = Sample {
        contigs: vec![reference.clone(), reference],
        reads: (0..2)
            .map(|contig| GenRead {
                contig,
                pos: 11,
                cigar: vec![(50, 'M')],
                seq: String::from_utf8(seq.clone()).expect("ACGT is ASCII"),
                qual: "I".repeat(50),
                flags: 0,
                mapq: 60,
            })
            .collect(),
        unmapped: Vec::new(),
    };
    let opts = CramOpts {
        version: "3.0",
        embed_ref: 2,
        seqs_per_slice: 10_000,
        slices_per_container: 1,
        multi_seq: Some(false),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    assert!(slice_headers(&cram).iter().any(|h| h.embedded_reference >= 0));
    assert_eq!(fetch_all(&cram, &fasta, &sample), fetch_all(&bam, &fasta, &sample));
}

/// A CRAM written with no reference at all, holding a read with an insertion.
///
/// `no_ref=1` makes samtools store every base literally instead of as edits
/// against a reference — and htslib falls into the same mode on its own when
/// `embed_ref` meets `multi_seq_per_slice`, which is how the option matrix
/// found this. seqair then reconstructs one base too many for any read with an
/// `I` op and refuses the record with `QualLenMismatch`; `D`, `S` and `H` are
/// fine, and so is the same read written against a reference. samtools reads
/// the file back without complaint.
// r[verify cram.scope.reference_required]
// r[verify cram.record.sequence]
#[test]
#[ignore = "seqair adds the insertion length twice when decoding a no-reference CRAM, \
            so any read with an I op fails with QualLenMismatch"]
fn no_ref_cram_decodes_reads_with_insertions() {
    let sample = Sample {
        contigs: vec!["ACGT".repeat(CONTIG_LEN as usize / 4)],
        reads: vec![GenRead {
            contig: 0,
            pos: 1,
            cigar: vec![(8, 'M'), (1, 'I'), (8, 'M')],
            seq: "ACGTACGTTACGTACGT".to_owned(),
            qual: "I".repeat(17),
            flags: 0,
            mapq: 60,
        }],
        unmapped: Vec::new(),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);

    let cram = dir.path().join("no_ref.cram");
    let mut cmd = Command::new("samtools");
    cmd.args([
        "view",
        "-C",
        "--output-fmt-option",
        "version=3.0",
        "--output-fmt-option",
        "no_ref=1",
        "-T",
    ])
    .arg(&fasta)
    .arg("-o")
    .arg(&cram)
    .arg(&bam);
    run(&mut cmd, "samtools view -C (no_ref)");
    run(Command::new("samtools").arg("index").arg(&cram), "samtools index (cram)");

    assert_eq!(fetch_all(&cram, &fasta, &sample), fetch_all(&bam, &fasta, &sample));
}

/// Two slices for the same reference inside one multi-ref container.
///
/// The reference window a container is decoded against is taken from *the
/// first* CRAI entry whose `container_offset` matches, but a multi-ref
/// container holds one entry per slice per reference — and htslib turns
/// multi-reference slices on by itself once a container would hold very few
/// records per reference, so this needs no option to ask for. When a later
/// slice
/// reaches further along the reference than the first one does, the window
/// stops short and the bases past its end are reconstructed as `N` — silently,
/// since falling off the end of the reference is only a warning
/// (`r[cram.slice.ref_bounds_warning]`). The window needs to be the union of
/// every matching entry, not the first one.
///
/// Both reads below are 8M on `c1`, one base apart, so the second slice needs
/// exactly one base more than the first entry's span reports.
// r[verify cram.index.multi_ref_slices]
// r[verify cram.record.sequence]
#[test]
#[ignore = "the reference window for a multi-ref container comes from the first matching \
            CRAI entry only, so later slices decode their tail as N"]
fn multi_ref_containers_use_every_slices_reference_range() {
    let reference = "ACGT".repeat(CONTIG_LEN as usize / 4);
    let read = |contig: usize, pos: u32| GenRead {
        contig,
        pos,
        cigar: vec![(8, 'M')],
        seq: reference[(pos as usize - 1)..(pos as usize + 7)].to_owned(),
        qual: "I".repeat(8),
        flags: 0,
        mapq: 60,
    };
    let sample = Sample {
        contigs: vec![reference.clone(), reference.clone()],
        reads: vec![read(0, 1), read(1, 1), read(1, 2)],
        unmapped: Vec::new(),
    };
    let opts = CramOpts {
        version: "3.0",
        embed_ref: 0,
        seqs_per_slice: 1,
        slices_per_container: 2,
        multi_seq: None,
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let headers = slice_headers(&cram);
    assert!(headers.iter().filter(|h| h.ref_seq_id == -2).count() >= 2, "{headers:?}");
    let entries = crai_entries(&cram);
    assert!(
        entries.iter().filter(|e| e.ref_id == 1).count() >= 2,
        "c1 needs two index entries in one container: {entries:?}"
    );
    assert_eq!(fetch_all(&cram, &fasta, &sample), fetch_all(&bam, &fasta, &sample));
}

// r[verify cram.index.query+2]
// r[verify cram.index.parse]
// r[verify cram.index.crai_per_slice]
// r[verify unified.fetch_equivalence]
/// A region query returns exactly the reads that overlap it.
///
/// The CRAI decides which containers and slices are even opened, so a query is
/// where an index bug shows up as missing records rather than as a parse
/// error. `seqs_per_slice` is kept small so the index has several entries to
/// route between, and the region is drawn to land anywhere on the contig,
/// including past the last read.
#[hegel::test(test_cases = 48)]
fn region_queries_return_exactly_the_overlapping_reads(tc: TestCase) {
    let sample = tc.draw(arb_sample().print_as_debug());
    let opts = tc.draw(arb_opts().print_as_debug());
    let contig =
        tc.draw(gs::integers::<usize>().max_value(sample.contigs.len() - 1).print_as_debug());
    let start0 = tc.draw(gs::integers::<u32>().max_value(CONTIG_LEN - 1).print_as_debug());
    let end0 =
        tc.draw(gs::integers::<u32>().min_value(start0).max_value(CONTIG_LEN - 1).print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let from_cram = fetch(&cram, &fasta, contig, start0, end0);
    let from_bam = fetch(&bam, &fasta, contig, start0, end0);
    let want = expected_qnames(&sample, contig, start0, end0);

    assert_eq!(qnames_of(&from_cram), want, "{}: c{contig}:{start0}-{end0}", opts.label());
    assert_eq!(from_cram, from_bam, "{}: CRAM and BAM disagree", opts.label());
    assert_matches_generator(&from_cram, &sample, &opts.label());

    // The index must not be more selective than the records are: every entry
    // the query keeps has to be one a decoded record could live in.
    let kept = CramIndex::from_path(&cram.with_extension("cram.crai"))
        .expect("read CRAI")
        .query(contig as i32, u64::from(start0), u64::from(end0))
        .len();
    assert!(kept > 0 || want.is_empty(), "index pruned every slice but reads overlap");

    tc.event_value("reads in region", want.len() as f64);
    if want.is_empty() {
        tc.event("empty region");
    }
    if want.len() == sample.reads_on(contig).count() {
        tc.event("whole contig");
    }
}

// r[verify cram.index.zero_span+2]
// r[verify cram.index.query+2]
/// A CRAI that reports no extent still returns every record.
///
/// Writers that do not compute a span emit `alignment_span == 0`, which means
/// "unknown extent", not "covers nothing". Treating it as zero-width would
/// prune the slice and lose its records silently, so this rewrites a real
/// index into that shape — the entries stay, only their spans go to zero — and
/// asks for the same regions again. Over-inclusive pruning is harmless;
/// `decode_slice` re-filters every record.
#[hegel::test(test_cases = 24)]
fn zero_span_index_entries_still_return_every_record(tc: TestCase) {
    let sample = tc.draw(arb_sample().print_as_debug());
    let opts = CramOpts {
        version: tc.draw(gs::sampled_from(&["3.0", "3.1"])),
        embed_ref: 0,
        seqs_per_slice: tc.draw(gs::sampled_from(&[1u32, 2, 10_000])),
        slices_per_container: tc.draw(gs::sampled_from(&[1u32, 2])),
        multi_seq: Some(false),
    };
    let contig =
        tc.draw(gs::integers::<usize>().max_value(sample.contigs.len() - 1).print_as_debug());
    let start0 = tc.draw(gs::integers::<u32>().max_value(CONTIG_LEN - 1).print_as_debug());
    let end0 =
        tc.draw(gs::integers::<u32>().min_value(start0).max_value(CONTIG_LEN - 1).print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let before = fetch(&cram, &fasta, contig, start0, end0);

    let rewritten = zero_the_spans(&cram);
    assert!(rewritten > 0, "no mapped CRAI entry to rewrite");
    let entries = crai_entries(&cram);
    assert!(
        entries.iter().any(|e| e.ref_id >= 0 && e.alignment_span == 0 && e.alignment_start > 0),
        "rewrite did not produce a span=0 entry with a start: {entries:?}"
    );

    let after = fetch(&cram, &fasta, contig, start0, end0);
    assert_eq!(after, before, "{}: span=0 index dropped records", opts.label());
    assert_eq!(qnames_of(&after), expected_qnames(&sample, contig, start0, end0));

    tc.event_value("entries rewritten", rewritten as f64);
}

/// Rewrite the CRAI so every mapped entry reports `alignment_span == 0`,
/// keeping the byte offsets intact. Returns how many entries changed.
fn zero_the_spans(cram: &Path) -> usize {
    use std::io::Read as _;

    let crai = cram.with_extension("cram.crai");
    let compressed = std::fs::read(&crai).expect("read CRAI");
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(&compressed[..])
        .read_to_string(&mut text)
        .expect("decompress CRAI");

    let mut rewritten = 0;
    let mut out = String::new();
    for line in text.lines() {
        let mut fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 6, "CRAI line has 6 fields: {line}");
        let ref_id: i32 = fields[0].parse().expect("CRAI ref id");
        if ref_id >= 0 && fields[2] != "0" {
            fields[2] = "0";
            rewritten += 1;
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }

    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(out.as_bytes()).expect("compress CRAI");
    std::fs::write(&crai, encoder.finish().expect("finish CRAI")).expect("write CRAI");
    rewritten
}

// r[verify cram.index.unmapped]
// r[verify cram.edge.unmapped_reads]
// r[verify cram.index.parse]
/// Unmapped reads get their own index entry, and it stays out of the way.
///
/// samtools writes reads with no reference into slices with `ref_seq_id == -1`
/// and indexes them as `-1 0 0` — start and span both zero. From CRAM v3.1 the
/// span really is zero (v3.0 rounds it up to one), which is the one place a
/// zero span means "nothing here" rather than "extent unknown": these entries
/// must be skipped, and the mapped reads must come back unaffected.
#[hegel::test(test_cases = 16)]
fn unmapped_reads_are_indexed_apart_from_the_mapped_ones(tc: TestCase) {
    let sample = tc.draw(arb_sample_with_unmapped().print_as_debug());
    let opts = CramOpts {
        version: "3.1",
        embed_ref: 0,
        seqs_per_slice: tc.draw(gs::sampled_from(&[2u32, 4])),
        slices_per_container: 1,
        multi_seq: Some(false),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);
    let cram = write_cram(dir.path(), &bam, &fasta, opts);

    let entries = crai_entries(&cram);
    let unmapped_entries: Vec<&CraiEntry> = entries.iter().filter(|e| e.ref_id < 0).collect();
    assert!(!unmapped_entries.is_empty(), "no unmapped CRAI entry: {entries:?}");
    assert!(
        unmapped_entries.iter().any(|e| e.alignment_span == 0 && e.alignment_start == 0),
        "v3.1 writes the unmapped slice as start=0 span=0: {unmapped_entries:?}"
    );

    let index = CramIndex::from_path(&cram.with_extension("cram.crai")).expect("read CRAI");
    for contig in 0..sample.contigs.len() {
        assert!(
            index
                .query(contig as i32, 0, u64::from(CONTIG_LEN) - 1)
                .iter()
                .all(|e| e.ref_id == contig as i32),
            "a whole-contig query pulled in an entry for another reference"
        );
        let from_cram = fetch(&cram, &fasta, contig, 0, CONTIG_LEN - 1);
        let from_bam = fetch(&bam, &fasta, contig, 0, CONTIG_LEN - 1);
        assert_eq!(from_cram.len(), sample.reads_on(contig).count(), "unmapped reads leaked in");
        assert_eq!(from_cram, from_bam, "contig {contig}: CRAM and BAM disagree");
    }

    tc.event_value("unmapped reads", sample.unmapped.len() as f64);
}

// r[verify cram.scope.versions]
/// CRAM 2.1 is out of scope, and says so instead of misreading the file.
///
/// samtools can still write it, so the version matrix above deliberately stops
/// at 3.0: the v2 container header has no CRC and a different field list, and
/// the reader rejects major version 2 up front rather than parsing it as v3.
/// This pins that — a v2 file must fail to open, not decode into nonsense.
#[test]
fn cram_v2_is_out_of_scope() {
    let sample = Sample {
        contigs: vec!["ACGT".repeat(CONTIG_LEN as usize / 4)],
        reads: vec![GenRead {
            contig: 0,
            pos: 1,
            cigar: vec![(20, 'M')],
            seq: "ACGT".repeat(5),
            qual: "I".repeat(20),
            flags: 0,
            mapq: 60,
        }],
        unmapped: Vec::new(),
    };

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &sample);
    let bam = write_bam(dir.path(), &sample);

    let cram = dir.path().join("v21.cram");
    let written = Command::new("samtools")
        .args(["view", "-C", "--output-fmt-option", "version=2.1", "-T"])
        .arg(&fasta)
        .arg("-o")
        .arg(&cram)
        .arg(&bam)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    if !written.success() {
        println!("skipping: this samtools cannot write CRAM 2.1");
        return;
    }

    let magic = std::fs::read(&cram).expect("read CRAM");
    assert_eq!(&magic[..4], b"CRAM");
    assert_eq!(magic[4], 2, "samtools wrote a major version other than 2");

    let err = Readers::open(&cram, &fasta).expect_err("v2.1 must not open");
    let text = format!("{err:?}");
    assert!(text.contains("UnsupportedVersion"), "expected UnsupportedVersion, got {text}");
}
