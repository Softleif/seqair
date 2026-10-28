//! The indexed readers (`IndexedReader::fetch_into`) against the streaming
//! readers (`BamStreamReader`, `SamStreamReader`, `AlignmentStream`) on one
//! coordinate-sorted BAM and bgzipped SAM of the same records, with samtools
//! as a third opinion on which records a region holds.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/sam_lines.rs"]
mod sam_lines;

use core::range::RangeInclusive;
use hegel::TestCase;
use hegel::generators as gs;
use sam_lines::{Dialect, SEQ_CHARS, aux_field, chars, cigar};
use seqair::bam::stream::BamStreamReader;
use seqair::bam::{BgzfStream, MateFields, RawBamRecord, RecordStore};
use seqair::reader::{AlignmentStream, IndexedReader, StreamFormat};
use seqair::sam::SamStreamReader;
use seqair_types::{Base, BaseQuality, Pos0};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::io::{Read, Write};
use std::num::NonZeroU8;
use std::path::{Path, PathBuf};
use std::process::Command;

// ── generated records ───────────────────────────────────────────────────

/// SAM QNAME characters (`[!-?A-~]`).
const QNAME_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

/// Contig lengths: one base, a few, either side of a 16 kbp BAI window, and
/// long enough for several.
#[hegel::composite]
fn contigs(tc: &TestCase) -> Vec<u32> {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    (0..n)
        .map(|_| tc.draw(gs::sampled_from(vec![1u32, 7, 300, 16_384, 16_385, 70_000, 250_000])))
        .collect()
}

/// Where a template record sits.
#[derive(Debug, Clone, Copy)]
enum Kind {
    Mapped,
    /// A reference and position, FLAG 0x4, CIGAR `*` or a stale one.
    PlacedUnmapped {
        cigar: bool,
    },
    /// RNAME `*`, POS 0.
    Unplaced,
}

/// One record, its position kept apart so copies can move it.
#[derive(Debug, Clone)]
struct Template {
    /// Contig index and the highest 1-based POS a copy may take.
    place: Option<(usize, u32)>,
    pos: u32,
    flag: u16,
    name: String,
    /// MAPQ onwards.
    tail: String,
}

#[hegel::composite]
fn template(tc: &TestCase, lens: Vec<u32>, dialect: Dialect) -> Template {
    let kind = tc.draw(gs::sampled_from(vec![
        Kind::Mapped,
        Kind::Mapped,
        Kind::Mapped,
        Kind::PlacedUnmapped { cigar: false },
        Kind::PlacedUnmapped { cigar: true },
        Kind::Unplaced,
    ]));
    let mut flag = tc.draw(gs::integers::<u16>().max_value(0xFFF)) & !0x4;
    if !matches!(kind, Kind::Mapped) {
        flag |= 0x4;
    }
    let contig = tc.draw(gs::integers::<usize>().max_value(lens.len() - 1));
    let len = lens[contig];
    // Every op at least 1 long; H, I, P and S alone span no reference.
    let (cigar_text, qlen, rlen) = match kind {
        Kind::Mapped | Kind::PlacedUnmapped { cigar: true } => {
            let (text, qlen, rlen) = tc.draw(cigar(Dialect { noodles: true }, 5));
            (text, qlen as usize, rlen as u32)
        }
        _ => ("*".to_owned(), tc.draw(gs::integers::<usize>().max_value(30)), 0),
    };
    let place = match kind {
        Kind::Unplaced => None,
        // The last base no further than the contig's end, where it fits.
        _ => Some((contig, len.saturating_sub(rlen.max(1)) + 1)),
    };
    let pos = match place {
        None => 0,
        Some((_, max)) => match tc.draw(gs::integers::<u8>().max_value(3)) {
            0 => 1,
            1 => max,
            _ => tc.draw(gs::integers::<u32>().min_value(1).max_value(max)),
        },
    };
    let (seq, qual) = if qlen == 0 || tc.draw(gs::weighted_booleans(0.1)) {
        ("*".to_owned(), "*".to_owned())
    } else {
        let seq = tc.draw(chars(SEQ_CHARS, qlen, qlen));
        let qual = if tc.draw(gs::weighted_booleans(0.2)) {
            "*".to_owned()
        } else {
            tc.draw(chars(QUAL_CHARS, qlen, qlen))
        };
        (seq, qual)
    };
    let (rnext, pnext) = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => ("*".to_owned(), 0),
        1 if place.is_some() => {
            ("=".to_owned(), tc.draw(gs::integers::<u32>().min_value(1).max_value(len)))
        }
        _ => {
            let mate = tc.draw(gs::integers::<usize>().max_value(lens.len() - 1));
            (format!("c{mate}"), tc.draw(gs::integers::<u32>().min_value(1).max_value(lens[mate])))
        }
    };
    let mut tail = format!(
        "{}\t{cigar_text}\t{rnext}\t{pnext}\t{}\t{seq}\t{qual}",
        tc.draw(gs::integers::<u8>()),
        tc.draw(gs::integers::<i32>()),
    );
    for i in 0..tc.draw(gs::integers::<usize>().max_value(4)) {
        let first = tc.draw(chars("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", 1, 1));
        write!(tail, "\t{}", tc.draw(aux_field(format!("{first}{i}"), dialect))).unwrap();
    }
    Template { place, pos, flag, name: tc.draw(chars(QNAME_CHARS, 0, 12)), tail }
}

/// Unsorted SAM: templates, each written `copies` times with its position
/// moved by `stride` per copy (wrapping inside its contig), so a few draws
/// give thousands of records, piled up or spread out.
#[hegel::composite]
fn sam_input(tc: &TestCase, lens: Vec<u32>) -> String {
    let dialect = Dialect { noodles: tc.draw(gs::booleans()) };
    let n = tc.draw(gs::integers::<usize>().max_value(30));
    let templates: Vec<Template> =
        (0..n).map(|_| tc.draw(template(lens.clone(), dialect))).collect();
    let copies = tc.draw(gs::sampled_from(vec![1u32, 1, 10, 100]));
    let stride = tc.draw(gs::integers::<u32>().max_value(3000));
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (i, len) in lens.iter().enumerate() {
        writeln!(text, "@SQ\tSN:c{i}\tLN:{len}").unwrap();
    }
    for (i, t) in templates.iter().enumerate() {
        for k in 0..copies {
            let (rname, pos) = match t.place {
                None => ("*".to_owned(), 0),
                Some((contig, max)) => (format!("c{contig}"), (t.pos - 1 + k * stride) % max + 1),
            };
            // Unique names: the template and the copy come first.
            writeln!(text, "{i}.{k}.{}\t{}\t{rname}\t{pos}\t{}", t.name, t.flag, t.tail).unwrap();
        }
    }
    text
}

// ── files ───────────────────────────────────────────────────────────────

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().unwrap();
    assert!(out.status.success(), "{cmd:?}: {}", out.stderr.escape_ascii());
    out.stdout
}

/// One BGZF block holding `data`; empty, it is the EOF marker's shape.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflate = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
    deflate.write_all(data).unwrap();
    let cdata = deflate.finish().unwrap();
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let mut out = vec![31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0];
    // BSIZE: the whole block's size minus one.
    out.extend_from_slice(&u16::try_from(cdata.len() + 25).unwrap().to_le_bytes());
    out.extend_from_slice(&cdata);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

/// How to cut a file's uncompressed bytes into BGZF blocks, wherever the
/// cuts fall: inside a record, a line, the header.
#[derive(Debug, Clone)]
struct Layout {
    /// Block sizes, used in turn until the data runs out.
    sizes: Vec<usize>,
    /// An empty block first and after every this many blocks.
    empty_every: Option<usize>,
}

/// `None`: the layout samtools or bgzip chose.
#[hegel::composite]
fn layout(tc: &TestCase) -> Option<Layout> {
    if tc.draw(gs::booleans()) {
        return None;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(5));
    let sizes = (0..n)
        .map(|_| {
            let max = tc.draw(gs::sampled_from(vec![40usize, 600, 60_000]));
            tc.draw(gs::integers::<usize>().min_value(1).max_value(max))
        })
        .collect();
    let empty_every = tc.draw(gs::optional(gs::integers::<usize>().min_value(1).max_value(4)));
    Some(Layout { sizes, empty_every })
}

/// `data` as BGZF in `layout`'s blocks, ending in the EOF marker.
fn reblock(data: &[u8], layout: &Layout) -> Vec<u8> {
    let mut out = Vec::new();
    if layout.empty_every.is_some() {
        out.extend(bgzf_block(&[]));
    }
    let mut rest = data;
    for (i, &size) in layout.sizes.iter().cycle().enumerate() {
        if rest.is_empty() {
            break;
        }
        let (block, tail) = rest.split_at(size.min(rest.len()));
        out.extend(bgzf_block(block));
        rest = tail;
        if layout.empty_every.is_some_and(|n| (i + 1) % n == 0) {
            out.extend(bgzf_block(&[]));
        }
    }
    out.extend(bgzf_block(&[]));
    out
}

fn inflate(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// The same records as a sorted, indexed BAM and bgzipped SAM: sorted by
/// samtools to SAM text once, then encoded to BAM by samtools and compressed
/// as it is, so both hold what that text says (a float the BAM held at full
/// precision would otherwise come back from `%g` text as another float).
struct Files {
    bam: PathBuf,
    sam: PathBuf,
}

fn files(tc: &TestCase, dir: &Path, input: &str) -> Files {
    let unsorted = dir.join("in.sam");
    std::fs::write(&unsorted, input).unwrap();
    let sorted = dir.join("sorted.sam");
    run(Command::new("samtools")
        .args(["sort", "--no-PG", "-O", "sam", "-o"])
        .arg(&sorted)
        .arg(&unsorted));

    let bam = dir.join("x.bam");
    run(Command::new("samtools").args(["view", "--no-PG", "-b", "-o"]).arg(&bam).arg(&sorted));
    if let Some(layout) = tc.draw(layout()) {
        std::fs::write(&bam, reblock(&inflate(&bam), &layout)).unwrap();
    }
    let mut index = Command::new("samtools");
    index.arg("index");
    if tc.draw(gs::booleans()) {
        index.arg("-c");
    }
    run(index.arg(&bam));

    let sam = dir.join("x.sam.gz");
    let bgzf = match tc.draw(layout()) {
        None => run(Command::new("bgzip").arg("-c").arg(&sorted)),
        Some(layout) => reblock(&std::fs::read(&sorted).unwrap(), &layout),
    };
    std::fs::write(&sam, bgzf).unwrap();
    // htslib's tabix gives a SAM record the span of its M, D and N ops only
    // (`tbx_parse1`, tbx.c), so a read with `=` or `X` ops is indexed short and
    // a query of its last bases cannot find it — through htslib as much as
    // seqair, which can only follow the index it is given. Such inputs get a
    // `samtools index -b` index, which spans `=`/`X` (`bam_cigar2rlen`).
    let text = std::fs::read(&sorted).unwrap();
    let eq_or_x = text
        .split(|&b| b == b'\n')
        .filter(|l| !l.starts_with(b"@"))
        .filter_map(|l| l.split(|&b| b == b'\t').nth(5))
        .any(|cigar| cigar.iter().any(|&b| b == b'=' || b == b'X'));
    if !eq_or_x && tc.draw(gs::booleans()) {
        run(Command::new("tabix").args(["-p", "sam"]).arg(&sam));
    } else {
        run(Command::new("samtools").args(["index", "-b"]).arg(&sam));
    }
    Files { bam, sam }
}

// ── the two paths ───────────────────────────────────────────────────────

/// What a placed record holds, as the store shows it.
#[derive(Debug, Clone, PartialEq)]
struct View {
    qname: Vec<u8>,
    flags: u16,
    tid: Option<u32>,
    pos: Pos0,
    end_pos: Pos0,
    mapq: u8,
    cigar: Vec<u32>,
    seq: Vec<Base>,
    qual: Vec<u8>,
    aux: Vec<u8>,
    mate: MateFields,
}

/// A streamed record as the store should hold it; `None` when unplaced.
/// `end_pos` by the spec, computed here: `pos` for an unmapped record
/// (`r[record_store.end_pos_htslib]`) or one whose CIGAR spans no reference
/// (`r[bam.record.zero_refspan]`), else `pos` + reference length − 1.
fn expected(record: &RawBamRecord) -> Option<View> {
    Some(placed_view(record, record.tid()?.as_u32(), record.pos()?))
}

fn placed_view(record: &RawBamRecord, tid: u32, pos: Pos0) -> View {
    let cigar: Vec<u32> = record.cigar().map(|op| op.to_bam_u32()).collect();
    let rlen: u32 =
        cigar.iter().filter(|&&op| matches!(op & 0xF, 0 | 2 | 3 | 7 | 8)).map(|op| op >> 4).sum();
    let end_pos = if record.flags().is_unmapped() || rlen == 0 {
        pos
    } else {
        Pos0::new(pos.as_u32() + rlen - 1).unwrap()
    };
    let mut ascii = Vec::new();
    record.decode_seq_into(&mut ascii);
    let seq = ascii
        .iter()
        .map(|b| match b {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            _ => Base::Unknown,
        })
        .collect();
    // A missing QUAL is 0xFF throughout (`r[sam.edge.missing_qual]`).
    let qual = record
        .qual()
        .map_or_else(|| vec![0xFF; record.seq_len()], |q| BaseQuality::slice_to_bytes(q).to_vec());
    View {
        qname: record.qname().to_vec(),
        flags: record.flags().raw(),
        tid: Some(tid),
        pos,
        end_pos,
        mapq: record.mapq().as_byte(),
        cigar,
        seq,
        qual,
        aux: record.aux().as_bytes().to_vec(),
        mate: record.mate_fields(),
    }
}

fn store_views(store: &RecordStore) -> Vec<View> {
    store
        .records()
        .map(|r| View {
            qname: r.qname(store).unwrap().to_vec(),
            flags: r.flags.raw(),
            tid: r.tid.map(|t| t.as_u32()),
            pos: r.pos,
            end_pos: r.end_pos,
            mapq: r.mapq.as_byte(),
            cigar: r.cigar(store).unwrap().iter().map(|op| op.to_bam_u32()).collect(),
            seq: r.seq(store).unwrap().to_vec(),
            qual: BaseQuality::slice_to_bytes(r.qual(store).unwrap()).to_vec(),
            aux: r.aux(store).unwrap().as_bytes().to_vec(),
            mate: r.mate_fields,
        })
        .collect()
}

/// The streamed placed records on `tid` that overlap `span` by
/// `r[interval.overlap_test]`, in file order; none for an empty span
/// (`r[interval.empty_span]`).
fn overlapping(placed: &[View], tid: u32, span: RangeInclusive<Pos0>) -> Vec<View> {
    if span.last < span.start {
        return Vec::new();
    }
    placed
        .iter()
        .filter(|v| v.tid == Some(tid) && v.pos <= span.last && v.end_pos >= span.start)
        .cloned()
        .collect()
}

fn stream_payloads(records: &[RawBamRecord]) -> Vec<&[u8]> {
    records.iter().map(RawBamRecord::as_payload).collect()
}

/// Every record of the BAM, the SAM, and both through `AlignmentStream`,
/// which must agree; then the placed ones as views.
fn stream_all(tc: &TestCase, files: &Files) -> Vec<View> {
    // `None`: inflate inline.
    let threads = tc
        .draw(gs::optional(gs::integers::<u8>().min_value(1).max_value(3)))
        .map(|n| NonZeroU8::new(n).unwrap());
    let open = |path: &Path| {
        let builder = BgzfStream::builder(std::fs::File::open(path).unwrap());
        match threads {
            Some(n) => builder.threads(n).build(),
            None => builder.build(),
        }
    };
    let mut bam = BamStreamReader::from_bgzf(open(&files.bam)).unwrap();
    let bam_records: Vec<RawBamRecord> = bam.records().collect::<Result<_, _>>().unwrap();
    let mut sam = SamStreamReader::new(open(&files.sam)).unwrap();
    let sam_records: Vec<RawBamRecord> = sam.records().collect::<Result<_, _>>().unwrap();
    assert_eq!(bam.header().header_text(), sam.header().header_text(), "header");
    assert!(
        stream_payloads(&bam_records) == stream_payloads(&sam_records),
        "BAM and SAM streams differ"
    );
    for (path, format) in [(&files.bam, StreamFormat::Bam), (&files.sam, StreamFormat::BgzfSam)] {
        let file = std::fs::File::open(path).unwrap();
        let builder = AlignmentStream::builder(file);
        let mut stream = match threads {
            Some(n) => builder.threads(n).build(),
            None => builder.build(),
        }
        .unwrap();
        assert_eq!(stream.format(), format);
        let mut records = Vec::new();
        let mut record = RawBamRecord::default();
        while stream.read_record(&mut record).unwrap() {
            records.push(record.as_payload().to_vec());
        }
        assert!(records == stream_payloads(&bam_records), "AlignmentStream on {format:?} differs");
    }
    bam_records.iter().filter_map(expected).collect()
}

struct Indexed {
    bam: IndexedReader,
    sam: IndexedReader,
    store: RecordStore,
}

impl Indexed {
    fn open(files: &Files) -> Self {
        Self {
            bam: IndexedReader::open(&files.bam).unwrap(),
            sam: IndexedReader::open(&files.sam).unwrap(),
            store: RecordStore::new(),
        }
    }

    /// The region from the BAM and from the SAM, which must both be exactly
    /// the streamed records that overlap it, field by field and in order;
    /// returns what the store held.
    fn check(&mut self, placed: &[View], contig: usize, span: RangeInclusive<Pos0>) -> Vec<View> {
        let want = overlapping(placed, contig as u32, span);
        let mut got = Vec::new();
        for (format, reader) in [("BAM", &mut self.bam), ("SAM", &mut self.sam)] {
            let tid = reader.header().tid(&format!("c{contig}")).unwrap();
            let n = reader.fetch_into(tid, span, &mut self.store).unwrap();
            got = store_views(&self.store);
            assert_eq!(n, got.len());
            assert!(
                got == want,
                "{format} c{contig}:{}-{} (0-based, inclusive): got {:#?}\nwant {:#?}",
                span.start.as_u32(),
                span.last.as_u32(),
                names(&got),
                names(&want),
            );
        }
        got
    }
}

fn names(views: &[View]) -> Vec<String> {
    views
        .iter()
        .map(|v| {
            format!(
                "{} {}-{} flag {:#x}",
                v.qname.escape_ascii(),
                v.pos.as_u32(),
                v.end_pos.as_u32(),
                v.flags
            )
        })
        .collect()
}

fn span(start: u32, last: u32) -> RangeInclusive<Pos0> {
    RangeInclusive { start: Pos0::new(start).unwrap(), last: Pos0::new(last).unwrap() }
}

/// A region of a contig `len` long: whole, at either end, anywhere, empty,
/// or with an end on a streamed record's first or last base.
#[hegel::composite]
fn region(tc: &TestCase, len: u32, ends: Vec<(u32, u32)>) -> RangeInclusive<Pos0> {
    let last = len - 1;
    let int = |lo: u32, hi: u32| tc.draw(gs::integers::<u32>().min_value(lo).max_value(hi));
    match tc.draw(gs::integers::<u8>().max_value(6)) {
        0 => span(0, last),
        1 => span(int(0, last), last),
        2 => span(0, int(0, last)),
        3 => {
            let start = int(0, last);
            let width = tc.draw(gs::sampled_from(vec![0u32, 50, 2000, len]));
            span(start, int(start, last.min(start.saturating_add(width))))
        }
        // `last` before `start`; `start` may sit one past the contig.
        4 => {
            let start = int(1, len);
            span(start, int(0, start - 1))
        }
        _ if ends.is_empty() => span(0, last),
        k => {
            let (pos, end) = ends[tc.draw(gs::integers::<usize>().max_value(ends.len() - 1))];
            let reach = int(0, 300);
            if k == 5 {
                span(end, last.min(end + reach))
            } else {
                span(pos.saturating_sub(reach), pos)
            }
        }
    }
}

/// The names `samtools view` gives for a region, in order.
fn samtools_names(bam: &Path, contig: usize, span: RangeInclusive<Pos0>) -> Vec<String> {
    let region = format!("c{contig}:{}-{}", span.start.as_u32() + 1, span.last.as_u32() + 1);
    let out = run(Command::new("samtools").arg("view").arg(bam).arg(region));
    String::from_utf8(out)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect()
}

// r[verify bam.reader.fetch_into+2]
// r[verify bam.reader.overlap_filter]
// r[verify bam.reader.unmapped_skipped]
// r[verify sam.reader.overlap_filter+4]
// r[verify sam.reader.overlap_inclusive]
// r[verify sam.reader.sorted_order]
// r[verify interval.overlap_test]
// r[verify interval.empty_span]
// r[verify record_store.end_pos_htslib]
// r[verify bam.record.zero_refspan]
// r[verify unified.fetch_equivalence]
// r[verify unified.stream.detect]
// r[verify unified.stream.records]
// r[verify bam.stream]
// r[verify sam.stream]
// r[verify region_buf.fetch_into+2]
/// Thousands of records over several contigs — mapped, placed but unmapped
/// (with and without a stale CIGAR), unplaced, zero-span, IUPAC, every aux
/// type — sorted by samtools into a BAM and a bgzipped SAM, each in
/// samtools'/bgzip's BGZF blocks or random ones (cut anywhere, empty blocks
/// between), indexed by `samtools index` (BAI or CSI) and `tabix -p sam` or
/// `samtools index -b`:
///
/// - the four streams (BAM, SAM, `AlignmentStream` on each) yield the same
///   records;
/// - `fetch_into` on any region — whole contig, at either end, empty,
///   ending on a read's first or last base — returns from the BAM and from
///   the SAM exactly the streamed records that overlap it, in file order,
///   field by field, and never an unplaced one;
/// - over a partition of each contig into adjacent regions, every placed
///   record turns up, in each region it overlaps;
/// - samtools names the same records for a sample of the regions.
#[hegel::test(test_cases = pinned::cases(40))]
fn indexed_regions_are_the_streamed_records_that_overlap_them(tc: TestCase) {
    let lens = tc.draw(contigs());
    let input = tc.draw(sam_input(lens.clone()));
    let dir = tempfile::tempdir().unwrap();
    let files = files(&tc, dir.path(), &input);
    let placed = stream_all(&tc, &files);
    let mut indexed = Indexed::open(&files);

    let mut samtools_checks = 0;
    for (contig, &len) in lens.iter().enumerate() {
        let on_contig: Vec<&View> =
            placed.iter().filter(|v| v.tid == Some(contig as u32)).collect();
        let ends: Vec<(u32, u32)> =
            on_contig.iter().map(|v| (v.pos.as_u32(), v.end_pos.as_u32())).collect();

        for _ in 0..tc.draw(gs::integers::<usize>().max_value(5)) {
            let span = tc.draw(region(len, ends.clone()));
            let want = indexed.check(&placed, contig, span);
            if samtools_checks < 4 && span.start <= span.last {
                samtools_checks += 1;
                let want: Vec<String> =
                    want.iter().map(|v| String::from_utf8(v.qname.clone()).unwrap()).collect();
                assert_eq!(samtools_names(&files.bam, contig, span), want, "samtools");
            }
        }

        let mut cuts: Vec<u32> = (0..tc.draw(gs::integers::<usize>().max_value(5)))
            .map(|_| tc.draw(gs::integers::<u32>().min_value(1).max_value(len)))
            .collect();
        cuts.push(0);
        cuts.push(len);
        cuts.sort_unstable();
        cuts.dedup();
        let mut seen = BTreeSet::new();
        for part in cuts.windows(2) {
            for v in indexed.check(&placed, contig, span(part[0], part[1] - 1)) {
                seen.insert(v.qname);
            }
        }
        let all: BTreeSet<Vec<u8>> = on_contig.iter().map(|v| v.qname.clone()).collect();
        assert!(seen == all, "c{contig}: a partition misses or invents records");
    }
}
