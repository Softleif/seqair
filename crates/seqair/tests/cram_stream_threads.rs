//! A CRAM `AlignmentStream` decoding on threads against the same stream
//! decoding inline: the same records and errors in the same order, however
//! the input arrives, on whichever pool, whatever fields are read.
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

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::RawBamRecord;
use seqair::cram::container::ContainerHeader;
use seqair::cram::stream::{CramContainerDecoder, CramStreamError, CramStreamReader};
use seqair::reader::{
    AlignmentStream, AlignmentStreamError, AuxTags, CramInput, PartialRecord, RecordFields,
    RecordKind,
};
use std::io::{Read, Write};
use std::num::NonZeroU8;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, OnceLock};

fn data(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data")).join(name)
}

/// The first 3000 records of `test.cram` rewritten by samtools with 100
/// records per container (30 containers), with names stored and without.
fn many_containers() -> &'static [PathBuf; 2] {
    static FILES: OnceLock<(tempfile::TempDir, [PathBuf; 2])> = OnceLock::new();
    let (_, files) = FILES.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let fasta = data("test.fasta.gz");
        let sam = Command::new("samtools")
            .args(["view", "-h", "--no-PG", "-T"])
            .arg(&fasta)
            .arg(data("test.cram"))
            .output()
            .unwrap();
        assert!(sam.status.success());
        let text = String::from_utf8(sam.stdout).unwrap();
        let (header, body): (Vec<&str>, Vec<&str>) = text.lines().partition(|l| l.starts_with('@'));
        let mut sam =
            header.iter().chain(body.iter().take(3000)).copied().collect::<Vec<_>>().join("\n");
        sam.push('\n');
        let write = |name: &str, opts: &str| {
            let out = dir.path().join(name);
            let mut child = Command::new("samtools")
                .args(["view", "--no-PG", "-C", "-T"])
                .arg(&fasta)
                .args(["--output-fmt-option", "seqs_per_slice=100"])
                .args(["--output-fmt-option", opts])
                .arg("-o")
                .arg(&out)
                .arg("-")
                .stdin(Stdio::piped())
                .spawn()
                .unwrap();
            child.stdin.take().unwrap().write_all(sam.as_bytes()).unwrap();
            assert!(child.wait().unwrap().success());
            out
        };
        let files = [write("named.cram", "lossy_names=0"), write("lossy.cram", "lossy_names=1")];
        (dir, files)
    });
    files
}

/// CRAM files with their reference, the many-container ones first.
fn inputs() -> Vec<(PathBuf, PathBuf)> {
    let fidelity = data("cram-stream/fidelity.fa");
    let mut out: Vec<_> =
        many_containers().iter().map(|f| (f.clone(), data("test.fasta.gz"))).collect();
    for f in ["fidelity_slices", "fidelity_multiref", "fidelity_lossy_names", "onlyu"] {
        out.push((data(&format!("cram-stream/{f}.cram")), fidelity.clone()));
    }
    out
}

/// Hands out the bytes in reads of the given sizes, in turn, like a pipe.
struct Trickle {
    bytes: Arc<[u8]>,
    at: usize,
    chunks: Vec<usize>,
    next: usize,
}

impl Read for Trickle {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let chunk = self.chunks[self.next % self.chunks.len()];
        self.next += 1;
        let n = buf.len().min(chunk).min(self.bytes.len() - self.at);
        buf[..n].copy_from_slice(&self.bytes[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

/// The byte offset where each container starts, then the end of the file.
fn container_bounds(bytes: &[u8]) -> Vec<usize> {
    let mut bounds = vec![26usize];
    while *bounds.last().unwrap() < bytes.len() {
        let at = *bounds.last().unwrap();
        let h = ContainerHeader::parse(&bytes[at..]).unwrap();
        bounds.push(at + h.header_size + h.length as usize);
    }
    bounds
}

#[derive(Debug, Clone, Copy)]
enum Damage {
    None,
    /// A flipped byte in a data container's body: that container fails.
    Body,
    /// A flipped byte in a container header: the stream ends there.
    Header,
}

#[derive(Debug, Clone)]
enum Threads {
    /// `threads(n)`, on the global pool.
    Global(u8),
    /// `threads_in(n, pool)` with a pool of `size` threads, read from one of
    /// the pool's own threads if `inside`.
    Pool { n: u8, size: usize, inside: bool },
}

/// What the stream gave, item by item: a record, or an error's `Debug`.
type Items<T> = Vec<Result<T, String>>;

fn is_poisoned(e: &AlignmentStreamError) -> bool {
    matches!(e, AlignmentStreamError::Cram { source: CramStreamError::Poisoned })
}

/// Every item of `stream`: through `records()`, or through `read_record`
/// into one reused record (until the end or `Poisoned`).
fn items<K: RecordKind, T>(
    mut stream: AlignmentStream<K>,
    reuse: bool,
    snap: impl Fn(&K::Record) -> T,
) -> Items<T> {
    if !reuse {
        return stream
            .records()
            .map(|r| r.map(|r| snap(&r)).map_err(|e| format!("{e:?}")))
            .collect();
    }
    let mut out = Vec::new();
    let mut record = K::Record::default();
    loop {
        match stream.read_record(&mut record) {
            Ok(true) => out.push(Ok(snap(&record))),
            Ok(false) => return out,
            Err(e) => {
                out.push(Err(format!("{e:?}")));
                if is_poisoned(&e) {
                    return out;
                }
            }
        }
    }
}

struct Case {
    bytes: Arc<[u8]>,
    fasta: PathBuf,
    chunks: Vec<usize>,
    reuse: bool,
}

impl Case {
    fn builder(&self) -> seqair::reader::AlignmentStreamBuilder<Trickle> {
        let input =
            Trickle { bytes: Arc::clone(&self.bytes), at: 0, chunks: self.chunks.clone(), next: 0 };
        AlignmentStream::builder(input).reference(&self.fasta)
    }

    /// The items of the stream `fields` and `threads` configure.
    fn read<T: Send>(
        &self,
        fields: Option<&RecordFields>,
        threads: Option<&Threads>,
        snap_all: fn(&RawBamRecord) -> T,
        snap_partial: fn(&PartialRecord) -> T,
    ) -> Items<T> {
        let run = |pool: Option<(u8, &Arc<rayon_core::ThreadPool>)>, global: Option<u8>| {
            let builder = self.builder();
            let n = |n: u8| NonZeroU8::new(n).unwrap();
            macro_rules! configure {
                ($b:expr) => {
                    match (pool, global) {
                        (Some((k, pool)), _) => $b.threads_in(n(k), Arc::clone(pool)),
                        (None, Some(k)) => $b.threads(n(k)),
                        (None, None) => $b,
                    }
                };
            }
            match fields {
                None => items(configure!(builder).build().unwrap(), self.reuse, snap_all),
                Some(f) => {
                    let b = builder.fields(f.clone());
                    items(configure!(b).build().unwrap(), self.reuse, snap_partial)
                }
            }
        };
        match threads {
            None => run(None, None),
            Some(Threads::Global(k)) => run(None, Some(*k)),
            Some(Threads::Pool { n, size, inside }) => {
                let pool = Arc::new(
                    rayon_core::ThreadPoolBuilder::new().num_threads(*size).build().unwrap(),
                );
                if *inside {
                    pool.install(|| run(Some((*n, &pool)), None))
                } else {
                    run(Some((*n, &pool)), None)
                }
            }
        }
    }
}

#[hegel::composite]
fn arb_threads(tc: &TestCase) -> Threads {
    let n = tc.draw(gs::integers::<u8>().min_value(1).max_value(8));
    if tc.draw(gs::booleans()) {
        Threads::Global(n)
    } else {
        let size = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
        Threads::Pool { n, size, inside: tc.draw(gs::booleans()) }
    }
}

#[hegel::composite]
fn arb_fields(tc: &TestCase) -> Option<RecordFields> {
    match tc.draw(gs::integers::<u8>().max_value(3)) {
        0 => None,
        1 => Some(RecordFields::all()),
        2 => Some(RecordFields::core()),
        _ => {
            Some(RecordFields::core().with_names().with_quality().with_aux(AuxTags::only([*b"RG"])))
        }
    }
}

// r[verify unified.stream.cram_threads]
// r[verify unified.stream.iter]
// r[verify bgzf.stream.current_pool]
// r[verify bgzf.stream.shared_pool]
/// Decoding on any number of threads, on the global pool or a given one
/// (read from inside it too), through a pipe that delivers any number of
/// bytes per read, gives exactly the inline stream's items: the same
/// records, `PartialRecord`s included, and for a damaged container the same
/// error at the same place, then the same continuation or end.
#[hegel::test(test_cases = pinned::cases(48))]
fn threaded_cram_streams_give_the_inline_items(tc: TestCase) {
    let (cram, fasta) = tc.draw(gs::sampled_from(inputs()));
    let mut bytes = std::fs::read(&cram).unwrap();
    let bounds = container_bounds(&bytes);
    // Data containers are 1..last: the first is the file header, the last EOF.
    let damage = tc.draw(gs::sampled_from(vec![Damage::None, Damage::Body, Damage::Header]));
    let data_containers = bounds.len() - 3;
    if !matches!(damage, Damage::None) && data_containers > 0 {
        let container = tc.draw(gs::integers::<usize>().min_value(1).max_value(data_containers));
        let start = bounds[container];
        let header_size = ContainerHeader::parse(&bytes[start..]).unwrap().header_size;
        let range = match damage {
            Damage::Body => start + header_size..bounds[container + 1],
            _ => start..start + header_size,
        };
        let at = tc.draw(gs::integers::<usize>().min_value(range.start).max_value(range.end - 1));
        bytes[at] ^= tc.draw(gs::integers::<u8>().min_value(1));
    }
    let case = Case {
        bytes: bytes.into(),
        fasta,
        chunks: tc.draw(
            gs::vecs(gs::integers::<usize>().min_value(1).max_value(1 << 16))
                .min_size(1)
                .max_size(4),
        ),
        reuse: tc.draw(gs::booleans()),
    };
    let fields = tc.draw(arb_fields());
    let threads = tc.draw(arb_threads());

    let payload = |r: &RawBamRecord| Ok(r.as_payload().to_vec());
    let partial = |r: &PartialRecord| Err(r.clone());
    type Snap = Result<Vec<u8>, PartialRecord>;
    let want: Items<Snap> = case.read(fields.as_ref(), None, payload, partial);
    let got: Items<Snap> = case.read(fields.as_ref(), Some(&threads), payload, partial);
    assert_eq!(got.len(), want.len(), "{damage:?} {threads:?}");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        assert!(g == w, "item {i} of {}: {damage:?} {threads:?}\n{g:?}\n{w:?}", want.len());
    }
    if want.iter().any(Result::is_err) {
        tc.event("an error");
    }
}

/// `into_cram`'s reader and decoder decode the rest: the records `read`
/// had not given.
fn decode_rest(
    (mut reader, decoder): (CramStreamReader<CramInput>, CramContainerDecoder),
) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut decoder = decoder.fork().unwrap();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut out).unwrap();
    }
    out.into_iter().map(RawBamRecord::into_payload).collect()
}

// r[verify unified.stream.into_cram+2]
// r[verify unified.stream.cram_threads]
/// A threaded stream hands over before its first read and after its end;
/// in between it holds containers read ahead and comes back unchanged, and
/// reading on gives every record. Dropping it mid-stream does not hang.
#[test]
fn a_threaded_stream_hands_over_only_when_nothing_is_read_ahead() {
    let cram = &many_containers()[0];
    let fasta = data("test.fasta.gz");
    let four = NonZeroU8::new(4).unwrap();
    let open = || {
        AlignmentStream::open_path(cram).unwrap().reference(&fasta).threads(four).build().unwrap()
    };
    let inline: Vec<Vec<u8>> = AlignmentStream::open_path(cram)
        .unwrap()
        .reference(&fasta)
        .build()
        .unwrap()
        .records()
        .map(|r| r.unwrap().into_payload())
        .collect();
    assert_eq!(inline.len(), 3000);

    assert_eq!(decode_rest(open().into_cram().unwrap()), inline);

    let mut stream = open();
    let mut got = Vec::new();
    let mut record = RawBamRecord::default();
    for _ in 0..150 {
        assert!(stream.read_record(&mut record).unwrap());
        got.push(record.as_payload().to_vec());
        stream = stream.into_cram().expect_err("containers are read ahead");
    }
    while stream.read_record(&mut record).unwrap() {
        got.push(record.as_payload().to_vec());
    }
    assert_eq!(got, inline);
    let (mut reader, _) = stream.into_cram().unwrap();
    assert!(reader.next_container().unwrap().is_none());

    let mut stream = open();
    for _ in 0..10 {
        assert!(stream.read_record(&mut record).unwrap());
    }
    drop(stream);
}

// r[verify unified.stream.cram_threads]
/// A pool of one thread, read from that thread, with more decoders asked
/// for than it has: every job still runs, on the reading thread meanwhile.
#[test]
fn one_pool_thread_read_from_inside_does_not_deadlock() {
    let cram = &many_containers()[1];
    let fasta = data("test.fasta.gz");
    let pool = Arc::new(rayon_core::ThreadPoolBuilder::new().num_threads(1).build().unwrap());
    let eight = NonZeroU8::new(8).unwrap();
    let count = pool.install(|| {
        let mut stream = AlignmentStream::open_path(cram)
            .unwrap()
            .reference(&fasta)
            .threads_in(eight, Arc::clone(&pool))
            .build()
            .unwrap();
        stream.records().map(Result::unwrap).count()
    });
    assert_eq!(count, 3000);
}
