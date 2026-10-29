//! The CRAM stream on damaged input, a wrong reference, a pipe that
//! delivers any number of bytes per read, and sizes past the allocation
//! limit, on the `tests/data` fixtures.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::RawBamRecord;
use seqair::cram::container::ContainerHeader;
use seqair::cram::reader::{CramError, MAX_ALLOC_SIZE};
use seqair::cram::stream::{CramContainerDecoder, CramStreamError, CramStreamReader};
use seqair::cram::varint::decode_itf8;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::io::Read;
use std::path::{Path, PathBuf};

/// Records the largest single allocation made on the current thread.
struct Tracking;

thread_local! {
    static LARGEST: Cell<usize> = const { Cell::new(0) };
}

fn track(size: usize) {
    // `try_with`: allocations during thread teardown find no slot.
    if let Err(error) = LARGEST.try_with(|l| l.set(l.get().max(size))) {
        std::hint::black_box(error);
    }
}

// SAFETY: every call forwards to `System` with the caller's arguments.
unsafe impl GlobalAlloc for Tracking {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        track(layout.size());
        // SAFETY: forwarded unchanged.
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        track(layout.size());
        // SAFETY: forwarded unchanged.
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: forwarded unchanged.
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        track(new_size);
        // SAFETY: forwarded unchanged.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Tracking = Tracking;

/// `f`'s result and the largest allocation it made.
fn largest_allocation<T>(f: impl FnOnce() -> T) -> (T, usize) {
    LARGEST.with(|l| l.set(0));
    let out = f();
    (out, LARGEST.with(Cell::get))
}

fn data(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data")).join(name)
}

fn fidelity_fa() -> PathBuf {
    data("cram-stream/fidelity.fa")
}

/// Fixtures with their reference.
fn fixtures() -> Vec<(PathBuf, PathBuf)> {
    let mut out = vec![
        (data("test.cram"), data("test.fasta.gz")),
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
        out.push((data(&format!("cram-stream/{f}.cram")), fidelity_fa()));
    }
    out
}

/// Fixtures whose every slice decodes against `fidelity.fa` and checks it
/// by MD5: one reference per slice, the reference external.
const MD5_CHECKED: [&str; 4] =
    ["fidelity_v31", "fidelity_v30", "fidelity_slices", "fidelity_lossy_names"];

/// Hands out the bytes in reads of the given sizes, in turn, like a pipe.
struct Trickle<'a> {
    bytes: &'a [u8],
    chunks: Vec<usize>,
    next: usize,
}

impl Read for Trickle<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let chunk = self.chunks[self.next % self.chunks.len()];
        self.next += 1;
        let n = buf.len().min(chunk);
        self.bytes.read(&mut buf[..n])
    }
}

/// Every record of `input`, decoded against `fasta`.
fn records(input: impl Read, fasta: &Path) -> Vec<RawBamRecord> {
    let mut reader = CramStreamReader::new(input).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap();
    let mut out = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut out).unwrap();
    }
    out
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

#[hegel::composite]
fn arb_fixture(tc: &TestCase) -> (PathBuf, PathBuf) {
    tc.draw(gs::sampled_from(fixtures()))
}

// r[verify cram.stream.container_header+2]
/// A reader that delivers 1 to 4096 bytes per call, in any sequence, gives
/// the records of reading the file in one go.
#[hegel::test(test_cases = pinned::cases(64))]
fn any_read_sizes_give_the_same_records(tc: TestCase) {
    let (cram, fasta) = tc.draw(arb_fixture());
    let chunks = tc.draw(
        gs::vecs(gs::integers::<usize>().min_value(1).max_value(4096)).min_size(1).max_size(8),
    );
    let bytes = std::fs::read(&cram).unwrap();
    let want = records(&bytes[..], &fasta);
    let got = records(Trickle { bytes: &bytes, chunks, next: 0 }, &fasta);
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        assert!(g.as_payload() == w.as_payload());
    }
}

// r[verify cram.stream.framing+2]
// r[verify cram.container.header]
/// A flipped byte in a container header — CRC-covered, the CRC included —
/// fails that container's read; every later read is `Poisoned`. Each
/// container before it still reads.
#[hegel::test]
fn a_damaged_container_header_fails_and_poisons(tc: TestCase) {
    let (cram, _) = tc.draw(arb_fixture());
    let mut bytes = std::fs::read(&cram).unwrap();
    let bounds = container_bounds(&bytes);
    let container = tc.draw(gs::integers::<usize>().max_value(bounds.len() - 2));
    let start = bounds[container];
    let header_size = ContainerHeader::parse(&bytes[start..]).unwrap().header_size;
    let at = start + tc.draw(gs::integers::<usize>().max_value(header_size - 1));
    bytes[at] ^= tc.draw(gs::integers::<u8>().min_value(1));

    let Ok(mut reader) = CramStreamReader::new(&bytes[..]) else {
        assert_eq!(container, 0, "only a damaged file header container fails the open");
        return;
    };
    assert_ne!(container, 0, "a damaged file header container opened");
    let mut read = 1;
    let err = loop {
        match reader.next_container() {
            Ok(Some(c)) => read = c.number() + 1,
            Ok(None) => panic!("the damage in container {container} went unnoticed"),
            Err(e) => break e,
        }
    };
    assert!(!matches!(err, CramStreamError::Poisoned), "{err:?}");
    assert_eq!(read, container as u64, "containers before the damage: {err:?}");
    assert!(matches!(reader.next_container(), Err(CramStreamError::Poisoned)));
    assert!(matches!(reader.next_container(), Err(CramStreamError::Poisoned)));
}

/// Decode every container of `bytes` against `fasta`, each into an output
/// that already holds `held` records: a failed decode must leave them as
/// they were. The number of failed decodes.
fn decode_keeping_output(bytes: &[u8], fasta: &Path, held: &[RawBamRecord]) -> usize {
    let mut reader = CramStreamReader::new(bytes).unwrap();
    let mut decoder: CramContainerDecoder = reader.decoder(fasta).unwrap();
    let mut failed = 0;
    while let Some(container) = reader.next_container().unwrap() {
        let mut out = held.to_vec();
        if decoder.decode(&container, &mut out).is_err() {
            failed += 1;
            assert_eq!(out.len(), held.len(), "container {}", container.number());
            assert!(out.iter().zip(held).all(|(o, h)| o.as_payload() == h.as_payload()));
        }
    }
    failed
}

// r[verify cram.stream.decoder]
// r[verify cram.edge.reference_mismatch]
/// Against a FASTA with some bases changed, every container either decodes
/// or fails leaving the output as it was, and one whose slice covers a
/// changed base with a mapped read fails (its MD5 does not match).
#[hegel::test]
fn a_wrong_reference_fails_without_touching_the_output(tc: TestCase) {
    let name = tc.draw(gs::sampled_from(MD5_CHECKED.to_vec()));
    let cram = data(&format!("cram-stream/{name}.cram"));
    let bytes = std::fs::read(&cram).unwrap();
    let want = records(&bytes[..], &fidelity_fa());

    // One sequence line per contig, so a byte changed in place keeps the
    // `.fai` valid.
    let mut fa = std::fs::read(fidelity_fa()).unwrap();
    let fai = std::fs::read_to_string(data("cram-stream/fidelity.fa.fai")).unwrap();
    let contigs: Vec<(String, usize, usize)> = fai
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_owned(), f[1].parse().unwrap(), f[2].parse().unwrap())
        })
        .collect();
    let contig = tc.draw(gs::integers::<usize>().max_value(contigs.len() - 1));
    let (contig_name, len, offset) = &contigs[contig];
    let start = tc.draw(gs::integers::<usize>().max_value(len - 1));
    let run = tc.draw(gs::integers::<usize>().min_value(1).max_value(20)).min(len - start);
    let shift = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    for b in &mut fa[offset + start..offset + start + run] {
        let i = b"ACGT".iter().position(|&a| a == b.to_ascii_uppercase()).unwrap_or(0);
        *b = b"ACGT"[(i + shift) % 4];
    }
    let changed = |pos: usize| (start..start + run).contains(&pos);
    let dir = tempfile::tempdir().unwrap();
    let wrong = dir.path().join("wrong.fa");
    std::fs::write(&wrong, &fa).unwrap();
    std::fs::write(dir.path().join("wrong.fa.fai"), &fai).unwrap();

    let header = CramStreamReader::new(&bytes[..]).unwrap().into_header();
    let covered = want.iter().any(|r| {
        let (Some(tid), Some(pos), Some(end)) = (r.tid(), r.pos(), r.end_pos()) else {
            return false;
        };
        !r.flags().is_unmapped()
            && header.target_name(tid) == Some(contig_name.as_str())
            && (pos.as_u64() as usize..=end.as_u64() as usize).any(changed)
    });
    let held = &want[..tc.draw(gs::integers::<usize>().max_value(3)).min(want.len())];
    let failed = decode_keeping_output(&bytes, &wrong, held);
    if covered {
        assert!(failed > 0, "a changed base under a mapped read went unnoticed");
    }
}

// r[verify cram.stream.decoder]
// r[verify cram.block.crc32]
/// A flipped byte in a container's body — its blocks, CRC-covered — fails
/// that container's decode and leaves the output as it was.
#[hegel::test]
fn a_damaged_container_body_fails_without_touching_the_output(tc: TestCase) {
    // `test.cram` takes a second to decode unoptimised.
    let small: Vec<_> = fixtures().into_iter().filter(|(c, _)| !c.ends_with("test.cram")).collect();
    let (cram, fasta) = tc.draw(gs::sampled_from(small));
    let mut bytes = std::fs::read(&cram).unwrap();
    let want = records(&bytes[..], &fasta);
    let bounds = container_bounds(&bytes);
    // Data containers only: the first is the file header, the last EOF.
    let container = tc.draw(gs::integers::<usize>().min_value(1).max_value(bounds.len() - 3));
    let start = bounds[container];
    let header_size = ContainerHeader::parse(&bytes[start..]).unwrap().header_size;
    let body = start + header_size..bounds[container + 1];
    let at = tc.draw(gs::integers::<usize>().min_value(body.start).max_value(body.end - 1));
    bytes[at] ^= tc.draw(gs::integers::<u8>().min_value(1));
    let held = &want[..tc.draw(gs::integers::<usize>().max_value(3)).min(want.len())];
    assert_eq!(decode_keeping_output(&bytes, &fasta, held), 1);
}

/// ITF8 in its five-byte form, which holds any `u32`.
fn itf8(value: u32) -> [u8; 5] {
    let out = [
        0xF0 | (value >> 28) as u8,
        (value >> 20) as u8,
        (value >> 12) as u8,
        (value >> 4) as u8,
        (value & 0x0F) as u8,
    ];
    assert_eq!(decode_itf8(&out), Some((value, 5)));
    out
}

/// A container header of `length` body bytes on reference 0 at 1 holding
/// `records` records in `blocks` blocks, its landmark count the given ITF8
/// bytes and no landmarks follow; CRC32 right.
fn container_header(length: i32, records: u8, blocks: u8, landmark_count: &[u8]) -> Vec<u8> {
    // reference 0, start 1, span 0, records, record counter 0, bases 0, blocks
    let mut bytes = length.to_le_bytes().to_vec();
    bytes.extend_from_slice(&[0, 1, 0, records, 0, 0, blocks]);
    bytes.extend_from_slice(landmark_count);
    let mut crc = libdeflater::Crc::new();
    crc.update(&bytes);
    bytes.extend_from_slice(&crc.sum().to_le_bytes());
    bytes
}

/// A raw compression header block of `data`, claiming `uncompressed` bytes;
/// CRC32 right.
fn compression_header_block(compressed: u32, uncompressed: u32, data: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0, 1, 0];
    bytes.extend_from_slice(&itf8(compressed));
    bytes.extend_from_slice(&itf8(uncompressed));
    bytes.extend_from_slice(data);
    let mut crc = libdeflater::Crc::new();
    crc.update(&bytes);
    bytes.extend_from_slice(&crc.sum().to_le_bytes());
    bytes
}

/// The file definition and header container of a fixture, then `rest`.
fn after_file_header(rest: &[u8]) -> Vec<u8> {
    let bytes = std::fs::read(data("cram-stream/fidelity_v31.cram")).unwrap();
    let mut file = bytes[..container_bounds(&bytes)[1]].to_vec();
    file.extend_from_slice(rest);
    file
}

/// Far below every size a test claims.
const SMALL: usize = 1 << 20;

#[derive(Debug, Clone, Copy)]
enum Oversized {
    /// A container length past the limit.
    ContainerLength(i32),
    /// A landmark count whose array is past the limit.
    LandmarkCount(u32),
    /// A container length within the limit, the input ending after `got`
    /// body bytes.
    ShortBody { length: i32, got: usize },
    /// A block that decodes to more than the limit.
    BlockSize(u32),
    /// A block that stores more bytes than its container holds.
    ShortBlock(u32),
}

#[hegel::composite]
fn arb_oversized(tc: &TestCase) -> Oversized {
    let limit = MAX_ALLOC_SIZE as u32;
    match tc.draw(gs::integers::<u8>().max_value(4)) {
        0 => Oversized::ContainerLength(
            tc.draw(gs::integers::<i32>().min_value(MAX_ALLOC_SIZE as i32 + 1)),
        ),
        1 => Oversized::LandmarkCount(
            tc.draw(gs::integers::<u32>().min_value(limit / 4 + 1).max_value(i32::MAX as u32)),
        ),
        2 => {
            let length =
                tc.draw(gs::integers::<i32>().min_value(1).max_value(MAX_ALLOC_SIZE as i32));
            let got = tc.draw(gs::integers::<usize>().max_value((length as usize - 1).min(4096)));
            Oversized::ShortBody { length, got }
        }
        3 => Oversized::BlockSize(tc.draw(gs::integers::<u32>().min_value(limit + 1))),
        _ => Oversized::ShortBlock(tc.draw(gs::integers::<u32>().min_value(64))),
    }
}

// r[verify io.fuzz.alloc_limits]
// r[verify cram.stream.framing+2]
/// A container or block size past the allocation limit, or past the bytes
/// the input holds, is an error, and nothing near that size is allocated.
#[hegel::test]
fn sizes_past_the_limit_are_errors_not_allocations(tc: TestCase) {
    let case = tc.draw(arb_oversized());
    let file = match case {
        Oversized::ContainerLength(length) => {
            after_file_header(&container_header(length, 1, 1, &[0]))
        }
        Oversized::LandmarkCount(count) => {
            after_file_header(&container_header(0, 1, 1, &itf8(count)))
        }
        Oversized::ShortBody { length, got } => {
            let mut rest = container_header(length, 1, 1, &[0]);
            rest.resize(rest.len() + got, 0);
            after_file_header(&rest)
        }
        Oversized::BlockSize(size) | Oversized::ShortBlock(size) => {
            let block = match case {
                Oversized::BlockSize(_) => compression_header_block(0, size, &[]),
                _ => compression_header_block(size, size, &[0; 16]),
            };
            let mut rest = container_header(block.len() as i32, 1, 1, &[0]);
            rest.extend_from_slice(&block);
            after_file_header(&rest)
        }
    };
    let mut reader = CramStreamReader::new(&file[..]).unwrap();
    let mut decoder = reader.decoder(&fidelity_fa()).unwrap();
    let (err, largest) = largest_allocation(|| match reader.next_container() {
        Ok(Some(container)) => decoder.decode(&container, &mut Vec::new()).unwrap_err(),
        Ok(None) => panic!("{case:?}: the stream ended"),
        Err(e) => e,
    });
    assert!(largest < SMALL, "{case:?}: allocated {largest} bytes for {err:?}");
    let alloc = |e: &CramStreamError| {
        matches!(e, CramStreamError::Cram { source: CramError::AllocationTooLarge { .. } })
    };
    match case {
        Oversized::ContainerLength(_) | Oversized::LandmarkCount(_) | Oversized::BlockSize(_) => {
            assert!(alloc(&err), "{case:?}: {err:?}");
        }
        Oversized::ShortBody { .. } => {
            assert!(matches!(err, CramStreamError::TruncatedContainer { container: 1 }), "{err:?}");
        }
        Oversized::ShortBlock(_) => {
            assert!(
                matches!(err, CramStreamError::Cram { source: CramError::Truncated { .. } }),
                "{case:?}: {err:?}"
            );
        }
    }
}
