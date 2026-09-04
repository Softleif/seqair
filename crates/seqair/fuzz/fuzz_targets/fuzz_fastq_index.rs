//! Fuzz indexed FASTA/FASTQ fetching against an adversarial index.
//!
//! `fuzz_fasta_index` covers FAI *parsing* only. This target goes a step
//! further and drives `fetch_seq_into`, because the interesting attack surface
//! is downstream of the parser: `byte_offset` multiplies caller-controlled
//! positions by index-supplied `linewidth`, and the difference between the
//! first and last byte offset sizes a `Vec` before a single byte is read. An
//! index claiming a huge `linewidth` therefore turns a small, in-bounds query
//! into a very large allocation.
//!
//! The index text is assembled from structured fields rather than taken as raw
//! bytes, so the fuzzer spends its budget on offset arithmetic instead of on
//! rediscovering the tab-separated grammar. `qual_offset` is generated too,
//! which is what makes this a FASTQ target: six columns exercise
//! `r[fastq.access.index_parse]` and `qual_byte_offset`.
#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::Pos0;
use seqair::fasta::{FastaIndex, IndexedFastaReader};

#[derive(Arbitrary, Debug)]
struct Entry {
    name: String,
    length: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
    /// `Some` makes this a 6-column FASTQ entry.
    qual_offset: Option<u64>,
}

#[derive(Arbitrary, Debug)]
struct Input {
    entries: Vec<Entry>,
    /// Becomes a GZI index: `(compressed_offset, uncompressed_offset)` pairs.
    gzi_entries: Vec<(u64, u64)>,
    /// Stands in for the BGZF payload. It need not be valid: the allocation
    /// under test happens before any block is read.
    bgzf_data: Vec<u8>,
    queries: Vec<(u8, u32, u32)>,
}

/// Tabs and newlines would just produce parse errors, so they are stripped —
/// the grammar is not what is being fuzzed here.
fn sanitize(name: &str) -> String {
    let cleaned: String =
        name.chars().filter(|c| *c != '\t' && *c != '\n' && *c != '\r').take(64).collect();
    if cleaned.is_empty() {
        "seq".to_owned()
    } else {
        cleaned
    }
}

fn build_fai(entries: &[Entry]) -> (String, Vec<String>) {
    let mut text = String::new();
    let mut names = Vec::new();
    for e in entries.iter().take(16) {
        let name = sanitize(&e.name);
        text.push_str(&name);
        text.push('\t');
        text.push_str(&format!("{}\t{}\t{}\t{}", e.length, e.offset, e.linebases, e.linewidth));
        if let Some(q) = e.qual_offset {
            text.push('\t');
            text.push_str(&q.to_string());
        }
        text.push('\n');
        names.push(name);
    }
    (text, names)
}

fn build_gzi(entries: &[(u64, u64)]) -> Vec<u8> {
    // Entries must be strictly ascending by uncompressed offset or the parser
    // rejects the index before any fetch can run.
    let mut sorted: Vec<(u64, u64)> = entries.iter().copied().take(64).collect();
    sorted.sort_unstable_by_key(|(_, u)| *u);
    sorted.dedup_by_key(|(_, u)| *u);

    let mut out = Vec::with_capacity(8 + sorted.len() * 16);
    out.extend_from_slice(&(sorted.len() as u64).to_le_bytes());
    for (c, u) in &sorted {
        out.extend_from_slice(&c.to_le_bytes());
        out.extend_from_slice(&u.to_le_bytes());
    }
    out
}

fuzz_target!(|input: Input| {
    let (fai_text, names) = build_fai(&input.entries);
    if names.is_empty() {
        return;
    }

    // The parser is exercised on its own too: a rejected index must be a clean
    // error, never a panic.
    let Ok(index) = FastaIndex::from_contents(&fai_text) else {
        return;
    };
    // Offset arithmetic must stay total for any in-range position.
    for name in &names {
        if let Some(entry) = index.get(name) {
            for pos in [0u64, 1, entry.length.saturating_sub(1)] {
                let _ = entry.byte_offset(pos);
                let _ = entry.qual_byte_offset(pos);
            }
        }
    }

    let gzi = build_gzi(&input.gzi_entries);
    let Ok(mut reader) = IndexedFastaReader::from_bytes(input.bgzf_data, &fai_text, &gzi) else {
        return;
    };

    let mut buf = Vec::new();
    for (which, start, stop) in input.queries.iter().take(32) {
        let Some(name) = names.get(usize::from(*which) % names.len()) else {
            continue;
        };
        let (Some(a), Some(b)) = (Pos0::new(*start), Pos0::new(*stop)) else {
            continue;
        };
        // Errors are fine and expected; a panic, an abort, or an unbounded
        // allocation is not.
        let _ = reader.fetch_seq_into(name, a, b, &mut buf);
    }
});
