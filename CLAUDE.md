# seqair

Pure-Rust BAM/SAM/CRAM/FASTA reader + pileup engine. I/O backend for [rastair](https://github.com/bsblabludwig/rastair).

Workspace: `crates/seqair` (readers, pileup, BGZF, CRAM) and `crates/seqair-types` (Base, Strand, Phred, Probability, RmsAccumulator, RegionString).

## Tracey specs

Specs in `docs/spec/*.md` with `r[rule.id]`. Code: `// r[impl rule.id]`, tests: `// r[verify rule.id]`. Use the `tracey` skill. Add spec rules before implementing; update specs when changing behavior.

## Coding style

Expert Rust. Modern idioms. Types are the primary abstraction.

* Correctness and clarity first. Comments explain "why" only.
* No `mod.rs` — use `src/some_module.rs`.
* No `unwrap()` — propagate with `?`.
* No indexing — use `.get()`. If indexing is unavoidable, `#[allow(clippy::indexing_slicing)]` + `debug_assert!`.
* No `let _ =` on fallible ops — propagate, log with `warn!`, or handle.
* No `from_utf8_lossy` — use `from_utf8()?` with typed errors.
* Error enums: one per module, typed fields only (never `String`), `#[from]` for wrapping. Never use `io::Error::other("message")` — add a typed variant instead. Hierarchy: `BgzfError` → `BamHeaderError`/`BaiError` → `BamError`; `FaiError`/`GziError` → `FastaError`; `FormatDetectionError` → `ReaderError`; `VcfHeaderError`/`VcfEncodeError`/`AllelesError` → `VcfError`.
* `color_eyre` for errors, `tracing` for logging.
* Sequence names are `SmolStr`.
* Tests: `cargo test`. Prefer `proptest` where applicable. Round-trip tests against noodles and bcftools are the strongest validation — always add them for new output formats.
* `SmallVec` in this project uses 2-arg form `SmallVec<T, N>` (not `SmallVec<[T; N]>`) due to debug-mode Vec alias in seqair-types. Use `vec![]` and `.to_vec()` instead of `SmallVec::from_buf`/`from_slice`.

## Architecture notes

**RecordStore**: 4 contiguous Vecs (records, names, bases, data). Zero per-record heap alloc.

**RegionBuf**: bulk-reads compressed bytes for a region, decompresses from memory. Uses `Vec<RangeMapping>` for disjoint chunks — never subtract offsets directly.

**ChunkCache**: `BamIndex::query_split()` separates nearby (L3–L5) from distant (L0–L2) BAI chunks. Distant chunks loaded once per tid per thread.

**CigarMapping**: `Linear` fast-path for clip+match (~90%), `Complex` with `SmallVec<6>`. Pre-extracted at construction.

**PileupAlignment**: base/qual/mapq/flags/strand pre-extracted. Hot loop reads flat fields only.

**PileupOp enum**: type-safe indel reporting — `Match`/`Insertion` carry `qpos`/`base`/`qual`, `Deletion` carries only `del_len`, `RefSkip` carries nothing. Compiler prevents reading a base from a deletion. Deletions and ref-skips are included in columns (not filtered out). `depth()` counts all alignments (matches htslib); `match_depth()` counts only those with a query base. Insertions attach to the last M/=/X position before the I op; `D I M` patterns produce orphaned insertions that are not reported. `del_len` is the total D op length at every position within the deletion (not remaining bases). `PileupOp` has a compile-time size guard (≤16 bytes).

**FASTA**: returns raw `Vec<u8>` (not `Vec<Base>`) — CRAM MD5 needs exact bytes. Conversion to `Base` at app boundary.

**Base::known_index()**: A/C/G/T → `Some(0..3)`, Unknown → `None`. Zero-depth pileups and Unknown ref bases are valid states.

**Forkable readers**: `Arc<BamShared>` (index + header) parsed once; `fork()` gives fresh File handle + ChunkCache.

**CRAM**: v3.0/v3.1. Multi-ref slices (ref_seq_id == -2), span=0 CRAI entries included in queries, embedded references, coordinate clamping to i64::MAX, rANS order-1 chunk-based interleaving, per-slice MD5 verification.

**Accumulator pattern**: `Default` struct, `add(&mut self, item)`, `finish(self) -> Result<T>`. Group by `Base::known_index()`, extract with `take`.

## I/O layers

1. `BgzfReader` (header only): BufReader 128KB → compressed 64KB → decompressed 64KB
2. `RegionBuf` (hot path): raw File → `data: Vec<u8>` (all compressed) → 64KB decompressed blocks
3. BAI: `fs::read()` into memory at open
4. `RecordStore`: ~900KB total for typical 30x region

## VCF/BCF writing

**Two encoding paths**: `BcfWriter::write_record(&VcfRecord)` (generic, allocates per record) and `BcfRecordEncoder` with typed handles (zero-alloc, pre-resolved dict indices). Both must produce BCF that noodles and bcftools parse identically.

**Type-safe Alleles**: `Reference`/`Snv`/`Insertion`/`Deletion`/`Complex` — enforces VCF structural invariants at construction. `write_ref_into`/`write_alts_into` for zero-alloc serialization. `begin_record()` on the encoder writes the BCF fixed header + alleles directly.

**Typed field handles**: `ScalarInfoHandle<T>`, `PerAlleleInfoHandle<T>`, `FlagInfoHandle`, `GtFormatHandle`, etc. Pre-resolved from header at setup. `handle.encode(&mut enc, value)` writes directly into BCF buffers. `BcfValue` trait: `scalar_type_code()` selects smallest int type for i32; arrays scan all values for uniform type.

**BCF format pitfalls** (caught by tests):
- BCF magic is `BCF\x02\x02` (v2.2), NOT `\x02\x01` (v2.1). bcftools rejects v2.1.
- BCF string dictionary order MUST match VCF header text emission order. Emit FILTER (PASS first) → INFO → FORMAT. If these don't match, noodles/htslib compute different dict indices and fields decode wrong.
- Integer arrays with missing values: scan only concrete values for `smallest_int_type`, then use per-type sentinel (int8=0x80, int16=0x8000, int32=0x80000000). Never use `i32::MIN` as a universal missing marker.
- `BgzfWriter::virtual_offset()`: cap `buf.len()` at `u16::MAX` to prevent truncation when buffer is exactly 65536 bytes.

**IndexBuilder**: single-pass TBI/CSI co-production during writing. Mirrors htslib's `hts_idx_push` state machine. Uses `BTreeMap` (not `HashMap`) for deterministic bin order.

## Profiling

`SEQAIR_PROFILE_JSON=/path/to.jsonl` → analyze with `python3 tools/analyze_profile.py`.
