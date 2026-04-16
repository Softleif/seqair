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
* Error enums: one per module, typed fields only (never `String`), `#[from]` for wrapping. Never use `io::Error::other("message")` — add a typed variant instead. Hierarchy: `BgzfError` → `BamHeaderError`/`BaiError` → `BamError`; `BamWriteError` (parallel to `BamError` for the write path); `FaiError`/`GziError` → `FastaError`; `FormatDetectionError` → `ReaderError`; `VcfHeaderError`/`VcfEncodeError`/`AllelesError` → `VcfError`.
* `color_eyre` for errors, `tracing` for logging.
* Sequence names are `SmolStr`.
* Tests: `cargo test`. Prefer `proptest` where applicable. Round-trip tests against noodles and bcftools are the strongest validation — always add them for new output formats. Avoid tautological tests — don't verify code using a copy of the same logic (e.g., compute expected bin with the same `reg2bin`). Use independent oracles. For tests that write temp files (e.g., bcftools round-trips), use the `tempfile` crate (`tempfile::NamedTempFile`) rather than `std::env::temp_dir().join(...)` to avoid filename collisions under parallel test execution.
* `SmallVec` in this project uses 2-arg form `SmallVec<T, N>` (not `SmallVec<[T; N]>`) due to debug-mode Vec alias in seqair-types. Use `vec![]` and `.to_vec()` instead of `SmallVec::from_buf`/`from_slice`.
* No silent `as i32`/`as u32` truncation at serialization boundaries — use `i32::try_from()` or equivalent with typed errors. At allocation boundaries (parsing untrusted counts), apply practical upper bounds matching real-world data, not format maximums (`i32::MAX` is not a limit). See `r[io.writer_limits]` and `r[io.fuzz.alloc_limits]` in the general spec.

## Architecture notes

**RecordStore**: 4 contiguous Vecs (records, names, bases, data). Zero per-record heap alloc.

**RegionBuf**: bulk-reads compressed bytes for a region, decompresses from memory. Uses `Vec<RangeMapping>` for disjoint chunks — never subtract offsets directly.

**ChunkCache**: `BamIndex::query_split()` separates nearby (L3–L5) from distant (L0–L2) BAI chunks. Distant chunks loaded once per tid per thread.

**CigarMapping**: `Linear` fast-path for clip+match (~90%), `Complex` with `SmallVec<6>`. Pre-extracted at construction.

**PileupAlignment**: base/qual/mapq/flags/strand pre-extracted. Hot loop reads flat fields only.

**PileupOp enum**: type-safe indel reporting — `Match`/`Insertion` carry `qpos`/`base`/`qual`, `Deletion` carries only `del_len`, `ComplexIndel` carries `del_len`+`insert_len` (deletion with following insertion, e.g. `D I M` in CIGAR), `RefSkip` carries nothing. Compiler prevents reading a base from a deletion. Deletions and ref-skips are included in columns (not filtered out). `depth()` counts all alignments (matches htslib); `match_depth()` counts only those with a query base. Insertions attach to the last M/=/X position before the I op; `D I M` patterns emit `ComplexIndel` at the last D position (matching htslib's `is_del=true, indel>0`). `del_len` is the total D op length at every position within the deletion (not remaining bases). `PileupOp` has a compile-time size guard (≤16 bytes).

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

**IndexBuilder**: single-pass TBI/CSI/BAI co-production during writing. Mirrors htslib's `hts_idx_push` state machine. Uses `BTreeMap` (not `HashMap`) for deterministic bin order. `push()` counts `n_mapped`; `push_unmapped()` counts `n_unmapped` for the pseudo-bin (placed-unmapped BAM reads). `write_bai()` requires `finish()` first (enforced by `finished` flag).

## BAM writing

**OwnedBamRecord**: mutable record with `Vec<CigarOp>`, `Vec<Base>`, `AuxData`. Separate from the read-path `BamRecord` (which uses `Box<[u8]>` for zero-copy decode). `from_raw_bam()` decodes all fields including mate info that the read-path drops (next_ref_id, next_pos, template_len). `to_bam_bytes()` appends into a caller-provided buffer; the caller clears it. `bin()` is BAI-only (u16, max 37449) and recomputed at serialize time.

**BamWriter**: wraps `BgzfWriter`, writes header eagerly at construction (no separate `write_header()`). Error poisoning: after any write failure, all subsequent writes return `Poisoned`. Index co-production validates records BEFORE writing to BGZF to avoid poisoning after partial writes.

**Index record dispatch** (three cases for BAI co-production):
1. Mapped (flags & 0x4 == 0, ref_id ≥ 0): `index.push(tid, pos, end_pos, voff)`
2. Placed unmapped (flags & 0x4 != 0, ref_id ≥ 0): `index.push_unmapped(tid, pos, pos+1, voff)`
3. Fully unmapped (ref_id == -1): not pushed to index
4. Mapped with ref_id == -1: rejected with `MappedWithoutReference` error

**AuxData**: raw BAM bytes with set/get/remove. `set_int()` validates range BEFORE writing tag bytes (orphaned bytes on error was a real bug). Unsigned-first type selection for non-negative values (C/S/I before c/s/i, matching htslib).

**Reader/writer limit parity** (`r[io.writer_limits]`): writers enforce the same field-size limits as readers. BAM header: l_text ≤ 256 MiB, n_ref ≤ 1M, l_name ≤ 256 KiB. BAI index: n_ref ≤ 100K, n_bin ≤ 100K, n_chunk ≤ 1M, n_intv ≤ 500K. Record size ≤ 2 MiB.

## Fuzzing

Fuzz targets live in `crates/seqair/fuzz/fuzz_targets/`. CI runs them nightly via `fuzz/run_all.sh`.

**Reproducing a crash locally** (requires Docker on macOS — cargo-fuzz needs Linux/ASAN):
```bash
# Copy the crash artifact to crates/seqair/fuzz/artifacts/<target>/
docker run --platform linux/amd64 --rm \
  -v "$PWD":/workspace -w /workspace/crates/seqair rust:1.93 \
  bash -c "rustup toolchain install nightly && cargo +nightly install cargo-fuzz && \
    cargo +nightly fuzz run <target> fuzz/artifacts/<target>/<crash-file> -- -runs=1 2>&1"
```

**Common crash patterns**: `debug_assert!` panics (cargo-fuzz enables `-Cdebug-assertions`). If a `debug_assert!` guards an `as i32`/`as u32` cast on untrusted input, replace it with `i32::try_from().trace_ok("msg")?` returning `None`/`Err`. Keep `debug_assert!` only for true internal invariants where the caller already validated the input.

## Profiling

`SEQAIR_PROFILE_JSON=/path/to.jsonl` → analyze with `python3 tools/analyze_profile.py`.
