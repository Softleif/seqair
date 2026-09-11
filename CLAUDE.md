# seqair

Pure-Rust BAM/SAM/CRAM/FASTA reader + pileup engine. I/O backend for [rastair](https://github.com/bsblabludwig/rastair).

Workspace: `crates/seqair` (readers, pileup, BGZF, CRAM) and `crates/seqair-types` (Base, Strand, Phred, Probability, RmsAccumulator, RegionString, Pos/QPos).

## Tracey specs

Specs in `docs/spec/*.md` with `r[rule.id]`. Code: `// r[impl rule.id]`, tests: `// r[verify rule.id]`. Use the `tracey` skill. Add spec rules before implementing; update specs when changing behavior.

## Coding style

Expert Rust. Modern idioms. Types are the primary abstraction.

- Correctness and clarity first. Comments explain "why" only.
- No `mod.rs` — use `src/some_module.rs`.
- No `unwrap()` — propagate with `?`.
- No indexing — use `.get()`. If indexing is unavoidable, `#[allow(clippy::indexing_slicing)]` + `debug_assert!`.
- No `let _ =` on fallible ops — propagate, log with `warn!`, or handle.
- No `from_utf8_lossy` — use `from_utf8()?` with typed errors.
- Error enums: one per module, typed fields only (never `String`), `#[from]` for wrapping. Never use `io::Error::other("message")` — add a typed variant instead. Hierarchy: `BgzfError` → `BamHeaderError`/`BaiError` → `BamError`; `BamWriteError` (parallel to `BamError` for the write path); `FaiError`/`GziError` → `FastaError`; `FormatDetectionError` → `ReaderError`; `VcfHeaderError`/`VcfEncodeError`/`AllelesError` → `VcfError`.
- `color_eyre` for errors, `tracing` for logging.
- Sequence names are `SmolStr`.
- Tests: prefer `cargo nextest run` (parallel test-binary execution; ~8× faster than `cargo test` on this repo, config in `.config/nextest.toml`). Doctests aren't run by nextest — append `cargo test --doc --quiet` when you need to cover them. Prefer property-based tests with `hegel` (`hegeltest`) where applicable: `#[hegel::test] fn f(tc: TestCase)` with `tc.draw(gs::..)`, `#[hegel::composite]` for custom generators, plain `assert!`/`assert_eq!`, `tc.assume(..)` to reject a case. Case counts live in `hegel.toml`; pin an expensive test with `#[hegel::test(test_cases = crate::pinned::cases(N))]` (integration tests add `#[path = "../src/pinned.rs"] mod pinned;`) so `SEQAIR_PINNED_CASES=N` can still override it — a bare `test_cases = N` beats `HEGEL_TEST_CASES`. Prefer constructing inputs in the shape you need over rejecting them — a rejection loop makes hegel's smallest generated input large and its shrinking useless. Round-trip tests against noodles and bcftools are the strongest validation — always add them for new output formats. Avoid tautological tests — don't verify code using a copy of the same logic (e.g., compute expected bin with the same `reg2bin`). Use independent oracles. For tests that write temp files (e.g., bcftools round-trips), use the `tempfile` crate (`tempfile::NamedTempFile`) rather than `std::env::temp_dir().join(...)` to avoid filename collisions under parallel test execution.
- `SmallVec` in this project uses 2-arg form `SmallVec<T, N>` (not `SmallVec<[T; N]>`).
- No silent `as i32`/`as u32` truncation at serialization boundaries — use `i32::try_from()` or equivalent with typed errors. At allocation boundaries (parsing untrusted counts), apply practical upper bounds matching real-world data, not format maximums (`i32::MAX` is not a limit). See `r[io.writer_limits]` and `r[io.fuzz.alloc_limits]` in the general spec.

## Pre-push CI parity

CI runs clippy and tests with `-D warnings`. **Before pushing, always run with the same flag:**

Frequent CI-only failures to watch for: `clippy::cast_possible_truncation`, `clippy::arithmetic_side_effects`, `clippy::doc_markdown`, `clippy::field_reassign_with_default`, `clippy::needless_update`, `clippy::useless_conversion`, `clippy::unnecessary_cast`

## Architecture notes

**RecordStore**: 4 contiguous Vecs (records, names, bases, data). Zero per-record heap alloc.

**RegionBuf**: streams a region's compressed bytes through a bounded sliding window (`WINDOW_BUDGET`, 64 MiB), refilling forward on demand — peak RAM is the window, not the whole region. `RegionBuf::new` borrows the reader and plans `Vec<MergedRange>` (lazy, reads nothing); `ensure_available(need)` refills the window (compacting the consumed prefix, keeping `window_file_start + cursor` = true file offset), `advance_range` steps to the next disjoint range. `block_offset` (for virtual offsets) MUST be set as `window_file_start + cursor` *after* the header `ensure_available` and before the cursor advances — a stale read there silently corrupts every virtual offset. `with_budget(reader, chunks, budget)` (doc-hidden) forces tiny windows in tests. There is no batching/`MAX_REGION_BYTES` anymore — one streaming `RegionBuf` spans all of a query's chunks (see `BamQuery`).

**ChunkCache**: `BamIndex::query_split()` separates nearby (L3–L5) from distant (L0–L2) BAI chunks. Distant chunks loaded once per tid per thread.

**CigarMapping**: `Linear` fast-path for clip+match (~90%), `Complex` with `SmallVec<6>`. Pre-extracted at construction.

**PileupAlignment**: base/qual/mapq/flags/strand pre-extracted. Hot loop reads flat fields only.

**PileupOp enum**: type-safe indel reporting — `Match`/`Insertion` carry `qpos` (a `QPos`, see below)/`base`/`qual`, `Deletion` carries only `del_len`, `ComplexIndel` carries `del_len`+`insert_len` (deletion with following insertion, e.g. `D I M` in CIGAR), `RefSkip` carries nothing. Compiler prevents reading a base from a deletion. Deletions and ref-skips are included in columns (not filtered out). `depth()` counts all alignments (matches htslib); `match_depth()` counts only those with a query base. Insertions attach to the last M/=/X position before the I op; `D I M` patterns emit `ComplexIndel` at the last D position (matching htslib's `is_del=true, indel>0`). `del_len` is the total D op length at every position within the deletion (not remaining bases). `PileupOp` has a compile-time size guard (≤12 bytes).

**QPos vs Pos**: query offsets (`QPos`, 0-based into a read's SEQ) and genomic positions (`Pos0`/`Pos1`) are distinct, non-interconvertible newtypes — no `From`/`Into` between either and integer types; extract explicitly (`as_u32()`/`get()`). `AlignedPair`/`PileupOp`/`CigarPosInfo`/`BaseModState` query-space APIs all carry `QPos`. NB: `PileupOp::Insertion.qpos` is the matched base *before* the insertion; `AlignedPair::Insertion.qpos` is the *first inserted* base.

**FASTA**: returns raw `Vec<u8>` (not `Vec<Base>`) — CRAM MD5 needs exact bytes. Conversion to `Base` at app boundary.

**Base::known_index()**: A/C/G/T → `Some(0..3)`, Unknown → `None`. Zero-depth pileups and Unknown ref bases are valid states.

**Forkable readers**: `Arc<BamShared>` (index + header) parsed once; `fork()` gives fresh File handle + ChunkCache.

**CRAM**: v3.0/v3.1. Multi-ref slices (ref_seq_id == -2), span=0 CRAI entries included in queries, embedded references, coordinate clamping to i64::MAX, rANS order-1 chunk-based interleaving, per-slice MD5 verification.

**Accumulator pattern**: `Default` struct, `add(&mut self, item)`, `finish(self) -> Result<T>`. Group by `Base::known_index()`, extract with `take`.

## I/O layers

1. `BgzfReader` (header only): BufReader 128KB → compressed 64KB → decompressed 64KB
2. `RegionBuf` (hot path): raw File → bounded `window: Vec<u8>` (≤64 MiB compressed, refilled forward) → 64KB decompressed blocks
3. BAI: `fs::read()` into memory at open
4. `RecordStore`: ~900KB total for typical 30x region

## VCF/BCF writing

**Two output paths, one encoder**: `Writer` takes an `OutputFormat` (`Vcf` text, `VcfGz`, `Bcf`) and records go through the same typestate `RecordEncoder` chain (`Begun` → `Filtered` → `WithSamples` → `emit()`), but VCF text and BCF binary are serialized by entirely separate code below that. The two must render the same logical record to the same field values — that is `r[record_encoder.vcf_bcf_equivalence]`, and the float printers in particular share nothing. (There is no `BcfWriter` and no `VcfRecord`; both were removed.)

**Type-safe Alleles**: `Reference`/`Snv`/`Insertion`/`Deletion`/`Complex` — enforces VCF structural invariants at construction. `write_ref_into`/`write_alts_into` for zero-alloc serialization. `begin_record()` on the encoder writes the BCF fixed header + alleles directly.

**Typed field handles**: `InfoKey<V>`/`FormatKey<V>` with the value kind as the parameter (`Scalar<i32>`, `Arr<f32>`, `OptArr<i32>`, `Flag`, `Str`, `Gt`), and an alias per combination — `InfoInt`, `InfoFloats`, `InfoIntOpts`, `FormatGt`, `FormatString`, etc. Pre-resolved from the header at setup; the `FieldId` inside is `pub(crate)`, so the `InfoEncoder`/`FormatEncoder` trait methods are not reachable from outside the crate and `handle.encode(&mut enc, value)` is the whole public surface. `BcfValue` trait: `scalar_type_code()` selects smallest int type for i32; arrays scan all values for uniform type.

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

## compair (pair-HMM)

**Five kernels, one recurrence**: `align_full` (`f64`, whole matrix, the oracle), `align_banded`/`align_banded_simd` (anti-diagonal traversal, `banded.rs`) and `align_strips`/`align_strips_simd` (row-strip traversal, `strips.rs`). Each banded pair is one generic function over `Lane` (`f32` or `f32x8` via `wide`); bit-parity is by construction. No `mul_add` in `Lane` — `wide` falls back to separate mul+add on targets without FMA, which would break parity. `vmax` is `fast_max`: `wide`'s `max` is NaN-preserving and three instructions on AVX2.

**Strip kernel** (GKL's traversal): lane `l` holds cell `(r0 + l, d - l)` at step `d`; `left` is the previous step, `up` is the previous step shifted one lane with the row above entering at lane 0 (from a row buffer), `diag` is the previous step's `up`. The I and D diagonals are carried as one sum (`indel_diag`) since one transition multiplies both. Renormalises every `STRIP_ROWS` (8) rows *whatever the lane count* — the scalar instance's strip is one row and must scale at the same points as the eight-lane one, or the flush threshold moves and parity breaks. The row buffer is one buffer for the row above and the row below (reads at column `d` precede the store at `d - LANES + 1`); column 0 is zeroed after strip 0 because the scalar instance's stores start at column 1 and would leave row 0's `init` there. Rows past the read need no mask: their plan tracks are zero. Its precision does not depend on band width (pinned by `the_strip_kernel_tracks_the_recurrence_at_any_width`). Profile (M4 Pro): ~134 instructions per 8 cells, SIMD-throughput bound; the prior (~30), flush (12) and three lane shifts (~20) are the floor under `wide` — the 4× estimate from "fewer loads" was wrong because loads were never the binding constraint. Ideas not taken: a `StandardEmission`-specialised prior (-20 instructions/step, TAPS gains nothing), lane permutes via intrinsics.

**Emission is hoisted, not called per cell**: `Plan::fill<E: Emission>` folds `site_weights` (per column, reversed) and `epsilon` (per row) into `f32` tracks; the hot loop is generic over `L` only and never sees `E`. `prior()` reconstructs `SiteWeights::probability` with lanewise selects; `emission_table_matches_the_trait` pins the two to 4 ulp.

**`Read` precomputes every `powf`**: base error probabilities and the five `Transition` values per base are computed in `Read::new`. An alignment does no transcendental arithmetic.

**`Workspace`** owns the `Plan` tracks (shared by both traversals), the diagonal kernel's `Ring` (three matrices × three anti-diagonals, phase-major so one diagonal's matrices are a contiguous region borrowed via `split_at_mut`) and the strip kernel's `RowBuffer`. Reuse it across haplotypes and across traversals; the free functions build a fresh one per call.

**`Band::DEFAULT_WIDTH` is 46, not 48**: a diagonal holds `width / 2 + 1` cells; 24 cells fill three 8-lane vectors exactly, 25 spill into a fourth (~8% slower). Widths `16n - 2` fill every vector.

**Scaling is lazy**: each anti-diagonal is computed and stored in one scale, chosen as the previous diagonal's scale plus the shift *that* diagonal's maximum asked for; the two diagonals a cell reads are lifted onto the current scale by powers of two on load (`lift1`, `lift_before`), applied as two factors, never their product: a collapsing diagonal maximum can ask for a shift over 100 and `2^(s1+s2)` overflows f32 (fuzzer-found; pinned by `a_collapsing_diagonal_maximum_does_not_overflow_the_lift`). No in-place rescale pass, no branch on the shift. Subnormals are flushed to zero on store (`Sink::store`), which is what makes "a band under-estimates, never over-estimates" a guarantee. `normalising_shift_f32` reports no shift for subnormals and `exp2_f32(-127)` is 0 — unreachable in the kernel, pinned by proptests, not handled.

**Hot-loop shape** (what the profile taught): the chunk loop loads 20 windows of 8 `f32`; per-window `get(..).first_chunk()` bounds checks were two thirds of its instructions, so `Sources`/`Sink` prove every track's span once per diagonal (`View::new`) and load unchecked (`View::window`, the crate's only `unsafe`). `Ring::phases`, `Sources::new`, `Sink::new` are `#[inline(always)]` because LLVM otherwise kept them out of line and spilled the twenty views to the stack (12.6 → 9.9 µs). Per-diagonal cost is instruction count, not the shift→lift latency chain (measured by breaking the chain: no gain). Width sweep (`examples/profile_150x48.rs` style) separates per-diagonal from per-chunk cost.

**Tests**: GATK's 104 vectors (`tests/gatk_vectors.rs`), a hand-computed 3×2 matrix pinning conventions the vectors cannot (`tests/conventions.rs`), an independent masked `f64` DP oracle in `tests/banded.rs`, and `fuzz_pair_hmm` in `crates/seqair/fuzz`. Proptest case generators (`tests/support/mod.rs`) include `N`; `derived_case` reads score near zero (use for value checks at any length), `arbitrary_case` reads score far below `-30` (geometry checks only past ~40 bases).

## Fuzzing

Fuzz targets live in `fuzz/fuzz_targets/`. CI runs them nightly via `fuzz/run_all.sh`. No nightly toolchain: `RUSTC_BOOTSTRAP=1` lets the pinned stable accept cargo-fuzz's `-Z` flags (`run_all.sh` sets it).

**Reproducing a crash locally** (requires Docker on macOS — cargo-fuzz needs Linux/ASAN):

```bash
# Copy the crash artifact to fuzz/artifacts/<target>/
docker run --platform linux/amd64 --rm \
  -v "$PWD":/workspace -w /workspace rust:1.93 \
  bash -c "cargo install cargo-fuzz --locked && \
    RUSTC_BOOTSTRAP=1 cargo fuzz run <target> fuzz/artifacts/<target>/<crash-file> -- -runs=1 2>&1"
```

**Common crash patterns**: `debug_assert!` panics (cargo-fuzz enables `-Cdebug-assertions`). If a `debug_assert!` guards an `as i32`/`as u32` cast on untrusted input, replace it with `i32::try_from().trace_ok("msg")?` returning `None`/`Err`. Keep `debug_assert!` only for true internal invariants where the caller already validated the input.
