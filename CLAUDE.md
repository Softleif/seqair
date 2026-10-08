# seqair

Pure-Rust BAM/SAM/CRAM/FASTA reader + pileup engine. I/O backend for [rastair](https://github.com/bsbludwig/rastair).

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

**Entry points**: `Workspace::align_reads` (many reads × few haplotypes, packs pairs eight to a vector), `Workspace::candidates` + `Candidates::align` (one read at a time × many haplotypes), `Workspace::align_candidates` (same dispatch, single read). They pick kernels via `BATCH_BREAK_EVEN`/`PAIRS_BREAK_EVEN`; no choice changes a score. Crate docs in `lib.rs` list every kernel.

**Kernels, one recurrence**: `align_full` (`f64`, whole matrix, the oracle), `align_banded`/`align_banded_simd` (anti-diagonal, `banded.rs`), `align_strips`/`align_strips_simd` (row strips, GKL's traversal, `strips.rs`), `align_batch` (eight haplotypes × one read) and `align_pairs` (eight unrelated pairs), which share one row sweep in `lanes.rs`. Every kernel is generic over `Lane` (`f32` or `fearless_simd`'s `f32x8<S>`, `simd.rs`), so scalar/SIMD parity, and batch/pairs-lane parity with the strip kernel, is bit-exact by construction. Keep it that way: no `mul_add` on `Lane` (FMA measured slower on Zen 2, and it would break parity), same ops in the same order everywhere. The `gpu` feature (wgpu) is the exception: its oracle is a tolerance.

**SIMD dispatch**: `dispatch!` compiles each kernel per level (NEON; SSE2/SSE4.2/AVX2/AVX-512) and picks one at run time; no build flags needed. **Kernels and every helper they call must be `#[inline(always)]`, and no vector op may sit in a closure** — an out-of-line helper or closure is compiled for the baseline, and on x86-64 each vector op becomes a call. aarch64 can't show this; check x86-64. `shift_in` is the one AVX2 escape hatch (`asm!`). `masked`/`masked_out` are an `and` of `SimdMask::to_vector` on AVX2/SSE4.2, where a `select` against zero is a sign-bit blend LLVM pads with a sign extension; `select` elsewhere.

**The loops are latency-bound, not port-bound** (`perf` on a 3950X): latency changes pay where instruction counts don't. **A/B by building two binaries and interleaving them under `perf stat`**, never by criterion's p-value across two runs — the box drifts several percent.

**Strip kernel**: lane `l` holds cell `(r0 + l, d - l)` at step `d`; the I and D diagonals are carried as one sum. Renormalises every `STRIP_ROWS` (8) rows *whatever the lane count* — the scalar instance must scale at the same points or parity breaks. The row buffer serves the row above and below (reads at column `d` precede the store at `d - LANES + 1`).

**Emission is hoisted**: `Plan::fill<E: Emission>` folds `site_weights` (per column) and `epsilon` (per row) into `f32` tracks; the hot loop never sees `E`. `prior()` reconstructs `SiteWeights::probability` lanewise; `emission_table_matches_the_trait` pins the two. `Read::new` precomputes error probabilities (256-entry table) and transitions — an alignment does no transcendental arithmetic.

**`Workspace`** owns every kernel's buffers; reuse it across haplotypes and reads (allocation-free). The free functions build a fresh one per call.

**`Band::DEFAULT_WIDTH` is 46, not 48**: a diagonal holds `width / 2 + 1` cells; 24 fill three 8-lane vectors, 25 spill into a fourth. Widths `16n - 2` fill every vector.

**Precision: every entry point returns the `f64` recurrence over its band, down to ≈ −600 (150 bp).** `align_banded_f64` keeps each row's largest cell at `2^F64_ROW_SCALE` (960, not `[1, 2)`: at `[1, 2)` a path 308 decades below its row's best rounded as a subnormal, and a pinned Q247/Q253 pair scored −430 instead of −352); below the floor it can be off by its bound either way (round-to-nearest, not a flush). The test oracle `align_masked` computes in `Wide` (mantissa + unbounded exponent per value) so it has no subnormals. The `f32` kernels flush stored cells below `2^-126`. The strip, batch, pairs and GPU kernels renormalise every eighth row so its largest cell sits in `[2^115, 2^116)` (`scaling::STRIP_SCALE`, `strip_shift_f32`, shift capped at 127). A cell is *not* a probability: GATK's recurrence takes a match's transitions from two rows (next row's `mm`/`mi`, own row's `md`) and a deletion's too (own `gc`, next row's `itm`), so where a read's gap qualities change a cell hands on more than it holds (alternating gap-continuation Q254/Q0 on a homopolymer scores +37.7). `transitions::Growth` bounds that per read: `gamma(i)` per row crossing (a sub-stochastic chain once divided out; deletion runs folded into match→match steps), `log10_total` over the read, `window` over one 8-row strip window, `spill` (deletion over its row's match max), `carry` (what a deletion hands on). `Read::growth` is one pass over the rows, so `Read::growth_bound` (cached at `Read::new`, three byte scans, only for a steady gap-continuation quality) is tried first — it only vouches where the exact factors would. `reference::trusted` rescores in `align_banded_f64` unless `strips_fit` (no value can overflow: cells ≤ `3·(1+spill)·window·2^(S+1)`, the total ≤ `(2+min(1,spill))·columns·window·2^(S+1)`; at 115 the widest band with steady qualities fits, at 116 it would not) and the score is above `trust_floor` = `6 + log10(3·r·columns·(2+carry)·total²·max(max(1,spill)·rho·2^(-126-S), 2^-253))` + rounding, with `rho` = (row-0 band columns feeding row 1)/h and the `3` covering a GPU that flushes subnormal intermediates (≈ −62 for 150 bp against 290 in a 64-wide band). The diagonal kernel keeps each anti-diagonal at the same scale, shifted by the larger of the two diagonals the next one reads (so a collapsing diagonal cannot lift the one two back past `f32`), and uses the same `trusted`; the floor takes `total` squared because a diagonal's scale comes from a later row than the cells it drops (free on steady reads). Every CPU entry point calls `trusted` on its kernel's raw score, exactly once (the strip-kernel fallback paths inside `align_candidates`/`route` are already rescued — don't wrap them again). `trusted` is crate-private: a GPU score comes back as `gpu::GpuScore` carrying its pair's floor (`reference::pair_floor`, taken at `GpuPairs::push_pair`, `min` of the bound's and the exact factors' floors so it keeps exactly what `trusted` keeps — `the_pair_floor_keeps_what_trusted_keeps`), and `GpuScore::trusted`/`or_rescore` apply it. The checks depend on the pair alone, so bit-parity between kernels survives. The proof needs every emission and every state's out-transitions ≤ 1: `emission::epsilon` caps `eps` at `MAX_EPSILON` = 3/4 (use it, never `Emission::epsilon` directly), and `Transition::from_qualities` raises gap-open qualities below `MIN_GAP_OPEN_QUALITY` (6) to 6 (GATK's `MIN_USABLE_Q_SCORE`). On whole chr12 0.97% of rastair's pairs are rescued (the replayed pair set is rescue-heavy: 2.2% there, by either kernel (strips 2.76% at `2^96`; diagonal 8.46% at `[1, 2)`)); a per-row flush flag was measured at +21% strip-kernel cycles and dropped. The value tests (`every_score_is_the_f64_recurrence_over_its_band`, the fuzz target) are unconditional — don't reintroduce `> -30` carve-outs; at default case counts they don't notice a floor ~12 decades too low, which `a_strip_flush_below_the_strip_floor_is_rescored` pins; `no_total_exceeds_the_free_start_times_the_growth` holds `Growth` against the `f64` recurrence and `the_widest_band_keeps_its_total_below_f32_max` pins the overflow bound (its raw total is `2^126.6`).

**The `f64` rescue (`rescue.rs`)** is the row recurrence two rows at a time: a vector prepass for row `a`'s match/insertion cells, then one sweep runs row `a`'s deletion chain and all of row `b` a vector behind, so the two chains (a multiply + add per cell, 7 cycles on M4) overlap. Row `b` reads row `a` at the scale of row `a`'s largest match/insertion cell; the chain's maximum is checked afterwards and row `b` redone if it differs (<1% of real rows). Bit-identical to `align_banded_f64_rows` (the old implementation, kept as the oracle; `the_rescue_is_the_row_recurrence_bit_for_bit`). Buffers are thread-local and never cleared: writers store a zero vector past each row and column 0 is zeroed; the test poisons them with NaN. Tried and dropped (M4): skipping the scale multiply on unscaled rows (+6%, the branch on a late-resolving scale mispredicts), making the next pair's prepass speculatively inside the sweep (+33%, the extra work fills the OoO window the two chains need).

**Scaling is lazy**: each anti-diagonal is stored in one scale (previous scale + the shift its maximum asked for); the two diagonals a cell reads are lifted by powers of two on load, as two factors, never their product (`2^(s1+s2)` overflows f32 — fuzzer-found, pinned by `a_collapsing_diagonal_maximum_does_not_overflow_the_lift`). Subnormals flush to zero on store, which is what makes "a band under-estimates, never over-estimates" a guarantee.

**Tests**: GATK's 104 vectors (`tests/gatk_vectors.rs`), a hand-computed 3×2 matrix (`tests/conventions.rs`), an independent masked `f64` DP oracle (`tests/banded.rs`), rust-bio as a second implementation (`tests/rust_bio.rs`), the 10s dataset against gkl on x86-64 (`tests/tenspeed.rs`), and `fuzz_pair_hmm` in `fuzz/`. Property tests are hegel, like the rest of the workspace (`crates/compair/src/pinned.rs` for pinned counts); generators in `tests/support/mod.rs` include `N` and draw lengths as integers, since hegel keeps a vector's own size short. `derived_case` reads score near zero (the kernels' own numbers, no rescue), `arbitrary_case` draws every quality Q0–Q254 and scores far below the trust floor (exercises the rescue).

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
