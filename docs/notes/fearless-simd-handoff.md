# fearless-simd: where we stopped (2026-09-24)

Handoff note for picking this work up again. **Drop this file before merging.**

## TL;DR

- Branch `fearless-simd` (worktree `../seqair-fearless-simd`, head `288df641`, 28 commits on `main`) replaces every hand-written `core::arch` kernel with one `fearless_simd` kernel. It also merges three rounds of profiling-driven work: CRAM, text/writers and BAM/pileup. Nothing is pushed.
- Tests: 1637 pass on the M4 and 1690 on the 3950X (samtools/bcftools round trips included). Every-level tests pass at 5000 cases. Clippy `-D warnings` and tracey are clean, and there are 0 out-of-line fearless calls in the x86 asm.
- **Blocker before a PR:** a real rastair run (chr12, `--gpu`, experimental-seqair) is **~10% more cycles** with compair + `fearless-simd@53caf6e0` than with compair + `main`:
  - 374.4 G vs 339.5 G cycles, with instructions flat (+0.8%) and byte-identical output.
  - The extra time lands in rastair's `ColumnDraft::accumulate` and seqair's per-alignment accessors, none of which changed.
  - Suspects: a code layout or inlining side effect.
  - Not bisected yet, and `288df641` hasn't been measured on rastair yet.
- The final consolidated 3950X criterion run (`bench-base` vs `fearless-simd`, aligned builds) was stopped after round 1 of the base side. Rerun it (see below) to get the PR table.

## Numbers we have

The 3950X figures are from interleaved criterion runs, taking the min of medians. They were measured **before** `-align-loops=64` was adopted, so treat differences under ~3% as noise. The M4 figures are only from runs with low spread.

| What | Before → after | Machine / workload |
|---|---|---|
| rANS Nx16 32-state order-0 kernel | 4.45 → 1.29 ms / MiB (0.29×) | 3950X, `simd_kernels` |
| ASCII→Base (150 B / 10 kB / 1 MiB) | 0.78× / 0.70× / 0.74× | 3950X, `simd_kernels` |
| `decode_bases_into` (150 B / 10 kB / 1 MiB) | 0.76× / 0.71× / 0.69× | 3950X |
| `decode_seq` (150 B / 1 MiB) | 0.98× / 0.79× | 3950X |
| FASTA fetch (plain, 1 kb → 1 Mb) | 0.85× → 0.80× | 3950X, `fasta_fetch_steady` |
| CRAM 3.1 end to end (chr20:10–30M, 813k reads) | 2.13 → 1.15 s (1.84×) · M4 1.078 → 0.541 s (2.0×) | c20_31.cram |
| CRAM 3.0 end to end | 1.86 → 1.34 s (1.39×) · M4 1.073 → 0.688 s (1.56×) | c20_30.cram |
| rANS Nx16 32× order-1 codec | 132 → 563 MB/s · M4 279 → 1027 | 1 MB QS stream |
| VCF text writing | 3.8× fewer cycles (M4 3.4×) | 3950X perf stat, `prof_textwrite vcf` |
| VCF.gz writing | 2.1× fewer cycles | 3950X |
| SAM region fetch | 2.7× fewer cycles (M4 2.9×) | 3950X, bgzipped SAM of test.bam |
| BAM serialize (encode) only | 1.9× fewer cycles; full BAM write −2% | 3950X |
| MM/ML base-mod resolution | 5.7–7.6× (sparse), 1.8–2.4× (dense) | 3950X, `simd_kernels/base_mod_parse` |
| `pileup_with_reference` / `pileup_tiled` | −13.5% / −6.5…−12.6% | 3950X, `bam` bench vs 75ca1a5b |
| rastair-shaped TAPS pileup (chr20:30–36M) | −16.5% cycles | M4, `examples/pileup_count` |
| **rastair chr12 end to end** | **+10% cycles (regression, see above)** | M4, fearless-simd@53caf6e0 |

## Next steps, in order

1. **Bisect the rastair regression.** Scratch setup, made by the profiling agent:
   - `scratch/compair-fearless` is compair + fearless.
   - `scratch/compair-main` is compair + main.
   - The rastair copies point at those. They were made from the `rastair2-secondary` worktree (c4a8bc6a), which predates rastair's own port. **rastair `main` already uses seqair 0.3**, so redo the scratch rastair from rastair `main` rather than the agent's hand port. The agent's port tripped on `fetch_base_seq`, whose end is now **inclusive**.
   - Workload, run from `rastair2-secondary`: `<bin> call -r tmp/taps/hg38.fa.gz -l chr12 -@ 8 --gpu --experimental-indels=ml -o X.bcf --bed X.bed.gz tmp/taps/NA12878_aa_chr12.bam`.
   - Compare cycles (`/usr/bin/time -l`), interleaved.
   - Measure `288df641` first; its pileup offset caching targets the hot accessors.
   - Then bisect `d4ca9536..53caf6e0`, starting with `0080e40a` (seqair-types) and `716ed125` (bam decode).
   - Also try building both sides with `-C llvm-args=-align-loops=64` / `-align-all-functions=6`.
2. **Rerun the final PR table on the 3950X.** Scripts are saved on the box in `~/seqair-bench-logs/box-final.sh` and `ab_table2.py`. The `bench-base` branch is `main` + the kernel benches (`eefd6910`) + the base-mod bench, so both sides run identical benches. Run it from the Mac. It bundles both branches and builds with `-align-loops=64`.
3. **Merge the leftovers:**
   - `explore/cram-2`, from the CRAM round 2 agent. Check which commits are measured wins and which are `wip:`.
   - compair: the branch **`compair-fearless`** (off `compair` 21ff717c, head `53cb5807`) is a clean series that makes fearless_simd the only SIMD lane. It drops `wide` and the hand intrinsics lane, and moves the diagonal kernel over too. `experiment/fearless-compair` is the exploratory history. This goes into `compair`, not `fearless-simd`, and needs a separate PR.
     - 3950X, aligned builds, cycles vs `compair` head: strips +0.8%, batch +0.1%, candidates +0.5%, diagonal −32%.
     - Criterion 10s: strips 17.22→17.41 ms, candidates 8.84→8.82 ms, diagonal 58.5→38.7 ms.
     - The M4 timings were inconclusive because of noise.
     - Costs: x86 `.text` +16%, and a compair lib release rebuild takes 0.69→1.85 s.
     - API: removed `align_strips_intrinsics`, `intrinsics_lane_available` and the `intrinsics` feature; added `simd_level()`. rastair only calls `align_strips_simd`, which is unchanged, so rastair needs no change.
     - Not done: the docs (README, lib docs, the CLAUDE.md compair section, and a benchmarking.md section on the Vec-header reloads, the codegen-units=1 caveat and `-align-loops=64`).
     - The dev-dependency `force_support_fallback` (needed for Fallback-level tests) makes test/bench/example rebuilds 1.2→5.0 s on the M4. Consider gating it behind a test-only feature.
4. **Follow-ups found but not done:**
   - Pileup `PileupAlignment` 44 → ~20 B (est. 5–10%) and a `ActiveRecord` CIGAR side arena (est. 3–5%).
   - CSI `query_split` still scans every bin.
   - `compute_end_pos` runs twice per record.
   - A base-mod decimal parser.
   - rastair-side accessor costs (19–23% of rastair CPU):
     - `extra_of` / `record` random loads: hoist one lookup per alignment.
     - A mate binary search per column: needs an O(1) per-column map.
   - CRAM: PACK unpack SIMD, RLE, eager `CramError` construction (6–10%), tok3 allocations, 4-state Nx16 x86 codegen, 4x8 o1 split table, MD5 caching.
   - SAM: `parse_aux_tags` (11% of SAM fetch).
   - VCF: a `SmolStr` per FORMAT field, and a float digit table.

## Rules learned (also in spec `io.simd_portable` and the Claude memory)

- Only the dispatched entry point is `#[simd]`; helpers are `#[inline(always)] fn f<S: Simd>`. Never put vector ops in a closure: on x86 each op becomes a `call …__fearless_simd_kernel`, 10–20× slower, and NEON hides it. Guard: grep the x86 asm for that symbol.
- `cargo asm` builds with codegen-units=1. Count instructions from `objdump` of the real build (compair's `tools/fearless_asm/objhot.py`).
- On Zen 2, build both A/B sides with `-C llvm-args=-align-loops=64`: a 16-byte function shift moved a kernel +9% with identical instructions.
- Overlapping-vector tails beat scalar tails. In-place kernels must load the tail before storing, or they stall on store forwarding.
- `swizzle_dyn` on the `Fallback` level wraps out-of-range indices instead of zeroing them.

## Box setup (3950X)

Connect with `ssh yulcal@ph-steamos.taile64a30.ts.net`, then:

- `export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"` (samtools/bcftools/bgzip)
- `export BINDGEN_EXTRA_CLANG_ARGS="-I$(ls -d /usr/lib/clang/*/include | tail -1)"`
- Use `~/.cargo/bin/cargo +1.98.1`.
- Wrap every timing in `flock ~/bench.lock`.
- Ship branches as a `git bundle`, since worktree `.git` files are pointers.
