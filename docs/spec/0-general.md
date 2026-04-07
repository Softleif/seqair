# Seqair — General Design

These are cross-cutting design rules that apply to all modules in the Seqair crate. They establish conventions for error handling, type safety, and platform portability.

> **Sources:** This file contains seqair-specific design rules with no direct upstream spec counterpart. Named constants (FLAG bits, CIGAR op codes) are derived from [SAM1] §1.4 and §4.2. See [references.md](references.md).

## Error handling

Seqair aims to provide helpful and complete error messages at all layers of the system.

r[io.errors]
Error types MUST implement `std::error::Error` using `thiserror` and capture necessary context and provide causes.

r[io.errors.typed_variants]
Every failure mode MUST be a distinct enum variant with typed fields — never a catch-all `String` or `reason: String`. Callers MUST be able to pattern-match on failure modes without inspecting display text. When a single enum variant is constructed at multiple call sites with different format strings, split it into separate variants. Paths, byte arrays, and counts MUST be stored in typed fields (e.g. `PathBuf`, `[u8; 4]`, `usize`), not interpolated into strings.

## Type safety

BAM files are full of magic numbers: flag bits (0x4 = unmapped, 0x10 = reverse strand), CIGAR operation codes (0 = M, 1 = I, ...), format constants. Using raw literals is error-prone and hard to read.

r[io.named_constants]
Magic numbers (bit flags, format constants, operation codes) MUST be defined as named constants, not used as raw literals. BAM flags MUST use named constants (e.g., `FLAG_UNMAPPED` instead of `0x4`). CIGAR operation types MUST use named constants (e.g., `CIGAR_M` instead of `0`).

r[io.typed_flags]
BAM flags MUST have a `BamFlags` newtype wrapping `u16` with named predicate methods (e.g., `is_unmapped()`, `is_reverse()`). The public API (`RecordRef`) MUST expose `BamFlags` rather than raw `u16` where appropriate.

r[io.typed_cigar_ops]
CIGAR operations MUST have a `CigarOpType` enum with variants for all 9 SAM-spec operations (M, I, D, N, S, H, P, =, X). The enum MUST provide `consumes_ref()` and `consumes_query()` methods. Invalid operation codes MUST be represented as `None` via `from_bam(u8) -> Option<CigarOpType>`.

## API design

r[io.non_exhaustive_enums]
All public enums (including error types) that may gain variants in future versions MUST be annotated with `#[non_exhaustive]`. Adding a variant to a non-exhaustive enum is not a semver-breaking change, allowing new variants and error conditions to be introduced in minor releases.

r[io.minimal_public_api]
Only intentionally public types should be exported. Internal implementation modules (compression codecs, container parsers, encoding details) MUST be `pub(crate)` or `#[doc(hidden)]` so that downstream crates do not depend on internal structure. The public API contract is defined by explicit `pub use` re-exports in each top-level module file.

## Arithmetic safety

r[io.checked_arithmetic]
All integer arithmetic on values derived from untrusted input (file data, parsed fields) MUST use checked or saturating operations. Default wrapping arithmetic (`+`, `-`, `*`) MUST NOT be used on untrusted values. The `clippy::arithmetic_side_effects` lint should be enabled.

r[io.writer_limits]
Writers MUST enforce the same field-size limits that the corresponding reader enforces. A writer that produces output the reader would reject is a bug — it creates files that appear valid but cannot be re-read. When a format stores a count or length as i32/u32, the writer MUST validate the value fits before casting, using `i32::try_from()` or equivalent with a typed error. Silent truncation via `as i32` is never acceptable at a serialization boundary. Where the format limit (e.g. i32::MAX ≈ 2.1 GB) is much larger than any practical value, writers SHOULD apply a tighter practical limit matching the reader (e.g. 256 MiB for header text, 1M references) to catch corruption early.

## Robustness — fuzzing

Seqair processes user-provided files that may be corrupt, truncated, or maliciously crafted. All parsers that read untrusted data SHOULD be fuzz-tested with `cargo-fuzz` (libFuzzer).

r[io.fuzz.readers]
All file format readers (BAM, SAM, CRAM, FASTA, BAI, FAI, GZI) and their sub-parsers (BGZF blocks, CIGAR ops, aux tags, rANS codecs, varints, Huffman tables) SHOULD have dedicated `cargo-fuzz` targets. Full-stack fuzz targets that exercise the complete read path (e.g., BGZF → header → records → pileup) SHOULD be maintained alongside unit-level targets.

r[io.fuzz.alloc_limits]
Parsers MUST NOT allocate memory proportional to an untrusted parsed size without a cap. All `Vec::with_capacity(n)` and `vec![0u8; n]` calls where `n` comes from parsed data MUST validate `n` against a reasonable upper bound and return an error if exceeded. This prevents OOM on corrupt files.

r[io.fuzz.no_panic]
Parsers MUST NOT panic on any input. Arithmetic overflow, out-of-bounds access, and shift overflow MUST be handled gracefully (via checked/saturating ops, `.get()` bounds checks, or early error returns). Fuzz targets SHOULD run with AddressSanitizer enabled to catch memory safety violations.

r[io.fuzz.seeds]
Fuzz targets SHOULD be seeded with structurally valid inputs derived from real test data files. Seeds dramatically improve coverage by letting the fuzzer start from valid file structures and mutate from there.

r[io.fuzz.simd]
Platform-specific SIMD code paths (NEON, SSSE3) SHOULD be fuzz-tested on their respective architectures to catch buffer overruns and alignment issues specific to vectorized implementations.

## Platform portability

r[io.platform_optimizations]
Every platform-specific optimization (SIMD, architecture-specific intrinsics) MUST have a scalar fallback that produces identical results on all platforms. The optimized and fallback paths MUST be tested with property-based tests that verify equivalence for arbitrary inputs.
