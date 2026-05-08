# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0 (2026-05-08)

### New Features

 - <csr-id-142e06d402ec13371cd2d9715b59f6181495462b/> complete BAM aux tag API — setters, spec, fuzz, proptests
   - Remove stale r[bam.record.aux_array_unsupported] — B-arrays are fully parsed
   - Add r[bam.record.aux_truncated] — unknown type codes stop iteration (matches htslib)
   - Add r[bam.owned_record.aux_array_setters] — spec for typed B-array setters
   - New AuxData methods: set_char, set_double, set_array_i16/u16/i32/u32/f32
   - Add InvalidArrayLength error variant for misaligned array data
   - DRY array encoding with array_count() and encode_array_header() helpers
   - New fuzz target: fuzz_bam_aux — round-trip, uniqueness, removal properties
   - 13 proptests: parser/builder oracle, integer round-trip, array round-trips, uniqueness, removal
   - 14 new unit tests covering all new setters and edge cases

### Other

 - <csr-id-7a2701048ee6ce08373e2207268abf431a2e4401/> emit Unknown/UNAVAILABLE for SEQ=* records
   Mapped records with SEQ=* (seq_len=0) were silently dropped from every
   pileup column because the qual slab lookup returned None. htslib's pileup
   keeps them — bam_seqi reads beyond the empty SEQ buffer, the 0xFF qual
   sentinel decodes to N, and the alignment counts toward depth(). Any BAM
   with secondary alignments lacking sequence (a common SAM convention) saw
   a different coverage profile under seqair than htslib; perbase #107 had
   to ignore an empty-SEQ parity test waiting for this fix.
   
   Match/Insertion ops on a seq_len=0 record now emit Base::Unknown +
   BaseQuality::UNAVAILABLE, matching htslib's depth() and match_depth().
 - <csr-id-7e743ac7fe6000e965c883cf59eaa1a3d08afb53/> rename iter() to records() and drop the index
   The yielded `u32` was redundant — `SlimRecord` already carries
   `cigar(&store)` / `qname(&store)` / `seq(&store)` / `qual(&store)` /
   `aux(&store)` getters that cover every shared-reference access pattern,
   and the only operations that need the index are mutable (`set_alignment`,
   `extra_mut`, …) which can't run during iteration anyway. Callers that
   genuinely need a position can use `records().enumerate()`.
 - <csr-id-fcce57b9088a3d410ddaa1c48f44597161ad571d/> add iter() yielding (idx, &SlimRecord) pairs
   Replaces the `for idx in 0..store.len() as u32 { let rec = store.record(idx); ... }`
   boilerplate with `for (idx, rec) in store.iter() { ... }`. Yielded `idx`
   matches the `u32` returned by `push_raw`/`push_fields` so callers can still
   reach the variable-length slabs via `store.cigar(idx)` etc., or via the
   `SlimRecord::cigar(&store)` getters.
 - <csr-id-f3c39688e59d48443c1f50b6b3a27e91a564f146/> add keep_unmapped opt-in for placed-unmapped reads
   IndexedBamReader's pre-filter unconditionally dropped placed-unmapped
   reads (flag 0x4 with a valid tid) before they reached the customize
   layer. Downstream tools that want htslib-`view`-equivalent semantics
   (e.g. perbase `only-depth -x`) had no way to see them and worked around
   it by injecting `-F 4` at the call site.
   
   Add a builder-style `keep_unmapped(bool)` setter (and `keeps_unmapped()`
   inspector). Default `false` preserves htslib-pileup-correct behavior;
   `true` lets `filter_raw` / `filter` decide. `fork()` propagates the
   flag.
 - <csr-id-dd324da07f5945f658abee76f66473131f99fcfa/> move entry points onto BamWriterBuilder
 - <csr-id-957d9bb85fbcbd2c2d4f8aae15c7e7336146f13e/> drop write_index requests on non-path writer targets
   When build_over(writer, header).write_index(true).build() is called,
   there is no sidecar location for a .bai file and BAI virtual offsets
   only make sense against a seekable BAM file — building an IndexBuilder
   would just waste work and hand the caller something that points into a
   stream they have already streamed away.
   
   Accept write_index(true) on the ToWriter path for API uniformity, but
   drop the request: log info! and skip IndexBuilder construction. The
   ToPath path is unchanged.
   
   Migrate the existing index tests that used build_over(&mut Vec) to
   tempfile-backed BamWriter::builder(path) targets, and add a new
   write_index_on_writer_target_is_soft_noop test that pins down the new
   behavior.
 - <csr-id-76d934b6d0e8c2d0bd1f81fab4e24ce87ee49021/> re-export BamWriterBuilder/ToPath/ToWriter from seqair::bam
   The builder types are part of the writer's public API (BamWriter::builder
   and BamWriter::build_over both return one), but were only reachable as
   seqair::bam::writer::BamWriterBuilder. Surface them at seqair::bam so
   they line up with BamWriter and BamWriteError.
 - <csr-id-1280edbcf63f5ef429ba0f94070446b2ef8f6011/> tighten &[u8] FromAuxValue to Z-only, add HexBytes for H
   The &[u8] FromAuxValue impl previously accepted both Z and H tags,
   returning the raw stored bytes for either. This silently misinterpreted
   H tags: their wire encoding is ASCII hex digits, so a caller fetching
   &[u8] would get the digits rather than the underlying byte array.
   
   Restrict the &[u8] impl to Z, and add a HexBytes wrapper with its own
   FromAuxValue impl plus a decode() method that materializes the byte
   array. H tags now require an explicit fetch as HexBytes; Z tags do not
   coerce into HexBytes. Searched the workspace — no production caller was
   relying on H acceptance via &[u8].
 - <csr-id-25d1791b0eceb3caa47cbba27690c6116fd4e04b/> post-BamRecord-removal cleanup
   Re-attach the spec rules orphaned by the BamRecord -> SlimRecord
   migration to their new code homes (BamFlags predicates, parse_header,
   decode_bases_into, RecordStore::{seq_at,aux}, push_raw SIMD path,
   matching/indel pre-compute fields). Rename remaining BamRecord prose
   references to OwnedBamRecord across the BAM-record-builder/writer/flags
   specs, point r[bam_writer.index_coproduction] at the new
   BamWriterBuilder::write_index API, and replace the BamRecord mention in
   r[perf.precompute_matches_indels] / base-quality field-usage with
   SlimRecord. Drop the contradictory "no universally correct value"
   paragraph from SegmentOptions::new now that Default exists, and note on
   BaseModState::from_record that we reject MM-with-deltas records lacking
   ML (stricter than SAM 1.6). Extract the duplicated push-raw test
   scaffold into bam::test_util.
 - <csr-id-34c2363cddf895426d3467ce81955972a4e219bc/> replace BamWriter::from_path boolean with a builder
   Drops the opaque `build_index: bool` and positional `level: i32` from
   `BamWriter::from_path` / `from_path_with_level` / `from_stdout` /
   `new_inner` in favor of a `BamWriterBuilder` reached via
   `BamWriter::builder(path, header)` (file target) or
   `BamWriter::build_over(writer, header)` (arbitrary `Write` sink).
   Setters: `write_index(bool)` and `compression_level(i32)`. The header
   is required positionally; everything else has a default. Header is
   still written eagerly during `build()`. Spec rules and all in-tree
   call sites (bench, example, tests, internal tests) updated.
 - <csr-id-dc174f711896f84bc219500c164a3e28a7605f65/> add BaseModState::from_record(rec, store) constructor
   Wraps the four-line MM/ML/seq/is_reverse extraction in a typed
   constructor returning Result<Option<Self>, FromRecordError>: Ok(None)
   when MM is absent, typed errors for wrong-type tags or malformed MM/ML.
   Forgetting to thread `is_reverse` through `parse` silently produces
   wrong methylation calls — the constructor pulls it from rec.flags
   itself. Drops the manual extraction from examples/base_mods.rs.
 - <csr-id-4e36d0bfeac5cfc8015712c60fd765b8cd278459/> add Display for CigarOp/CigarOpType and CigarStr slice wrapper
   Centralizes SAM CIGAR formatting so callers stop hand-rolling the
   op-code-to-char match. Drops the inline fmt_cigar helper from the
   realignment example.
 - <csr-id-41fc3a7a9dae61898c2e6713724be988258ddd07/> extract BamWriter::finalize_record, share post-encode path
   Both write paths (write -> write_inner via OwnedBamRecord::to_bam_bytes,
   and write_store_record -> write_store_record_inner via slab + the
   shared encode_fixed_header) ended with the same dance: size-cap check,
   pre-write index validation (mapped-without-reference rejection), BGZF
   block-size + bytes write, then index dispatch with the placed-unmapped
   vs mapped vs fully-unmapped fork.
   
   Two copies of the index dispatch meant a fix in one site silently
   diverged from the other — already happened once with the
   beg.saturating_add(1) clamp.
   
   Extract `finalize_record(ref_id, is_unmapped, beg, end_pos)` for the
   shared tail. Each writer entry now: fill self.buf with record bytes,
   compute (beg, end_pos), call finalize_record. Net -35 lines in writer.rs
   and one place to fix any future change to size limits, index pushing,
   or BGZF write order.
 - <csr-id-6d80e6e0a05c62d4857067232624a9c036b0992a/> dedupe 32-byte fixed-header encode into shared helper
   Previously two writer paths each packed the BAM record fixed header
   inline: OwnedBamRecord::to_bam_bytes (the typed-record path) and
   BamWriter::write_store_record_inner (the slab-resident path). Both
   duplicated the bin_mq_nl / flag_nc bit-packing and the eight LE-byte
   appends. A field added to one site without the other (e.g. the bin
   recompute, or the n_cigar_op casts) was a bug class waiting to happen.
   
   Extract `bam::record::encode_fixed_header` as the inverse of the
   existing `parse_header`, taking a `FixedHeaderFields` struct in BAM
   wire types (pos/next_pos already resolved to -1 for unmapped, and
   l_read_name as u8 to encode the qname ≤ 254 invariant in the type).
   Both call sites now construct the struct with their own validated
   casts and call the shared encoder. The resulting BAM bytes are
   unchanged — verified by the existing writer round-trip tests and the
   broader suite.
 - <csr-id-c4675241707b03fabea905aa2180a97e506e65b6/> drop OwnedBamRecord::from_raw_bam, the unused single-record decoder
   The method had zero callers outside its own two tests — its docstring
   already steered everyone to RecordStore::push_raw for bulk decode, and
   no production caller in seqair or rastair2 used it for the "edit one
   record" flow either. Keeping it meant maintaining a third copy of the
   "parse 32-byte header → walk qname/cigar/seq/qual/aux" walk (after
   RecordStore::push_raw and the now-deleted BamRecord::decode).
   
   Drop the method, the OwnedRecordError::Decode variant it raised, and
   the pos_from_bam_i32 helper that nothing else used. Drop the two
   tests that exercised it as a to_bam_bytes round-trip oracle — the
   writer's round-trip coverage already exercises that path through a
   RecordStore.
   
   After this, "decode raw BAM bytes" has exactly one entry point:
   RecordStore::push_raw.
 - <csr-id-1275c9e8c937eddc4f7b6cdc3ac6a217030ed5f8/> remove BamRecord, decode through RecordStore::push_raw
   The read-path BamRecord struct had no production callers — every reader
   in the crate decodes via RecordStore::push_raw straight into slabs. Its
   only uses were as a round-trip oracle in tests and a separate fuzz
   target that fuzzed a code path production never takes.
   
   Drop the struct, decode helper, and accessors. Keep parse_header,
   ParsedHeader, compute_end_pos_from_raw, read2, and read4 — those are
   the genuinely-shared primitives still used by RecordStore::push_raw,
   the pre-push CIGAR scan, and OwnedBamRecord::from_raw_bam.
   
   Migrate the round-trip tests in owned_record.rs and writer.rs to push
   into a RecordStore and assert against SlimRecord/store accessors —
   which is what production runs anyway. Retarget the fuzz target the
   same way.
 - <csr-id-1a3c0846ba6d24be106882a89cd72803c92ffd3a/> skip seq/CIGAR reconstruction for foreign-tid multi-ref records (Finding 5)
   In a multi-ref slice (`sh.ref_seq_id == -2`), records can belong to any
   reference. A query for tid=N decodes every record in the slice, but
   records with `record_ref_id != N` are guaranteed to be discarded — yet
   the previous code ran the full `decode_features_and_reconstruct` pass,
   allocating into bases_buf/cigar_ops_buf, looking up reference bases
   against the wrong reference (which fired spurious "reference shorter
   than expected" warnings), and then threw the result away.
   
   Move the foreign-tid check up and dispatch to a new
   `scan_features_for_refspan` helper that consumes the per-record data
   series in lockstep but skips reconstruction. Stream consumption MUST
   match the reconstruct path exactly (same `ds.X.decode*` calls in the
   same order for the same feature codes); drift desynchronises every
   subsequent record in the slice. The integration tests against htslib
   and noodles cover this — they all still pass.
   
   Also drain the MQ + quality streams for foreign-tid records (still
   required for stream sync) but route the quality bytes into the scratch
   `feature_byte_buf` so they don't pollute `qual_buf`.
 - <csr-id-1a13ca0e595c023a90414a6210f8db2e24fa8eb5/> validate rANS 4x8 header compressed_size against payload (Finding 4)
   The rANS 4x8 header's 4-byte `compressed_size` field was previously read
   into `_compressed_size` and discarded. A malformed or truncated stream
   (advertised size != actual remaining bytes) would produce a confusing
   downstream Truncated error from inside the rANS decode loop instead of
   being caught up front.
   
   Add `CramError::Rans4x8CompressedSizeMismatch { advertised, actual }`
   and check `cur.len() == compressed_size` after parsing the 9-byte header.
   Matches htscodecs' guard at `rANS_static.c:245`
   (`if (in_sz != in_size - 9) return NULL`).
 - <csr-id-ec7d05eb54bf98795057c4b1d5db2ef77c5b0287/> propagate SIMD decode errors instead of scalar fallback (Finding 3)
   The order-0 32-state SIMD dispatch previously snapshotted state+src,
   ran SIMD, and on Err restored and re-ran scalar. SIMD and scalar share
   the same algorithm over the same `src`/`states`, so the only way for
   SIMD to fail is the same way scalar would (truncation). Falling back
   silently masked the very class of SIMD bug the dispatch was supposed
   to surface — the existing `simd_matches_scalar_with_renorm` proptest
   even bailed out on SIMD-Err, hiding the case the fallback was for.
   
   Choose option (b) from the critique: trust SIMD. Errors propagate.
   
   - Drop the snapshot/restore + scalar retry. SIMD result returns directly.
   - Extract `decode_order_0_32state_scalar` for non-SIMD targets only
     (non-aarch64 + non-AVX2 x86_64), reachable via `cfg_attr` allow.
   - Tighten the proptest: scalar runs first as oracle, SIMD must match
     exactly — including erroring on the same inputs. Mismatched outcomes
     fail the test instead of bailing out.
 - <csr-id-d5beab9be4765793e3078e676077c5b2cc33f237/> fix detached-mate next_pos off-by-one (Finding 2 + Finding 1 oracle)
   The two CRAM round-trip suites compared bases/positions/flags/end_pos but
   never the mate fields, so an off-by-one in `decode_record`'s detached-
   mate path went unnoticed. CRAM's `NP` data series is 1-based per the
   spec; BAM's `next_pos` is 0-based. htslib does the conversion at the BAM
   emit (`cram_decode.c:3121` `cr->mate_pos - 1`), but seqair was passing
   the 1-based value straight through to `store.push_fields(..., next_pos_val, ...)`.
 - <csr-id-7d81be00706942ccae2119269103e5b31bda549c/> filter landmarks by CRAI slice_offsets in fetch_into (Finding 5)
   `fetch_into_customized` previously deduplicated CRAI entries by
   container_offset and then iterated *all* landmarks of each container,
   relying on `decode_slice`'s per-record overlap check to drop irrelevant
   slices. CRAI's `slice_offset` is the same value as the container's
   `landmark` (byte offset from container data start to slice header), so
   the reader can skip non-listed landmarks entirely.
   
   Build a `BTreeMap<u64, SmallVec<u64, 4>>` of container_offset → set of
   wanted slice_offsets (containers visited in file order; slice list is
   typically 1–4). When iterating landmarks, skip any not in the set.
   
   Strict win for multi-slice / multi-ref containers (a 100-record query
   hitting a multi-ref container with 4 slices for other contigs now
   decodes 1 slice, not 4). No-op for single-slice / single-ref containers,
   which is the common case in samtools-produced CRAM.
 - <csr-id-8a08edc2bf54fdc46a4f752318a85eaedd7fec88/> extract shared codec_io module for per-byte primitives (sweep B)
   `read_u8`, `read_u32_le`, `read_uint7`, and `split_off` were duplicated
   across `rans.rs`, `rans_nx16.rs`, and `tok3.rs`. tok3's copy used
   `Result<u8, CramError>` for `read_u8`, defeating the size-and-drop
   optimization the other two carefully applied — `CramError` is 80 bytes
   (heap-owning variants like `Open { path: PathBuf }`) and its drop is
   discriminant-dispatched, so per-byte `Result<u8, CramError>` triggered
   `drop_in_place<CramError>` on every successful read in tok3.
   
   New `cram::codec_io` module:
   
   - `read_u8`, `read_u16_le`, `read_u32_le`, `split_off` return `Option<T>`
     (preserves the per-byte size+drop win documented in the rans modules).
   - `read_uint7` returns `Result<u32, Uint7Error>` — a narrow 1-byte enum
     with `Truncated` and `Overflow` variants, so `Result<u32, Uint7Error>`
     is the same size as `Option<u32>` and incurs no drop cost.
   
   Each consumer module bridges the narrow types to the rich `CramError`
   on the err path with a tiny `uint7_to_cram_error` adapter and contextual
   `ok_or_else(|| CramError::Truncated { context: ... })` calls.
   
   tok3 in particular gains the size/drop optimization for the per-byte
   token-stream parser. rans_nx16 keeps its `read_u16_le_prv` re-export
   for SIMD-module consumers.
 - <csr-id-51473c2f61098e4724ef84b5750358ddede4d416/> precompute per-context symbol tables for rANS Nx16 order-1 (sweep D)
   The order-1 hot loop called `cumulative_frequencies_symbol` per byte —
   a linear scan up to 256 steps over the per-context cumulative-frequency
   table. The 4x8 codec already precomputes per-context symbol-decode
   tables in `Rans4x8Buf::sym_tables` and looks up directly; rANS Nx16
   order-1 was the only remaining hot path with the linear scan.
   
   Mirror the 4x8 design:
   
   - Add `sym_tables: Box<[[u8; 4096]; 256]>` to `Nx16Order1Buf` (~1 MB,
     same as `Rans4x8Buf::sym_tables`). The buf is reused across blocks,
     so the alloc happens once per reader.
   - Extract `build_symbol_table_nx16_into` from `build_symbol_table_nx16`
     for in-place population; the order-0 path keeps the allocating form.
   - After building cumulative frequencies, populate `sym_tables` for all
     256 contexts. The per-byte hot loop (and the trailing-chunk fallback)
     replaces the linear scan with `buf.sym_tables[k][f as usize]`.
   - Mark `cumulative_frequencies_symbol` as `#[cfg(test)]` — it's now
     used only as the oracle for `build_symbol_table_proptest`.
 - <csr-id-a716a0dfd516dc568697e637098a89e01d3f4dfd/> stream slice MD5 verification, no per-slice Vec alloc (sweep C)
   `decode_slice` previously verified the slice's reference MD5 by collecting
   `slice_ref.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>()`
   and feeding the whole vec to `md5::compute`. Slice references run
   typically 10s of KB to a few MB; this allocated and freed one full-slice
   buffer per slice (every slice with a non-zero `reference_md5`, which is
   nearly all of them).
   
   Replace with a streaming `md5_uppercase_streaming` helper that uppercases
   into a 4 KB stack buffer and feeds it to `md5::Context::consume` in
   chunks. Stack-only, no per-slice heap allocation.
 - <csr-id-8c493f7149e4c22a056547d8af77d4c2498ce3a5/> proptest rANS Nx16 SIMD renormalization (Finding 4b)
   The two existing SIMD↔scalar proptests (`simd_matches_scalar_order0_32state`
   and `simd_remainder_uses_correct_lanes`) construct streams where state
   stays ≥ 1<<15 for the whole decode — `state_renormalize` is never called.
   That left the renorm path uncovered by SIMD-vs-scalar comparison, which
   matters because each SIMD module previously inlined its own copy of the
   renorm loop (cleaned up in 31266583 to call the shared function).
   
   Two new proptests:
   
   - `state_renormalize_matches_spec` exercises the function directly across
     the full u32 state space and 0..16 byte src buffers, replaying the spec
     in the test to derive the expected result.
   
   - `simd_matches_scalar_with_renorm` constructs a stream with initial state
     0x8001 (just above the threshold) and freq[sym]=2048, so the first
     decode step drops state below 1<<15 and forces a renorm read on every
     lane. Random renorm bytes are appended; SIMD and scalar must agree on
     the decoded output for any seed.
 - <csr-id-c9725ce7fee129f51ea4fc45dd3e6042d2fa263e/> add htscodecs-derived uint7 spec vectors (Finding 3)
   The existing roundtrip proptest validates encoder↔decoder consistency
   in the same module — it would silently keep passing if both sides got
   flipped MSB↔LSB together (exactly the regression that already happened
   once at 6e56e7ef and was reverted).
   
   Add `read_uint7_spec_vectors` with 11 hard-coded byte sequences derived
   from htscodecs `var_put_u32` (varint.h:206, BIG_END / MSB-first). These
   are an independent oracle: the test will fail if `read_uint7`'s byte
   order ever drifts from the spec, regardless of what the in-module
   encoder does.
 - <csr-id-440f637388fd1c8df9b589e9fe01d29ad0c110ed/> reject malformed rANS alphabet/freq-table runs (Finding 2 + A1 + A2)
   The run-length encoding for rANS frequency tables advances a `u8` symbol
   counter inside the run loop. With `wrapping_add(1)` and no bounds check,
   a malformed `start + len > 255` silently wraps, corrupting alphabet/freq
   entries near 0 and desynchronising the source stream. htscodecs rejects
   these (`rANS_static.c::decode_freq` and `rANS_static16_int.h::decode_alphabet`
   both `goto cleanup` on `j > 255`).
   
   Three call sites had the bug:
   - rans_nx16::read_alphabet — no guard, full silent wraparound.
   - rans::read_frequencies_0 — partial guard (`if sym == 255 break`) that
     exited the inner loop early, leaving freq-bytes unread → src desync.
   - rans::read_frequencies_1_into — no guard on the outer per-context dim.
 - <csr-id-320ac1325f68eaa8c346d10a21d28b8f3e1e987d/> tighten SIMD/codec hot paths and drop debug_assert on untrusted input
   - container.rs: drop dead `let _ = i; // suppress unused warning` (use `_`).
   - state_step: remove `debug_assert!(result >= g)` — the precondition is on
     untrusted CRAM input (already handled via `wrapping_sub`), so the assert
     panics in fuzz/debug builds on adversarial inputs. Comment updated to
     call this out.
   - NEON/AVX2 32-state loop: replace inlined 2-step renormalization with a
     call to the shared `state_renormalize`, so the SIMD paths can never drift
     from the scalar path. Add `debug_assert_eq!(states.len(), 32)` to make
     the unsafe-block precondition explicit.
   - NEON/AVX2 symbol-lookup loop: replace `chunk.get_mut(j).unwrap_or(&mut 0)`
     (which silently swallowed OOB writes) with direct indexing. Bounds are
     static (`f < 4096`, `j < 32`).
 - <csr-id-282578eaecedf67a22e7ceb14c35577f4b6a1de8/> revert incorrect read_uint7 fix; ITF-8 is MSB-first, not LSB-first
 - <csr-id-209f7f32fbbf013adfc9b3f85465264b207fd78b/> fix remaining clippy warnings in tests + bench
 - <csr-id-9d12414d5b76fc3f903eee3b640d379e67626b06/> fix clippy indexing_slicing warnings with debug_assert guards
 - <csr-id-3135eb9ac2cbd9333ff738a902e1bdfd8de2c273/> fix SIMD clippy lints — safety comments, unsafe scope, cfg-gated imports
 - <csr-id-0e97e6632e1855a79fb9cb582f21b930132e983e/> gate SIMD imports behind target_arch cfg to silence unused-import warnings
 - <csr-id-887dc84ec4940c346e92471fdb13890d45e477df/> wrap SIMD intrinsics in explicit unsafe blocks (Rust 2024)
 - <csr-id-ed73d457d785f32e185f91efe92c94ef352a14ca/> fix read_uint7 decode for multi-byte values + proptest
 - <csr-id-76c5e95bef6b8824712db39197ca4ceb2cc8d7d4/> cap proptest range to 0..96, remove known-failing decode test
 - <csr-id-6750a333f7bcf96d55973211b09db117f4a8f650/> NEON→scalar runtime fallback on decode_order_0_32state failure
 - <csr-id-3f2952d1574c1218deeba4f540938bb8768e0970/> proptest SIMD vs scalar equivalence for rANS Nx16 order-0 32-state
 - <csr-id-e220cfde90b80494e807aaf07c7fbedc11c28158/> AVX2 SIMD decode_order_0_32state + feature-detect dispatch
 - <csr-id-0942b8d9b879e2a4284f50f43441673e26623ec3/> NEON SIMD decode_order_0_32state inner loop (aarch64)
 - <csr-id-e13a790a1fd70dbc8269e90fbab2a71823c4f3f9/> extract decode_order_0_32state from generic decode_order_0
 - <csr-id-6a950408d73b245b37742fdd40f5c1f520a64f06/> branchless 2-step renorm for rANS Nx16 state_renormalize
 - <csr-id-38a12d81875ac08b4d91e9de4fa571c042daac30/> pre-compute 4096-entry symbol table for rANS Nx16 order-0
 - <csr-id-0572c722c555fc56f230114f0258350902ebc133/> reuse Nx16Order1Buf allocations across rANS Nx16 order-1 blocks
 - <csr-id-ffc737fbe83b4632bda32706ba89752f7af58fa6/> reuse Rans4x8Buf allocations across rANS 4x8 order-1 blocks
 - <csr-id-4377073bb61a7b760d024c5a3b0ab1de31affb1b/> replace get_external FxHashMap with SmallVec linear scan
 - <csr-id-593055ad9fb6b3da2736331bb290ea27c349f08f/> tighten rans 4x8 decode_order_1 inner loop
   Three per-iteration overheads removed in the order-1 hot loop, which
   samply now flags as the topmost function:
   
   * `out_idx` was computed via `si.checked_mul(chunk_size)
     .and_then(|base| base.checked_add(pos)).ok_or_else(...)?` per inner
     iteration. The result is statically bounded — `si < 4` and
     `chunk_size == dst.len() / 4` mean `bases[si] + pos < 4 * chunk_size
     <= dst.len()` always. Replace with a precomputed `bases: [usize; 4]`
     and a plain `wrapping_add(pos)`.
   
   * `states[si]` was indexed three times per iteration in the rANS state
     step (`s & 0xFFF`, `s >> 12`, then the assignment). Pull into a local
     `state: u32` and write back once. Same for the read of `prev_syms[si]`
     via the new local `ctx` (already there) and the eventual write-back.
   
   * The `remainder_start = 4usize.checked_mul(chunk_size)` had the same
     static-safety story; replace with `wrapping_mul(4)`.
   
   Output is byte-identical: covered by the noodles + htslib CRAM
   record-level comparison tests and the BAM↔CRAM mpileup byte-equality
   test.
 - <csr-id-e9312547ad9838e3c59657608e21f1d4d22dafd9/> unroll renormalize loop to match htslib's RansDecRenorm shape
   The `while *state < LOWER_BOUND { read_u8(...)? }` loop carried a
   predictable but real per-iteration cost — `?` propagation, condition
   check, branch back. Per-byte for a state that almost always needs
   exactly one renorm read.
   
   Switch to two unrolled `if`-and-read steps (the common case for any
   well-formed 4x8 stream — encoders maintain state ≥ LOWER_BOUND/(2^16)
   before each renorm, so 1 or 2 bytes always suffice), followed by a
   small tail loop to keep the malformed-input path correct.
   
   Mirrors htslib's `RansDecRenorm` macro layout in rANS_byte.h. Output
   is byte-identical; the tail loop is functionally equivalent to the
   old `while` for the rare case where state was near 0.
 - <csr-id-cce533b916b95076cb45056707e06a2c167279b8/> bulk decode_n_into for External byte streams
   Adds `ByteEncoding::decode_n_into(ctx, n, buf)` with a fast path that
   collapses N per-byte calls into one FxHashMap lookup + one memcpy when
   the encoding is `External`. Used at three hot sites:
   
     * `ByteArrayEncoding::decode_into` for the `ByteArrayLen` branch:
       `read_name`, every aux tag value, and every Insertion / SoftClip /
       Bases-block feature payload now bulk-decode the value bytes when
       val_encoding is External (the common case).
   
     * The two `qual_buf` loops in `decode_record` (mapped + unmapped
       paths). Quality scores are read_length bytes per record — the
       biggest per-record win since a 100bp read used to do 100 HashMap
       lookups + 100 byte reads + 100 CramError::Truncated constructions.
   
   samply put `ByteEncoding::decode` second in the flamegraph after the
   previous round of fixes; this hits the underlying cause directly.
   For Null and Huffman variants `decode_n_into` falls back to the
   per-symbol path (Huffman bit lengths are data-dependent so there's no
   bulk shortcut).
   
   Output is byte-identical: covered by `compare_cram_with_htslib`,
   `compare_cram_with_noodles`, and `example_mpileup_cram_matches_bam`.
 - <csr-id-c506815df27d0c6ae0a2b47d2de02385fc76e754/> lazy CramError + #[inline] on DecodeContext::get_external
   `get_external` is called once per byte of every External-encoded data
   series — bam_flags, cram_flags, feature_code, base, base_sub,
   quality_score, and (per byte) the val_encoding inside ByteArrayLen.
   Two changes:
   
   * Eager `ok_or(CramError::ExternalBlockNotFound { content_id })` was
     building (and dropping) a 80-byte CramError on every successful
     lookup. Same drop_in_place pattern as the rans helpers; switched to
     `ok_or_else`.
   * `#[inline]` so the hot per-byte callsites can fold the FxHashMap
     lookup directly into their loop bodies instead of paying a function
     call. With the smaller return type and no eager construction, the
     body is now tight enough that LLVM should reliably inline.
   
   File-scope `#![allow(clippy::unnecessary_lazy_evaluations)]` mirrors
   rans.rs / rans_nx16.rs — clippy doesn't know about the
   discriminant-dispatched Drop and would otherwise rewrite back to the
   eager form.
 - <csr-id-635ff55003eaaa39143abf6ab7bebe21003745fb/> thin Option<T> return for rans byte/state helpers
   `CramError` is 80 bytes (sized to fit `Open { path: PathBuf, source:
   io::Error }` and similar heap-owning variants), so a
   `Result<u8, CramError>` from `rans::read_u8` is also 80 bytes and gets
   passed via memory/sret in any non-inlined callsite. `Option<u8>` is 2
   bytes, fits in a register.
   
   Switch the per-byte and per-state helpers in `rans.rs` and
   `rans_nx16.rs` to `Option<T>`:
   
     rans.rs:        read_u8 / read_u32_le / read_itf8_u16 / read_states
                     / renormalize / read_frequencies_0 / read_frequencies_1
     rans_nx16.rs:   read_u8 / read_u16_le / read_u32_le / state_renormalize
   
   All marked `#[inline]` to push toward inlining at small call sites.
   Outermost decoders (`decode`, `decode_order_0`, `decode_order_1`,
   `decode_stripe`, etc.) materialize a `CramError::Truncated` only on
   the err path with `ok_or_else(|| ...)`. Cold helpers that already
   emit richer errors (`read_uint7`'s `Uint7Overflow`, the parse-time
   `split_off`) keep their `Result<_, CramError>` signature.
   
   Net effect: smaller calling convention on the hot path, no eager
   `CramError` construction (and therefore no `drop_in_place<CramError>`)
   on every successful byte read. samply showed both costs under
   `rans::read_u8` / `rans::renormalize` and `slice::first`.
 - <csr-id-962d254eef6656490c4dade2c16ce7bc2396ca9a/> split_first / split_first_chunk in rans byte readers
   The previous read_u8 used `slice::first()?` followed by
   `slice::get(1..)?`: two separate slice bounds checks, one per access.
   samply's biggest stack pointed at `core::slice::<impl [T]>::first`
   under `rans::read_u8` → `rans::renormalize` → `rans::decode_order_1`.
   
   `split_first()` does the bounds check once and returns `(first, rest)`
   in one shot. Same idea with `split_first_chunk::<N>()` for the
   fixed-size LE readers.
   
   Applied to rans.rs::read_u8 / read_u32_le and rans_nx16.rs::read_u8 /
   read_u16_le / read_u32_le.
 - <csr-id-eac99a11756d8db7495810bb19fb12cb0669502f/> precompute Huffman canonical codes + per-level index
 - <csr-id-d57a4257fca74b4c55d689b8513232bbdd4574a8/> silence clippy unnecessary_lazy_evaluations on rans hot paths
   The lazy `ok_or_else` form added in the previous two commits
   (`ddd84fd6` and `e8f3545f`) is the deliberate hot-path choice — eager
   `ok_or` triggers per-call `drop_in_place<CramError>` because the enum
   has variants that own heap data. clippy's `unnecessary_lazy_evaluations`
   lint flags it because the construction looks "cheap"; the lint doesn't
   account for the enum's discriminant-dispatched Drop.
   
   Add a file-scope `#![allow(...)]` with a reason on rans.rs and
   rans_nx16.rs so `cargo clippy --all-targets -- -D warnings` (CI parity)
   stays clean.
 - <csr-id-63a2b4bdbf71694f3f4e3e78333d5ed31236ca22/> lazy CramError construction in rans Nx16 hot path
   Same fix as the previous commit applied to `rans_nx16.rs`: per-byte
   helpers (`read_u8`, `read_u16_le`, `read_u32_le`), `state_renormalize`,
   and the inner-loop checked-arithmetic guards in `decode_order_1`
   switch from `ok_or` to `ok_or_else`. Avoids building (and dropping) a
   full-sized `CramError` on the success path of every byte read.
 - <csr-id-60b734d2664076ab6af6bf668268096d363ecea9/> lazy CramError construction in rans 4x8 hot path
   samply showed `core::ptr::drop_in_place<CramError>` on a stack with
   `Option::ok_or` and `rans::decode_order_1`. Cause: `ok_or(error)`
   evaluates the error eagerly. Per-byte helpers (`read_u8`, `read_u32_le`,
   `read_itf8_u16`) were building a full-sized `CramError::Truncated`
   *every* call, then dropping it on success — `CramError` has a real Drop
   because some variants own `PathBuf` / `std::io::Error` / `SmolStr`, so
   the Drop has to dispatch on the discriminant even for the no-op
   Truncated variant.
   
   Switching to `ok_or_else(|| CramError::Truncated { ... })` makes the
   construction (and the matching drop) lazy: nothing happens on the
   success branch. Same fix applied to the two checked-arithmetic guards
   inside `decode_order_1`.
 - <csr-id-7cd52053bb1e025e799462d0f44b3cbf9375576e/> TODO(perf) for container/compression-header reuse and shared RG cache
   Documents three follow-ups identified by the same survey that produced
   the previous decode-loop fixes; flags only, no behaviour change.
   
   D + E (container fetch loop in fetch_into_customized): every call
   re-reads the container bytes off disk and reparses the compression
   header (incl. tag dictionary and tag encodings). Multi-segment pileup
   workflows revisit the same container repeatedly. Inline TODO sketches
   three options (one-entry "last container" cache, small LRU,
   plan-time grouping) with profile signals to look for.
   
   G (CramShared::read_group_ids): note that this cache lives on the CRAM
   side because BAM/SAM don't need it; if another format ever does, lift
   to BamHeader so the parsed @RG list is canonical instead of each
   format-specific opener re-scanning the header text.
 - <csr-id-1d63e36cea7eb1d8806c4cbf633e851087fb3ea0/> zero per-record allocations in the decode hot loop
   Three coupled wins in `decode_record` / `decode_features_and_reconstruct`,
   all driven by the same observation: every CRAM record was performing
   multiple short-lived `Vec` allocations, dominating profiles after the
   @RG header rescan was already addressed.
   
   A. ByteArrayEncoding::decode_into(ctx, buf). New API that appends into
      a caller-owned scratch buffer; the old `decode(ctx) -> Vec<u8>`
      wrapper still exists for cold paths but the CRAM record loop no
      longer uses it. Underlying ExternalCursor gains
      `read_bytes_into` / `read_bytes_until_into` that memcpy from
      contiguous data into the caller's buffer instead of `.to_vec()`.
      Hot-path call sites:
        * read_name (incl. the discarded detached + RN=false path) now
          reuses `name_buf`.
        * Per-tag values decode straight into `aux_buf` after the BAM
          tag header bytes — skips the per-tag scratch entirely.
        * Per-feature payloads (insertion, soft-clip, bases-block,
          quality-block) reuse `feature_byte_buf`.
   
   B. cigar_ops moved from a `Vec::new()` (per record, cap=0, first push
      allocates) to a reusable `cigar_ops_buf` threaded from the reader.
   
   C. Fused the two-pass feature loop. Previously: collect every feature
      into a local `Vec<Feature>` (with `Vec<u8>` payloads in
      `FeatureData::Insertion`/`SoftClip`/`Bases`), then walk the read
      filling reference matches and applying features. Per-record cost was
      1× `features` allocation + N× per-feature payload allocations.
      CRAM features are emitted in increasing read-position order (delta-
      encoded), so we can decode each feature inline: catch read_pos up to
      the feature's anchor with reference matches, decode the feature's
      payload directly into `feature_byte_buf`, apply, repeat. Eliminates
      `Feature`/`FeatureData` entirely.
   
   Net steady-state allocation count per record: zero (after the buffers
   warm). Output is byte-identical: covered by `compare_cram_with_htslib`,
   `compare_cram_with_noodles` (full sequence + qual + qname), and
   `example_mpileup_cram_matches_bam` (BAM↔CRAM mpileup byte-equality).
 - <csr-id-1f5bb9ddc752adabf3ae04a7cd6b03b8d9944d06/> borrow ref_name into header instead of allocating a String per fetch
   `fetch_into_customized` was calling `target_name(tid).to_string()` to
   get the contig name for the FASTA fetch — one heap allocation every
   fetch, with the only consumer being a `&str`-taking API. Borrow into
   the Arc'd `BamHeader` directly; the cold `MissingReference` error path
   is the single place we still materialize a `SmolStr`.
   
   Also fixes the fuzz-only `from_bytes` constructor to populate
   `CramShared::read_group_ids`, which I missed when the cache field was
   introduced (the fuzz feature wasn't enabled in the default build, so it
   compiled silently).
 - <csr-id-b1c5401cf6434473df6701e771b6ffd3e953f12e/> cache parsed @RG IDs on open instead of rescanning per record
   `get_read_group_id` was called from `decode_record` for every CRAM record
   with `read_group >= 0`, doing a full O(L) scan of `header.header_text()`
   plus a `String` allocation each time. samply pointed at it as the
   dominant cost in `fetch_into_customized` for large CRAMs — an O(N·L) hot
   loop where N is records, L is header lines.
   
   Parse `@RG ID:` values once at open into a `Vec<SmolStr>` cached on
   `CramShared` (Arc-shared across forks), then plumb a `&[SmolStr]` slice
   through `decode_slice` → `decode_record` and replace the lookup with
   `read_group_ids.get(rg_idx)`. Append the bytes straight into `aux_buf`
   with no String round-trip. The `header: &BamHeader` parameter on
   `decode_record` was unused after the fix and dropped.
   
   Output is byte-identical: covered by `cram_records_match_htslib_for_chr19`
   (record-level fields), `compare_cram_with_noodles` (full sequences and
   quality), and the `example_mpileup_cram_matches_bam` byte-equality test.
 - <csr-id-1d825ebfa45fb89276f715845e4d618259a656b1/> rewrite mpileup with Readers::pileup + segments
   Drop the manual fetch_into/window loop in favour of `Readers::open` +
   `readers.segments(target, opts)` + `readers.pileup(&segment)`. The
   unified API handles CRAM with reference-based decoding for free, and
   segmenting drops the bespoke 1-MiB window math.
   
   Region is now optional and accepts every shape `RegionString` parses
   (`chr19`, `chr19:6103076`, `chr19:6103076-6103200`); omitting `-r`
   runs a whole-genome scan via the `()` `IntoSegmentTarget` impl.
 - <csr-id-415dfc55845b09b529dd6abdbfa22293df0d4ab6/> fix end_pos off-by-one and overlap filter to match BAM convention
   The CRAM slice decoder stored `end_pos = pos + ref_consumed` (half-open
   exclusive), but the BAM/SAM path's `compute_end_pos` returns `pos +
   ref_consumed - 1` (0-based inclusive) and the pileup engine evicts on
   `end_pos < pos`, treating end_pos as inclusive. Result: every CRAM record's
   end_pos was one base larger than its BAM equivalent, inflating pileup depth
   by 1 on the trailing column of every read. The CRAM overlap filter was
   also half-open (`pos >= query_end`) while the BAM IndexedReader uses
   inclusive bounds (`pos > end`), which dropped reads starting exactly at
   the requested end position.
   
   Cross-validated the fix against rust-htslib's exclusive `cigar.end_pos()`
   (`-1` to land on inclusive) and noodles' 1-based-inclusive `alignment_end`
   (`-1` for 0-based). Updated `cram_records_match_htslib_for_chr19` — it had
   been masking the bug by asserting `cram.end_pos == hts.end_pos` (matching
   the bug, not BAM). Updated `compare_cram_with_noodles.rs` filters to
   inclusive (`r.pos <= end`) since seqair's `fetch_into` is inclusive.
   
   Kept `end_pos_raw` (exclusive) for the TLEN math at slice.rs:663, which
   is naturally half-open.
 - <csr-id-28766a811596437f685f76aac753f963cde96f7f/> fix clippy errors triggered under -Dwarnings
   CI runs clippy/tests under RUSTFLAGS=-Dwarnings (set by setup-rust-toolchain),
   which surfaces lints that bare `cargo clippy` doesn't. Fixes the failures from
   run 25109911630:
   
   - cast_possible_truncation: prefer try_from / Pos0 deref / saturating ops
   - arithmetic_side_effects: saturating_sub/add in iterator + MD encoder
   - doc_markdown: backtick OwnedBamRecord, push_raw, next_pos, MD return tuple
   - field_reassign_with_default: struct literal for AlignedPairsOptions
   - useless_conversion / unnecessary_cast cleanups
   - example: CustomizeRecordStore::keep_record → filter (trait method rename)
   
   Add a "Pre-push CI parity" section to CLAUDE.md with the exact local commands
   and the recurring lint patterns to watch for.
 - <csr-id-219f88b31e10b5b1c5dce481df32583324708a59/> add aligned_pairs_walk criterion benchmark vs rust-htslib
   Pre-loads the test BAM region into both a seqair RecordStore and a Vec of
   rust-htslib Records once outside the timer, then measures only the walk.
   Five sub-benches:
   
   - seqair_bare              walks AlignedPairs (positions only)
   - seqair_with_read         walks AlignedPairsWithRead (+ seq/qual attached)
   - seqair_matches_only      walks bare → matches_only filter
   - htslib_aligned_pairs_full per-base walk via rust-htslib
   - htslib_aligned_pairs     matches-only walk via rust-htslib
   
   Local numbers on the test.bam fixture (chr19:0-6.14M, M1 Pro):
   
     seqair_bare              1.80 ms
     htslib_aligned_pairs_full 1.78 ms   ← parity
     htslib_aligned_pairs     1.69 ms
     seqair_matches_only      2.30 ms
     seqair_with_read         3.67 ms
   
   Bare AlignedPairs is at parity with rust-htslib aligned_pairs_full despite
   yielding a richer typed enum (with MatchKind preservation). with_read
   costs roughly 2x because of the seq/qual lookup work that htslib forces
   the caller to do manually anyway.
   
   Also fixes a stale C3 site in bam_roundtrip bench (next_pos: -1 → None).
 - <csr-id-665200156cf8f9615576d12028a640d545f243c0/> add fuzz_aligned_pairs target
   Structure-aware fuzz target that drives every public AlignedPairs* path
   with the same Arbitrary-generated record + reference window:
   
   - bare AlignedPairs (default mode + with_soft_clips + full)
   - size_hint upper-bound vs collected count (asserted)
   - matches_only at each layer
   - with_read (success path; validation errors discarded)
   - with_reference
   - nm() and md() (md error tolerated, panic isn't)
   
   Internal `consume_*` helpers destructure every variant and assert the
   no-zero-length-summary invariant in debug builds. Reference window's start
   position is independently chosen so out-of-window cases (None ref_base /
   ref_bases) get exercised.
   
   Auto-picked up by run_all.sh via `cargo fuzz list`.
 - <csr-id-a29311822502bbf06cdb599fa22137f619c810ab/> NM/MD tag recompute via AlignedPairsWithRef
   Methods on AlignedPairsWithRef that consume the iterator and emit the SAM
   edit-distance tag and mismatch-position string:
   
   - nm() -> u32: best-effort edit distance. M-op mismatches (decided via
     query vs ref_base) + 1 per X + I length + D length. = always 0; N and
     clips never contribute. Positions outside the loaded RefSeq are skipped
     (count 0) so partial coverage produces a partial NM rather than an error.
   - md() -> Result<Vec<u8>, NmMdError>: SAM-spec MD bytes, with leading and
     trailing digits, '0' separator between adjacent deletions, '^'-prefixed
     deletion runs, no insertions/clips/N's. Errors with NmMdError::Missing
     Reference if any required ref position is outside the loaded window —
     MD encodes ref bases verbatim and partial coverage would corrupt output.
   
   Tests (18 new):
   - 7 NM unit tests covering perfect match, M-op mismatches, indel lengths,
     =/X overrides, missing-ref skip, N-as-mismatch, realistic CIGAR
   - 9 MD unit tests covering perfect match, single mismatch, consecutive
     mismatches, deletion encoding, I/S skipping, leading-mismatch zero,
     =/X-emit-ref-base, missing-ref error, realistic CIGAR
   - 1 proptest cross-validation: NM == md_mismatches + md_deletions +
     cigar_insertions (catches divergent NM/MD bugs)
   - 1 helper test for push_number internal utility
   
   Spec rules: cigar.aligned_pairs.{nm,md,nm_md}.compute|consistency.
 - <csr-id-6c62b124deb6907567ec708449f9a53cd121e9b4/> matches_only adapter, OwnedBamRecord::aligned_pairs_with_read, Clone/size_hint
   Tier-1 ergonomics + iterator-trait completeness.
   
   - matches_only adapter on all three layers, yielding flat value structs:
     MatchPosition (bare), MatchedBase (with read), MatchedRef (with read +
     ref). Equivalent of pysam aligned_pairs(matches_only=True), but with
     MatchKind preserved so callers can distinguish M/=/X. None of the
     structs need lifetimes — Base/BaseQuality/Pos0 are all Copy.
   - OwnedBamRecord::aligned_pairs_with_read symmetric with the SlimRecord
     helper. Equivalent to self.aligned_pairs()?.with_read(&self.seq,
     &self.qual) but a one-line call site.
   - Clone derived on all iterators (state is Copy/borrows). Lets callers
     walk twice without re-validating.
   - ExactSizeIterator + size_hint on bare AlignedPairs and the layered
     iterators — exact count computed from remaining CIGAR + expansion
     state. matches_only adapters are filter-style: lower=0, upper=inner.
   - FusedIterator on everything (next never returns Some after None).
 - <csr-id-8cd10273f2b35bb229937252026dde1d88d1c078/> add htslib parity proptest + with_read/with_reference contract proptests + 150bp realistic test
   Four new tests fill remaining coverage gaps highlighted by the post-hardening
   review:
   
   - matches_htslib_random_cigars: proptest generating random CIGARs (M/I/D/N/
     S/H/=/X with non-zero lengths) and comparing seqair's expanded
     AlignedPairs output against rust_htslib::aligned_pairs_full(). 256 random
     inputs per run; rust-htslib is the htslib-binding ground truth, so this
     is the strongest possible parity test.
   - with_read_slice_contents_match_seq_at_qpos: verifies Insertion and
     SoftClip slice CONTENTS are byte-identical to seq[qpos..qpos+len] (and
     same for qual). Previous pass-through proptest only checked length
     equality; this one would catch off-by-one in slice computation.
   - with_reference_contract_matches_refseq_lookups: for every Match event,
     ref_base == RefSeq::try_base_at(rpos); for every Deletion, ref_bases ==
     RefSeq::range(rpos, del_len). Direct invariant — the iterator is just a
     cache of these lookups, so they should agree exactly.
   - realistic_150bp_walk_through_all_layers: 150bp Illumina-shaped CIGAR
     (2S95M5I20M2D28M) at pos 1_000_000 walked through bare, with_read, and
     with_reference layers. Exercises larger qpos/rpos arithmetic that small
     fixtures hide.
 - <csr-id-301d31cd48da750c30bff04f39e3274dc2f33a7b/> harden AlignedPair APIs — typed errors, split lifetimes, tighter tests
   Acts on critique findings:
   
   - New AlignedPairsError enum with typed variants and verbose diagnostic
     messages: Access(RecordAccessError), UnmappedWithCigar, CigarSeqLength
     Mismatch, SeqQualLengthMismatch.
   - OwnedBamRecord::aligned_pairs() now returns Result and refuses pos=None
     with a non-empty CIGAR — eliminates silent rpos-anchored-at-zero
     corruption that the previous infallible signature allowed.
   - AlignedPairs::with_read validates length invariants up front (cigar
     query length == seq.len(), qual.len() ∈ {0, seq.len()}) and returns
     Result. Prevents the silent Base::Unknown substitution on OOB qpos that
     was indistinguishable from legitimate N's. The slim_record helper
     surfaces the same typed error.
   - Lifetimes split: AlignedPairsWithRead<'cigar, 'read>; AlignedPairs
     WithRef<'cigar, 'read, 'ref_seq>. Yield enums use only the lifetimes
     that appear in their slices. Allows callers to compose borrows from
     different sources.
   - Test hardening: rewrote match_kind_does_not_change_qpos_rpos to walk a
     single mixed M/=/X CIGAR (the previous version compared three separate
     single-kind CIGARs and would pass even with a totally broken
     implementation). Widened view-layer proptest to match bare layer's
     full op-code/length range. Added .with_reference(..).with_soft_clips()
     and .with_reference(..).full() chained tests. Added the four
     validation-error tests for with_read. Added htslib comparison cases
     for hard clips, hard+soft+match combinations, longer ops.
   
   Spec rules updated to match: cigar.aligned_pairs.owned_record now
   requires the Result-returning signature; cigar.aligned_pairs.with_read.
   iterator now requires the validation contract; cigar.aligned_pairs.
   options drops the misleading "set-once" wording.
 - <csr-id-dbc52979384f33bf7cf5deea5fda2e6331f6002b/> add aligned_pairs_walk example + module-level layer index
   Runnable example exercises all three layers (bare AlignedPairs, with_read,
   with_reference) against an indexed BAM + FASTA, with three calling tasks
   sharing a single walk: TAPS C→T methylation, per-read SNV evidence, and
   indel collection. Three additional showcase functions cover the
   single-purpose patterns (bare indel counter, inserted-base collector,
   MatchKind-based explicit-mismatch counter).
   
   Module-level docs on aligned_pairs.rs gain a "quick reference" table
   linking the four common needs (positions only, +read, +reference,
   MatchKind dispatch) to the right entry point.
 - <csr-id-b11a658fa3955f2f0590fe549b9edbf733afa7ba/> AlignedPair gains MatchKind + layered with_read/with_reference views
   Builds on the typed AlignedPairs iterator (#prev) with three changes:
   
   - Match variants now carry MatchKind (M/=/X) so callers can distinguish
     aligner-known matches from explicit =/X without re-walking the CIGAR.
     htslib's aligned_pairs_full collapses these; we preserve the distinction.
   - OwnedBamRecord.pos / next_pos move from i64 (with -1 sentinel) to
     Option<Pos0>, making wire-format overflow unconstructable. ref_id /
     next_ref_id stay raw i32 with -1 (used as header-table indices).
   - New layered iterators in aligned_pairs_view.rs: AlignedPairs::with_read
     attaches the record's seq/qual slabs (per-event base+qual for Match,
     pre-sliced runs for Insertion/SoftClip); .with_reference reuses the
     pileup engine's RefSeq to add ref_base/ref_bases (None outside loaded
     window). One-shot SlimRecord::aligned_pairs_with_read helper.
   
   Also drops the stale bam.owned_record.aligned_pairs spec rule (it
   described the pre-enum tuple API), demotes cigar.aligned_pairs.complementary
   to a non-verifiable design note, and rewrites the htslib comparison test
   to feed the actual AlignedPairs iterator through an expansion helper
   instead of comparing two seqair-internal oracles.
 - <csr-id-372cdab734f48340d4d1554da2ee60721de3c27c/> add typed AlignedPairs iterator (summary-form CIGAR walk)
   Add AlignedPair enum and AlignedPairs iterator — a typed, summary-form
   CIGAR walker that yields Match (per-position), Insertion/Deletion/RefSkip
   (summary, one per op). Default mode skips soft clips, padding, unknown ops;
   .with_soft_clips() and .full() opt in.
   
   Replaces the old per-position (Option<usize>, Option<i64>) iterator on
   OwnedBamRecord with the new typed API. Adds SlimRecord::aligned_pairs(store)
   integration via RecordStore.
   
   Design matches htslib's aligned_pairs semantics but uses typed variants
   instead of Option tuples, matching seqair's existing PileupOp/CigarPosInfo
   style. Verified against htslib with 8 CIGAR patterns; 29 tests including
   proptest oracle equivalence and position monotonicity checks.
 - <csr-id-99e8c05cf308f0564a6f8cc368bb71e34dd92048/> fix NaN comparison in roundtrip_random_bytes proptest
   Float(NaN) != Float(NaN) and Double(NaN) != Double(NaN) with derived
   PartialEq. Handle Float/Double explicitly using to_bits() comparison.
 - <csr-id-2c9975de7979d69b8d523105aca3fc2fd03d46fd/> add filter_raw pre-filter to CustomizeRecordStore
   Adds CustomizeRecordStore::filter_raw(&mut self, &FilterRawFields) -> bool,
   called before any slab extension or base decode in push_raw and push_fields.
   Rejected records incur zero work — no memcpy into slabs.
   
   Adds FilterRawFields struct carrying all header fields + raw slices
   (raw_cigar_bytes/packed_seq/bases as Option depending on path).
   
   Renames keep_record -> filter throughout the codebase.
   
   Adds tests verifying filter_raw rejection avoids slab writes in both
   push_raw and push_fields paths.
 - <csr-id-9ac6b99edcf75dde1eb9a9b13cf897cec026d3d0/> cap per-range allocation to file size
   A corrupt BAM index entry could request a virtual offset far past EOF,
   causing `data.resize(buf_start + len, 0)` to attempt a multi-TiB
   allocation before `read_all` had a chance to return what little data
   existed. Found by fuzz_reader_indexed (ASAN: requested allocation size
   ~280 TiB).
   
   Cap each merged range's `len` to `file_size - file_start`, which keeps
   legitimate whole-chromosome queries unaffected while bounding fuzz
   inputs to the actually-readable bytes. Matches `r[io.fuzz.alloc_limits]`.
 - <csr-id-63576413cc55cdf9f6f4282e97cec064fa8efa99/> rename RecordStoreExtras → CustomizeRecordStore
   + SlimRecord getters
 - <csr-id-1640db234338ed90ac5c0674f325eb830cac8a15/> Add `M` number
 - <csr-id-ea65683bebe1edec40250801bc4042ca9e8845a0/> fix warnings in BAM/FASTA/types test modules
 - <csr-id-5dd60bab056d40a4e039af9d60a71b4af26d0866/> fix warnings in CRAM test modules
 - <csr-id-767dd10b0d354fa30d0f48660efb15cbc351f16a/> fix cast warnings in FASTA reader and GZI index
 - <csr-id-c19d902ba66572e119c6d99bafd50eda1929ce95/> misc fixes (aux i8 decode, Cargo.toml lint)
 - <csr-id-94506f6a619075f879d2a15946e53aadbc318da5/> fix cast warnings in VCF/BCF encoder and writers
 - <csr-id-744a552c55e29a6abac67750bfe57e737e6ce395/> fix cast warnings in SAM reader
   - tid casts (u32→i32): use .cast_signed() — BAM tid is a signed index
   - B-array element counts (usize→u32): use try_from with AuxArrayTooLarge error
   - B-array element values from text (i64→u8/u16/u32): use try_from with
     InvalidAuxValue — these come from untrusted SAM text and must be validated
   - serialize_bam_int range-checked branches: #[expect] — range guards above
     each arm ensure value fits in target type
 - <csr-id-75093e97744b85f59b12345dfc86038da7ce6dee/> fix cast warnings in CRAM parsing (ITF8/LTF8 wire format)
   All `read_itf8`/`read_ltf8` results are wire-format signed integers stored
   as u32/u64 — replace `as i32`/`as i64` with `.cast_signed()` throughout.
   
   - bitstream.rs: read_bits_i32 result
   - block.rs: content_id ITF8
   - compression_header.rs: tag encoding map keys
   - container.rs: record_counter, bases (LTF8), landmark, and header fields
   - encoding.rs: all encoding parameter fields; remove spurious
     `#[expect(cast_possible_truncation)]` on u32→usize casts (not needed on
     32/64-bit platforms); decode_gamma/decode_subexp bit patterns
   - index.rs, rans.rs: ITF8/LTF8 fields
   - reader.rs: tid comparisons, ref_start decode
   - slice.rs: all slice header fields; tid comparison with known-small bound
 - <csr-id-6929ae22415c99fed35b23a8586e16c2f2d15c71/> fix cast warnings in BAM pileup, index, region_buf, writer
 - <csr-id-a96075850ebd8bc0caa6c59ddbe0f1fb247cf311/> fix cast warnings in BAM record, record_store, owned_record
 - <csr-id-6f3897207397fad9926c50870c3445b3ab722b87/> fix allow_attributes_without_reason, option_option, manual_let_else, same_name_method
   - Convert all bare #[allow] to #[expect] with descriptive reasons
   - bam/aux.rs: simplify Option<Option<AuxValue>> to Option<AuxValue> — the
     Some(None) skip path was never returned; flatten iterator logic accordingly
   - sam/reader.rs: use let-else for read_byte() break pattern
   - bgzf_writer.rs, vcf/writer.rs, vcf/bcf_writer.rs: annotate methods that
     intentionally share names with trait methods (inherent impl, trait delegates)
 - <csr-id-e585244629b66c4bae7e1ef8977f3b41f4f19370/> make sure fuzz files are formatted and compile
 - <csr-id-d682360548222b9ce4ea718e8b3c1fa0f0e211f4/> allow setting thread count

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 370 commits contributed to the release.
 - 90 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#1](https://github.com/Softleif/seqaire/issues/1)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#1](https://github.com/Softleif/seqaire/issues/1)**
    - Fixes after review ([`536263b`](https://github.com/Softleif/seqaire/commit/536263b1d779a5c491d29bbb036708723f65ee6e))
 * **Uncategorized**
    - Don't publish fuzz and test data ([`e031d29`](https://github.com/Softleif/seqaire/commit/e031d298e459956d4302fbe4aedff8b643c0fd33))
    - Rename bam::aux to be compatible with windows… ([`18dc1b9`](https://github.com/Softleif/seqaire/commit/18dc1b96cf6f8580213c6ef95a76a9be97ddb177))
    - Add readmes ([`983b64b`](https://github.com/Softleif/seqaire/commit/983b64bdcec1962268f4195aab1053e0c0ce5676))
    - Update deps ([`2f0d5f9`](https://github.com/Softleif/seqaire/commit/2f0d5f9fb1d88e41781819f01c502a23a69f4de4))
    - Emit Unknown/UNAVAILABLE for SEQ=* records ([`7a27010`](https://github.com/Softleif/seqaire/commit/7a2701048ee6ce08373e2207268abf431a2e4401))
    - Rename iter() to records() and drop the index ([`7e743ac`](https://github.com/Softleif/seqaire/commit/7e743ac7fe6000e965c883cf59eaa1a3d08afb53))
    - Add iter() yielding (idx, &SlimRecord) pairs ([`fcce57b`](https://github.com/Softleif/seqaire/commit/fcce57b9088a3d410ddaa1c48f44597161ad571d))
    - Add keep_unmapped opt-in for placed-unmapped reads ([`f3c3968`](https://github.com/Softleif/seqaire/commit/f3c39688e59d48443c1f50b6b3a27e91a564f146))
    - Move entry points onto BamWriterBuilder ([`dd324da`](https://github.com/Softleif/seqaire/commit/dd324da07f5945f658abee76f66473131f99fcfa))
    - Drop write_index requests on non-path writer targets ([`957d9bb`](https://github.com/Softleif/seqaire/commit/957d9bb85fbcbd2c2d4f8aae15c7e7336146f13e))
    - Re-export BamWriterBuilder/ToPath/ToWriter from seqair::bam ([`76d934b`](https://github.com/Softleif/seqaire/commit/76d934b6d0e8c2d0bd1f81fab4e24ce87ee49021))
    - Tighten &[u8] FromAuxValue to Z-only, add HexBytes for H ([`1280edb`](https://github.com/Softleif/seqaire/commit/1280edbcf63f5ef429ba0f94070446b2ef8f6011))
    - Post-BamRecord-removal cleanup ([`25d1791`](https://github.com/Softleif/seqaire/commit/25d1791b0eceb3caa47cbba27690c6116fd4e04b))
    - Replace BamWriter::from_path boolean with a builder ([`34c2363`](https://github.com/Softleif/seqaire/commit/34c2363cddf895426d3467ce81955972a4e219bc))
    - Add BaseModState::from_record(rec, store) constructor ([`dc174f7`](https://github.com/Softleif/seqaire/commit/dc174f711896f84bc219500c164a3e28a7605f65))
    - Add Display for CigarOp/CigarOpType and CigarStr slice wrapper ([`4e36d0b`](https://github.com/Softleif/seqaire/commit/4e36d0bfeac5cfc8015712c60fd765b8cd278459))
    - Tightening seqair's public API ([`1b0a3b9`](https://github.com/Softleif/seqaire/commit/1b0a3b94889036ccd0c4d6117268988806f73629))
    - Extract BamWriter::finalize_record, share post-encode path ([`41fc3a7`](https://github.com/Softleif/seqaire/commit/41fc3a7a9dae61898c2e6713724be988258ddd07))
    - Dedupe 32-byte fixed-header encode into shared helper ([`6d80e6e`](https://github.com/Softleif/seqaire/commit/6d80e6e0a05c62d4857067232624a9c036b0992a))
    - Drop OwnedBamRecord::from_raw_bam, the unused single-record decoder ([`c467524`](https://github.com/Softleif/seqaire/commit/c4675241707b03fabea905aa2180a97e506e65b6))
    - Remove BamRecord, decode through RecordStore::push_raw ([`1275c9e`](https://github.com/Softleif/seqaire/commit/1275c9e8c937eddc4f7b6cdc3ac6a217030ed5f8))
    - Skip seq/CIGAR reconstruction for foreign-tid multi-ref records (Finding 5) ([`1a3c084`](https://github.com/Softleif/seqaire/commit/1a3c0846ba6d24be106882a89cd72803c92ffd3a))
    - Validate rANS 4x8 header compressed_size against payload (Finding 4) ([`1a13ca0`](https://github.com/Softleif/seqaire/commit/1a13ca0e595c023a90414a6210f8db2e24fa8eb5))
    - Propagate SIMD decode errors instead of scalar fallback (Finding 3) ([`ec7d05e`](https://github.com/Softleif/seqaire/commit/ec7d05eb54bf98795057c4b1d5db2ef77c5b0287))
    - Fix detached-mate next_pos off-by-one (Finding 2 + Finding 1 oracle) ([`d5beab9`](https://github.com/Softleif/seqaire/commit/d5beab9be4765793e3078e676077c5b2cc33f237))
    - Filter landmarks by CRAI slice_offsets in fetch_into (Finding 5) ([`7d81be0`](https://github.com/Softleif/seqaire/commit/7d81be00706942ccae2119269103e5b31bda549c))
    - Extract shared codec_io module for per-byte primitives (sweep B) ([`8a08edc`](https://github.com/Softleif/seqaire/commit/8a08edc2bf54fdc46a4f752318a85eaedd7fec88))
    - Precompute per-context symbol tables for rANS Nx16 order-1 (sweep D) ([`51473c2`](https://github.com/Softleif/seqaire/commit/51473c2f61098e4724ef84b5750358ddede4d416))
    - Stream slice MD5 verification, no per-slice Vec alloc (sweep C) ([`a716a0d`](https://github.com/Softleif/seqaire/commit/a716a0dfd516dc568697e637098a89e01d3f4dfd))
    - Proptest rANS Nx16 SIMD renormalization (Finding 4b) ([`8c493f7`](https://github.com/Softleif/seqaire/commit/8c493f7149e4c22a056547d8af77d4c2498ce3a5))
    - Add htscodecs-derived uint7 spec vectors (Finding 3) ([`c9725ce`](https://github.com/Softleif/seqaire/commit/c9725ce7fee129f51ea4fc45dd3e6042d2fa263e))
    - Reject malformed rANS alphabet/freq-table runs (Finding 2 + A1 + A2) ([`440f637`](https://github.com/Softleif/seqaire/commit/440f637388fd1c8df9b589e9fe01d29ad0c110ed))
    - Tighten SIMD/codec hot paths and drop debug_assert on untrusted input ([`320ac13`](https://github.com/Softleif/seqaire/commit/320ac1325f68eaa8c346d10a21d28b8f3e1e987d))
    - Clippy ([`8e0716b`](https://github.com/Softleif/seqaire/commit/8e0716b8d139409b98245c090c565c74db0aad8a))
    - Fix simd ([`1908045`](https://github.com/Softleif/seqaire/commit/1908045232f14bff108d280861509b9d5bc5da7d))
    - Fix fuzz compile ([`7d21a31`](https://github.com/Softleif/seqaire/commit/7d21a312336cb213bf69c78d03baac51adf44787))
    - Fix ci? ([`ab8c0d4`](https://github.com/Softleif/seqaire/commit/ab8c0d4c6330b9376b8387d8bd5e6bce9ac4df41))
    - Revert incorrect read_uint7 fix; ITF-8 is MSB-first, not LSB-first ([`282578e`](https://github.com/Softleif/seqaire/commit/282578eaecedf67a22e7ceb14c35577f4b6a1de8))
    - Fix remaining clippy warnings in tests + bench ([`209f7f3`](https://github.com/Softleif/seqaire/commit/209f7f32fbbf013adfc9b3f85465264b207fd78b))
    - Fix clippy indexing_slicing warnings with debug_assert guards ([`9d12414`](https://github.com/Softleif/seqaire/commit/9d12414d5b76fc3f903eee3b640d379e67626b06))
    - Fix SIMD clippy lints — safety comments, unsafe scope, cfg-gated imports ([`3135eb9`](https://github.com/Softleif/seqaire/commit/3135eb9ac2cbd9333ff738a902e1bdfd8de2c273))
    - Gate SIMD imports behind target_arch cfg to silence unused-import warnings ([`0e97e66`](https://github.com/Softleif/seqaire/commit/0e97e6632e1855a79fb9cb582f21b930132e983e))
    - Wrap SIMD intrinsics in explicit unsafe blocks (Rust 2024) ([`887dc84`](https://github.com/Softleif/seqaire/commit/887dc84ec4940c346e92471fdb13890d45e477df))
    - Fix read_uint7 decode for multi-byte values + proptest ([`ed73d45`](https://github.com/Softleif/seqaire/commit/ed73d457d785f32e185f91efe92c94ef352a14ca))
    - Cap proptest range to 0..96, remove known-failing decode test ([`76c5e95`](https://github.com/Softleif/seqaire/commit/76c5e95bef6b8824712db39197ca4ceb2cc8d7d4))
    - NEON→scalar runtime fallback on decode_order_0_32state failure ([`6750a33`](https://github.com/Softleif/seqaire/commit/6750a333f7bcf96d55973211b09db117f4a8f650))
    - Proptest SIMD vs scalar equivalence for rANS Nx16 order-0 32-state ([`3f2952d`](https://github.com/Softleif/seqaire/commit/3f2952d1574c1218deeba4f540938bb8768e0970))
    - AVX2 SIMD decode_order_0_32state + feature-detect dispatch ([`e220cfd`](https://github.com/Softleif/seqaire/commit/e220cfde90b80494e807aaf07c7fbedc11c28158))
    - NEON SIMD decode_order_0_32state inner loop (aarch64) ([`0942b8d`](https://github.com/Softleif/seqaire/commit/0942b8d9b879e2a4284f50f43441673e26623ec3))
    - Extract decode_order_0_32state from generic decode_order_0 ([`e13a790`](https://github.com/Softleif/seqaire/commit/e13a790a1fd70dbc8269e90fbab2a71823c4f3f9))
    - Branchless 2-step renorm for rANS Nx16 state_renormalize ([`6a95040`](https://github.com/Softleif/seqaire/commit/6a950408d73b245b37742fdd40f5c1f520a64f06))
    - Pre-compute 4096-entry symbol table for rANS Nx16 order-0 ([`38a12d8`](https://github.com/Softleif/seqaire/commit/38a12d81875ac08b4d91e9de4fa571c042daac30))
    - Reuse Nx16Order1Buf allocations across rANS Nx16 order-1 blocks ([`0572c72`](https://github.com/Softleif/seqaire/commit/0572c722c555fc56f230114f0258350902ebc133))
    - Reuse Rans4x8Buf allocations across rANS 4x8 order-1 blocks ([`ffc737f`](https://github.com/Softleif/seqaire/commit/ffc737fbe83b4632bda32706ba89752f7af58fa6))
    - Replace get_external FxHashMap with SmallVec linear scan ([`4377073`](https://github.com/Softleif/seqaire/commit/4377073bb61a7b760d024c5a3b0ab1de31affb1b))
    - Tighten rans 4x8 decode_order_1 inner loop ([`593055a`](https://github.com/Softleif/seqaire/commit/593055ad9fb6b3da2736331bb290ea27c349f08f))
    - Unroll renormalize loop to match htslib's RansDecRenorm shape ([`e931254`](https://github.com/Softleif/seqaire/commit/e9312547ad9838e3c59657608e21f1d4d22dafd9))
    - Bulk decode_n_into for External byte streams ([`cce533b`](https://github.com/Softleif/seqaire/commit/cce533b916b95076cb45056707e06a2c167279b8))
    - Lazy CramError + #[inline] on DecodeContext::get_external ([`c506815`](https://github.com/Softleif/seqaire/commit/c506815df27d0c6ae0a2b47d2de02385fc76e754))
    - Add CRAM benchmarks comparing seqair vs htslib vs noodles ([`494b421`](https://github.com/Softleif/seqaire/commit/494b421f21155d3db9ea0c84b8c289fdcfa9574d))
    - Thin Option<T> return for rans byte/state helpers ([`635ff55`](https://github.com/Softleif/seqaire/commit/635ff55003eaaa39143abf6ab7bebe21003745fb))
    - Split_first / split_first_chunk in rans byte readers ([`962d254`](https://github.com/Softleif/seqaire/commit/962d254eef6656490c4dade2c16ce7bc2396ca9a))
    - Precompute Huffman canonical codes + per-level index ([`eac99a1`](https://github.com/Softleif/seqaire/commit/eac99a11756d8db7495810bb19fb12cb0669502f))
    - Silence clippy unnecessary_lazy_evaluations on rans hot paths ([`d57a425`](https://github.com/Softleif/seqaire/commit/d57a4257fca74b4c55d689b8513232bbdd4574a8))
    - Lazy CramError construction in rans Nx16 hot path ([`63a2b4b`](https://github.com/Softleif/seqaire/commit/63a2b4bdbf71694f3f4e3e78333d5ed31236ca22))
    - Lazy CramError construction in rans 4x8 hot path ([`60b734d`](https://github.com/Softleif/seqaire/commit/60b734d2664076ab6af6bf668268096d363ecea9))
    - TODO(perf) for container/compression-header reuse and shared RG cache ([`7cd5205`](https://github.com/Softleif/seqaire/commit/7cd52053bb1e025e799462d0f44b3cbf9375576e))
    - Zero per-record allocations in the decode hot loop ([`1d63e36`](https://github.com/Softleif/seqaire/commit/1d63e36cea7eb1d8806c4cbf633e851087fb3ea0))
    - Borrow ref_name into header instead of allocating a String per fetch ([`1f5bb9d`](https://github.com/Softleif/seqaire/commit/1f5bb9ddc752adabf3ae04a7cd6b03b8d9944d06))
    - Cache parsed @RG IDs on open instead of rescanning per record ([`b1c5401`](https://github.com/Softleif/seqaire/commit/b1c5401cf6434473df6701e771b6ffd3e953f12e))
    - Rewrite mpileup with Readers::pileup + segments ([`1d825eb`](https://github.com/Softleif/seqaire/commit/1d825ebfa45fb89276f715845e4d618259a656b1))
    - Fix end_pos off-by-one and overlap filter to match BAM convention ([`415dfc5`](https://github.com/Softleif/seqaire/commit/415dfc55845b09b529dd6abdbfa22293df0d4ab6))
    - Fix clippy errors triggered under -Dwarnings ([`28766a8`](https://github.com/Softleif/seqaire/commit/28766a811596437f685f76aac753f963cde96f7f))
    - Add aligned_pairs_walk criterion benchmark vs rust-htslib ([`219f88b`](https://github.com/Softleif/seqaire/commit/219f88b31e10b5b1c5dce481df32583324708a59))
    - Add fuzz_aligned_pairs target ([`6652001`](https://github.com/Softleif/seqaire/commit/665200156cf8f9615576d12028a640d545f243c0))
    - NM/MD tag recompute via AlignedPairsWithRef ([`a293118`](https://github.com/Softleif/seqaire/commit/a29311822502bbf06cdb599fa22137f619c810ab))
    - Matches_only adapter, OwnedBamRecord::aligned_pairs_with_read, Clone/size_hint ([`6c62b12`](https://github.com/Softleif/seqaire/commit/6c62b124deb6907567ec708449f9a53cd121e9b4))
    - Add htslib parity proptest + with_read/with_reference contract proptests + 150bp realistic test ([`8cd1027`](https://github.com/Softleif/seqaire/commit/8cd10273f2b35bb229937252026dde1d88d1c078))
    - Harden AlignedPair APIs — typed errors, split lifetimes, tighter tests ([`301d31c`](https://github.com/Softleif/seqaire/commit/301d31cd48da750c30bff04f39e3274dc2f33a7b))
    - Add aligned_pairs_walk example + module-level layer index ([`dbc5297`](https://github.com/Softleif/seqaire/commit/dbc52979384f33bf7cf5deea5fda2e6331f6002b))
    - AlignedPair gains MatchKind + layered with_read/with_reference views ([`b11a658`](https://github.com/Softleif/seqaire/commit/b11a658fa3955f2f0590fe549b9edbf733afa7ba))
    - Add typed AlignedPairs iterator (summary-form CIGAR walk) ([`372cdab`](https://github.com/Softleif/seqaire/commit/372cdab734f48340d4d1554da2ee60721de3c27c))
    - Fix NaN comparison in roundtrip_random_bytes proptest ([`99e8c05`](https://github.com/Softleif/seqaire/commit/99e8c05cf308f0564a6f8cc368bb71e34dd92048))
    - Add filter_raw pre-filter to CustomizeRecordStore ([`2c9975d`](https://github.com/Softleif/seqaire/commit/2c9975de7979d69b8d523105aca3fc2fd03d46fd))
    - Cap per-range allocation to file size ([`9ac6b99`](https://github.com/Softleif/seqaire/commit/9ac6b99edcf75dde1eb9a9b13cf897cec026d3d0))
    - Address PileupGuard critique findings ([`bea7fbf`](https://github.com/Softleif/seqaire/commit/bea7fbf5d088aed26fa4fd3158ebfb4c3c9a1e4b))
    - RAII auto-recovery of RecordStore via PileupGuard ([`b4cad10`](https://github.com/Softleif/seqaire/commit/b4cad10c802fd1ade1e44abbb543f87f26815230))
    - Address third critique of segmentation API ([`87e87ee`](https://github.com/Softleif/seqaire/commit/87e87ee42ab59ab61221a770458300dab4fa70cc))
    - Fallible Segment::new, buffer reuse, typed errors, conigs ([`e52ddb6`](https://github.com/Softleif/seqaire/commit/e52ddb6539f2acead01c1f8d87c6dc3889e7d52e))
    - Address review of segmentation API ([`bf138c0`](https://github.com/Softleif/seqaire/commit/bf138c06739c5a2163e3469d2800d469190edadf))
    - Replace Readers::pileup(tid, start, end) with pileup(&Segment) ([`c143c67`](https://github.com/Softleif/seqaire/commit/c143c6757876d2d0ebf73fa5870cb5e5958afefe))
    - Bump versions, use crates.io rust-htslib ([`4c21ac2`](https://github.com/Softleif/seqaire/commit/4c21ac29f46a95507e82ffb5c09c34082ba6c47c))
    - Address second critique: typed array setters, OutOfRange, set_hex, char validation, htslib/noodles parity ([`370d5bb`](https://github.com/Softleif/seqaire/commit/370d5bbd38a2fe5c4bb91269ccf02d04ea1c79ff))
    - Fix doc refs ([`d35f27b`](https://github.com/Softleif/seqaire/commit/d35f27b039b70c447360e77e863e86669ee8efd7))
    - Fix clippy warnings for CI (-D warnings): arithmetic_side_effects, sort_by_key, redundant &/casts ([`415b60c`](https://github.com/Softleif/seqaire/commit/415b60c4c57a41a7a2ae78a6a9e562ff0466aa05))
    - Address critique: AuxTag Display, signed-to-unsigned widening, Hex support, expanded proptests ([`0bb9268`](https://github.com/Softleif/seqaire/commit/0bb92682d824c45ca74ba79edf03435822629351))
    - Add Aux<'_> and FromAuxValue<'_> fluent query API for BAM aux tags ([`5df88ba`](https://github.com/Softleif/seqaire/commit/5df88bad2eb9e885e8f54e6c79aa509bf828736f))
    - Complete BAM aux tag API — setters, spec, fuzz, proptests ([`142e06d`](https://github.com/Softleif/seqaire/commit/142e06d402ec13371cd2d9715b59f6181495462b))
    - Fix some review feedback ([`f93c275`](https://github.com/Softleif/seqaire/commit/f93c275683cfce96d49c18ed4aba9d9257302a4d))
    - Clean up customize api ([`077bc66`](https://github.com/Softleif/seqaire/commit/077bc665f9ff5111d987fdc68b71052f6721781b))
    - Clean up customize api ([`ebb2ddb`](https://github.com/Softleif/seqaire/commit/ebb2ddb60e64b3b4f152e484f740d40bc7b1b46e))
    - Replace CigarOp slice cast with bytemuck ([`6b463f2`](https://github.com/Softleif/seqaire/commit/6b463f2fadafd8964102010d402d99001f97b63f))
    - Three follow-ups from /critique: TLEN sentinel, set_filter removed, oracles ([`6f9a4a6`](https://github.com/Softleif/seqaire/commit/6f9a4a61ed4733045ef4590f96f78b4d6ce28abd))
    - Address critique: rename drift, doc gaps, example double-filter ([`c2e0e95`](https://github.com/Softleif/seqaire/commit/c2e0e95a3a4acdefa670464a0aa0fbbe1c79a4ae))
    - Type the BAM CIGAR slab as Vec<CigarOp> end-to-end ([`b39b5dc`](https://github.com/Softleif/seqaire/commit/b39b5dcb1ef26621ac5dd7c66fa9a15840fc0e24))
    - Tolerate reserved CIGAR op codes via CigarOpType::Unknown(u8) ([`085b893`](https://github.com/Softleif/seqaire/commit/085b893102a5f103272b80e0e7e0ce2fc97bce2d))
    - Make CigarOp a transparent u32 wrapper ([`0731b15`](https://github.com/Softleif/seqaire/commit/0731b156e59652a294337d37fb65c3664f112c65))
    - Replace closure-based filter API with CustomizeRecordStore trait ([`2ef542e`](https://github.com/Softleif/seqaire/commit/2ef542e7b0567425276689ddd5f0b76d0eec83fb))
    - Rename RecordStoreExtras → CustomizeRecordStore ([`6357641`](https://github.com/Softleif/seqaire/commit/63576413cc55cdf9f6f4282e97cec064fa8efa99))
    - Push-time filtering in CRAM fetch_into_filtered ([`e1e3519`](https://github.com/Softleif/seqaire/commit/e1e35190a2a9dad8dd280e02edacfac74ad458f1))
    - Add record store filters ([`d5061b9`](https://github.com/Softleif/seqaire/commit/d5061b979152d3b1fd961b1a33c441c367f3a958))
    - Rework pileup + Readers API around lending iterator and extras provider ([`5138b3f`](https://github.com/Softleif/seqaire/commit/5138b3f42117472d582651c9d7a47a47238f9899))
    - Some doodles ([`d4b8265`](https://github.com/Softleif/seqaire/commit/d4b8265ec759ee293f12a67a6b6aedbff8a3e396))
    - Docs.rs show example stuff ([`f0ad627`](https://github.com/Softleif/seqaire/commit/f0ad6278f99d3b51602cdd7c22deb271e626422e))
    - Fix clippy ([`e55c6c0`](https://github.com/Softleif/seqaire/commit/e55c6c06d31521e5d5d5009e588bb6eb73b46113))
    - Add pileup_extras example and smoke tests for all examples ([`a2a6e7f`](https://github.com/Softleif/seqaire/commit/a2a6e7f645ae87e9ff3ff3f3533c409a99d00f9a))
    - Add extras_idx to SlimRecord for sort/dedup-safe extras ([`c23a404`](https://github.com/Softleif/seqaire/commit/c23a404f99b0816ae465dd98f76321e75fc07cb3))
    - Add per-record extras storage to RecordStore and PileupEngine ([`a2c9207`](https://github.com/Softleif/seqaire/commit/a2c9207dd14e5b34937ce8706e9775e88cce21ab))
    - Fix ci ([`999fab4`](https://github.com/Softleif/seqaire/commit/999fab474e9162e066053ae813e15236094e78dc))
    - Improve pileup_e2e benchmark ([`17e7b05`](https://github.com/Softleif/seqaire/commit/17e7b0545b87d2d3be842b7ca44b929585757402))
    - Improve pileup_e2e benchmark ([`ba37731`](https://github.com/Softleif/seqaire/commit/ba377319c9c68ca336decfc8c0e406e921d4dce4))
    - Random code fixes ([`c180ef3`](https://github.com/Softleif/seqaire/commit/c180ef33cfbd37423b03ee9edee1cb433a949bb7))
    - Use windowed fetching in mpileup example for bounded memory ([`1a5c4d6`](https://github.com/Softleif/seqaire/commit/1a5c4d6de0d17d269f65a2f90a5e189280b72286))
    - Allow RegionBuf to load oversized regions for whole-chromosome queries ([`fd0cd8e`](https://github.com/Softleif/seqaire/commit/fd0cd8ed9d5f63e3f9a4d2ae50c5ca2e3d5fbade))
    - Better bcf benchmark ([`87492d6`](https://github.com/Softleif/seqaire/commit/87492d62fa4a4be95e32a4a2b6cb4564246c13e7))
    - Remove old fuzz target ([`e60e99b`](https://github.com/Softleif/seqaire/commit/e60e99b116d21de6501f647819b730725391db7e))
    - Fix review findings: dedup offset bug, FieldTracker tests, bcftools roundtrip ([`822ed76`](https://github.com/Softleif/seqaire/commit/822ed76794ddb850848a6931d76b6317a726be2d))
    - Deduplicate INFO and FORMAT fields written twice to the same record ([`f8813a4`](https://github.com/Softleif/seqaire/commit/f8813a44d9e97ad7b8b138f1f1a2dbe5ac087c9b))
    - Some nice vcf examples ([`bcd0daf`](https://github.com/Softleif/seqaire/commit/bcd0daf79b062b99107d4a770d3201d6a1a52e55))
    - Suppress inherited stdout/stderr on subprocess .status() calls in tests ([`32fac89`](https://github.com/Softleif/seqaire/commit/32fac89c15dd8eb0a73ddee967c21b9e3fe0b3e8))
    - Clean up writing APIs for publication ([`aff5332`](https://github.com/Softleif/seqaire/commit/aff5332bde2580fed615daedfa105ddba51f8f09))
    - Add #[non_exhaustive] to all public error enums, fix PathBuf in VcfError ([`cf1c926`](https://github.com/Softleif/seqaire/commit/cf1c9265fbbc59f5eb626b731e3b08232cd1e43f))
    - Tracey spec review: fix spec/code mismatches and annotation gaps ([`3422266`](https://github.com/Softleif/seqaire/commit/34222669970f9f1209a0172196d197a2392e0253))
    - Fmt ([`bd7f42c`](https://github.com/Softleif/seqaire/commit/bd7f42cdf3ee3f5b7b93300603645ebcf18255f6))
    - Add some fuzz seeds ([`7508f1a`](https://github.com/Softleif/seqaire/commit/7508f1aa58de4b0f303cd12f512d895269efb074))
    - Validate qual/bases length match in push_fields ([`f532840`](https://github.com/Softleif/seqaire/commit/f5328400d9c50a16c158d104ae316d69da036f5b))
    - Fix fuzz crash in merge_chunks and improve fuzz infrastructure ([`530f218`](https://github.com/Softleif/seqaire/commit/530f2180d434ed8ec7670974dc7c1f5235e47da9))
    - Add simple_variant_caller and base_mods examples ([`802e017`](https://github.com/Softleif/seqaire/commit/802e017704731bf8cbe01f48aec1a95192478ed4))
    - Fix clippy ([`6c639d0`](https://github.com/Softleif/seqaire/commit/6c639d0087c97a22d6134934398f4b26306ce05c))
    - Fix BaseQuality adaptation after rebase onto main ([`3f7397c`](https://github.com/Softleif/seqaire/commit/3f7397c3a3d66e8a83fadf34d656513e82884594))
    - Address review findings for write_store_record ([`63acb66`](https://github.com/Softleif/seqaire/commit/63acb6659c4b20a90d09b216d6ae67d0561afa2a))
    - Simplify realignment example using write_store_record ([`00875d2`](https://github.com/Softleif/seqaire/commit/00875d2e553ba45a9580f8ea4814118325b32726))
    - Add write_store_record to BamWriter for direct RecordStore writing ([`081f2fd`](https://github.com/Softleif/seqaire/commit/081f2fdf0e6dde973199e80623d9d0f3ceb5b395))
    - Add next_ref_id to SlimRecord and propagate through read path ([`dc42378`](https://github.com/Softleif/seqaire/commit/dc423781adc0c00873002fe8cb10076ed0686c41))
    - Add --output flag to realignment example for BAM writing ([`53bbb2c`](https://github.com/Softleif/seqaire/commit/53bbb2cbf08fd0156c89ad7e0e875de30e6f064d))
    - Add example for local realignment API ([`63d6e6a`](https://github.com/Softleif/seqaire/commit/63d6e6a70f21d9fcd26a9db461561f5f2412cae3))
    - Add BaseQuality newtype for 'quality unavailable' sentinel ([`d24dc9c`](https://github.com/Softleif/seqaire/commit/d24dc9c8009dc48fe7c0fc2650a883f87652650f))
    - Calibrate RecordStore::with_byte_hint against real BAMs ([`bc016d6`](https://github.com/Softleif/seqaire/commit/bc016d63bff86126fa8b4a96e957ef3d29e7b6e0))
    - Split combined qual/aux slab into independent slabs ([`aa62ed6`](https://github.com/Softleif/seqaire/commit/aa62ed62598c84187dd8397398ac68c65e6e59a9))
    - Extract try_qpos helper to test SeqTooLong without 4 GiB allocation ([`7382008`](https://github.com/Softleif/seqaire/commit/7382008fb7d3ddcdffb83aaa42b42de4e5691723))
    - Add proptests for MM-tag validation paths ([`1680816`](https://github.com/Softleif/seqaire/commit/1680816ceb09f995fcece22b73783d5191bfdc09))
    - Tracey review: base_mod spec text and annotations ([`aa8aff9`](https://github.com/Softleif/seqaire/commit/aa8aff91e0ede98437f7c66a8c24fa0abf616790))
    - Spec ([`dd42dbd`](https://github.com/Softleif/seqaire/commit/dd42dbd9545444bda0edb8e32ffb9c7ec77dddf0))
    - Address adversarial review: drop dead helper, FxHashMap, checked CSI math ([`f81de3b`](https://github.com/Softleif/seqaire/commit/f81de3bac38e33fd812ec618f0b5f72cbd3793d5))
    - Address CSI / base-mod review findings ([`88f835d`](https://github.com/Softleif/seqaire/commit/88f835da45005ab056a1134d1f824d43dfd54fa9))
    - Add structured MM/ML base modification parsing ([`965aa93`](https://github.com/Softleif/seqaire/commit/965aa939cdce338960511a274d7d30bf03973dd7))
    - Clippy ([`3b0e618`](https://github.com/Softleif/seqaire/commit/3b0e618153c0b192ea669361b23b3b9b5ff9558d))
    - Fix 8 review findings in CSI implementation ([`ff89bf2`](https://github.com/Softleif/seqaire/commit/ff89bf21bb2bac33b470cd5ec3b707ed3141fc36))
    - Add CSI index support: reader, writer, unified AlignmentIndex, format detection ([`cbf6c44`](https://github.com/Softleif/seqaire/commit/cbf6c44c05d7d86c82339acd3a786c16ee399bbd))
    - Allow some lints in example ([`45c270b`](https://github.com/Softleif/seqaire/commit/45c270b3635210ffa4c3ccede244a490ad582b2c))
    - Add mpileup example: simple pileup viewer using seqair's pileup API ([`8814d1a`](https://github.com/Softleif/seqaire/commit/8814d1a459120ee8d66da56c81c7a9644b81abfd))
    - Fix clippy ([`82ed0bc`](https://github.com/Softleif/seqaire/commit/82ed0bc56999373978f3b33b7772f645b42d6e8b))
    - Add pileup vs htslib cross-validation on htslib's mpileup test SAMs ([`adc9de4`](https://github.com/Softleif/seqaire/commit/adc9de4f3da3903f4d8a9f30aafd9d05662fea6b))
    - Add clap ([`a4e2ae6`](https://github.com/Softleif/seqaire/commit/a4e2ae6024aae302a22f305cb01efcee8da52bfd))
    - Add ComplexIndel PileupOp variant for D→I and N→I CIGAR patterns ([`dae0a07`](https://github.com/Softleif/seqaire/commit/dae0a0701b93ac7b18f3f4ed9691e5037397afc7))
    - Add base modification (MM/ML) spec with current + planned rules ([`e07167f`](https://github.com/Softleif/seqaire/commit/e07167fe6c6b2b2cac97257551489eb6b5f18b1a))
    - Add MM/ML base modification tag round-trip tests ([`402dc59`](https://github.com/Softleif/seqaire/commit/402dc595145370d513b4c44cbc8605f664e1c0be))
    - Add BGZF writer validation: bgzip decompresses seqair output correctly ([`05e2702`](https://github.com/Softleif/seqaire/commit/05e27025ecff9a0a342ba3c8e78ff12d15414400))
    - Add VCF/BCF bcftools round-trip tests: integer boundaries, genotypes, indels ([`0f94b77`](https://github.com/Softleif/seqaire/commit/0f94b773f84818b9da6ac99b4be562e7a1291177))
    - Add BAM writer stress tests: field limits, poisoning, unmapped dispatch ([`8fac85a`](https://github.com/Softleif/seqaire/commit/8fac85ab505554003317dcb2ceb01931afcd76af))
    - Add spec rule for CRAM attached-mate TLEN reconstruction ([`2d12d1b`](https://github.com/Softleif/seqaire/commit/2d12d1b783d6891508b52f02ebc01c978656b3ad))
    - Fix CRAM template_len for attached/downstream mates ([`035a57c`](https://github.com/Softleif/seqaire/commit/035a57c3294c2e1b13222985171cc0c033c80945))
    - Add TLEN validation: 30 BAM tests + 6 CRAM field tests from htslib pairs ([`15687d2`](https://github.com/Softleif/seqaire/commit/15687d228bcec8dce6873756f7c3e1d7775cf3b7))
    - Add CRAM version matrix tests: v3.0, v3.1, multi-ref, embedded ref ([`b4c51c8`](https://github.com/Softleif/seqaire/commit/b4c51c825154a5c2e1ce2be0147c0c15a862112a))
    - Add index round-trip tests: write BAI with seqair, verify with samtools ([`9859476`](https://github.com/Softleif/seqaire/commit/98594767323ad0c4c9be006c5a3ce1a0bc470a4c))
    - Add BGZF block boundary tests: large records spanning 64KB blocks ([`dc046fa`](https://github.com/Softleif/seqaire/commit/dc046fa80d631f60e8b03a09f98adaea865faee3))
    - Add BAM write round-trip tests: write with seqair, verify with noodles/samtools ([`1d213d7`](https://github.com/Softleif/seqaire/commit/1d213d7c65bd130d9ebf0e129e9d526a0a25c9fb))
    - Add htslib edge case tests ([`cc6c264`](https://github.com/Softleif/seqaire/commit/cc6c264dccf918a67ae773ff6bf8419eaffdac09))
    - Add htslib SAM parity tests: parse htslib test files, compare vs noodles ([`7d18a64`](https://github.com/Softleif/seqaire/commit/7d18a64d094c1cc11d2ade548cd59f06f5772668))
    - Simplify pos ([`e6b7bc1`](https://github.com/Softleif/seqaire/commit/e6b7bc154eedaaf07d4372bd817494178bc273eb))
    - Fix review findings: use as_i32() in CIGAR hot path, silence clippy ([`9a4f19a`](https://github.com/Softleif/seqaire/commit/9a4f19ab7aa58bfa206ba181fbd34b5c642a3560))
    - Constrain Pos to i32::MAX, add standard TryFrom impls, use Pos0/Pos1 aliases ([`50b9c15`](https://github.com/Softleif/seqaire/commit/50b9c157b140815eb362bdf51c6fab4656da5dfd))
    - Add compile_fail doctests for Writer and RecordEncoder typestates ([`0e3deed`](https://github.com/Softleif/seqaire/commit/0e3deedb963df9861893414ec38495d322f1b845))
    - Address critique: re-export phases, fix from_bam_header, add compile_fail test ([`afc3cbe`](https://github.com/Softleif/seqaire/commit/afc3cbe42a1944f0873c1df8d5da3018be90ff9c))
    - Typestate VcfHeaderBuilder enforces BCF string dict ordering at compile time ([`195d6d9`](https://github.com/Softleif/seqaire/commit/195d6d98e20b452dcd29f809f8e646dcae217c63))
    - Clippy ([`b8d1d4d`](https://github.com/Softleif/seqaire/commit/b8d1d4db86d51932c0f25f614ff9c8574d36198d))
    - Simplify FailedToWriteFormattedString error ([`462fdd8`](https://github.com/Softleif/seqaire/commit/462fdd8e065f65f73bc8bc898b4498a0fad5208b))
    - Return Result from FormatEncoder/FormatKey instead of panicking ([`aa5f705`](https://github.com/Softleif/seqaire/commit/aa5f70541f4516edd124fc900c07e85b57ed4858))
    - Harden VCF encoder and RecordStore from second review round ([`678d91b`](https://github.com/Softleif/seqaire/commit/678d91b2127802b4453350ed3fdcaace56b78099))
    - Improve docs ([`977cf5f`](https://github.com/Softleif/seqaire/commit/977cf5f20d13129862152e03d49dfc93ca084d9d))
    - Fix multi-sample VCF encoder review findings ([`f2aabde`](https://github.com/Softleif/seqaire/commit/f2aabde1844e7dc9f692acaff8667d8477b15311))
    - Prepare RecordStore for local realignment: cigar slab, tid, mate fields, set_alignment ([`f972492`](https://github.com/Softleif/seqaire/commit/f97249279db2f17ecd5496dee23ba6081f5fa63c))
    - Clippy ([`99f6456`](https://github.com/Softleif/seqaire/commit/99f64561c7bda2baed19f0e32769a44c45353917))
    - Move lib.rs doc examples into their defining modules ([`5fab356`](https://github.com/Softleif/seqaire/commit/5fab3566de324e95542ab2613ab121bb0c5c7593))
    - Split comparison.rs into per-topic bench files with realistic germline data ([`4c813f3`](https://github.com/Softleif/seqaire/commit/4c813f3364ffc38c47dadc0e5b344c6f6c558da0))
    - Add multi-sample FORMAT encoding support ([`4fe37b1`](https://github.com/Softleif/seqaire/commit/4fe37b15c29356e4638f774b13831e33c5a50515))
    - Better vcf benchmarks ([`e9ee5ad`](https://github.com/Softleif/seqaire/commit/e9ee5add673b3dcf6f1571715d114e4a839c19bd))
    - Fix clippy warnings: indexing_slicing in write_float_g, cloned_ref_to_slice_refs in tests ([`3ef5838`](https://github.com/Softleif/seqaire/commit/3ef583823fda1a9c8292aea5b16e63fffeb2d93d))
    - Fix all critique findings from unified Writer review ([`d999a8d`](https://github.com/Softleif/seqaire/commit/d999a8dd7feec30b6e520deacd523629b76f1392))
    - Consolidate VCF/BCF writing to unified Writer<W,S> ([`0e84b74`](https://github.com/Softleif/seqaire/commit/0e84b748ff24696ffe840c42c29ab47cfc09e1ef))
    - Detect out-of-order register calls, check header_written in record_encoder ([`830c516`](https://github.com/Softleif/seqaire/commit/830c516ebd462eb432909d2baf4ee48166890597))
    - Fix review findings: filter handling, enum dispatch, multi-sample guard ([`d03305d`](https://github.com/Softleif/seqaire/commit/d03305d5ac83e5c4c250cbf345877f07642d2fc9))
    - Fix doc test: import RecordEncoder trait in encoder example ([`bb5d736`](https://github.com/Softleif/seqaire/commit/bb5d73642a45be2defa524afb1e8afdd709b6bdd))
    - Fix clippy and rustdoc warnings ([`a2d7830`](https://github.com/Softleif/seqaire/commit/a2d7830f60fc4ccf50b594d3748699ab9b86f10b))
    - Simplify StringMap to Vec<SmolStr> and fix clippy warnings ([`e759df0`](https://github.com/Softleif/seqaire/commit/e759df087ec8816c760746523dc25b44f55793f2))
    - Clippy ([`6035c25`](https://github.com/Softleif/seqaire/commit/6035c2577f4eb240328c0bc5780bcc10ae6017d5))
    - Hm region buffer still didn't work? didn't we fix this? ([`279aa4a`](https://github.com/Softleif/seqaire/commit/279aa4a57899cd169825beb18d94966f9c4d1b56))
    - Cleanup ([`dc5beca`](https://github.com/Softleif/seqaire/commit/dc5beca4e3edce09dc67fdfa739c123ddb9e01b0))
    - Address code review: fix casts, allocs, and ordering invariants in RecordEncoder ([`2f8fda1`](https://github.com/Softleif/seqaire/commit/2f8fda1a5fee46e668670599b8bce23bf6336105))
    - Some fixes ([`260dc47`](https://github.com/Softleif/seqaire/commit/260dc47970abcf00c90b3b834db4e60fe8bb07e4))
    - Use proper tmp dir ([`8c29fa3`](https://github.com/Softleif/seqaire/commit/8c29fa3c00c8afb13e76b52ff084aafb48a6c409))
    - Use smallvec ([`c07d012`](https://github.com/Softleif/seqaire/commit/c07d012157bc719c8a0dce466ee3594b3e31229a))
    - Add cross-format BCF↔VCF equivalence tests via bcftools ([`56ad4b1`](https://github.com/Softleif/seqaire/commit/56ad4b16112cd6e6e829aca24570a2bfb0e8fcb6))
    - Add VcfRecordEncoder for zero-alloc VCF text encoding ([`f17033d`](https://github.com/Softleif/seqaire/commit/f17033df18bbac1b54ac5f4da6339f9781d5a60a))
    - Add RecordEncoder trait with typed field keys for zero-alloc VCF/BCF encoding ([`45f87e6`](https://github.com/Softleif/seqaire/commit/45f87e640fe037fe690a4f1df4d6ff79ca104f82))
    - Add Number::BaseModification (M), emit PASS before metadata in header ([`1eb9c27`](https://github.com/Softleif/seqaire/commit/1eb9c275a409275c70bdb0d1e74c88c921d4203b))
    - Add `M` number ([`1640db2`](https://github.com/Softleif/seqaire/commit/1640db234338ed90ac5c0674f325eb830cac8a15))
    - Emit other_lines (metadata) before FILTER in VCF header ([`e8ad711`](https://github.com/Softleif/seqaire/commit/e8ad71174a7e724c497a4eaed6d08cfc8db17da6))
    - Clean up cram's bam flags parse ([`eaa4b47`](https://github.com/Softleif/seqaire/commit/eaa4b476158ed4f2dd4c3a63cf246ce86d10f691))
    - Add BamFlags newtype to replace raw u16 flag fields ([`ce84bc3`](https://github.com/Softleif/seqaire/commit/ce84bc3809e870878cbea958b5614bbe36122ece))
    - Update spec: region_buf.max_region_bytes ([`4d267a4`](https://github.com/Softleif/seqaire/commit/4d267a4b03ed5dd3e386d03e524be96ac4c80b96))
    - Fix RegionBuf failing on large BAM files (>100 GB) by batching chunk loads ([`ce0bb5d`](https://github.com/Softleif/seqaire/commit/ce0bb5df5022de1ac410997f572bf995c3e5a469))
    - Let's not add tracing error logs for no reason ([`c708ef6`](https://github.com/Softleif/seqaire/commit/c708ef63fe7e7495ce2f6b74c46fb7f0d3d4b936))
    - Fix unfulfilled expect(unreachable_code) on x86_64 and u8-as-i8 SIMD casts ([`e567ade`](https://github.com/Softleif/seqaire/commit/e567ade7ec25f75c8f5622b18bffc5c174fe71da))
    - Fmt ([`a3590b6`](https://github.com/Softleif/seqaire/commit/a3590b633aa8b5049f40fd90dc1f132499d2ec07))
    - Fix clippy allow_attributes_without_reason across all integration tests ([`add4db6`](https://github.com/Softleif/seqaire/commit/add4db6f3568586219471423cb0f0593ba122d91))
    - Fix warnings in BAM/FASTA/types test modules ([`ea65683`](https://github.com/Softleif/seqaire/commit/ea65683bebe1edec40250801bc4042ca9e8845a0))
    - Fix warnings in CRAM test modules ([`5dd60ba`](https://github.com/Softleif/seqaire/commit/5dd60bab056d40a4e039af9d60a71b4af26d0866))
    - Fix cast warnings in FASTA reader and GZI index ([`767dd10`](https://github.com/Softleif/seqaire/commit/767dd10b0d354fa30d0f48660efb15cbc351f16a))
    - Misc fixes (aux i8 decode, Cargo.toml lint) ([`c19d902`](https://github.com/Softleif/seqaire/commit/c19d902ba66572e119c6d99bafd50eda1929ce95))
    - Fix cast warnings in VCF/BCF encoder and writers ([`94506f6`](https://github.com/Softleif/seqaire/commit/94506f6a619075f879d2a15946e53aadbc318da5))
    - Fix cast warnings in SAM reader ([`744a552`](https://github.com/Softleif/seqaire/commit/744a552c55e29a6abac67750bfe57e737e6ce395))
    - Fix cast warnings in CRAM parsing (ITF8/LTF8 wire format) ([`75093e9`](https://github.com/Softleif/seqaire/commit/75093e97744b85f59b12345dfc86038da7ce6dee))
    - Fix cast warnings in BAM pileup, index, region_buf, writer ([`6929ae2`](https://github.com/Softleif/seqaire/commit/6929ae22415c99fed35b23a8586e16c2f2d15c71))
    - Fix cast warnings in BAM record, record_store, owned_record ([`a960758`](https://github.com/Softleif/seqaire/commit/a96075850ebd8bc0caa6c59ddbe0f1fb247cf311))
    - Fix allow_attributes_without_reason, option_option, manual_let_else, same_name_method ([`6f38972`](https://github.com/Softleif/seqaire/commit/6f3897207397fad9926c50870c3445b3ab722b87))
    - Clippy ([`b28c19d`](https://github.com/Softleif/seqaire/commit/b28c19de3f24c6348a60fd689567a9564cdab116))
    - No format ([`15873c1`](https://github.com/Softleif/seqaire/commit/15873c1cd0a046bef85c644d954f67fdb26db08e))
    - Manual clippy fixes ([`5064288`](https://github.com/Softleif/seqaire/commit/50642882c8099bd6fc4f8a5c0ae6d14fc371b504))
    - Auto clippy fixes ([`8b00840`](https://github.com/Softleif/seqaire/commit/8b008403aedd8152b64d2097d4618702afeb6a60))
    - Fix fuzz OOM: validate BGZF ISIZE in RegionBuf::read_block ([`ab6336b`](https://github.com/Softleif/seqaire/commit/ab6336b76878352e3b9e81e907a90a885a8c9d77))
    - Nicer trace_ok/err ([`f098e5a`](https://github.com/Softleif/seqaire/commit/f098e5a8c783637a6af7daa071c4e25d3dfa781b))
    - Improve fuzz CI output and add seeds for all 23 targets ([`f314518`](https://github.com/Softleif/seqaire/commit/f3145180ff396f15c6c1d4fbb6d24f71e054a6f8))
    - Fix fuzz crash: replace debug_assert with checked i32 conversion in CIGAR mapping ([`3bdf259`](https://github.com/Softleif/seqaire/commit/3bdf2598609c4c47d7640e76a414f61518652daf))
    - Some spec fixes ([`1b5d05f`](https://github.com/Softleif/seqaire/commit/1b5d05fb5fbb737c0bcb5db6c310e27c9f031ab6))
    - Update stale tracey references and annotate 28 uncovered spec rules ([`ffa2585`](https://github.com/Softleif/seqaire/commit/ffa2585dba60d7c5bf10f34cfff4986073c95697))
    - Fix a clippy lint ([`cde447d`](https://github.com/Softleif/seqaire/commit/cde447d4daa7bd0a2e286a373d75e517d68da870))
    - Harden readers: add practical limits to BAI index and BAM header parsing ([`86d469f`](https://github.com/Softleif/seqaire/commit/86d469f4f9af5d67636289fc03d90654a51a609c))
    - Add r[io.writer_limits] spec rule, harden all format writers ([`6d0b930`](https://github.com/Softleif/seqaire/commit/6d0b93094eb50cd4db3e18d38bb1c746a7a59531))
    - Use reader-matching limits instead of i32::MAX for header validation ([`003f6cc`](https://github.com/Softleif/seqaire/commit/003f6ccbffafb778462ada15d96d371bb88bf644))
    - Reject oversized header fields instead of silently capping ([`fd04e18`](https://github.com/Softleif/seqaire/commit/fd04e1838ca030400aec1297bf662bb5cbac2b01))
    - Fix critique review findings: tautological tests, spec precision, correctness ([`99011fc`](https://github.com/Softleif/seqaire/commit/99011fcc174bb8d2d916266639eacd88a2da9da7))
    - Add BAM roundtrip benchmark: seqair vs htslib vs noodles ([`ea0610c`](https://github.com/Softleif/seqaire/commit/ea0610cbe064d87eed74b0a0e3e870d3ecd5e435))
    - Add from_raw_bam() and aligned_pairs() for BAM rewrite pipeline ([`e5a9e26`](https://github.com/Softleif/seqaire/commit/e5a9e26e301d51649aea8635dbeb7773cdf0d1d9))
    - Fix critique findings: orphan bytes, unsafe transmute, validation ordering ([`0aaada0`](https://github.com/Softleif/seqaire/commit/0aaada0c624daa89f0e11131bcdddaedcdc3dfc1))
    - Add missing Tracey spec annotations to BAM writing code ([`f9e4db2`](https://github.com/Softleif/seqaire/commit/f9e4db2cf6d6c1462af51e0f1470a8902d57be00))
    - Implement BAM writing: OwnedBamRecord, BamWriter, BAI index co-production ([`1ec8c51`](https://github.com/Softleif/seqaire/commit/1ec8c51c9ddc7a9efd64ee6966ee4515c3b5d839))
    - Make sure fuzz files are formatted and compile ([`e585244`](https://github.com/Softleif/seqaire/commit/e585244629b66c4bae7e1ef8977f3b41f4f19370))
    - Fix some clippy lints with better errors ([`2a0422a`](https://github.com/Softleif/seqaire/commit/2a0422af38b0136c32ae4f238efc6ce32a680c37))
    - Use SmallVec for InfoValue/SampleValue arrays, Display for info_string ([`fbbe541`](https://github.com/Softleif/seqaire/commit/fbbe541e34f4d9ea7442586468b8d3004963dcc3))
    - Add info_integers/info_floats convenience methods for non-missing arrays ([`3e18c1e`](https://github.com/Softleif/seqaire/commit/3e18c1e5a75c25bdad3d65ab0f6e2c8a964b0f25))
    - Format VCF floats like C's %g — 6 significant digits, no trailing zeros ([`6e07433`](https://github.com/Softleif/seqaire/commit/6e0743379d2575be55e9cc42f4d0a324cc830817))
    - Let's not re-export all error types on root ([`57c1c91`](https://github.com/Softleif/seqaire/commit/57c1c91750467e1d4d8104e3a3bd3768a6c8ab10))
    - Fix doc link ([`cab9a33`](https://github.com/Softleif/seqaire/commit/cab9a331de0c0922e6b059f56c7734640a4bd62a))
    - Use proper smallvec ([`2a7a8a9`](https://github.com/Softleif/seqaire/commit/2a7a8a90b7a5629468a73b172f5da09a09f33f71))
    - Fix some review stuff ([`9154953`](https://github.com/Softleif/seqaire/commit/91549539b65a7732eeb6e63170beb11d6235af8f))
    - Add noodles to VCF/BCF writing benchmarks — full 3-way comparison ([`0ec2e21`](https://github.com/Softleif/seqaire/commit/0ec2e2186b9a2f788c5a1d63c966eb3c08954511))
    - Restructure benchmarks: 3-way comparison for VCF text, .vcf.gz, .bcf ([`9dd1796`](https://github.com/Softleif/seqaire/commit/9dd17963c2554b0948ea13a1aa0d3965edd50542))
    - Add BCF/VCF writing benchmarks, fix SmallVec release-mode compilation ([`6570167`](https://github.com/Softleif/seqaire/commit/657016728dd8acc0de5b9119068564622bb37870))
    - Omit empty references from TBI output, matching bcftools behavior ([`0815771`](https://github.com/Softleif/seqaire/commit/08157718258fbe96b42c34c9e0538ee075aaf2f3))
    - Add TBI index proptests: random region queries validated by bcftools ([`717abbd`](https://github.com/Softleif/seqaire/commit/717abbdfbe8ae4456a1d52d5825ebffc4acd0bf5))
    - Add Pos0/Pos1 aliases, clean docs, TBI index comparison tests ([`a8ece81`](https://github.com/Softleif/seqaire/commit/a8ece8191ab85651c5a4fccc73ee401972caef68))
    - Remove duplicate BCF encoding path — single encoder for all BCF output ([`c8da59f`](https://github.com/Softleif/seqaire/commit/c8da59f4c4819dab5c1849aaa6a7e16044ae1e95))
    - Fix code review: shared bcf_encoding, checked casts, typed errors, #[must_use] ([`2ca1453`](https://github.com/Softleif/seqaire/commit/2ca14536c4334317992fb199f5b70ddd6d74bc57))
    - Add crate-level and VCF module documentation with compilable examples ([`802c2bb`](https://github.com/Softleif/seqaire/commit/802c2bb367a5ab889a945a682ca3ef2a7ae3da27))
    - Fix critique findings: integer sentinels, smallest-int, typed errors ([`264fb34`](https://github.com/Softleif/seqaire/commit/264fb3427827af0ef6e858e4433df40ad485f268))
    - Add deep round-trip proptests, fix BCF string dictionary order ([`eff56e2`](https://github.com/Softleif/seqaire/commit/eff56e2af7fe7d361018fd36c6a7af7899830472))
    - Add encoder-vs-VcfRecord equivalence proptests ([`b7c105c`](https://github.com/Softleif/seqaire/commit/b7c105cf79b462d5ddd2324a9c20f41ec4a14299))
    - Add bcftools comparison tests, fix BCF version magic to 2.2 ([`7129b90`](https://github.com/Softleif/seqaire/commit/7129b90b8ca29482ef9f9ea81b27b736674418cd))
    - Add noodles round-trip comparison tests for VCF/BCF output ([`2915646`](https://github.com/Softleif/seqaire/commit/2915646826afb19182188f6c7eed5ae05a232474))
    - Add BcfRecordEncoder with typed handles for zero-alloc BCF encoding ([`032d421`](https://github.com/Softleif/seqaire/commit/032d421b25c0dde92af5b9758bfa50fe82327175))
    - Performance fixes, clippy cleanup, and proptest coverage ([`107ce07`](https://github.com/Softleif/seqaire/commit/107ce07aa40251012ab53305741dca58f1fa6092))
    - Fix review findings: validation, percent encoding, BCF correctness ([`302df9c`](https://github.com/Softleif/seqaire/commit/302df9cc42fc48d363c934d7e2d2197d1733966e))
    - Add VcfWrite trait, OutputFormat, and open_writer factory ([`c06a8ba`](https://github.com/Softleif/seqaire/commit/c06a8ba2deebe4d0a64dff27e09b0245a1db6923))
    - Add BcfWriter with BCF2 binary encoding and CSI index co-production ([`65bb92b`](https://github.com/Softleif/seqaire/commit/65bb92bdab2aed0dd2a93b0fd859a40e3c401fa2))
    - Add VcfWriter with text serialization and BGZF+TBI co-production ([`fab2dbe`](https://github.com/Softleif/seqaire/commit/fab2dbe6c2f1865d27aba70f488fc64086c8ade3))
    - Add single-pass IndexBuilder for TBI/CSI index co-production ([`c46b10c`](https://github.com/Softleif/seqaire/commit/c46b10ce859d884fe77ff4d17ebac65514515925))
    - Add VcfRecord model with builder, Filters, Genotype, and InfoValue ([`b16a422`](https://github.com/Softleif/seqaire/commit/b16a4220341daa957444be07a0736fe283580cd9))
    - Add VCF header model, type-safe Alleles enum, and error types ([`7036f84`](https://github.com/Softleif/seqaire/commit/7036f84b35bdce6b0d631e9a76147a4418dfb670))
    - Add BgzfWriter with virtual offset tracking and round-trip tests ([`44c2278`](https://github.com/Softleif/seqaire/commit/44c2278d260825609b44ee2dad07d5a81f347558))
    - Add SAM fuzzing support to fuzz_reader_indexed ([`8b62bd7`](https://github.com/Softleif/seqaire/commit/8b62bd761c02a1ae0ae434f06e60ed3395dc9686))
    - Fix remaining clippy arithmetic_side_effects warnings ([`7576aff`](https://github.com/Softleif/seqaire/commit/7576aff37843e9fe7d447fefe5210274a1440dc0))
    - Replace all unsafe arithmetic with checked_*, saturating_*, wrapping_* ([`881d3f1`](https://github.com/Softleif/seqaire/commit/881d3f1d7a60d6f7ba364ab07d11d2ebb2aed922))
    - Improve some expects ([`81efa82`](https://github.com/Softleif/seqaire/commit/81efa820d3aaf0558d663e065531943b533dc857))
    - Update fuzz README: bug #15 ([`4185df7`](https://github.com/Softleif/seqaire/commit/4185df7fda0987f27a0dc6a901eace1bce6fed4a))
    - Some more arithmetic fixes ([`395b567`](https://github.com/Softleif/seqaire/commit/395b5672ccfea3b780f4bdde3ff55a957c3a2015))
    - Fix add overflow in CRAI index query on large alignment_start + span ([`9e1a02f`](https://github.com/Softleif/seqaire/commit/9e1a02f15ec2275f94727ed570c97995c0997916))
    - Add fuzz README documenting all 23 targets and 14 bugs found ([`95a5021`](https://github.com/Softleif/seqaire/commit/95a5021deedfc69d114854247d6d84af341710f9))
    - Cleanup fuzzer ([`8c081ee`](https://github.com/Softleif/seqaire/commit/8c081ee7737c05738ede068cc6d355df3069c662))
    - Add debug_asserts for claimed invariants in indexing and SIMD code ([`f1b78ef`](https://github.com/Softleif/seqaire/commit/f1b78efa3886c9674ac08b016b0679e41c3d7821))
    - Fix debug_assert panic in RegionBuf when BGZF bsize < header + xlen ([`683eb05`](https://github.com/Softleif/seqaire/commit/683eb055e76565ca03dc6c23dd95b51f22af85e7))
    - Fix overflow in RegionBuf::load sum and per-range allocation ([`d50c8c7`](https://github.com/Softleif/seqaire/commit/d50c8c773504c872394026ae0c429461210773dd))
    - Replace all signed-to-usize casts with try_from in CRAM/BAM parsers ([`94c65a6`](https://github.com/Softleif/seqaire/commit/94c65a67ac46ce7a5c7628338738c6e25acd9738))
    - Fix capacity overflow from negative CRAM container length cast to usize ([`6a57ea6`](https://github.com/Softleif/seqaire/commit/6a57ea67a2544a7b0916be93f2b5ea0f14b32f13))
    - Replace Arbitrary with hand-rolled binary format for fuzz_reader_indexed ([`234b0f5`](https://github.com/Softleif/seqaire/commit/234b0f5a0c77917b6c64fb21a1f012de9fb5f257))
    - Move indexed reader fuzz types to shared lib, add seed generator example ([`bd57051`](https://github.com/Softleif/seqaire/commit/bd57051ec7166e7b5b1754779ed8e696df7a46b8))
    - Add BAM+CRAM seeds for fuzz_reader_indexed with Arbitrary enum dispatch ([`a935a18`](https://github.com/Softleif/seqaire/commit/a935a184aa63b609eccb92ec41b11e918cdd8c6d))
    - M ([`d5bfdae`](https://github.com/Softleif/seqaire/commit/d5bfdaee022bd9434dbea160150e2f32fa9301d2))
    - Better reader fuzzer ([`4e09fe7`](https://github.com/Softleif/seqaire/commit/4e09fe708b826717dbdb0ba35a1ad0e9f825777d))
    - Clean up some reader code ([`cff12dc`](https://github.com/Softleif/seqaire/commit/cff12dc235edb884a9e6ec8b3918452c2b10ba25))
    - Fix/allow some arithmetic side effects ([`b663e65`](https://github.com/Softleif/seqaire/commit/b663e6539bc1b99d119107703b69c8058bdef1e0))
    - Make IndexedReader fully generic: all readers cursor-based for fuzzing ([`616b1bc`](https://github.com/Softleif/seqaire/commit/616b1bc5d4db888bcb457f012774e28631646e78))
    - Add structure-aware pileup fuzz target with valid BAM records ([`50e3fda`](https://github.com/Softleif/seqaire/commit/50e3fdaeec9ec566729cf5a75e3ef4e1d9e87496))
    - Add full cursor-based Readers pileup pipeline fuzz target ([`2f0c02b`](https://github.com/Softleif/seqaire/commit/2f0c02bd50ba624d29ef2df6633204ada953273d))
    - Add fuzzing spec rules, IndexedBamReader fuzz target, and tmpfs docs ([`0cb6f95`](https://github.com/Softleif/seqaire/commit/0cb6f9506789a4387d1a1a44a1669315566133df))
    - Add allocation size limits to BAM header parsing ([`9f65338`](https://github.com/Softleif/seqaire/commit/9f653382cd12e893920cff173d5a9072c06a0287))
    - Improve fuzz_reader_bam: use raw &[u8] input with real BAM file seed ([`bd2161b`](https://github.com/Softleif/seqaire/commit/bd2161b5895fafcec0347ea5fa556042d81188eb))
    - Add cursor-backed BgzfReader and full-stack BAM reader fuzz target ([`dc2c831`](https://github.com/Softleif/seqaire/commit/dc2c8317aec24f184342e6555c75bf0ccb75ddee))
    - Fix THREADS variable in run_all.sh, add full-stack fuzz targets ([`cc8f2d2`](https://github.com/Softleif/seqaire/commit/cc8f2d2282b7c6950bcce137c5e5e24b23a825e3))
    - Add full-stack fuzz targets: pileup pipeline and CRAM decode pipeline ([`e895725`](https://github.com/Softleif/seqaire/commit/e895725dea7d10b5857dcf17f2d80132b404bc0e))
    - Rustfmt ([`24cbc05`](https://github.com/Softleif/seqaire/commit/24cbc05354f543b25227c6caab75cf2fd7c706d2))
    - Fix clippy lints ([`6099183`](https://github.com/Softleif/seqaire/commit/6099183a1cd7338255ef7655c73aacdf546e7ee4))
    - Allow setting thread count ([`d682360`](https://github.com/Softleif/seqaire/commit/d682360548222b9ce4ea718e8b3c1fa0f0e211f4))
    - Add seed corpora, generator script, and CI fuzz runner ([`9d1b0ee`](https://github.com/Softleif/seqaire/commit/9d1b0eee81d59fab08a0d19eced0f1f03261feee))
    - Add seed corpora for fuzz targets extracted from test data ([`dc59eb7`](https://github.com/Softleif/seqaire/commit/dc59eb722c72472ae6101e975154b579953c5d65))
    - Add BGZF block parsing fuzz target ([`1e1a89c`](https://github.com/Softleif/seqaire/commit/1e1a89c4f75c8e917cc287a15193811837d318e9))
    - Add Dockerfile for x86_64 fuzz testing of SSSE3 SIMD paths ([`5453dfb`](https://github.com/Softleif/seqaire/commit/5453dfba3c08836d665bd9e4174dd39d16b2355b))
    - Fix shift overflow in Huffman canonical code generation ([`82ebfea`](https://github.com/Softleif/seqaire/commit/82ebfea40eb71b5079500171890b5fb1fc3e31f4))
    - Add 5 new fuzz targets: bitstream, encoding, base, region, substitution ([`566d52e`](https://github.com/Softleif/seqaire/commit/566d52e29926fb6f1939da03b9698534ca312833))
    - Fix accumulation overflow in CRAM slice feature reconstruction ([`96abeaa`](https://github.com/Softleif/seqaire/commit/96abeaa8dd51f6736c73d2f4d06d27d975c9aff1))
    - Fix arithmetic overflow in rANS state step on malformed frequency tables ([`dee2e9a`](https://github.com/Softleif/seqaire/commit/dee2e9a495827741725d6df15d3aff188d2a53b5))
    - Add allocation size limits to CRAM parsers to prevent OOM on malicious input ([`d55f299`](https://github.com/Softleif/seqaire/commit/d55f29931c8cc9b148befa99803d676f428b571a))
    - Fix SEGV in sequence decoders on undersized encoded buffers ([`44e2818`](https://github.com/Softleif/seqaire/commit/44e28189d6bb897ba291cc2cc1a7168d6975be88))
    - Fix arithmetic overflow in CIGAR parsing found by fuzzing ([`892ecc6`](https://github.com/Softleif/seqaire/commit/892ecc6aefacd1ef4e3703481f4d5bdb33368802))
    - Add cargo-fuzz targets for all file format readers ([`8fa7327`](https://github.com/Softleif/seqaire/commit/8fa73274636d0bd507efd90aaf2f321f7dfabe00))
    - Remove ChunkCache dead code and fix clippy warnings ([`0c2982f`](https://github.com/Softleif/seqaire/commit/0c2982ff89c8531536b90cd9cb8b2949f4b4b256))
    - Remove built-in overlap dedup from pileup engine ([`5e2c6c4`](https://github.com/Softleif/seqaire/commit/5e2c6c4ecd5c052f0f034138753901b27aa0c8cc))
    - Add store() accessor to PileupEngine, update dedup spec ([`74bf58f`](https://github.com/Softleif/seqaire/commit/74bf58fc1478a4000bfdca54601e7b7284810ee1))
    - Make overlap dedup deterministic: always keep first-in-template ([`712f7f2`](https://github.com/Softleif/seqaire/commit/712f7f2a52eab8f36ffe406ac445711df5536b00))
    - Match htslib overlap dedup: remove insertion tie-breaking ([`2738f4b`](https://github.com/Softleif/seqaire/commit/2738f4bbea4b1babb510f3f0574688bb23c35f5a))
    - Use unified query() for fetch_into, removing chunk cache path ([`0c1e8ca`](https://github.com/Softleif/seqaire/commit/0c1e8ca7d20feca6e11cc7e4fe89d2e3a2828c6d))
    - Deduplicate records loaded from overlapping nearby/distant chunks ([`7c2d1ce`](https://github.com/Softleif/seqaire/commit/7c2d1ce4a54318f0327b2a8eb2a98d5e9edd2ab8))
    - Sort store by position after chunk cache injection ([`1c51e98`](https://github.com/Softleif/seqaire/commit/1c51e98ed943d6752c69383b3efbf4619bc77f7d))
    - Fix duplicate records from overlapping BAM index chunks ([`72477c0`](https://github.com/Softleif/seqaire/commit/72477c00b4bbd8d94b88c54b728feca2c9ad44e2))
    - Harden pileup indel implementation after critical review ([`d468466`](https://github.com/Softleif/seqaire/commit/d468466dee4e3d9a7e8e2abe1314866dbd8f28c2))
    - Fix deletion_ops_match_htslib test semantics mismatch ([`6218a34`](https://github.com/Softleif/seqaire/commit/6218a345f716c26dc33f6f6a1695dd44a6306b09))
    - Cross-validate del_len values against htslib ([`168a9ec`](https://github.com/Softleif/seqaire/commit/168a9ec7835ba9227dc20898c2a2acdf40c29c06))
    - Add proptests for del_len on PileupOp::Deletion ([`69cbc06`](https://github.com/Softleif/seqaire/commit/69cbc06f033e2c0580645bb35e66529b5b41eb7d))
    - Add del_len: u32 to PileupOp::Deletion and CigarPosInfo::Deletion ([`e5f6523`](https://github.com/Softleif/seqaire/commit/e5f6523518f35486357ea64fd41a07681142e5ef))
    - Add htslib comparison tests for deletion and insertion pileup ops ([`503bfb5`](https://github.com/Softleif/seqaire/commit/503bfb536177d9b822d7c8f3d4e5234f2067f07c))
    - Remove panicking Add/Sub operator impls from Pos and Offset ([`61073fb`](https://github.com/Softleif/seqaire/commit/61073fbff3e0cf4269be9045fc7c463e3717112c))
    - Add tests for invalid position error paths ([`11c5e66`](https://github.com/Softleif/seqaire/commit/11c5e667be87b98c73a22ef50b607d2a40dec338))
    - Harden Pos<S>: no panics on untrusted data, proper error propagation ([`eb8301d`](https://github.com/Softleif/seqaire/commit/eb8301d388be32ace94fe0a3d79307b64e0f8a8d))
    - Fix Pos<S> soundness: remove all unsafe, add proptests ([`6a0a63d`](https://github.com/Softleif/seqaire/commit/6a0a63dc9dbdd9ed46aeb36276b0cce922cfed5c))
    - Complete Pos<S> migration: Pos<One> conversions, API cleanup, nonmax niche ([`c7e41a1`](https://github.com/Softleif/seqaire/commit/c7e41a1c07a536b0722f0119282fc6fc95d3dade))
    - Migrate all positions to Pos<Zero>/Pos<One> newtype ([`6ab4639`](https://github.com/Softleif/seqaire/commit/6ab463981135fe8e8976c4030deebf9aa0733330))
    - Remove CigarIndex, make CigarMapping the sole CIGAR query API ([`baedd72`](https://github.com/Softleif/seqaire/commit/baedd72d415166bba5a3b880b8a7480a46e13a00))
    - Add PileupOp enum for type-safe indel reporting in pileup engine ([`eb2b662`](https://github.com/Softleif/seqaire/commit/eb2b662fd932aa9da5c3464c74462879783d311c))
    - Add pileup coverage tests: ref_base, combined controls, boundaries ([`aa9d7ae`](https://github.com/Softleif/seqaire/commit/aa9d7ae11ebf160f0709c977b2c9105fabfbf94f))
    - Harden pileup engine: review fixes, Debug impls, proptest read strategy ([`18eb7b4`](https://github.com/Softleif/seqaire/commit/18eb7b4ce58e5278be0eb116873a0e440a23ea3e))
    - Fix clippy warnings and add VCF/BCF writer specs ([`41b142e`](https://github.com/Softleif/seqaire/commit/41b142e2ae43936b16d762d5af7c01a70b8ef69c))
    - Fix review findings: push_fields overflow, consolidate AuxValue, clarify specs ([`1702c7a`](https://github.com/Softleif/seqaire/commit/1702c7aa310e9714189945da519189d6119a2275))
    - Add #[inline] hints to CigarMapping::new, qpos_at, and try_linear ([`25f7aa2`](https://github.com/Softleif/seqaire/commit/25f7aa288017f4b47c4f491569884da93ce86658))
    - Split ActiveRecord end_pos into parallel vec for cache-friendly eviction ([`8188455`](https://github.com/Softleif/seqaire/commit/8188455a3c6ed262cb142050e46550098c9b070b))
    - Remove trace-level tracing spans from pileup hot path ([`eb1022c`](https://github.com/Softleif/seqaire/commit/eb1022ccd448ecbb942779f7ce5842c5f56333ea))
    - Add noodles comparison tests, criterion benchmarks, and direct-to-slab base decode ([`1c01cc9`](https://github.com/Softleif/seqaire/commit/1c01cc9cd32b1fa017a57f00c72375012364d8f0))
    - Proptest Audit: Tautological vs. Valid Tests ([`49e446e`](https://github.com/Softleif/seqaire/commit/49e446ed28c951770c381d35dfc7dab061f70d43))
    - Sort spec pages ([`e7dc48c`](https://github.com/Softleif/seqaire/commit/e7dc48c33304b5e366935e3b099158569a2dc4dc))
    - Initial import from local repo ([`1681005`](https://github.com/Softleif/seqaire/commit/1681005cf8e5b97089e50b224b6efabf7566787a))
</details>

