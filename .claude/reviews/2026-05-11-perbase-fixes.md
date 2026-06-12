Now I have sufficient context. Let me compose the critique.

**7 findings: 0 showstoppers, 3 gaps, 1 inconsistency, 2 underspecified, 1 suggestion.**

---

[Gap] [1] `FilterRawFields.end_pos` is a semantic trap in the BAM `push_raw` path

The `FilterRawFields` struct exposes `end_pos` as a public field, but in the `push_raw` (BAM) path it's always set to `pos` because the CIGAR hasn't been decoded yet. In `push_fields` (SAM/CRAM), it carries the proper computed value. The docs on `FilterRawFields` mention this, but the field is visibly present and typed `Pos0` — a user writing `filter_raw` will naturally try `f.end_pos > some_threshold` and get silently wrong results on BAM input.

**Scenario:** User implements `RejectUnmapped`-style `filter_raw` but adds an `end_pos` filter: `f.end_pos > start_of_region`. Works perfectly in tests against SAM/CRAM. Deploys with BAM files, gets zero retained records because `end_pos == pos` for every read, and the filter rejects all of them. The root cause is invisible — no error, no panic, just silently empty results.

Sources: `src/bam/record_store.rs:348-396` (push_raw), `src/bam/record_store.rs:662-725` (FilterRawFields docs), `src/bam/record_store.rs:376-378` (end_pos set to h.pos)

[Gap] [2] `RejectUnmapped` duplicated across 9 test files with no public export

The most common customization — "drop unmapped reads like htslib does" — requires implementing `CustomizeRecordStore` by hand. The pattern is identical in 9 test files (5 lines each), plus examples. Yet there's no `pub struct RejectUnmapped` in the library. Users coming from `Readers::open()` who get unmapped reads in their pileup must discover the `CustomizeRecordStore` mechanism, understand `filter_raw`, and write the struct themselves.

**Scenario:** User opens a BAM file with `Readers::open()`, runs the pileup, sees unmapped reads polluting columns at depth 1. They check the spec (`r[bam.reader.unmapped_skipped]`), see "filter_raw is the single gate", look for `RejectUnmapped` or a comparable built-in type, and don't find one. They then need to read the `CustomizeRecordStore` trait docs, understand `FilterRawFields`, write 5 lines of boilerplate, and switch from `Readers::open` to `Readers::open_customized`.

Sources: `tests/compare_bam_with_htslib.rs:16-29`, same in 8 other test files, `src/bam/record_store.rs:726-836`

[Gap] [3] `CustomizeRecordStore::compute` runs before `filter`, making expensive computes wasted on filtered records

The sequencing is `filter_raw → compute → filter`. `filter_raw` is documented as "cheap checks", `filter` as "post-compute validation". But if a user puts anything non-trivial in `compute` (aux parsing, base counting, read group extraction) and then has `filter` reject some records (e.g., based on aux tags, CIGAR complexity, read group membership), the `compute` work is wasted. There's no post-filter compute hook.

**Scenario:** User implements a customizer that extracts read group from aux in `compute` and filters by read group in `filter`. Every record gets its aux parsed (expensive, allocates), then 80% get rolled back because they're the wrong read group. The `filter` logic that determines rejection could have been folded into `filter_raw` only if the user also re-implements aux parsing there (no aux access in `filter_raw` because the aux slab hasn't been extended yet). The user must either accept wasted work, or do aux parsing in both `filter_raw` and `compute`, which is a worse pattern.

Sources: `src/bam/record_store.rs:489-493` (compute before filter), `src/bam/record_store.rs:494-503` (filter after compute), `src/bam/record_store.rs:727-767` (trait docs showing the sequencing)

[Inconsistency] [4] `FilterRawFields` has 5 path-dependent fields with no type-level distinction

The struct carries `end_pos`, `matching_bases`, `indel_bases` that are accurate only from `push_fields`; and `raw_cigar_bytes`, `packed_seq`, `cigar_ops`, `bases` that are mutually exclusive between paths. A customizer that checks `f.cigar_ops.is_some()` to branch on the path is guessing — there's no enum or type-level guarantee. The trait docs describe what's available per-path in long prose, but the type system doesn't enforce it. If a future reader change adds a third path (e.g., streaming decode), the invariants get more complex.

**Scenario:** User writes `match f.cigar_ops { Some(ops) => { ... }, None => { ... } }` in `filter_raw`, assuming the `None` branch means BAM and trying to use `f.raw_cigar_bytes`. Works today. In six months, a new code path leaves both as `None` (e.g., header-only records), and the user's branch fires incorrectly.

Sources: `src/bam/record_store.rs:662-725` (FilterRawFields definition)

[Underspecified] [5] `column.store()` exposes raw `RecordStore<U>` access — unclear boundary between `AlignmentView` and direct access

`AlignmentView` offers `extra()`, `qname()`, `aux()`, `store()`, and `alignment()`. But `<PileupColumn as store()>` returns the full `&RecordStore<U>`, enabling `column.store().seq(aln.record_idx())`, `column.store().qual(aln.record_idx())`, `column.store().cigar(aln.record_idx())`, and `column.store().record(aln.record_idx()).seq(store)` (with Result). `AlignmentView` doesn't expose `seq()` or `qual()` — those must be fetched through `store()`. Is this intentional? Users need to know whether `AlignmentView` is the blessed API or just a convenience wrapper.

**Scenario:** User wants to access sequence bases at a pileup position. They see `AlignmentView` has `extra()`, `qname()`, `aux()` but not `seq()`. They search docs, discover `column.store()`, try `store.seq(aln.record_idx())`, then `store.seq_at(aln.record_idx(), qpos)`. This works but feels like they've gone around the public API. They wonder if `seq_at` is stable or internal.

Sources: `src/bam/pileup.rs:308-353` (AlignmentView), `src/bam/pileup.rs:226-229` (PileupColumn::store), `src/bam/record_store.rs:955-970` (seq, seq_at)

[Underspecified] [6] `PileupGuard::take_store()` is accessible through `Deref` — documented footgun with no compile-time prevention

`PileupGuard` derefs to `PileupEngine`, so `guard.take_store()` is reachable. Calling it disables store recovery — the guard's `Drop` finds an empty engine and leaves the slot as a default (empty, zero-capacity) `RecordStore`, causing the next `pileup()` call to allocate fresh ~39 MB. This is documented and tested, but a footgun that `#[deprecated]` or a clippy lint could prevent at compile time. Users calling `guard.store()` for inspect-only access are one dot away from `guard.take_store()`.

**Scenario:** User iterates pileup columns, wants to inspect the store mid-iteration, types `pileup.take_store()` by muscle memory (this was the API before `PileupGuard`), sees no error. The code compiles and runs. On the next segment, pileup allocates a fresh 39 MB store. The performance regression is silent — the old code iterated just fine, but the allocation overhead doubles. Only visible in profiles or memory monitors.

Sources: `src/bam/pileup.rs:571-671` (PileupGuard), `src/reader/readers.rs:356` (mem::take), `tests/readers.rs:506-534` (test pinning this behavior)

[Suggestion] [7] Consider a `PileupColumnLike` trait for `PileupColumn` / `PileupColumnPinned`

Both types have identical API surfaces (pos, depth, reference_base, alignments, raw_alignments, match_depth, store) but are separate types. This prevents writing functions generic over "any pileup column." For a library this focused, it's acceptable — but if more column types appear (e.g., async-compatible), a trait would reduce API surface duplication. Not blocking.

Sources: `src/bam/pileup.rs:178-240` (PileupColumn), `src/bam/pileup.rs:252-301` (PileupColumnPinned)

---

**What looks good:**

- `PileupGuard` RAII store recovery: correct on `break`, `?`, early return, and normal scope exit. Well-tested.
- `filter_raw` rejecting before slab extension: zero wasted memcpy or base decode for rejected reads. Right design.
- `pileups_into(&mut buf)`: caller-owned buffer reuse avoids per-column Vec allocation. Benchmarked alongside the lending iterator.
- `SlimRecord::end_pos_htslib()`: cleanly separates htslib-compatible unmapped-read span from the pileup engine's own `end_pos`, without breaking existing code.
- `AlignmentView`'s `Deref<Target=PileupAlignment>`: callers who don't need slab data write `aln.base()` and get the same API as before.
