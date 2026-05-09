# Plan: pileup column summaries and custom aggregation

## Why this exists

This is a high-level but implementation-oriented plan for making seqair's pileup engine more directly useful to downstream tools such as perbase.

perbase currently uses pileup columns as an intermediate representation. For every output row it walks the pileup alignments and recomputes command-specific summary state:

- output depth,
- A/C/G/T/N counts,
- insertion count,
- deletion count,
- reference-skip count,
- failed-read count,
- optionally mate-overlap resolution.

seqair already walks the active records for each reference position to construct each `PileupAlignment`. The opportunity is to let seqair either:

1. compute common column summaries while doing that work, or
2. let callers install a custom per-column accumulator and avoid both a downstream second pass and, for simple use cases, the internal per-column alignment `Vec` materialization.

This is the same broad motivation as the older noodles pileup draft: make the pileup engine provide richer per-column information directly, instead of requiring every downstream consumer to re-analyze the stack of reads.

## Current seqair shape

Relevant current code is in `crates/seqair/src/bam/pileup.rs`.

Today, the public API is:

```rust
pub fn pileups(&mut self) -> Option<PileupColumn<'_, U>>
```

`PileupColumn` internally owns a per-column alignment vector and exposes iterators over it:

```rust
pub struct PileupColumn<'store, U = ()> {
    pos: Pos0,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    store: &'store RecordStore<U>,
}
```

Important wording: seqair does **not** publicly return a `Vec<PileupAlignment>` directly. It returns a `PileupColumn`. The `Vec` is an internal implementation detail that matters for performance/API design.

The hot path currently does roughly:

```rust
let mut alignments = Vec::with_capacity(self.active.len());

for active in &self.active {
    let Some(info) = active.cigar.pos_info_at(pos) else { continue };
    let op = match info { ... };
    alignments.push(PileupAlignment { ... });
}

if let Some(max) = self.max_depth {
    alignments.truncate(max as usize);
}

return Some((pos, reference_base, alignments));
```

Downstream perbase then iterates `column.raw_alignments()` or `column.alignments()` to build its own row.

## Perbase requirements to keep in mind

### Non-mate `base-depth`

For each alignment in a column, perbase currently applies these semantics:

1. If the read fails the user read filter, it contributes to `FAIL` and to no other output count.
2. If the read is a refskip, it contributes to `REF_SKIP` and not to output depth.
3. If the read is a deletion, it contributes to `DEL` and does contribute to output depth.
4. If the read has a base, it contributes to output depth and one of A/C/G/T/N.
5. If there is a base-quality cutoff and the observed base quality is below that cutoff or unavailable, the base is counted as N.
6. If the read has a simple insertion anchored at this position, it contributes to `INS`.

A direct accumulator can compute this from zero instead of starting from raw pileup depth and subtracting fail/refskip counts.

### Mate-aware `base-depth -m`

Mate-aware mode is harder. It groups alignments by QNAME within each column and picks a representative based on perbase's mate resolution strategy. That requires comparing multiple alignments from the same read name and sometimes using base qualities and recommended IUPAC bases.

Do not make mate-aware mode the first target. The first useful milestone is faster non-mate `base-depth`. Later, the same custom aggregation hook may support mate-aware aggregation, but it will still need some per-column grouping storage.

### `only-depth`

`only-depth` is mostly a separate optimization path. It wants reference intervals/CIGAR spans, not full pileup columns. Do not block this plan on `only-depth`, but keep the same design principle in mind: downstream consumers should be able to request less work when they need less information.

## Goals

1. Add cheap built-in column summaries computed while building each pileup column.
2. Add an opt-in custom per-column aggregation/fold API.
3. Let perbase classify records into count/fail/drop using existing seqair extension points first.
4. Avoid breaking or changing `PileupEngine::pileups()` behavior.
5. Validate everything against the current htslib parity tests and perbase output parity.

## Non-goals for the first prototype

- Do not remove `PileupColumn::alignments()` or `raw_alignments()`.
- Do not change default read filtering semantics.
- Do not make perbase-specific fields part of the default `PileupColumn` API.
- Do not require all users to adopt the accumulator API.
- Do not solve mate-aware aggregation before the non-mate path proves useful.
- Do not change max-depth behavior except where an internal refactor is provably equivalent.

---

# Phase 0: specs, terminology, and measurement

## Add Tracey spec rules first

Before implementation, add spec rules in `docs/spec/4-pileup.md`.

Potential rules:

- `r[pileup.summary]`: `PileupColumn` exposes a column summary.
- `r[pileup.summary_matches_alignments]`: `PileupColumn::summary()` reports the same post-max-depth counts that a caller would compute by iterating `raw_alignments()`.
- `r[pileup.summary_base_counts]`: summary base counts include only alignments with an actual query base.
- `r[pileup.summary_indel_counts]`: summary insertion/deletion/refskip counts are defined in terms of `PileupOp` variants.
- `r[pileup.custom_aggregation]`: an opt-in custom aggregation API can observe each emitted alignment in a column.
- `r[pileup.custom_aggregation_no_materialization]`: custom aggregation does not need to allocate a per-column alignment vector unless the accumulator chooses to store alignments.
- `r[pileup.custom_aggregation_max_depth]`: custom aggregation observes the same alignments that `pileups()` would expose after max-depth truncation.
- `r[pileup.custom_aggregation_existing_api_unchanged]`: the existing lending `pileups()` API remains unchanged.
- `r[pileup.read_disposition_via_extras]`: downstream read disposition can be represented as `RecordStore` extras without changing default binary filtering.

## Define terms carefully

Avoid ambiguous names around insertions and complex indels.

Current `PileupOp` variants:

```rust
Match { qpos, base, qual }
Insertion { qpos, base, qual, insert_len }
Deletion { del_len }
ComplexIndel { del_len, insert_len, is_refskip }
RefSkip
```

Potential summary terminology:

- `depth`: number of emitted alignments in the column; same as current `column.depth()`.
- `match_depth`: number of emitted alignments with a query base; `Match` + `Insertion`.
- `base_counts`: counts bases from `Match` + `Insertion` only.
- `op_counts.match_ops`: count of `PileupOp::Match`.
- `op_counts.insertion_ops`: count of `PileupOp::Insertion`.
- `op_counts.deletion_ops`: count of `PileupOp::Deletion`.
- `op_counts.complex_indel_ops`: count of `PileupOp::ComplexIndel`.
- `op_counts.refskip_ops`: count of `PileupOp::RefSkip`.
- `insertion_count`: decide explicitly whether this means only `Insertion` or any op with `insert_len > 0` (`Insertion` + `ComplexIndel`). For a generic seqair summary, `insertion_count` should probably mean any insertion signal. If this is too ambiguous, prefer `insertion_signal_count` or keep only `op_counts` initially.
- `deletion_count`: `Deletion` + `ComplexIndel { is_refskip: false, .. }`.
- `ref_skip_count`: `RefSkip` + `ComplexIndel { is_refskip: true, .. }`.

perbase's current `INS` semantics may not equal a generic `insert_len > 0` count for `ComplexIndel`. Keep perbase-specific counting in the custom accumulator, not in the built-in summary.

## Baseline before implementation

Before changing code, record current benchmark/profile numbers in seqair and perbase.

perbase baseline from the HG00157 10 Mb BAM subset after seqair v0.1.0 wiring:

- `base-depth`: htslib `4.992 ± 0.207 s`, seqair `5.323 ± 0.127 s`, exact parity.
- `base-depth -m`: htslib `49.455 s`, seqair `32.429 s`, seqair `1.53x` faster, with known sparse default `-F 0` mate-order diffs.
- `only-depth`: htslib `914.5 ± 38.7 ms`, seqair `1.024 ± 0.072 s`, exact parity.
- `only-depth -x`: htslib `757.2 ± 6.5 ms`, seqair `1.012 ± 0.088 s`, exact parity.

Expected first win: non-mate `base-depth`, because it currently materializes/iterates pileup alignments and then performs perbase's counting pass.

---

# Phase 1: built-in `PileupSummary`

## API sketch

Add small summary/count structs to `pileup.rs`.

```rust
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub unknown: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PileupOpCounts {
    pub match_ops: usize,
    pub insertion_ops: usize,
    pub deletion_ops: usize,
    pub complex_indel_ops: usize,
    pub refskip_ops: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PileupSummary {
    pub depth: usize,
    pub match_depth: usize,
    pub insertion_count: usize,
    pub deletion_count: usize,
    pub ref_skip_count: usize,
    pub op_counts: PileupOpCounts,
    pub base_counts: BaseCounts,
}
```

Expose it from `PileupColumn`:

```rust
pub struct PileupColumn<'store, U = ()> {
    pos: Pos0,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    summary: PileupSummary,
    store: &'store RecordStore<U>,
}

impl<'store, U> PileupColumn<'store, U> {
    pub fn summary(&self) -> &PileupSummary;

    // Optional convenience accessors backed by summary:
    pub fn match_depth(&self) -> usize;
    pub fn insertion_count(&self) -> usize;
    pub fn deletion_count(&self) -> usize;
    pub fn ref_skip_count(&self) -> usize;
    pub fn base_counts(&self) -> BaseCounts;
}
```

`match_depth()` currently iterates alignments. After this change it can return `self.summary.match_depth`.

## Summary semantics

`summary()` should describe the emitted column exactly. That means:

- post-max-depth truncation,
- same alignments as `raw_alignments()` exposes,
- same `depth()` as current `self.alignments.len()`,
- no hidden perbase read filtering,
- no base-quality threshold.

If pre-truncation information becomes useful later, add separate fields such as:

```rust
pub raw_depth_before_max_depth: usize,
pub truncated_by_max_depth: bool,
```

Do not mix pre- and post-truncation counts in the first API.

## Implementation approach

Refactor the current column-building loop so summary is updated only for alignments that will be emitted.

Instead of building every alignment and then truncating, apply the max-depth limit while pushing:

```rust
let mut alignments = Vec::with_capacity(self.active.len());
let mut summary = PileupSummary::default();
let max = self.max_depth.map(|m| m as usize);

for active in &self.active {
    if max.is_some_and(|max| alignments.len() >= max) {
        break;
    }

    let Some(info) = active.cigar.pos_info_at(pos) else { continue };
    let alignment = build_alignment(active, info, &self.store);
    summary.observe(&alignment);
    alignments.push(alignment);
}
```

This should be equivalent to the current `truncate(max)` behavior because only the first `max` emitted alignments survive. It also avoids constructing discarded `PileupAlignment`s when max depth is hit.

Add helper methods:

```rust
impl BaseCounts {
    fn add_base(&mut self, base: Base);
}

impl PileupSummary {
    fn observe(&mut self, alignment: &PileupAlignment);
}
```

Keep these helpers private initially unless there is a clear downstream reason to expose them.

## Tests

Add unit tests in `pileup.rs` or an existing pileup test module.

For each test, compute expected summary independently from the emitted alignments in the test body, not by calling the same `PileupSummary::observe` helper.

Cases:

1. Only matches: base counts and match depth.
2. Empty `SEQ=*`: `Base::Unknown`, unavailable quality, still counts in `match_depth` and `base_counts.unknown`.
3. Deletion: `depth` includes deletion, `match_depth` does not, `deletion_count` increments.
4. Refskip: `depth` includes refskip, `match_depth` does not, `ref_skip_count` increments.
5. Simple insertion: `match_depth` and base count increment; insertion op/count increments.
6. Complex indel with deletion: complex op increments; deletion-like count increments; insertion signal behavior matches documented choice.
7. Complex indel with refskip: complex op increments; refskip-like count increments; insertion signal behavior matches documented choice.
8. Max depth: summary equals only the emitted/truncated alignments.
9. `match_depth()` accessor no longer scans but returns the same value as manual iteration.

## Why this helps perbase

This is not enough to replace perbase's counting path because perbase has read filters, fail counts, and optional base-quality filtering. However, it is still useful because:

- it removes repeated scans for generic counts,
- it gives a correctness oracle for the later custom accumulator,
- it establishes clear summary semantics around complex indels and max depth,
- it is a small safe first PR before the more experimental fold API.

---

# Phase 2: internal column observer refactor

Before exposing a public custom aggregation API, refactor the hot loop around an internal observer/sink. The goal is to have one piece of code that:

1. advances positions,
2. maintains the active set,
3. builds `PileupOp`/`PileupAlignment` values,
4. applies max-depth semantics,
5. either materializes alignments or feeds them to a custom accumulator.

## Internal helper shape

Possible internal helper:

```rust
struct ColumnStart {
    pos: Pos0,
    reference_base: Base,
}

struct ColumnEnd {
    depth: usize,
    truncated_by_max_depth: bool,
}

trait ColumnSink<U> {
    type Output;

    fn begin_column(&mut self, start: ColumnStart);

    fn observe_alignment(
        &mut self,
        alignment: PileupAlignment,
        store: &RecordStore<U>,
    );

    fn finish_column(&mut self, end: ColumnEnd) -> Option<Self::Output>;
}
```

Use by-value `PileupAlignment` in `observe_alignment`. The engine has already constructed the value. A materializing sink can push it into a `Vec`; a counting sink can inspect and drop it. This avoids giving the accumulator references that might escape the hot loop.

The existing `pileups()` path can use a private materializing sink that returns:

```rust
struct MaterializedColumn {
    pos: Pos0,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    summary: PileupSummary,
}
```

Then `pileups()` attaches `store: &self.store` and returns `PileupColumn<'_, U>` as it does today.

## Important max-depth rule

The observer must see the same alignments that `pileups()` would expose.

Recommended behavior:

- Maintain an emitted alignment count for the current column.
- If `max_depth` is set and emitted count reaches the limit, stop observing/building additional alignments for that column.
- This matches current `alignments.truncate(max)` output while avoiding work after the limit.

Do not let the custom API observe alignments that `pileups()` would later truncate unless the API explicitly grows a pre-max-depth mode.

## Expected code movement

Likely split current `advance()` into smaller helpers:

```rust
fn advance_active_set_to(&mut self, pos: Pos0) -> Option<()>;
fn next_reference_base(&self, pos: Pos0) -> Base;
fn build_alignment(&self, active: &ActiveRecord, pos: Pos0) -> Option<PileupAlignment>;
fn advance_with_sink<S>(&mut self, sink: &mut S) -> Option<S::Output>
where
    S: ColumnSink<U>;
```

Keep helper visibility private at first.

## Tests

At this phase, existing public behavior should be unchanged. Good tests:

- Run existing pileup parity tests unchanged.
- Add a test-only sink that collects simple `(pos, depth, match_depth)` tuples and compare to `pileups()`.
- Add a max-depth test proving sink-observed output is exactly the materialized output.

---

# Phase 3: public custom aggregation API

After the internal sink refactor is working, expose an opt-in accumulator API.

## API sketch

A minimal infallible API:

```rust
pub struct PileupColumnContext<'store, U> {
    pub pos: Pos0,
    pub reference_base: Base,
    pub store: &'store RecordStore<U>,
    pub depth: usize,
    pub truncated_by_max_depth: bool,
}

pub trait PileupColumnAccumulator<U> {
    type Output;

    fn begin_column(&mut self, pos: Pos0, reference_base: Base);

    fn observe_alignment(
        &mut self,
        alignment: PileupAlignment,
        store: &RecordStore<U>,
    );

    fn finish_column(&mut self, context: PileupColumnContext<'_, U>) -> Option<Self::Output>;
}

impl<U> PileupEngine<U> {
    pub fn pileup_with<A>(&mut self, accumulator: &mut A) -> Option<A::Output>
    where
        A: PileupColumnAccumulator<U>;
}
```

Notes:

- `Output` should be owned. Do not let outputs borrow from the engine/store in v1; that keeps lifetimes tractable and fits perbase rows.
- `observe_alignment` gets `PileupAlignment` by value. Accumulators that need to retain it can store it. Accumulators that only count can avoid per-column allocation.
- Store access is passed to support `RecordStore<U>::extra(record_idx)`, QNAME lookup, aux lookup, and future downstream grouping.
- `finish_column` can return `None` if the accumulator wants to skip output for a column.

If fallibility is needed later, add a second API rather than complicating v1:

```rust
pub trait TryPileupColumnAccumulator<U> {
    type Output;
    type Error;
    ...
}

pub fn try_pileup_with<A>(&mut self, accumulator: &mut A) -> Option<Result<A::Output, A::Error>>;
```

perbase's first use case should not need fallibility.

## Alternative closure API

A closure API may be ergonomic for simple summaries:

```rust
engine.pileups_with(|column, alignment| {
    // update caller state
});
```

However, a trait is probably better for the first serious prototype because it gives explicit lifecycle hooks:

- initialize per-column state,
- observe alignments,
- finish/emit output,
- reuse allocation inside the accumulator across columns.

A closure helper can be layered on later.

## Example: reproduce `PileupSummary`

A seqair test/example accumulator:

```rust
#[derive(Default)]
struct SummaryAccumulator {
    summary: PileupSummary,
}

impl<U> PileupColumnAccumulator<U> for SummaryAccumulator {
    type Output = PileupSummary;

    fn begin_column(&mut self, _pos: Pos0, _reference_base: Base) {
        self.summary = PileupSummary::default();
    }

    fn observe_alignment(&mut self, alignment: PileupAlignment, _store: &RecordStore<U>) {
        self.summary.observe(&alignment);
    }

    fn finish_column(&mut self, _context: PileupColumnContext<'_, U>) -> Option<Self::Output> {
        Some(self.summary)
    }
}
```

Test it by comparing `engine.pileup_with(&mut SummaryAccumulator)` against `engine.pileups().summary()` over the same synthetic store.

## Example: perbase-like non-mate accumulator

Pseudo-code for the shape perbase wants:

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PerbaseDisposition {
    Count,
    FailOnly,
}

struct PerbaseColumnAccumulator {
    min_base_quality: Option<u8>,
    row: PerbaseRow,
}

impl PileupColumnAccumulator<PerbaseDisposition> for PerbaseColumnAccumulator {
    type Output = PerbaseRow;

    fn begin_column(&mut self, pos: Pos0, reference_base: Base) {
        self.row = PerbaseRow::new(pos, reference_base);
    }

    fn observe_alignment(
        &mut self,
        alignment: PileupAlignment,
        store: &RecordStore<PerbaseDisposition>,
    ) {
        match store.extra(alignment.record_idx()) {
            PerbaseDisposition::FailOnly => {
                self.row.fail += 1;
                return;
            }
            PerbaseDisposition::Count => {}
        }

        match alignment.op() {
            PileupOp::RefSkip => {
                self.row.ref_skip += 1;
            }
            PileupOp::ComplexIndel { is_refskip: true, .. } => {
                self.row.ref_skip += 1;
            }
            PileupOp::Deletion { .. } | PileupOp::ComplexIndel { is_refskip: false, .. } => {
                self.row.depth += 1;
                self.row.del += 1;
            }
            PileupOp::Match { base, qual, .. } => {
                self.row.depth += 1;
                self.row.add_base_with_quality(*base, *qual, self.min_base_quality);
            }
            PileupOp::Insertion { base, qual, .. } => {
                self.row.depth += 1;
                self.row.add_base_with_quality(*base, *qual, self.min_base_quality);
                self.row.ins += 1;
            }
        }
    }

    fn finish_column(&mut self, _context: PileupColumnContext<'_, PerbaseDisposition>) -> Option<Self::Output> {
        Some(self.row.clone())
    }
}
```

This mirrors perbase's current non-mate semantics but avoids:

- constructing a public `PileupColumn`,
- storing an internal per-column `Vec<PileupAlignment>` for simple counting,
- a downstream second pass over `column.raw_alignments()`.

## Tests

1. Summary accumulator equals `PileupColumn::summary()`.
2. A toy accumulator that counts only depth equals `column.depth()`.
3. A perbase-like accumulator equals manual iteration over `raw_alignments()` for:
   - passing reads,
   - fail-only reads,
   - base-quality demotion to N,
   - deletion,
   - refskip,
   - insertion,
   - complex indel.
4. Max-depth gives identical observed alignments/counts to current `pileups()`.

---

# Phase 4: read disposition for pass/fail/drop

## Problem

`CustomizeRecordStore::filter_raw` and `filter` are binary: keep or reject.

That is ideal for use cases where rejected reads truly disappear, but perbase `base-depth` has a `FAIL` output column. Reads that fail user filters should still contribute to `FAIL` if they cover the position.

So perbase needs three conceptual states:

```rust
enum ReadDisposition {
    Count,
    FailOnly,
    Drop,
}
```

Meaning:

- `Count`: read contributes normally.
- `FailOnly`: read remains in the pileup for coverage overlap but contributes only to `FAIL`.
- `Drop`: read is omitted completely.

## Prototype using existing seqair APIs

Do **not** start by changing `CustomizeRecordStore`.

Use existing extras first:

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PerbaseDisposition {
    Count,
    FailOnly,
}
```

Then implement a customizer in perbase:

```rust
impl CustomizeRecordStore for PerbaseCustomize {
    type Extra = PerbaseDisposition;

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        // Return false only for true Drop cases.
        // For base-depth, normal user read-filter failures are NOT dropped,
        // because they need to be counted as FAIL per covered position.
        true
    }

    fn compute(
        &mut self,
        rec: &SlimRecord,
        _store: &RecordStore<PerbaseDisposition>,
    ) -> PerbaseDisposition {
        if self.read_filter_passes(rec.flags, rec.mapq) {
            PerbaseDisposition::Count
        } else {
            PerbaseDisposition::FailOnly
        }
    }
}
```

Then the custom pileup accumulator reads:

```rust
store.extra(alignment.record_idx())
```

and decides whether an observed alignment contributes to normal counts or fail-only counts.

## Why extras-first is the right first step

- No public seqair filtering API churn.
- Perbase can test the semantics without waiting on a new record-store trait.
- It uses the existing `RecordStore<U>` design exactly as intended.
- The custom accumulator API must support store extras anyway.

## Potential future seqair API

If many users need this pattern, consider a first-class classification trait later:

```rust
pub enum ReadDisposition {
    Count,
    KeepAsExtraOnly,
    Drop,
}

pub trait ClassifyRecordStore: Clone {
    type Extra;

    fn classify_raw(&mut self, fields: &FilterRawFields<'_>) -> ReadDisposition;

    fn compute(
        &mut self,
        rec: &SlimRecord,
        store: &RecordStore<Self::Extra>,
    ) -> Self::Extra;
}
```

But this should be driven by evidence from the perbase prototype. It may turn out that `CustomizeRecordStore::compute` + extras is sufficient.

## Optional future decode-profile optimization

For perbase fail-only reads, the accumulator only needs to know that the read covers a position. It does not need base/quality values because fail-only reads do not contribute to base counts.

A future, larger optimization could combine disposition with decode profiles:

- `Count`: decode/store CIGAR + SEQ + QUAL + needed metadata.
- `FailOnly`: decode/store CIGAR + minimal metadata, maybe skip SEQ/QUAL.
- `Drop`: skip entirely.

That is outside the first pileup aggregation prototype, but it is worth keeping in mind if fail-only reads are a meaningful cost.

---

# Phase 5: perbase integration experiment

Once seqair has the custom aggregation API on a branch, wire perbase's experimental seqair backend to it.

## Initial perbase target

Start with non-mate `base-depth --seqair-pileup` only.

Suggested perbase flow:

1. Define `PerbaseDisposition` as a seqair `RecordStore` extra.
2. Fetch records into `RecordStore<PerbaseDisposition>`.
3. Use `PileupEngine<PerbaseDisposition>::pileup_with(&mut PerbaseColumnAccumulator)`.
4. Produce `PileupPosition` rows directly from accumulator output.
5. Keep `base-depth -m` on the current materialized `PileupColumn` path at first.

## Perbase parity checklist

Compare htslib and seqair custom aggregation for:

- default filters,
- `-F` exclude flags,
- `-f` include flags if applicable,
- `-q` min MAPQ,
- `-Q` min base quality,
- insertions,
- deletions,
- refskips,
- empty `SEQ=*`,
- max depth,
- zero-position filling.

Expected row semantics for non-mate custom aggregation:

- `DEPTH`: count of passing reads with match/insertion/deletion observations; excludes refskips and fail-only reads.
- `A/C/G/T/N`: passing reads with base observations after base-quality demotion.
- `INS`: passing reads with simple insertion observations matching current perbase behavior.
- `DEL`: passing reads with deletion observations.
- `REF_SKIP`: passing reads with refskip observations.
- `FAIL`: fail-only reads with any pileup observation at the position.

## Benchmark plan

Run three versions:

1. htslib baseline,
2. current seqair materialized pileup path,
3. seqair custom accumulator path.

Commands should match the latest perbase benchmark setup:

- release build,
- `--features seqair-pileup`,
- HG00157 chr1 10 Mb BAM subset,
- same BED/window split,
- TSV output to disk or same output sink for all runs.

Benchmarks:

- `base-depth` non-mate default filters: primary target.
- `base-depth -Q <threshold>`: validates base-quality path and overhead.
- `base-depth -F 2304`: useful exact-parity control for secondary/supplementary behavior.
- `base-depth -m`: control, should remain unchanged until a mate-aware accumulator exists.

Success criteria:

- Byte-for-byte parity with current htslib/reference outputs for non-mate mode.
- No regressions in existing seqair htslib pileup parity tests.
- Custom accumulator path is measurably faster than current seqair non-mate `base-depth`, or profiling clearly identifies a different bottleneck.

---

# Phase 6: mate-aware follow-up, only if worthwhile

If non-mate custom aggregation works, explore mate-aware aggregation.

Challenges:

- Need group-by-QNAME per column.
- Need compare multiple observations from the same QNAME using perbase's mate strategy.
- Need base quality and possibly recommended IUPAC base output.
- Current `AlignmentView::qname()` borrows from the store; storing borrowed QNAMEs in an accumulator is lifetime-sensitive.

Possible approaches:

1. Accumulator stores owned QNAME bytes or hashes for grouping. Simple, but allocates.
2. Accumulator stores `record_idx` and sorts/group by QNAME at finish using store lookups. This may resemble current materialized behavior, but could still avoid some wrapper work.
3. seqair provides a specialized qname-grouped column helper later. This is probably too downstream-specific for v1.

Do not optimize this until the non-mate path shows the API is viable.

---

# Implementation order

Recommended PR sequence in the fork:

1. **Spec-only PR or first commit**: add Tracey spec rules for summaries/custom aggregation.
2. **Summary PR**:
   - add `BaseCounts`, `PileupOpCounts`, `PileupSummary`,
   - add `PileupColumn::summary()`,
   - switch `match_depth()` to summary,
   - test summary against manual iteration.
3. **Internal sink refactor PR**:
   - introduce private `ColumnSink`,
   - keep `pileups()` output unchanged,
   - prove with tests that behavior matches.
4. **Public accumulator PR**:
   - expose `PileupColumnAccumulator`, context type, and `pileup_with`,
   - add examples/tests.
5. **perbase branch**:
   - use `RecordStore<PerbaseDisposition>`,
   - add non-mate accumulator,
   - benchmark and compare outputs.
6. **API cleanup PR**:
   - rename types/methods based on real perbase usage,
   - decide whether fallible accumulator or first-class disposition is needed.

---

# Open questions

## Summary semantics

- Should `insertion_count` include `ComplexIndel`, or should the API avoid a single insertion count and expose only precise `op_counts` initially?
- Should summary expose `truncated_by_max_depth` and pre-truncation depth now, or wait until needed?
- Should `BaseCounts` use `usize`, `u32`, or a small fixed array `[usize; 5]` internally?

## Accumulator API

- Should `observe_alignment` receive `PileupAlignment` by value, by reference, or an `AlignmentView`-like wrapper?
- Should output be required to be owned in v1?
- Should the first API be infallible only?
- Should there be both `pileup_with` for one column and `fold_columns` for full-region consumption?
- Should a custom accumulator be able to request pre-max-depth observations?

## Perbase integration

- How should perbase represent `PerbaseDisposition` without making seqair depend on perbase concepts?
- Can perbase's read filter be evaluated fully from `SlimRecord` fields (`flags`, `mapq`) for all current modes?
- How should complex indel insertion counts be handled to preserve current htslib parity?
- How much does avoiding the internal alignment `Vec` help compared with other costs such as CIGAR mapping, base decode, and output formatting?

## Future performance work

- If custom aggregation is still slower than htslib for non-mate `base-depth`, profile whether the bottleneck is:
  - CIGAR position lookup,
  - base/quality lookup,
  - active set maintenance,
  - record-store decode,
  - output construction/serialization.
- If fail-only reads are common, consider decode profiles to skip SEQ/QUAL for fail-only records.
- For `only-depth`, pursue a separate CIGAR/interval-only API rather than forcing it through pileup.
