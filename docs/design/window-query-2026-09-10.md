# Window query over a prepared store; load-time mutator sees the reference

2026-09-10, branch `feature/window-query`. Two small features the rastair 3
design asks for (`docs/plans/rastair3-design-2026-09-10.md` §7 row C3, §11.2,
and §12 E6 in the rastair repo).

## 1. `records_overlapping`

The rastair 3 slow path re-solves an active region — a few hundred bases
around a cluster of candidate indels — against the reads that span it, using
the store the pileup engine already holds. It needs the record indices
overlapping a span without walking the store.

### API

```rust
impl<U> PileupInput<U> {
    pub fn store(&self) -> &RecordStore<U>;
    pub fn records_overlapping(&self, start: Pos0, end: Pos0) -> impl Iterator<Item = u32> + '_;
}

impl<U> PileupEngine<U> {
    pub fn records_overlapping(&self, start: Pos0, end: Pos0) -> impl Iterator<Item = u32> + '_;
}
```

- `[start, end]` is **inclusive at both ends**, like every other interval in
  the API (`r[interval.inclusive_ends]`); overlap is
  `r[interval.overlap_test]`: `pos <= end && end_pos >= start`, `end_pos`
  being inclusive. `end < start` is the empty interval and yields nothing.
- Indices come out ascending and are the store's record indices, valid for
  `RecordStore::record` and for `PileupColumn::find_record` on any column of
  an engine built from the input.
- There is deliberately no such method on `RecordStore`. The query
  binary-searches on `pos`, which means nothing on an unsorted store, and a
  raw store cannot prove it is sorted — the same silent failure `PileupInput`
  was introduced to close. The crate-private `records_overlapping_sorted` on
  the store carries a `debug_assert` on the order flag; the two public
  entry points are the types that can vouch for it. The engine never
  reorders its store, so its answer is the input's; after
  `reclaim_allocation` the store is empty and the query yields nothing.

### Lower bound: tracked `max_ref_span`, not a backward scan

A record starting before `start` can still overlap the window, so the first
record with `pos >= start` is not the first candidate. The store now tracks
`max_ref_span`, the widest `end_pos - pos + 1` it has kept, widened on every
kept push and on `set_alignment`, reset by `clear`, carried by
`take_contents`. The query starts at `partition_point(pos < start - (max_span - 1))`
and scans forward while `pos <= end`, keeping records with `end_pos >= start`.

This was a decision by reasoning, not a measurement (the instruction changed
mid-task to run no timings on the machine):

- A backward scan from `partition_point(pos < start)` has **no stopping
  point without a span bound** — a record arbitrarily far back overlaps if
  it is long enough. So the bound is needed for correctness either way, and
  the only question is whether to walk to it or jump to it.
- Given the bound, the backward scan and the binary search stop at the same
  first record. At ~4k records per 10 kb region with 150 bp reads, the
  records inside one span's width are ~60 either way; the binary search is
  `O(log n)` regardless of span and is one `partition_point` call, which is
  the simpler thing to reason about.
- The bound may over-estimate (a rolled-back push does not widen it — it is
  updated after the keep decision — and a record `dedup` removes does not
  lower it) but never under-estimates. Over-estimating only widens the scan
  prefix; under-estimating would skip records with no symptom. One long
  record (a spliced read with a large `N`) widens the prefix for the whole
  store; acceptable, and the alternative is a per-query `max` over the
  records.
- Cost: one `max` per kept push, a `u32` on the store.

### Spec rules

- `r[record_store.window_query]` — the contract, and why it is not on
  `RecordStore`.
- `r[record_store.window_query.max_ref_span]` — the bound, its monotonicity,
  and the decision above.
- `r[pileup.records_overlapping]` — the engine forwards to the same query.

### Tests

- `record_store::tests::window_query::matches_brute_force` (proptest): random
  stores of 0–80 records with `pos` in `0..5000`, 1–16 query bases and, one
  time in five, a trailing deletion of 1–400 so spans vary independently of
  query length; pushed in random order (`prepare_for_pileup` sorts); 1–8
  spans per case with both ends drawn from `0..7000`, so spans are often
  inverted, before the first record, or past the last. Property: the index
  set equals a brute-force filter of every record by the overlap test (with
  `end < start` defined as empty).
- `span_bound_is_the_widest_record` (proptest): the tracked bound equals the
  widest generated record and resets on `clear`.
- Unit tests: a 300 bp-deletion read starting before the window is found; the
  first/last covered base; empty, before-first and past-last spans; an empty
  store; `set_alignment` lengthening a record widens the bound; a rejected
  push leaves it alone.
- `pileup::tests::window_query::engine_query_agrees_with_columns` (proptest):
  the engine's answer equals the brute force before iteration and from a
  column mid-iteration; every alignment a column inside the span reports is
  in the set, and every record in the set is reported by some column inside
  the span (records consume reference by construction); after
  `reclaim_allocation` the query is empty.

## 2. The load-time mutator sees the reference

rastair's E6 prototype found that the `Readers::pileup_with` hook ran before
seqair fetched the engine's `RefSeq`, so a hook that wanted to score against
the reference had to load a second copy.

### API

```rust
impl<E: CustomizeRecordStore> Readers<E> {
    // unchanged signature; now runs after the reference fetch
    pub fn pileup_with<F>(&mut self, segment: &Segment, depth: DepthLimit, mutate: F)
        -> Result<PileupGuard<'_, E::Extra>, ReaderError>
    where F: FnMut(&mut RecordStore<E::Extra>);

    // new
    pub fn pileup_mutate<F>(&mut self, segment: &Segment, depth: DepthLimit, mutate: F)
        -> Result<PileupGuard<'_, E::Extra>, ReaderError>
    where F: FnMut(&mut RecordStore<E::Extra>, &RefSeq);
}
```

`pileup_impl` now fetches (or validates the supplied) reference first and
runs the mutator with `&ref_seq` before `prepare_for_pileup`; `pileup_with`
is a wrapper closure that drops the reference. The `RefSeq` handed to the
mutator is the very value attached to the engine, so
`ref_seq.base_at(pos)` in the hook equals `column.reference_base()` later.

One limitation worth knowing: the reference covers **exactly the segment**
(`r[unified.readers_pileup]` step 3), not the reads. A record reaching past
either end has bases the `RefSeq` does not hold; read those with
`try_base_at` and treat `None` as unavailable, because `base_at` returns
`Base::Unknown` there, indistinguishable from an `N`. Widening the fetch to
the store's extent is a possible follow-up if the slow path needs it.

### Spec rules

- `r[unified.readers_pileup_store_mutation]` bumped to `+2`: the hook runs
  after step 3 now, with a paragraph on why the order moved.
- `r[unified.readers_pileup_mutate]` — the new entry point, the identity of
  the reference, the segment-only coverage, and that `pileup_with` stays.

### Test

`reader::readers::tests::pileup_mutate_sees_the_engines_reference`: the
mutator asserts the reference starts and ends at the segment (`try_base_at`
is `None` one base outside either end), records `base_at` for every segment
position, and the pileup loop asserts `reference_base()` at every column
equals the recorded base. The two existing `pileup_with` tests (no-op is
transparent; a soft-clipping realignment is observed) are unchanged apart
from their annotation version.

## Verification

- `cargo nextest run`: 1490 passed, 1 skipped.
- `cargo test --doc --quiet`: 9 passed.
- `cargo clippy --all-targets -- -D warnings`: clean.
- rastair build against this branch: see below.

### rastair against this branch

Scratch rastair worktree at the `feature/phasing` tip (`44c09c14`), with
`seqair`/`seqair-types` pointed at this worktree by path and
`rust-version = "1.98.1"` (cargo `-v` confirms both crates resolve from
`.claude/worktrees/window-query`):

- `cargo build --release --features experimental-seqair`: **builds**.
- `cargo test --features experimental-seqair --test call_cli`: **37 passed,
  0 failed**, with `Using experimental seqair backend` in the log, so the
  seqair reading path was the one exercised.

Both features are additive — `pileup_with` keeps its signature and the
store's new `max_ref_span` field is private — so rastair compiles without a
source change. Nothing in rastair calls the new entry points yet.

Tracey (`tracey query status` / `stale`): every new rule has impl and verify
references; no stale references against the `+2` bump.

## Done / Half-done / Next step

**Done.** Both features, their spec rules, unit tests and proptests, the
rastair build check, this note. Four commits on `feature/window-query`, not
pushed; the worktree is clean.

**Half-done.** Nothing.

**Next step.** In rastair (the rastair 3 slow path, design doc §7 row C3 /
§11.2 / §12 E6): call `engine.records_overlapping(start, end)` to gather the
reads spanning an active region and `Readers::pileup_mutate` where the
realignment hook needs the reference. Two things to decide there, both
recorded above rather than solved here:

- Whether the slow path needs reference bases beyond the segment for reads
  reaching past its ends. If so, widen the fetch in `pileup_impl` to the
  store's extent (`min pos`, `max end_pos`) — a seqair change, spec rule
  `r[unified.readers_pileup_mutate]` would need its coverage paragraph
  amended.
- Before merging into seqair `main`: rebase onto the current `main` (other
  sessions commit there concurrently), re-run `cargo nextest run`,
  `cargo test --doc --quiet` and `cargo clippy --all-targets -- -D warnings`,
  then push and move rastair's `Cargo.lock` seqair pin.
