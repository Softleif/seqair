# Window query over a prepared store; the load-time hook sees the reference

2026-09-10, revised 2026-09-11 after review and again after the second
review round, branch `feature/window-query`.
Two small features the rastair 3 design asks for
(`docs/plans/rastair3-design-2026-09-10.md` §7 row C3, §11.2, and §12 E6 in
the rastair repo).

## 1. `records_overlapping`

The rastair 3 slow path re-solves an active region — a few hundred bases
around a cluster of candidate indels — against the reads that span it, using
the store the pileup engine already holds. It needs the record indices
overlapping a span without walking the store.

### API

```rust
impl<U> PileupInput<U> {
    pub fn store(&self) -> &RecordStore<U>;
    pub fn records_overlapping(&self, start: Pos0, end: Pos0) -> impl Iterator<Item = RecordIdx> + '_;
}

impl<U> PileupEngine<U> {
    pub fn records_overlapping(&self, start: Pos0, end: Pos0) -> impl Iterator<Item = RecordIdx> + '_;
}

impl<'eng, U> PileupColumn<'eng, U> {
    pub fn records_overlapping(&self, start: Pos0, end: Pos0) -> impl Iterator<Item = RecordIdx> + 'eng;
}
```

- `[start, end]` is **inclusive at both ends**, like every other interval in
  the API (`r[interval.inclusive_ends]`); overlap is
  `r[interval.overlap_test]`: `pos <= end && end_pos >= start`, `end_pos`
  being inclusive. `end < start` is the empty interval and yields nothing.
- Only **mapped** records are yielded. A placed-unmapped read sits at its
  mate's position with `end_pos == pos` and no alignment; the pileup never
  reports one (`r[pileup.unmapped_excluded]`), and the query answers for
  what the pileup reports. Review caught this: the first version yielded
  them, and nothing in the tests could tell, because the generator never
  made one.
- Indices come out ascending and are the store's record indices (`RecordIdx`
  since the second review round), valid for `RecordStore::record` /
  `try_record` and for `PileupColumn::find_record` on any column of an engine
  built from the input. They mean nothing once the engine's store has been
  reclaimed (guard drop): resolve them before that.
- There is deliberately no such method on `RecordStore`. The query
  binary-searches an index that only exists on a store
  `prepare_for_pileup` has ordered, and a raw store cannot prove it has been
  — the same silent failure `PileupInput` was introduced to close. The
  crate-private `records_overlapping_sorted` on the store carries a
  `debug_assert` on the order flag; the three public entry points are the
  types that can vouch for it.
- The **column** entry point exists because `pileups(&mut self)` hands out a
  column that borrows the engine for as long as it lives, so
  `engine.records_overlapping` does not borrow-check while a column is held
  — and a caller who has just seen a trigger at a column is exactly the
  caller who wants the reads around it. Review found the first version
  promised "at any point of the iteration" while only an in-crate test could
  reach it.
- What a column reports and what the query yields agree on alignment overlap
  and differ in three documented cases: soft-clip overhang (column reports
  more), depth cap and engine region (column shows less).

### Lower bound: a running maximum of `end_pos`, built at `prepare_for_pileup`

A record starting before `start` can still overlap the window, so the first
record with `pos >= start` is not the first candidate, and with no bound at
all a backward scan has no stopping point — a long enough record arbitrarily
far back still overlaps. The bound is needed for correctness; the only
question was its shape.

The first version tracked one `u32` per store, `max_ref_span`, the widest
`end_pos - pos + 1` ever kept, widened on every kept push and on
`set_alignment`, reset by `clear`, carried by `take_contents`, and searched
`pos >= start - (max_span - 1)`. Review (adversarial probing, mutation
testing) confirmed it was a sound upper bound on every path — and found two
things wrong with it as a design:

- **One spliced read collapses it.** `N` consumes reference, so a spliced
  read's span is its intron; after one 100 kb intron, `start - max_span`
  saturates to 0 for every realistic query and the binary search lands on
  record 0 — a full scan of the store, silently, forever (the bound never
  shrinks). RNA-seq is not rastair's input today, but a design that
  degrades to O(n) on one read is the wrong one.
- **It is an invariant maintained by hand** across five sites, and mutation
  testing showed one of them (`push_raw`, the real BAM path) was never
  exercised by a query in the test suite.

The revision keeps `reach: Vec<Pos0>` on the store, `reach[i]` = the largest
`end_pos` among the mapped records at indices `..= i`, built in one pass by
`prepare_for_pileup` after the sort. It never decreases, so
`partition_point(reach[i] < start)` is exactly the first record that can
overlap — the record at that index is itself a hit or is already past `end`,
nothing before it reaches `start`, nothing after it is skipped. The forward
scan then runs to `pos > end` as before.

- Exact, not an over-estimate: a long read costs only the queries it
  genuinely overlaps (where it is a hit, and the records between it and the
  window must be inspected by any list-shaped structure anyway).
- Nothing to maintain: derived at the one point the order is established,
  cleared with the records, moved with them. `set_alignment` and pushes
  after that point already flip the order flag, which the query asserts.
- Cost: 4 bytes per record, one pass per prepared store. The alternative was
  an interval tree, which is the right answer if the slow path ever queries
  thousands of windows per tile and the forward scan shows up in a profile.

### Spec rules

- `r[record_store.window_query]` — the contract, mapped-only, and why it is
  not on `RecordStore`.
- `r[record_store.window_query.reach]` — the index, its exactness, and why
  it replaced the tracked span.
- `r[pileup.records_overlapping]` — engine and column entry points, the
  three ways a column may differ, and index validity.

### Tests

- `record_store::tests::window_query::matches_brute_force` (proptest):
  random stores of 0–80 records with `pos` in `0..5000`, 1–16 query bases,
  deletions bimodal (none, 1–40, or one long 200–1500 so a single record
  sets how far back a window must look), one in ten placed-unmapped; pushed
  in random order with unique qnames; 0–6 `set_alignment` moves before
  `prepare_for_pileup`; 1–8 spans per case, absolute (both ends from
  `0..7000`, so often inverted or past the last record) or anchored a base
  or two off a record's own start or end. Property: the index set equals a
  brute-force filter of every record by the definition.
- `first_reaching_is_where_the_linear_scan_stops` (proptest): the binary
  search lands where a linear scan for the first mapped record with
  `end_pos >= start` stops.
- Unit tests: a 300 bp-deletion read starting before the window is found;
  first/last covered base; empty, before-first, past-last spans; empty
  store; a placed-unmapped record is not yielded; `set_alignment`
  lengthening a record is seen; a store cleared and refilled answers for its
  new records, not a stale index.
- `pileup::tests::window_query::engine_query_agrees_with_columns` (proptest):
  the engine's answer equals the brute force before iteration, from a column
  mid-iteration (through `PileupColumn::records_overlapping`), and after the
  last column; every alignment a column inside the span reports is in the
  set, and every record in the set is reported by some column inside the
  span — the direction that keeps unmapped records out without a second copy
  of the rule; each yielded index resolves with `find_record` exactly at the
  columns it covers; after `reclaim_allocation` the query is empty.
- `tests/window_query.rs` (integration, public API only, real reads through
  `push_raw`): sampled columns of the fixture segment match brute force from
  the column and from the engine; and the two features composed — a hook
  that soft-clips leading bases, then the query names the moved reads where
  they now lie, and the hook's reference bases equal the columns'.
- `examples/active_region.rs`, smoke-tested: triggers on indel columns,
  gathers the spanning reads from the column, prints one line per region.

## 2. The load-time hook sees the reference

rastair's E6 prototype found that the `Readers::pileup_with` hook ran before
seqair fetched the engine's `RefSeq`, so a hook that wanted to score against
the reference had to load a second copy.

### API

```rust
impl<E: CustomizeRecordStore> Readers<E> {
    pub fn pileup<'a>(&'a mut self, segment: &'a Segment, depth: DepthLimit)
        -> Pileup<'a, E>;
}

impl<'a, E: CustomizeRecordStore, F> Pileup<'a, E, F>
where F: FnMut(&mut RecordStore<E::Extra>, &RefSeq) {
    pub fn with_reference(self, ref_seq: RefSeq) -> Self;
    pub fn mutate<G>(self, f: G) -> Pileup<'a, E, G>;
    pub fn run(self) -> Result<PileupGuard<'a, E::Extra>, ReaderError>;
}
```

`run` fetches (or validates the supplied) reference first and runs the
mutator with `&ref_seq` before `prepare_for_pileup`. The `RefSeq` handed to
the mutator is the very value attached to the engine, so
`ref_seq.base_at(pos)` in the hook equals `column.reference_base()` later.

The first version kept `pileup_with` with its old one-argument signature and
added `pileup_mutate` for the two-argument one. seqair has no stable API to
protect, so the revision has one hook with one signature; a mutator that
does not need the reference ignores the argument.

The second review round replaced the three entry points with the plan above.
`pileup`, `pileup_with` and `pileup_with_reference` covered three of the four
combinations of the two options; the fourth — hold the region's bases *and*
realign against them, which is what a tiled slow path does — had no spelling,
though `pileup_impl` already took both arguments. The options are
independent, so a fourth method would have been the wrong fix. rastair's E6
hook becomes `.pileup(seg, depth).mutate(|store, _| …).run()`.

Two consequences worth knowing:

- The reference covers **exactly the segment** (`r[unified.readers_pileup]`
  step 3), not the reads. A record reaching past either end has bases the
  `RefSeq` does not hold; read those with `try_base_at` and treat `None` as
  unavailable, because `base_at` returns `Base::Unknown` there,
  indistinguishable from an `N`. Widening the fetch to the store's extent is
  a possible follow-up if the slow path needs it.
- Because the hook runs after the fetch, it does not run when the fetch
  fails: the caller gets the error and an untouched store. The next
  `pileup` clears the store before refilling, so nothing leaks across.

### Spec rules

- `r[unified.readers_pileup_store_mutation+3]`: the hook runs after step 3,
  receives the reference, the coverage caveat, the error-path behaviour, and
  the version history.

### Tests

`reader::readers::tests::pileup_with_sees_the_engines_reference`: the
mutator asserts the reference starts and ends at the segment (`try_base_at`
is `None` one base outside either end), records `base_at` for every segment
position, and the pileup loop asserts `reference_base()` at every column
equals the recorded base. The two older `pileup_with` tests (no-op is
transparent; a soft-clipping realignment is observed) are unchanged apart
from the extra closure argument.

## 3. Done in the second review round

- **The plan** (§2 above) replaced the three entry points and made the
  supplied-reference-plus-hook combination expressible.
- **`RecordIdx`,** in its own module with **`RecordRef`**. The `u32` store
  indices this query, `find_record`, `record` and `mate_idx` shared are one
  newtype now, crate-wide; `u32::MAX` is unrepresentable, so
  `Option<RecordIdx>` is free and both mate-link sentinels are gone.

  `RecordStore::record(idx) -> Option<RecordRef<'_>>` is the only way in, and
  it does not panic. The nine by-index readers that each dereferenced the same
  untrusted index — and each panicked — are methods on the handle instead. The
  handle holds the store's borrow, not just a lifetime, which is what makes
  the guarantee real rather than documented: `clear`, a push, `sort_by_pos`
  and `dedup` are rejected at compile time while one is alive, and it cannot
  outlive its store. Four `compile_fail` doctests pin that to E0502/E0515, so
  the proof fails loudly if the shape ever changes.

  True *uniqueness* branding — a handle from store A rejected by store B — was
  not attempted: it needs invariant generative lifetimes and a closure-scoped
  API (`with_store(|s| …)`), which the forked-reader and `PileupGuard` shapes
  cannot absorb. The handle gets the same effect a different way: the
  accessors are methods on it, so there is no call that could name a second
  store. What remains is that two *live* borrows of different stores share a
  region, which no call site can exploit because none takes both.
- **`cargo doc -D warnings` in CI**, over the library and the examples. The
  22 + 7 warnings it was hiding were mostly `r[`rule.id`]` and `[SAM1]`
  being read as intra-doc links; one was a real rotted link in
  `simple_variant_caller`'s tutorial.
- **Reference beyond the segment for the hook**, as
  `Pileup::reference_covers_reads()`. Measured on `tests/data/test.bam`
  (80 bp reads, chr19:6.1–6.2 Mb), the share of reads reaching past their
  tile — and so with bases the hook's `RefSeq` did not hold — was 27.2 % at
  a 500 bp tile, 13.9 % at 1 kb, 3.0 % at 5 kb and 0 % at 50 kb, with the
  overhang bounded by the longest read (79 bases) at both ends at every tile
  size.

  The widening is exact rather than a padding constant, because step 2 runs
  before step 3: the span is the records' own `min(pos)`..`max(end_pos)`,
  unioned with the segment and clamped to the contig, so a tile whose reads
  all stop inside it fetches exactly what it fetched before. Columns cannot
  change — `RefSeq` resolves absolute positions — which is its own test.

  Opt-in, because the policy is the caller's: a hook that only scores within
  the tile wants the segment and the smaller fetch, and one spliced read can
  make the wider one large. Combined with `with_reference` there is nothing
  to widen, so it becomes a requirement on the supplied reference
  (`ReaderError::SuppliedReferenceMissesReads`) rather than silently doing
  nothing.

- **The fuzz crate builds again**, and its CI step gates. `PileupInput`, the
  `NonZeroU32` max-depth and the `AlignedPair::Insertion` field rename all
  landed before this branch and left five targets unbuildable; nobody was
  told because the "Fuzz compiles" step is `continue-on-error` and was not
  named in the job's failure condition. Both are fixed here.

## Next step

In rastair (the rastair 3 slow path, design doc §7 row C3 / §11.2 / §12 E6):
`col.records_overlapping(start, end)` at a triggering column to gather the
reads spanning the active region, resolved with `col.store().record(idx)`
or `find_record` on later columns; `.pileup(seg, depth).mutate(f).run()` for
the hook, with the reference. Before merging into seqair `main`: rebase, run
`cargo nextest run`, `cargo test --doc --quiet`,
`cargo clippy --all-targets -- -D warnings` and
`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps --all-features --examples`,
push, and move rastair's `Cargo.lock` seqair pin (its E6 hook needs the
extra closure argument, `.run()`, and `RecordIdx` wherever it holds a store
index).
