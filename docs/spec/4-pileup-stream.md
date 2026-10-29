# Streaming Pileup

The indexed pileup (`r[unified.readers_pileup+2]`) answers one segment at a
time: it needs an index and a seekable file, and it loads every record that
overlaps the segment before the first column. A tool that reads a whole
coordinate-sorted alignment from `stdin` or a pipe — `samtools mpileup` on
unindexed input — has neither, and wants the columns of every contig, in
order, as the records go past.

`StreamPileup` is that pileup. It reads an `AlignmentStream`
(`r[unified.stream.records]`) front to back and yields the same columns the
indexed pileup yields, driving the existing `RecordStore` and `PileupEngine`
over a sliding window of positions.

> **Sources:** htslib `bam_plp_push` (`sam.c`) is the behavioural reference
> for which records a pileup takes and for the sort check
> ("The input is not sorted"); its `overlap_push`/`tweak_overlap_quality`
> is what `samtools mpileup` does with overlapping mates. See
> [References](./99-references.md).

## What it yields

r[pileup_stream.columns]
`StreamPileup::pileups()` MUST yield, contig by contig in input order, exactly
the columns the indexed pileup yields for one segment spanning the whole contig
(`Readers::pileup` over a `Segment` from `Readers::segments(contig, ..)` with
`max_len` at least the contig length), given the same customize value and the
same `DepthLimit`: the same positions, and at each the same depth, reference
base, `matches_reference`, and the same alignments in the same order, each with
the same record, `PileupOp`, `indel_after`, `in_mate_overlap` and mate in the
column (`mate_of`). It MUST yield nothing at a position past the contig's last
base, as the segment does not reach there. Each column comes with its `Tid`,
since a `PileupColumn` carries only a position.

Store indices are the exception: `record_idx`, `mate_idx` and `fragment_idx`
index the window's store (`r[pileup_stream.window+2]`), not a whole-contig one,
and `mate_idx` is narrower (`r[pileup_stream.mates+1]`).

r[pileup_stream.lending]
Iteration MUST be lending, as the engine's is (`r[pileup.lending_iterator]`):
`pileups(&mut self) -> Result<Option<(Tid, PileupColumn<'_, U>)>, StreamPileupError>`,
the column borrowing the pileup until the next call. Its store is the window's,
so `column.store()`, `records_overlapping` and `find_record` work as on an
indexed column, over the records of the window.

r[pileup_stream.records+1]
Every record with a reference and a position — mapped or placed-unmapped — MUST
be pushed into a store through the customize value exactly once, in input
order, so `filter_raw`, `compute` and `filter` see the sequence of records an
indexed whole-contig fetch hands them (`r[bam.reader.unmapped_skipped]`). A
record without a reference or without a position MUST be skipped: no pileup
column can contain it, and the indexed reader skips it too.

Only the records the sort check covers (`r[pileup_stream.sorted+1]`) — mapped
ones the customize value kept — open, close or bound a window. A placed-unmapped
record is pushed into whichever store is being filled, whatever its position:
the engine never puts it in a column (`r[pileup.unmapped_excluded]`), so where
it sits changes nothing a column shows, and letting it move a window boundary
would let a record the sort check never saw decide which columns are final.

## Order

r[pileup_stream.sorted+1]
The sort check MUST cover exactly the records htslib's `bam_plp_push` checks,
translated to seqair's hooks. `bam_plp_push` returns before its check for a
record with no reference (`tid < 0`), for an unmapped one (`BAM_FUNMAP`), and
for one past its per-position cap (`maxcnt`); `samtools mpileup`'s read
callback (`bam_plp_auto`'s `func`) has already dropped what its filters
reject, so those never reach it. The rest it checks against the largest
`(tid, pos)` pushed so far and fails with "The input is not sorted" when the
reference id is smaller, or equal with a smaller position. So the check MUST
cover the records that are mapped (flag `0x4` clear), on a reference, and kept
by the customize value — which plays the read callback's part, and whose
`DepthCap` (`r[pileup_stream.depth_limit+1]`) plays `maxcnt`'s — and MUST skip
unplaced records, unmapped ones (placed or not) and the ones the customize
value rejects. A covered record whose `(tid, pos)` is less than the last covered
record's MUST fail with `StreamPileupError::Unsorted`, naming its 1-based record
number in the input and both positions. Skipping it is not an option: windows
already emitted may need it, so every later column could be wrong.

A mapped record on a reference but with no position (`-1`) is checked at
position `-1`, as htslib checks it, so it fails after any covered record on its
reference; when it passes (the first of its reference) it is skipped, since a
store holds placed records only — htslib would pile it up from `-1`. It never
reaches the customize value, which can only be asked through a push.

r[pileup_stream.errors]
Any error — from the stream, a record the store refuses, the sort check, the
store's hard limit (`r[pileup_stream.budget]`), or the reference — MUST end the
pileup: every later call returns
`StreamPileupError::Poisoned`. A pileup that dropped a record and carried on
would yield wrong columns without saying so.

## Windows

The engine needs every record that overlaps a column before it builds that
column. Records arrive sorted by start, so once a record starting past position
`p` has arrived, every record that can reach `p` is already in hand — however
far it extends. Nothing has to be guessed about read length: the next record's
start is the whole of the lookahead.

r[pileup_stream.window+2]
Records are collected into a window: an interval of positions `[start, end]`
and a store. In this section "record" means one the sort check covers
(`r[pileup_stream.records+1]`). A window is opened by a record at `pos`, with
`end = pos + len - 1` (`len` is the window length,
`r[pileup_stream.window_len]`). Each record whose position is at or before
`end`, on the same contig, is pushed into the window's store. The window closes
when a record arrives that starts after `end` or on another contig, when the
store outgrows its budget (`r[pileup_stream.budget]`), or at the end of the
input. Its store then holds every record
that overlaps `[start, end]`, and the engine is run over that interval and
yields its columns. The last window of a contig runs instead to the furthest
`end_pos` of its records, capped at the contig's last position: nothing else
will reach those positions, and none past that furthest end can have a column.
That bound, not the contig's end, MUST be what the window fetches reference
bases for (`r[pileup_stream.reference]`) — the contig's end would make the last
window of a 20 Mb prefix of chr12 read the other 113 Mb of it.

r[pileup_stream.carry]
When a window closes on its contig (not at a contig change or the end), the
records of its store that extend past `end` (`end_pos > end`) MUST be carried
into the next window's store, in store order and before any record of the next
window: copied with `RecordStore::copy_record_from`
(`r[record_store.copy_from]`), not decoded again and not handed to the customize
value again. The next window starts at `end + 1`, so consecutive windows cover
disjoint, adjacent intervals: no column is yielded twice, and a record seen by
several windows is still one alignment per column. When nothing is carried, the
next window starts at the record that opened it, so a gap in coverage costs no
empty windows. Either way the next window's `end` is that record's position
plus `len - 1`. (The record that closes a window is pushed through the customize
value before it is known to close it, since only a kept record can; it is then
moved into the next store, after the carried reads.)

By induction every record overlapping a window is in its store: it either
starts inside the window, or started earlier and overlapped the previous
window's end, so was carried. Its store order is input order, as the
whole-contig store's is, so the column order (`r[pileup.column_record_order]`)
— and with it which reads a depth cap keeps — is the same.

r[pileup_stream.window_len]
The window length MUST default to 64 Ki positions and MAY be set to any
`NonZeroU32`. No column may depend on it. It trades memory against per-window
work: a smaller window holds fewer records, but carries the reads crossing each
boundary once more and builds one engine more.

r[pileup_stream.budget]
A store cannot grow without bound: record indices and slab offsets are `u32`
(`r[record_store.checked_offsets]`), so a deep region — an amplicon, a 64 Ki
window at tens of thousands of reads per position — would otherwise fail with
`DecodeError::SlabOverflow`. A window MUST therefore also close early when,
after a record at `pos > start` is pushed, its store holds more records than
the budget's `records`, or its longest slab (in its own elements: bytes for
names, bases, qualities and aux, operations for CIGARs) is longer than the
budget's `slab_len`. The window then ends at `pos - 1` and the record opens the
next one, after the carried reads. The columns up to `pos - 1` are final:
every later record starts at or after `pos`. No column may depend on the
budget.

The budget MUST default to half the hard limits — 2^31 records and a 2^31
slab length — and MAY be set to any non-zero values. Half, because the reads
carried out of a closing window are at most its whole store, so the next store
starts at no more than the budget and keeps the other half for its own
records. A budget above half the limits gives up that headroom.

A record at the window's own `start` cannot close it, because there is no
column boundary inside one position: every record in the window then starts
at `start` or was carried into it. Such a pile MUST be let grow past the budget,
which is a soft limit, and only when a push or a carry would exceed the hard
limit MUST the pileup fail, with `StreamPileupError::PileTooDeep` naming the
reference and the position. Refusing the pile at the budget instead would
reject inputs the store can hold — the default budget's whole point is
headroom — and a caller who wants a smaller ceiling has `DepthLimit`.

## Mates

r[pileup_stream.mates+1]
Mates MUST be linked within each window's store
(`RecordStore::prepare_for_pileup`, `r[record_store.link_mates+2]`). Both reads
of a pair that overlap each other overlap every column of that overlap, so both
are in the store of every window containing such a column, and
`in_mate_overlap` and `mate_of` agree with the whole-contig pileup. A pair whose
reads do not overlap is linked only when both happen to share a window store:
its `mate_idx` MAY be `None` where the whole-contig pileup has `Some`, and
`fragment_idx` then names the read itself rather than the pair. Whether such a
pair is linked depends on the window length and budget, so a consumer MUST NOT
rely on `mate_idx` or `fragment_idx` for mates that do not overlap; for mates
that overlap they are always set, and that is what deduplicating overlapping
mates needs. `StreamPileup`'s documentation MUST say so.

htslib's mpileup is no different: `overlap_push` pairs a read with its mate only
while both are in its buffer, which is exactly the reads that overlap the
current position.

Input with repeated primary qnames (`r[record_store.link_mates.qname_uniqueness]`)
MAY pair differently from the whole-contig pileup, since "first free pair wins"
then depends on which reads share a store.

## Depth and reference

r[pileup_stream.depth_limit+1]
`DepthLimit::PerColumn(cap)` MUST wrap the customize value in a `DepthCap`, as
the indexed pileup does. The cap sees each contig's records once, in input
order, as a whole-contig fetch hands them, so it keeps the same records as an
indexed pileup over one segment per contig. It is not the same as an indexed
pileup tiled into several segments: there each segment's fetch restarts the
cap (`r[record_store.depth_cap.reset]`) and first sees the reads reaching into
the segment from before it, so near tile boundaries a tiled pileup can keep
different reads. Window boundaries do not restart it — windows are not
fetches. `StreamPileup`'s documentation MUST say so.

r[pileup_stream.reference]
Without a reference every column's reference base MUST be `Base::Unknown`. With
one (`with_reference(IndexedFastaReader)`), each window's engine MUST be given
the bases of exactly the interval it yields columns for, read from the FASTA by
the contig's name; a FASTA error MUST fail with
`StreamPileupError::FastaFetch`.

## Memory

r[pileup_stream.memory]
The pileup MUST hold at most two stores — the window being iterated and the
one being filled — and reuse their capacity and the engine's scratch from one
window to the next (`r[record_store.clear+2]`). A store holds the records
starting in its window plus the ones carried into it, so memory is bounded by
the window length times the coverage, plus the reads spanning a boundary; it
does not grow with the contig.

## Rationale

**Why windows over the unchanged engine.** The engine walks a sorted, linked
store and never looks ahead past the next record's start; it needs no change to
pileup a window, only a region that stops where the records known so far stop
mattering. Growing one engine's store as records arrive instead would need an
append path that keeps the store sorted and its mates and window-query index
valid while columns are borrowed, and a compaction that renumbers the active
set — a second engine, essentially, for what a window boundary does with one
copy per crossing read. Two `pub(crate)` seams were added (the engine's column
step split from building the `PileupColumn`, and the scratch round trip it
already had); neither changes a column.

**Why copy and not re-push.** Re-pushing a crossing read from its raw bytes
would decode it again and run the customize value on it again: a stateful one
(`DepthCap`, a counter) would see it twice, and a whole-contig fetch shows it
once. Copying the decoded record and its extras keeps "each record is pushed
once" exact. The extras are cloned, so the stream needs `U: Clone`.

**Why the window ends at a fixed distance from its first record.** It is the
simplest rule that bounds a store and needs no knowledge of read length. A read
of any length is correct: it is carried as long as it extends past a window's
end. What a long read costs is one copy per window it crosses: a 100 kb
ONT read crosses two 64 Ki windows; an RNA-seq read spanning a 1 Mb intron
(`N`), sixteen.

**Future work (not a rule): grow one store instead of copying.** A read is
copied into every window it crosses, which for a spliced read over a long
intron, or an ultra-long read, is many copies. An engine whose store accepts
appended records while it iterates, and drops records that end before the
current column (renumbering its active set), would decode and store each read
once. It needs the append path to keep the store sorted, its mates linked and
its window-query index valid while columns are borrowed — see the first
paragraph above — so it waits until long reads are shown to cost enough.

**Not supported yet.** The soft-clip overhang (`r[pileup.soft_clip_overhang.activation]`)
reaches `n` positions before a read's start, so it needs `n` positions of
lookahead (close at `pos > end + n`, carry at `end_pos + n > end`); the
column cap `set_max_depth` and a `mutate` hook have no spelling.
