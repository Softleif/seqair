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
index the window's store (`r[pileup_stream.window]`), not a whole-contig one,
and `mate_idx` is narrower (`r[pileup_stream.mates]`).

r[pileup_stream.lending]
Iteration MUST be lending, as the engine's is (`r[pileup.lending_iterator]`):
`pileups(&mut self) -> Result<Option<(Tid, PileupColumn<'_, U>)>, StreamPileupError>`,
the column borrowing the pileup until the next call. Its store is the window's,
so `column.store()`, `records_overlapping` and `find_record` work as on an
indexed column, over the records of the window.

r[pileup_stream.records]
Every record with a reference and a position — mapped or placed-unmapped — MUST
be pushed into a store through the customize value exactly once, in input
order, so `filter_raw`, `compute` and `filter` see the sequence of records an
indexed whole-contig fetch hands them (`r[bam.reader.unmapped_skipped]`). A
record without a reference or without a position MUST be skipped: no pileup
column can contain it, and the indexed reader skips it too.

## Order

r[pileup_stream.sorted]
Records with a reference and a position MUST arrive in ascending
`(tid, pos)` order — contigs in header order, positions ascending within one.
A record that precedes the previous such record MUST fail with
`StreamPileupError::Unsorted`, naming its 1-based record number in the input and
both positions, as htslib fails with "The input is not sorted". Skipping it is
not an option: windows already emitted may need it, so every later column could
be wrong. Unplaced records are exempt (they sort last but take no part).

This is stricter than htslib in one respect: htslib checks only mapped records,
seqair also checks placed-unmapped ones, because they enter the store and the
customize value (`r[pileup_stream.records]`). A coordinate-sorted file places
them in order either way.

r[pileup_stream.errors]
Any error — from the stream, a record the store refuses, the sort check, or the
reference — MUST end the pileup: every later call returns
`StreamPileupError::Poisoned`. A pileup that dropped a record and carried on
would yield wrong columns without saying so.

## Windows

The engine needs every record that overlaps a column before it builds that
column. Records arrive sorted by start, so once a record starting past position
`p` has arrived, every record that can reach `p` is already in hand — however
far it extends. Nothing has to be guessed about read length: the next record's
start is the whole of the lookahead.

r[pileup_stream.window]
Records are collected into a window: an interval of positions `[start, end]`
and a store. A window is opened by a record at `pos`, with `end = pos + len - 1`
(`len` is the window length, `r[pileup_stream.window_len]`). Each record whose
position is at or before `end`, on the same contig, is pushed into the window's
store. The window closes when a record arrives that starts after `end` or on
another contig, or at the end of the input. Its store then holds every record
that overlaps `[start, end]`, and the engine is run over that interval and
yields its columns — over `[start, last]` instead, `last` the contig's last
position, when the window is the last of its contig, since nothing else will
reach those positions.

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
plus `len - 1`.

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

## Mates

r[pileup_stream.mates]
Mates MUST be linked within each window's store
(`RecordStore::prepare_for_pileup`, `r[record_store.link_mates+2]`). Both reads
of a pair that overlap each other overlap every column of that overlap, so both
are in the store of every window containing such a column, and
`in_mate_overlap` and `mate_of` agree with the whole-contig pileup. A pair whose
reads do not overlap is linked only when both happen to share a window store:
its `mate_idx` MAY be `None` where the whole-contig pileup has `Some`. No column
shows the difference except through `mate_idx` and `fragment_idx`, and a
consumer deduplicating overlapping mates never needs the far mate.

htslib's mpileup is no different: `overlap_push` pairs a read with its mate only
while both are in its buffer, which is exactly the reads that overlap the
current position.

Input with repeated primary qnames (`r[record_store.link_mates.qname_uniqueness]`)
MAY pair differently from the whole-contig pileup, since "first free pair wins"
then depends on which reads share a store.

## Depth and reference

r[pileup_stream.depth_limit]
`DepthLimit::PerColumn(cap)` MUST wrap the customize value in a `DepthCap`, as
the indexed pileup does. The cap sees each contig's records once, in input
order, as a whole-contig fetch hands them, so it keeps the same records. It is
not the same as an indexed pileup tiled into several segments, whose cap
restarts at each segment.

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

**Not supported yet.** The soft-clip overhang (`r[pileup.soft_clip_overhang.activation]`)
reaches `n` positions before a read's start, so it needs `n` positions of
lookahead (close at `pos > end + n`, carry at `end_pos + n > end`); the
column cap `set_max_depth` and a `mutate` hook have no spelling. CRAM is not
streamable (`r[unified.stream.detect]`).
