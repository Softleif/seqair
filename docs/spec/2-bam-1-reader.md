# BAM Reader

The indexed BAM reader is the main entry point for loading aligned reads from a BAM file. It ties together three layers: BGZF decompression (reading compressed blocks), BAI index lookup (finding which blocks contain reads for a region), and record decoding (parsing the binary BAM records into usable structs).

The reader's primary operation is `fetch_into`: given a genomic region like `chr19:6,105,700–6,105,800`, it queries the index, bulk-reads the relevant compressed data into memory (see `region_buf.md`), decompresses it, decodes the records, and stores them in a `RecordStore` for downstream processing by the pileup engine.

> **Sources:** [SAM1] §4 "The BAM Format Specification" generally — BGZF, BAM header, record layout, and BAI index together define the reader's behaviour. The forking design and ChunkCache are seqair-specific optimisations with no upstream spec counterpart. See [references.md](references.md).

## Opening

r[bam.reader.open]
The reader MUST open a BAM file by path, parse its header, and locate the corresponding `.bai` index file (trying both `<path>.bai` and `<path_without_ext>.bai`).

r[bam.reader.header_access]
The reader MUST expose the parsed `BamHeader` for contig name/tid/length lookups.

## Region fetching

r[bam.reader.fetch_into+2]
`fetch_into(region, store)` MUST clear the store, query the BAI index for chunks overlapping the region, seek through those chunks via BGZF, decode all records, and add records that actually overlap the region to the store (RecordStore).

r[bam.reader.chunk_batching]
When the merged compressed byte range of all chunks exceeds `MAX_REGION_BYTES` (256 MiB), `fetch_into` MUST partition chunks into batches where each batch fits within the limit, loading each batch with its own `RegionBuf`. This occurs on very large BAM files (>100 GB) where BAI bin structure produces huge level-5 bins spanning hundreds of megabytes of compressed data. Chunks MUST be added greedily in order; records from all batches MUST be collected into the same store. Peak memory MUST NOT exceed `MAX_REGION_BYTES` per batch (the previous `RegionBuf` is dropped before loading the next).

r[bam.reader.overlap_filter]
A record overlaps a region if `record.pos <= region.end` AND `record.end_pos >= region.start`. Records outside the region that appear in index chunks MUST be filtered out. This is necessary because the BAM index is approximate — chunks may contain records slightly outside the queried region.

r[bam.reader.sorted_order+2]
Records in the store MUST be in coordinate-sorted order (by pos, then by end_pos as tiebreaker) after `fetch_into`, matching the order a sorted BAM file provides.

## Edge cases

r[bam.reader.unmapped_skipped+2]
Unmapped reads (flag 0x4) MUST be excluded from the store during `fetch_into`. htslib's pileup engine rejects them in `bam_plp_push`, so including them would inflate depth counts.

r[bam.reader.secondary_supplementary_included+2]
Secondary (0x100) and supplementary (0x800) reads MUST be included in the store by default. The caller's pileup filter is responsible for excluding them if desired, matching htslib's behavior.

## Forking (thread-safe shared state)

Rastair processes regions in parallel via rayon. Each worker thread needs its own `IndexedBamReader` because file handles have mutable seek state. Without forking, every thread independently opens the BAM file, re-reads the entire BAI index from disk, and re-parses the BAM header — wasting both I/O bandwidth and memory. On HPC clusters with NFS/Lustre, redundant index reads are especially costly.

The fork pattern separates the reader into shared immutable data (`BamIndex`, `BamHeader`) and per-thread mutable state (file handles, decompression buffers). A prototype reader is opened once on the main thread; worker threads call `fork()` to get a lightweight copy that shares the index and header via `Arc` but owns fresh file handles.

### Shared state

r[bam.reader.shared_state]
The BAM index and header MUST be stored behind `Arc` so that forked readers share a single copy. These structures are read-only after initial parsing and are accessed concurrently without synchronization.

### Fork operation

r[bam.reader.fork]
`fork()` MUST produce a new reader that:
1. Shares the same `Arc` holding `BamIndex` and `BamHeader` as the source reader (no re-parsing, no additional memory).
2. Opens a fresh raw `File` handle for bulk `RegionBuf` reads. No `BgzfReader` is needed — `fetch_into` only uses `RegionBuf`, and the `BgzfReader` used during `open()` for header parsing is not retained.
3. Allocates its own `ChunkCache` (starts empty, populated on first region fetch for a tid).

r[bam.reader.fork_equivalence]
A forked reader MUST produce identical `fetch_into` results as independently opening the same BAM file. The fork is purely an optimization — it MUST NOT change observable behavior.

r[bam.reader.fork_independence]
Forked readers MUST be fully independent for mutable operations. Seeking or fetching on one reader MUST NOT affect any other reader (original or forked). Each reader owns its own file descriptors and buffer state.

### Edge cases

r[bam.reader.fork_arc_identity]
`Arc::ptr_eq` on the shared state of the original and forked reader MUST return `true`, confirming they share the same allocation rather than having made a deep copy.

r[bam.reader.fork_concurrent]
Multiple forks from the same prototype MUST be safe to use concurrently from different threads. Since the shared state is immutable and each fork owns its file handles, no synchronization is needed beyond the `Arc` reference count.

## Coordinate validation

r[bam.reader.coordinate_overflow]
`fetch_into` MUST validate that `tid` fits in `i32` and that `start`/`end` fit in `i64` before using them. If any value exceeds the signed maximum, `fetch_into` MUST return a `BamError::CoordinateOverflow` rather than silently truncating via `as` casts.

## Error propagation

r[bam.reader.propagate_errors]
`fetch_into` (and `ChunkCache::load`) MUST propagate I/O and decompression errors from the BGZF layer. Only `BgzfError::UnexpectedEof` (signalling the end of a chunk's data) should silently break the read loop. All other errors — CRC failures, truncated blocks, decompression failures, invalid magic, et.c — MUST be returned to the caller.
