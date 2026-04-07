# BAM Index (BAI/CSI)

A BAM file can be hundreds of gigabytes. Without an index, finding all reads overlapping a small genomic region (e.g. chr19:6,105,700–6,105,800) would require decompressing and scanning the entire file. The BAM index (`.bai` file) solves this by recording, for each genomic region, which compressed byte ranges in the BAM file contain relevant reads.

The index uses a **binning scheme**: the genome is divided into hierarchical bins of decreasing size (the whole reference → 512 Mbp → 64 Mbp → ... → 16 KiB leaf bins). Each bin stores a list of **chunks** — pairs of virtual offsets (see `bgzf.md`) marking the start and end of compressed byte ranges in the BAM file that contain reads falling into that bin. A **linear index** provides an additional optimization: for each 16 KiB window, it stores the minimum virtual offset of any read starting in that window, allowing the reader to skip chunks that are entirely before the query start.

> **Sources:** [SAM1] §5 "Indexing BAM" — binning algorithm; §5.1.1 "Basic binning index" — bin hierarchy and sizes; §5.2 "The BAI index format for BAM files" — binary format; §5.3 "C source code for computing bin number and overlapping bins" — `reg2bin`/`reg2bins`. [CSI] — for CSI format. See [References](./99-references.md).

## BAI format

> _[SAM1] §5.2 "The BAI index format for BAM files" — magic, n_ref, n_bin, chunk virtual offsets, n_intv linear index_

r[bam.index.bai_magic]
A BAI file MUST begin with the magic bytes `BAI\1` (0x42, 0x41, 0x49, 0x01). The reader MUST reject files that do not match.

r[bam.index.bai_parse]
The reader MUST parse the BAI binary format: n_ref references, each containing bins with chunks (pairs of virtual offsets) and a linear index of 16 KiB window offsets.

r[bam.index.bgzf_block_validation]
When decompressing BGZF blocks (e.g. for tabix), the reader MUST validate that bsize is large enough for the block header and footer before indexing into the block. A bsize that is too small MUST return an error rather than panicking.

r[bam.index.region_query]
Given a region (tid, start, end), the index MUST return a list of BGZF chunk pairs `(begin, end)` that contain all reads overlapping the region. The query MUST use the binning scheme to identify candidate bins and the linear index to skip chunks that cannot contain reads starting before the query start.

r[bam.index.bin_calculation]
The BAI binning scheme uses `reg2bin(beg, end)` where bin 0 spans the whole reference, bins 1–8 span 64 Mbp each, and leaf bins span 16 KiB. The implementation MUST correctly compute bins for a region and enumerate all overlapping bins at every level.

## Bin 0

> _[SAM1] §5.1.1 "Basic binning index" — bin 0 spans 512 Mbp, bins 1–8 span 64 Mbp_

Bin 0 is the root of the binning hierarchy, covering the entire reference sequence (0 to 2^29). Each alignment is placed in the smallest bin that fully contains its `[pos, end)` interval. An alignment is placed in bin 0 when its interval straddles a level-1 boundary (multiples of 2^26 = 64 Mbp). This can happen to short reads (e.g. 150 bp) that cross positions like 67,108,864 or 134,217,728. See [SAM spec §5.3](https://samtools.github.io/hts-specs/SAMv1.pdf).

In practice, bin 0 contains very few reads (only those near the ~3–4 level-1 boundaries per chromosome), but its chunks may be located at a very different file offset from the rest of the query's chunks. In a 43 GB BAM, this can mean an extra seek of several GB — a significant cost on network storage.

r[bam.index.bin0_required]
The query MUST always include bin 0 in the candidate bins for any region query. Omitting it would silently drop reads that straddle 64 Mbp boundaries.

r[bam.index.bin0_separate+2]
The index SHOULD provide a `query_split` method that separates distant (level 0–2) chunks from nearby (level 3–5) chunks in the result, returning a `QueryChunks` struct with `nearby` and `distant` fields. This allows callers to handle them differently (e.g. caching, separate I/O). The original `query` method MUST remain available for callers that don't need the separation.

## Higher-level bins

The same distant-chunk problem affects bins at levels 1 (64 Mbp) and 2 (8 Mbp). For a query in chr5:100K–200K, the index includes bin 1 (covering chr5:0–64M) whose chunks span ~2 GB of the BAM file, and bin 9 (covering chr5:0–8M) whose chunks span ~0.27 GB. These chunks are **identical for every region query within the same reference** and contain many records that are outside the query region (~15% of decoded records are discarded by the position filter).

r[bam.index.chunk_separation+2]
The `query_split` method MUST separate chunks into `nearby` (levels 3–5, loaded via RegionBuf per query) and `distant` (levels 0–2, suitable for caching). The threshold (level ≤ 2) is based on profiling data showing that levels 0–2 contribute >90% of the file span in typical region queries.

The reader uses the unified `query()` path, which merges nearby and distant chunks and feeds them all to `RegionBuf` in a single pass. A per-tid `ChunkCache` was previously used to pre-load distant-bin records once per chromosome, but profiling showed the unified path is simpler and sufficient.

## CSI format

> _[CSI] — min_shift, depth, l_aux header; per-bin loffset; parameterised reg2bin/reg2bins_

CSI is an alternative index format that supports larger genomes and arbitrary bin sizes. It uses the same conceptual approach (bins + chunks + linear index) but with a configurable minimum shift and depth.

r[bam.index.csi_support]
CSI index support MAY be added as a follow-up. The same `region_query` interface MUST be used for both BAI and CSI indices.
