# Region Buffer

> **Sources:** The RegionBuf design is seqair-specific (bulk I/O optimisation for cluster storage). BGZF block decompression follows [SAM1] §4.1. Virtual offsets are defined in [SAM1] §4.1 "Random access". Index chunks come from [SAM1] §5.2 (BAI) or [TABIX] (TBI). See [References](./99-references.md).

## Background

When reading a BGZF-compressed file (BAM or bgzf-compressed SAM), the reader doesn't decompress the entire file. Instead, it uses an index (`.bai` for BAM, `.tbi` for SAM.gz) to find which compressed byte ranges — called **chunks** — contain reads overlapping the requested genomic region. Each chunk is a pair of virtual offsets (see `bgzf.md`) that point to a start and end position in the compressed file.

A typical region query returns a handful of chunks. The standard approach is to seek to each chunk's start, then read and decompress BGZF blocks one at a time until reaching the chunk's end. Each block read involves a header read (~18 bytes) followed by a compressed data read (~20–40 KB), then decompression.

On local SSDs this works fine — seek + read latency is microseconds. On **cluster storage** (NFS, Lustre, GPFS), each I/O operation can take milliseconds due to network round-trips. A region with 200 BGZF blocks means 200+ network round-trips, which dominates processing time.

## Approach

`RegionBuf` solves this by **pre-fetching all compressed bytes for a region into RAM in one shot**, then decompressing from the in-memory buffer. The I/O pattern changes from "hundreds of small reads" to "one large sequential read" per region. Since Rastair processes BAM segments in parallel (one thread per segment), each thread creates its own `RegionBuf` with no contention.

The flow:

1. Query the BAM index → get a list of chunks (compressed byte ranges)
2. Merge overlapping/adjacent chunks into larger contiguous ranges
3. Read each merged range with one `seek` + `read` into a `Vec<u8>`
4. Decompress BGZF blocks from the in-memory buffer (zero network I/O)

## Chunk merging

The BAM index may return multiple chunks whose compressed byte ranges overlap or are adjacent. Reading them separately would cause redundant I/O and seeking.

r[region_buf.merge_chunks]
Before reading, index chunks MUST be merged into non-overlapping byte ranges sorted by file offset. Overlapping or adjacent chunks MUST be coalesced into a single range. The end of each range MUST extend past the last chunk's end block offset by at least one maximum BGZF block size (64 KiB) to ensure the final block is fully included in the buffer.

r[region_buf.no_bin0]
`RegionBuf` MUST NOT include bin 0 chunks in its bulk read. Bin 0 chunks are at distant file offsets and contain very few reads; including them would cause an expensive additional seek. Bin 0 records are handled separately via the bin 0 cache (see `bam_index.md`).

## Loading

r[region_buf.load]
`RegionBuf::load` MUST accept a seekable reader and a slice of index chunks. It MUST merge chunks, then perform one `seek` + `read_exact` per merged range to load all compressed bytes into a contiguous in-memory buffer. On typical BAM files this results in a single read of a few hundred KB to a few MB per region.

r[region_buf.max_region_bytes]
`RegionBuf::load` MUST reject chunk sets whose total merged compressed byte range exceeds 256 MiB (`MAX_REGION_BYTES`), returning a `RegionTooLarge` error. This guards against OOM from corrupt indexes or extremely large query regions. The caller (`fetch_into`) is responsible for partitioning chunks into batches that each fit within this limit (see `r[bam.reader.chunk_batching]`).

r[region_buf.empty]
When given zero chunks, `load` MUST return an empty buffer that immediately signals EOF on any read.

## Seeking

After loading, the caller iterates through the original (unmerged) chunks. For each chunk, it seeks to the chunk's start virtual offset within the pre-loaded buffer and reads records until reaching the chunk's end.

r[region_buf.seek_virtual]
The region buffer MUST support seeking to a virtual offset within the pre-loaded data. Since disjoint ranges are concatenated, the buffer MUST maintain a range map that translates file offsets to buffer positions. A seek to a file offset that falls within any loaded range MUST succeed; a seek to an offset in a gap between ranges or before/after the loaded data MUST return an error.

## Block decompression

The in-memory buffer contains raw BGZF-compressed bytes, exactly as they appear on disk. Decompression follows the same block-by-block process as the file-based reader, but reads from RAM instead of the filesystem.

r[region_buf.decompress]
Block decompression from the in-memory buffer MUST follow the same BGZF format rules as the file-based reader: magic validation, BSIZE extraction, DEFLATE decompression, and CRC32 verification.

r[region_buf.fast_header]
The fast-path header parsing (XLEN=6, BC at fixed offset) MUST be used when applicable, with fallback to searching extra fields.

## Reading

r[region_buf.read_exact]
The region buffer MUST support reading an exact number of bytes, transparently crossing block boundaries, with the same API as `BgzfReader::read_exact_into`.

r[region_buf.drop_no_panic]
The `Drop` implementation MUST never panic. Gap-size calculations between ranges MUST use saturating arithmetic to handle overlapping or malformed range boundaries safely.

r[region_buf.virtual_offset]
The region buffer MUST track virtual offsets correctly: the block offset corresponds to the file position of the current block (not the buffer position), so that virtual offset comparisons with the BAM index chunk boundaries work correctly.

## Integration

r[region_buf.fetch_into+2]
Both `IndexedBamReader::fetch_into` and `IndexedSamReader::fetch_into` MUST use `RegionBuf` to bulk-read compressed data before record decoding. The file-based `BgzfReader` is retained for header parsing and other sequential access at file open time.

r[region_buf.not_cram]
CRAM does NOT use `RegionBuf`. CRAM containers are not BGZF-compressed — the reader seeks to container byte offsets directly via the CRAI index and reads container data with plain `File::read_exact`. Block-level decompression (gzip, rANS, etc.) happens after reading the container into memory.
