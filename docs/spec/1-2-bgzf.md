# BGZF Block Reader

BGZF (Blocked Gzip Format) is the compression layer used by BAM and other HTS (High-Throughput Sequencing) formats. Unlike plain gzip, which compresses an entire file as one stream, BGZF splits the data into independent blocks of up to 64 KiB uncompressed. Each block is a self-contained gzip member, meaning any block can be decompressed without reading the ones before it. This enables random access: a BAM index can point directly to a specific block offset in the compressed file, and the reader can decompress just that block to find the data it needs.

A BGZF file is simply a concatenation of these gzip blocks, ending with a special empty EOF marker block. Each block has a standard gzip header with an extra field that encodes the block's total compressed size (BSIZE), which the reader uses to locate where the next block begins.

> **Sources:** All rules in this file derive from [SAM1] §4.1 "The BGZF compression format", §4.1 "Random access", and §4.1 "End-of-file marker". See [references.md](references.md).

## Block structure

> *[SAM1] §4.1 "The BGZF compression format" — block layout, extra field, BSIZE*

r[bgzf.magic]
A BGZF block MUST begin with the gzip magic bytes `1f 8b 08 04` (gzip, DEFLATE, FEXTRA flag set).

r[bgzf.bsize]
The extra field MUST contain a `BC` subfield (SI1=0x42, SI2=0x43) with SLEN=2 whose value is the total block size minus one (BSIZE). The total compressed block size is `BSIZE + 1`.

r[bgzf.decompression]
Each block's DEFLATE payload MUST be decompressed independently. The last 8 bytes of the block are the gzip footer: a CRC32 checksum (4 bytes) followed by the uncompressed size ISIZE (4 bytes), both little-endian. Decompression MUST produce exactly ISIZE bytes.

r[bgzf.max_block_size]
BGZF blocks MUST NOT exceed 65536 bytes uncompressed. If the ISIZE footer field claims a larger value, the reader MUST return an error rather than allocating an unbounded buffer.

r[bgzf.eof]
An EOF marker block has ISIZE=0. When encountered, the reader MUST signal end-of-stream.

## Virtual offsets

> *[SAM1] §4.1 "Random access" — virtual file offset definition*

Because BGZF blocks are at known compressed file offsets and have known uncompressed sizes, any byte in the uncompressed stream can be addressed with a **virtual offset**: a packed 64-bit value that combines "which block" and "where within that block." BAM index files (`.bai`) store virtual offsets to point at specific records.

r[bgzf.virtual_offset]
A virtual offset is a 64-bit value where the upper 48 bits encode the compressed block offset in the file and the lower 16 bits encode the byte offset within the uncompressed block.

r[bgzf.seek]
The reader MUST support seeking to an arbitrary virtual offset by seeking the underlying file to the block offset, decompressing that block, and advancing to the within-block offset.

## Reading

BAM records are variable-length and can span block boundaries (a record may start near the end of one block and continue into the next). The reader must handle this transparently.

r[bgzf.read_exact]
The reader MUST support reading an exact number of bytes, transparently crossing block boundaries when the requested data spans multiple blocks.

r[bgzf.read_partial]
The reader MUST support partial reads (up to N bytes) returning the actual count, returning 0 at EOF.

r[bgzf.libdeflate]
Decompression SHOULD use the `libdeflater` crate (Rust bindings to the libdeflate C library) for performance parity with htslib's libdeflate usage. libdeflate provides hardware-accelerated DEFLATE decompression and CRC32 computation.

## Integrity

> *[SAM1] §4.1 "The BGZF compression format" — gzip footer CRC32 and ISIZE fields*

Each gzip block includes a CRC32 checksum of the uncompressed data. Verifying this catches silent data corruption from disk errors, network glitches on cluster storage, or truncated writes.

r[bgzf.crc32]
After decompression, the reader MUST verify the CRC32 checksum of the decompressed data against the expected CRC32 stored in the gzip footer (4 bytes before ISIZE). A mismatch MUST return a `ChecksumMismatch` error.

## Performance

r[bgzf.fast_header]
The reader SHOULD fast-path standard BGZF headers where XLEN=6 and the BC subfield is at the fixed offset (bytes 12–17 of the 18-byte header). All BAM files produced by samtools, htslib, and Picard use this layout. This avoids allocating an extra-fields buffer for the common case. Non-standard layouts MUST fall back to searching the extra fields for the BC subfield.

r[bgzf.resize_uninit]
When resizing buffers that will be immediately and fully overwritten (by `read_exact` or decompression), the reader SHOULD use an uninitialized resize to avoid redundant zero-filling. This applies to the decompressed block buffer (~64 KB per block) and the compressed data buffer.

r[bgzf.block_offset_tracking]
The reader MUST track the current block's compressed file offset for virtual offset calculation. The offset MUST be queried from the stream position before reading each new block, and set directly on seeks.

## Writing

> *[SAM1] §4.1 "The BGZF compression format" — block structure, gzip member format, BC extra field, EOF marker block*

r[bgzf.writer]
The writer MUST accept arbitrary byte sequences and emit valid BGZF blocks. Each block MUST contain a complete gzip member with the `BC` extra subfield, DEFLATE-compressed payload, CRC32 checksum, and ISIZE footer.

r[bgzf.writer.buffer]
The writer MUST accumulate uncompressed data in an internal buffer (up to 64 KB). When the buffer is full or `flush()` is called, the buffer MUST be compressed into a BGZF block and written to the underlying stream.

r[bgzf.writer.compression]
Compression MUST use the `libdeflater` crate (matching the reader's decompression backend) with configurable compression level. The default compression level SHOULD be 6 (matching htslib's default).

r[bgzf.writer.eof_marker]
`finish()` MUST write the standard 28-byte BGZF EOF marker block after flushing any remaining buffered data. The EOF marker is a valid gzip member with ISIZE=0.

r[bgzf.writer.virtual_offset]
The writer MUST track virtual offsets. After each block is written, the writer MUST record the compressed file offset of that block. A `virtual_offset()` method MUST return the current write position as a `VirtualOffset` (block offset + within-block offset). The within-block offset MUST be strictly less than 65536; converting buffer length to u16 MUST use checked conversion to prevent silent truncation when the buffer is exactly full.

r[bgzf.writer.flush_if_needed]
The writer MUST provide a `flush_if_needed(upcoming_bytes)` method that flushes the current block if the upcoming data would exceed the 64 KB uncompressed block limit. This allows callers (e.g., VCF/BCF writers) to keep records from spanning block boundaries when possible, improving seek granularity for index-based random access.

r[bgzf.writer.finish]
`finish()` MUST consume the writer and return the inner `io::Write` stream, allowing the caller to perform additional operations (e.g., syncing, closing). Dropping the writer without calling `finish()` SHOULD NOT silently discard buffered data; the writer SHOULD flush on drop (best-effort, ignoring errors).
