# BGZF Block Reader

BGZF (Blocked Gzip Format) is the compression layer used by BAM and other HTS (High-Throughput Sequencing) formats. Unlike plain gzip, which compresses an entire file as one stream, BGZF splits the data into independent blocks of up to 64 KiB uncompressed. Each block is a self-contained gzip member, meaning any block can be decompressed without reading the ones before it. This enables random access: a BAM index can point directly to a specific block offset in the compressed file, and the reader can decompress just that block to find the data it needs.

A BGZF file is simply a concatenation of these gzip blocks, ending with a special empty EOF marker block. Each block has a standard gzip header with an extra field that encodes the block's total compressed size (BSIZE), which the reader uses to locate where the next block begins.

> **Sources:** All rules in this file derive from [SAM1] §4.1 "The BGZF compression format", §4.1 "Random access", and §4.1 "End-of-file marker". See [References](./99-references.md).

## Block structure

> _[SAM1] §4.1 "The BGZF compression format" — block layout, extra field, BSIZE_

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

> _[SAM1] §4.1 "Random access" — virtual file offset definition_

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

> _[SAM1] §4.1 "The BGZF compression format" — gzip footer CRC32 and ISIZE fields_

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

> _[SAM1] §4.1 "The BGZF compression format" — block structure, gzip member format, BC extra field, EOF marker block_

r[bgzf.writer]
The writer MUST accept arbitrary byte sequences and emit valid BGZF blocks. Each block MUST contain a complete gzip member with the `BC` extra subfield, DEFLATE-compressed payload, CRC32 checksum, and ISIZE footer.

r[bgzf.writer.buffer]
The writer MUST accumulate uncompressed data in an internal buffer (up to the block size of `r[bgzf.writer.block_size]`). When the buffer is full or `flush()` is called, the buffer MUST be compressed into a BGZF block and written to the underlying stream. A write that fills the buffer exactly MUST flush before it returns, so the buffer is never observably full: 65536 is not a within-block offset, and the stream position after a block's last byte is the *next* block's `(offset, 0)`.

r[bgzf.writer.block_size]
A block MUST hold at most 65280 (`0xff00`) uncompressed bytes — htslib's `BGZF_BLOCK_SIZE`. The whole block, gzip header and footer included, must fit in 65536 bytes (BSIZE is a u16 of size − 1), and data that does not compress is *stored*: 65536 bytes would need 65546 bytes of DEFLATE plus 26 of framing, which does not fit, so every level-0 file and every block of random bytes would fail to write. 65280 leaves room for the worst case; it is the limit htslib writes with. Should a block still not fit, the writer MUST fail with `BgzfError::BlockTooLarge` carrying the block's size — a writer-side error, never a reader-side one such as `CorruptHeader`.

r[bgzf.writer.compression]
Compression MUST use the `libdeflater` crate (matching the reader's decompression backend) with configurable compression level. The default compression level SHOULD be 6 (matching htslib's default).

r[bgzf.writer.single_write]
Each block SHOULD be assembled contiguously — header, DEFLATE payload, footer — in one buffer that is allocated once, at the largest size a block can compress to, and never zero-filled again, and it SHOULD reach the inner stream as a single `write_all`. The payload is compressed straight into its place after the header. An unbuffered sink (rastair hands the VCF writer a `Box<dyn Write>` over a file) otherwise pays four `write` calls per block, and a per-block zero fill of the compressed buffer is 64 KiB of memset that the compressor overwrites anyway.

r[bgzf.writer.parallel]
The BAM and VCF/BCF writers MAY compress BGZF blocks on a pool of worker threads (`compression_threads(n)`, `n = 0` meaning the calling thread). The calling thread keeps producing uncompressed data, hands each full block to a worker, and writes finished blocks to the inner stream strictly in block order. At most a small multiple of `n` blocks MAY be outstanding, so memory stays bounded when the sink is slower than the workers. A block that fails to compress or write poisons the stream: every later write, flush or finish MUST return an error rather than wait for that block, and dropping the writer MUST NOT block on it. Worker threads MUST be joined by `finish` and by drop.

r[bgzf.writer.parallel.identical_output]
Parallel compression MUST produce the same bytes as the single-threaded writer at the same level. Block boundaries depend only on the uncompressed data and the `flush_if_needed` calls, never on compressed sizes or timing, and every worker compresses at the writer's level with the same libdeflate settings.

r[bgzf.writer.parallel.index_offsets]
A block's file offset is only known once every earlier block is compressed, so index co-production MUST NOT wait for it per record. The parallel writer hands out *index offsets* instead — the current block's sequence number where a virtual offset has the block's file offset, `(block_number << 16) | within_block` — which order exactly as the virtual offsets do. The writer MUST record every written block's file offset, and `finish` MUST, once all blocks are written, translate every offset the index builder holds to its virtual offset (the pseudo-bin's mapped/unmapped counts are not offsets and MUST be left alone). The resulting index MUST be byte-identical to the one the single-threaded writer co-produces.

r[bgzf.writer.eof_marker]
`finish()` MUST write the standard 28-byte BGZF EOF marker block after flushing any remaining buffered data. The EOF marker is a valid gzip member with ISIZE=0.

r[bgzf.writer.virtual_offset]
The writer MUST track virtual offsets. After each block is written, the writer MUST record the compressed file offset of that block. A `virtual_offset()` method MUST return the current write position as a `VirtualOffset` (block offset + within-block offset). The within-block offset MUST be strictly less than 65536; converting buffer length to u16 MUST use checked conversion, never a clamp, because clamping a full buffer to 65535 names a byte *inside* the record that just ended, and an index built from such an offset seeks readers into the middle of a record (see `r[bgzf.writer.buffer]`).

r[bgzf.writer.flush_if_needed]
The writer MUST provide a `flush_if_needed(upcoming_bytes)` method that flushes the current block if the upcoming data would exceed the uncompressed block limit (`r[bgzf.writer.block_size]`). This allows callers (e.g., VCF/BCF writers) to keep records from spanning block boundaries when possible, improving seek granularity for index-based random access.

r[bgzf.writer.finish]
`finish()` MUST consume the writer and return the inner `io::Write` stream, allowing the caller to perform additional operations (e.g., syncing, closing). Calling `finish()` or writing to a writer that has already been finished MUST return `BgzfError::AlreadyFinished` (a typed error variant, never `io::Error::other`). Dropping the writer without calling `finish()` SHOULD flush on drop (best-effort, logging failures with `warn!`).
