# FASTA Reader

The indexed FASTA reader provides random access to reference genome sequences. It reads `.fai` index files to locate sequence data within plain or bgzip-compressed FASTA files, enabling efficient region queries without scanning the entire file.

Rastair fetches reference sequences once per segment (~100 KB) for each worker thread. The forkable design shares the parsed index across threads while giving each thread its own file handle, matching the `IndexedBamReader` pattern.

> **Sources:** The FASTA format and FAI index are not covered by any hts-specs document. The FAI format is defined by samtools (`samtools faidx`); the byte offset formula and field layout are described in the `samtools faidx` man page and implemented in [htslib]. The GZI index format is a htslib convention for bgzip-compressed random-access files, implemented in [htslib] `bgzf.c`. BGZF decompression follows [SAM1] §4.1. The forking and buffer-reuse design are seqair-specific. See [References](./99-references.md).

## FAI Index

The FAI index (`.fai` file) is a tab-delimited text file with one line per sequence. It enables O(1) lookup of any base position within the FASTA file by recording where each sequence's data starts and how lines are formatted.

### Parsing

r[fasta.index.parse]
The FAI parser MUST read a `.fai` file and produce an index with one entry per sequence. Each entry MUST contain: sequence name (`String`), length (number of bases, `u64`), byte offset of the first base in the file (`u64`), bases per line (`u64`), and bytes per line including line terminator (`u64`). The parser MUST NOT trim whitespace from field values — tabs are the only delimiter, and sequence names may contain spaces.

r[fasta.index.fields]
Each line MUST have exactly 5 tab-separated fields: NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH. Lines with fewer or more fields MUST produce an error identifying the malformed line. Empty lines MUST be skipped.

r[fasta.index.validation]
The parser MUST reject entries where `linewidth < linebases` (line terminator would be negative), `linebases == 0` (would cause division by zero in offset calculations), or `length == 0` (empty sequences are not useful for region queries). Note: `linewidth == linebases` is valid — it means lines have no trailing newline character, which occurs when a sequence is stored on a single line without a mid-sequence newline.

r[fasta.index.name_lookup]
The index MUST support O(1) lookup of a sequence entry by name. Duplicate sequence names MUST produce an error at parse time.

### Location

r[fasta.index.location]
For a FASTA file at path `P`, the reader MUST look for the index at `P.fai` (i.e., appending `.fai` to the full path including any `.gz` extension). This matches samtools/htslib convention: `ref.fa` → `ref.fa.fai`, `ref.fa.gz` → `ref.fa.gz.fai`.

r[fasta.index.missing]
When the FAI index file is not found, the reader MUST produce a clear error message that names the expected index path and suggests running `samtools faidx <fasta_path>` to create it. The reader MUST NOT attempt to build the index itself. Rationale: auto-indexing hides user errors (wrong FASTA passed silently), fails on read-only shared reference directories (common on HPC clusters), introduces non-deterministic file creation and race conditions, and takes 10–30s for a human genome which makes the tool appear to hang with no feedback.

### Byte offset calculation

r[fasta.index.offset_calculation]
To find the byte offset of base at 0-based position `pos` within a sequence, the reader MUST compute: `offset + (pos / linebases) * linewidth + (pos % linebases)`, where `offset`, `linebases`, and `linewidth` come from the FAI entry. This accounts for newline characters inserted at line boundaries.

The offset formula is correct for all valid positions including those on the last line of a sequence, even when the last line is shorter than `linebases`. This works because integer division gives the correct line number regardless of last-line length. The raw read range `[byte_offset(start), byte_offset(last_pos) + 1)` may include newline bytes from line boundaries; these are stripped after reading.

## Sequence fetching

r[fasta.fetch.coordinates]
`fetch_seq(name, start, stop)` MUST use 0-based half-open coordinates `[start, stop)`. A request for `(name, 0, 4)` returns the first 4 bases of the named sequence.

r[fasta.fetch.bounds_check]
The reader MUST verify that `start < stop` (zero-length and reversed ranges are invalid) and `stop <= sequence_length`. Out-of-bounds requests MUST produce an error naming the sequence, the requested range, and the actual sequence length.

r[fasta.fetch.uppercase]
Returned sequence bytes MUST be uppercased. FASTA files may contain lowercase bases (indicating soft-masking), but Seqair always works with uppercase bases.

r[fasta.fetch.newline_stripping]
When reading from the file, the reader MUST strip line terminators (LF `\n` and CR `\r`) from the raw bytes. The returned sequence MUST contain only base characters.

r[fasta.fetch.unknown_sequence]
When the requested sequence name is not in the index, the reader MUST produce an error listing the requested name. It SHOULD include available sequence names if the index is small (fewer than 20 entries) to help diagnose typos.

### Raw byte return type

r[fasta.fetch.raw_bytes]
`fetch_seq` and `fetch_seq_into` MUST return raw uppercase ASCII bytes (`Vec<u8>` / `&mut Vec<u8>`), not `Base` values. CRAM reference MD5 verification computes checksums on the exact byte sequence from the FASTA file — converting to `Base` and back is lossy because IUPAC ambiguity codes (M, R, S, etc.) collapse to `Unknown`. Any u8→Base conversion MUST happen downstream, after CRAM verification is complete.

### Buffer reuse

r[fasta.fetch.buffer_reuse]
`fetch_seq_into(name, start, stop, buf)` MUST accept a caller-provided `&mut Vec<u8>` buffer, clear it, and write the fetched sequence into it. This avoids per-call allocation when fetching many segments. The `fetch_seq` convenience method MAY allocate a new `Vec<u8>` internally.

## Plain FASTA reading

r[fasta.plain.read]
For uncompressed FASTA files, the reader MUST seek to the computed byte offset and read the required byte range in a single I/O operation, then strip newlines and uppercase the result.

r[fasta.plain.read_size]
The raw read size MUST account for embedded newlines. For a fetch of `n` bases starting at position `pos`, the number of raw bytes to read is: `end_byte_offset - start_byte_offset` where both offsets are computed via `fasta.index.offset_calculation`.

## Bgzip FASTA reading

r[fasta.bgzf.detect]
The reader MUST auto-detect bgzip compression by reading the first 18 bytes of the file and checking for: (1) gzip magic bytes `1f 8b`, (2) DEFLATE compression method `08`, (3) FEXTRA flag `04`, and (4) the `BC` subfield identifier at bytes 12–13. All four checks MUST pass to identify a file as BGZF. Plain gzip files with the FEXTRA flag but without the BC subfield MUST NOT be misidentified as BGZF.

r[fasta.bgzf.gzi_required]
For bgzip FASTA files, a GZI index (`.gzi`) MUST be present alongside the FASTA. The GZI index maps compressed block offsets to uncompressed byte offsets, enabling random access into the compressed data.

r[fasta.bgzf.gzi_missing]
When the GZI index is not found for a bgzip FASTA, the reader MUST produce a clear error naming the expected GZI path and suggesting `samtools faidx <fasta_path>` (which creates both `.fai` and `.gzi`).

### GZI index format

r[fasta.gzi.parse]
The GZI parser MUST read a binary file containing: an 8-byte little-endian entry count (`u64`), followed by that many pairs of 8-byte little-endian values: `(compressed_offset: u64, uncompressed_offset: u64)`. Each pair records where a BGZF block starts in the compressed file and what uncompressed byte offset that corresponds to. The entry count multiplication MUST use checked arithmetic to prevent overflow on corrupted files.

r[fasta.gzi.sorted]
GZI entries MUST be sorted in strictly increasing order by `uncompressed_offset`. The parser MUST reject files with unsorted or duplicate entries, as `translate()` relies on binary search.

r[fasta.gzi.translate]
Given an uncompressed byte offset, the GZI index MUST be able to find the BGZF block that contains it: find the last entry whose `uncompressed_offset <= target`, then the compressed block starts at that entry's `compressed_offset`, and the within-block offset is `target - uncompressed_offset`.

r[fasta.gzi.within_block_bounds]
The within-block offset MUST fit in `u16` (0–65535), matching the BGZF specification's maximum uncompressed block size of 65,536 bytes. If the computed within-block offset exceeds `u16::MAX`, the reader MUST return an error indicating the GZI index may be corrupt or missing entries. The reader MUST NOT silently truncate the offset.

### Decompression

r[fasta.bgzf.decompress]
The reader MUST decompress BGZF blocks using the same libdeflater-based decompression as the BAM reader, including CRC32 verification.

r[fasta.bgzf.sequential_read]
For a fetch spanning multiple BGZF blocks, the reader MUST read and decompress blocks sequentially from the starting block until all requested bytes are obtained. Since FASTA fetches are contiguous byte ranges (unlike BAM which has disjoint index chunks), a simple sequential scan suffices.

## Forking

r[fasta.fork.shared_state]
The FAI index and GZI index (if applicable) MUST be stored behind `Arc` so forked readers share a single parsed copy. These structures are read-only after initial parsing.

r[fasta.fork.operation]
`fork()` MUST produce a new reader that shares the same `Arc` as the source (no re-parsing) and opens a fresh `File` handle. Each fork MUST own its own read buffer.

r[fasta.fork.equivalence]
A forked reader MUST produce byte-identical results to independently opening the same FASTA file.

r[fasta.fork.independence]
Forked readers MUST be fully independent for mutable operations. Seeking on one MUST NOT affect any other.

## Error handling

r[fasta.errors]
All FASTA-specific errors MUST contain context (file paths, sequence names, positions) sufficient to diagnose the problem without a debugger. The implementation MUST NOT use `panic!`, `unwrap()`, `expect()`, or `unreachable!()` — all error paths MUST propagate errors via `Result`.
