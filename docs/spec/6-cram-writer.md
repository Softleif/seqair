# CRAM Writer

The CRAM writer serializes alignment records into CRAM 3.0 against a reference FASTA. It is the inverse of seqair's CRAM reader ([CRAM Reader](./2-cram-1-reader.md), [CRAM Stream](./2-cram-2-stream.md)): every structure it writes is one the reader parses, and every record it writes decodes, through seqair's reader and through htslib's, to the record it was given, except for the normalizations `r[cram_writer.fidelity]` lists. It sits next to the BAM writer ([BAM Writer](./6-bam-writer.md)) and takes the same inputs.

The purpose is smaller output for rastair/postalt than BAM at close to `BamWriter` speed. The writer is reference-based: it stores each aligned read as its differences from the reference, or, for a reference the FASTA lacks, from one built from the reads and stored in the slice.

> **Sources:** [CRAM3] §6 "File definition", §7 "Container structure", §8 "Block structure", §8.4 "Compression header block", §8.5 "Slice header block", §9 "End of file container", §10 "Record structure", §10.6 "Read feature records", §12 "Indexing"; [CRAMcodecs] §2 "rANS 4×8". htslib's `cram/cram_encode.c` (`process_one_read`, `cram_encode_aux`) for the choices the format leaves open. See [References](./99-references.md).

## Scope

r[cram_writer.scope+2]
The writer MUST write CRAM version 3.1, or 3.0 when asked (`r[cram_writer.version+2]`, with its codec profile). Its input MUST be coordinate-sorted (`r[cram_writer.sort_order]`). A slice holds the records of one reference, only unplaced records, or, where references hold few records each, the records of several references (`r[cram_writer.multi_ref]`). A slice on a reference the FASTA lacks embeds a reference built from its records (`r[cram_writer.embedded_reference]`) or is written without one (`r[cram_writer.reference_free+2]`). Read names are always stored (`RN` true). Mates in one slice are attached where the decode rebuilds their mate fields exactly (`r[cram_writer.mates]`); every other record's are stored as they are (`r[cram_writer.record.mate]`).

## Writer

r[cram_writer.builder]
`CramWriterBuilder::to_path(path, header, reference)` and `CramWriterBuilder::to_writer(writer, header, reference)` MUST configure a writer of the records of `header` (a `BamHeader`) against `reference` (an `IndexedFastaReader` the writer owns). `build()` MUST write the file definition (`r[cram_writer.file_definition]`) and the header container (`r[cram_writer.header_container]`) before it returns, so a header the writer cannot write fails the build, not a later record. `to_path` creates (truncating) the file at `path`.

r[cram_writer.inputs]
`write(&OwnedBamRecord)`, `write_raw(&RawBamRecord)` and `write_store_record(&RecordStore, RecordIdx)` MUST encode the same record to the same bytes, whichever form it arrives in. A raw record whose CIGAR is held in a `CG` tag (`r[bam.record.long_cigar]`) is written with its real CIGAR and without the `CG` tag; CRAM stores CIGARs of any length.

r[cram_writer.store]
`write_store(&RecordStore)` MUST write every record of the store, in store order, to the bytes `write_store_record` would write for each in turn: the slices, codec choices, file and index are the same. The store holds each record's fields in contiguous slabs and outlives the call, so the writer need not copy a record into its slice as it arrives: a slice that the store's records complete is encoded (its read features against the reference) and compressed whole, from the store, on one of `compression_threads` threads (on the calling thread for `0`), while the calling thread only checks records and decides where slices end; the slice the store's last records leave open is encoded into the writer as `write_store_record` would, for later records to complete. A record the writer refuses ends the call with its error, the records before it written and it and the ones after it not, as if `write_store_record` had been called for each until the refusal.

r[cram_writer.validate_before_write]
A record the writer refuses MUST leave the writer as it was: it is not part of the output and the next record is written as if it had not been offered. Records are buffered until their slice is complete, so a refusal costs nothing. An error writing a container to the sink (I/O, or a codec failure) MUST poison the writer: every later call returns `CramWriteError::Poisoned`, as `r[bam_writer.error_poisoning]`.

r[cram_writer.finish]
`finish()` MUST write the pending slice, then the EOF container (`r[cram_writer.eof]`), flush, and return the inner writer and, when `write_index(true)` was set, the `CramIndex` of what it wrote (`r[cram_writer.crai]`). The writer does not write the index file; `CramIndex::write_crai` does, where the caller chooses. Dropping the writer without `finish()` leaves a file without its last slice and without the EOF container, which htslib reports as truncated.

## File structure

r[cram_writer.file_definition+2]
The file MUST begin with `CRAM`, major version 3, the minor version `r[cram_writer.version+2]` chose (0 or 1), and a 20-byte file id: the file name of the path, cut to 20 bytes and padded with NUL, as htslib writes it; all NUL for `to_writer`.

r[cram_writer.header_container]
The second container MUST hold one block, of content type `FILE_HEADER`, whose data is the header text's length (i32, little-endian) followed by the text. The container's reference id, alignment start and span, record count, record counter and base count are 0, and it has one landmark, 0.

r[cram_writer.header_text]
CRAM has no binary reference list: the header text is the only place the reader learns the references. If the header's text has no `@SQ` lines, the writer MUST write one per target (`SN`, `LN`) after the `@HD` line (or first, if there is none). If it has `@SQ` lines, their names MUST be the targets', in order, or the build fails with `CramWriteError::HeaderTargetsMismatch`.

r[cram_writer.header_m5]
Every `@SQ` line of a reference the FASTA holds MUST carry `M5`, the lowercase hex MD5 of the reference's bases uppercased, which is what htslib checks a reference against. The writer MUST add `M5` where the line has none. An `LN` that is not the FASTA's length MUST fail the build (`ReferenceLengthMismatch`). An `M5` that is not the FASTA's MUST refuse the first record on that reference and every later one (`ReferenceMd5Mismatch`) without reading the reference again: the file would claim a reference it was not written against. It is checked when the writer first reads the reference rather than at build, so a header that carries every `M5` (one copied from a CRAM) costs no pass over the FASTA. A reference the FASTA lacks keeps its line as it is, and its records are written against a reference built from them (`r[cram_writer.embedded_reference]`) or without one (`r[cram_writer.reference_free+2]`). Computing the missing MD5s reads those references once, at `build()`.

r[cram_writer.header_m5.reuse]
The bases of the first reference whose `M5` `build()` computes MUST be kept for the first slice on that reference, rather than read from the FASTA a second time. They are dropped, or their buffer reused, when a slice on another reference needs its bases.

r[cram_writer.header_m5.threads]
With `compression_threads(n)`, `build()` MUST compute the missing MD5s on `n` threads, each reading the FASTA through a handle of its own (`IndexedFastaReader::fork`): hashing a whole genome's references is most of the time a CRAM of one chromosome takes to write, and htslib does it on one thread. The header is the same whatever the thread count.

r[cram_writer.eof]
The file MUST end with the CRAM 3 EOF container, the 38 bytes htslib's `cram_check_EOF` matches.

## Containers and slices

r[cram_writer.sort_order]
Records MUST arrive ordered by reference id, then position, with every unplaced record (no reference) after all placed ones; unplaced records are in no order among themselves, though they keep the `POS` they have (htslib writes 1 for many). A record that breaks the order MUST be refused with `CramWriteError::Unsorted`. A record with a reference but no position, and a mapped record without a reference, MUST be refused (`InvalidPlacement`).

r[cram_writer.slicing+2]
Each container MUST hold one slice. A slice MUST end when it holds 10,000 records or 5,000,000 bases, htslib's defaults, or under the archive profile (`r[cram_writer.profile]`) 100,000 records or 50,000,000 bases, as htslib's archive profile: fewer, larger slices compress better, at the cost of coarser random access and more memory per slice. It MUST also end when the next record cannot share it: a single-reference slice ends before a record on another reference (or the first unplaced one), a multi-reference slice before a record on a reference the FASTA lacks or the first unplaced one (`r[cram_writer.multi_ref]`). Unplaced records go in slices whose reference id is `-1`.

r[cram_writer.multi_ref]
The writer MUST switch to multi-reference slices as htslib's `cram_put_bam_seq` does: when a slice ends (before a record on another reference, or at its records or bases cap, as htslib), and it and the slice before it each held fewer than `(m + 8) / 4 - 8` records, `m` the records a slice may hold (2,494 by default), the next slice is a multi-reference one. A multi-reference slice holds the records of every reference the FASTA has, in their order; a record on a reference the FASTA lacks ends it and is written in a slice of its own reference, after which the next slice is multi-reference again. Unplaced records end it too and are not written in one, unlike htslib: they keep their `POS` in no order (`r[cram_writer.sort_order]`), and `samtools index` refuses a multi-reference slice whose records of one reference id, unplaced ones included, go backwards. A multi-reference slice that ends holding records of one reference only switches the writer back to one reference per slice. A reference the FASTA lacks never shares a slice: in a slice whose records are stored against the FASTA (`RR` true), htslib's, seqair's and noodles' decoders fetch the reference of every record and fail where the FASTA lacks it, and an embedded reference serves one reference only.

r[cram_writer.multi_ref.records]
In a multi-reference slice each record MUST store its reference id in `RI`, read after `CF`, and be stored against its own reference as in a single-reference slice. `AP` is the 1-based position minus the previous record's in the slice, from 0 for the first, and is negative where the reference changes. The compression header MUST list `RI` only for a multi-reference slice, `EXTERNAL` in a block of its own.

r[cram_writer.container_header+2]
A data container's header MUST give its slice's reference id (`-2` multi-reference), alignment start and span, its record count, the number of records before it in the file (the record counter), its total read length (`bases`), its block count, one landmark (the offset of the slice header block from the container's data start, which is the size of the compression header block), and the CRC32 of the header's bytes before it.

r[cram_writer.slice_extent+2]
A slice's alignment start MUST be the 1-based position of its first record, and its span MUST reach the end of every record's reference span (an unmapped record spans its own position), cut at the reference's end. A slice of unplaced records, and a multi-reference slice, has start 0 and span 0, as the specification recommends.

r[cram_writer.slice_header+2]
The slice header MUST give the slice's reference id (`-2` multi-reference), alignment start and span, record count, record counter, block count (the core block and every external block), the content ids of its external blocks, the content id of its embedded reference's block or `-1`, and the MD5 of the uppercased reference bases over its extent: of the embedded bases for a slice that embeds them, as htslib writes and checks it, and sixteen zero bytes for a slice of unplaced records, a multi-reference slice, or one written without a reference.

r[cram_writer.blocks+2]
Every block MUST be written as `r[cram.block.structure]` reads it: method, content type, content id, compressed and uncompressed size, data, CRC32 of all the bytes before it. Each slice MUST have a core data block (content id 0), empty, as every series is external. The compression header MUST list an encoding for every data series, and the slice MUST have an external block for every `EXTERNAL` series that holds data, and one for each tag it holds: htslib's decoder refuses a slice that lacks the encoding of a series it reads, and noodles one that lacks the block of an `EXTERNAL` series it reads, even when the read is of nothing (`BA` and `QS` for a record without bases), which is why a series with no values is a constant code instead. An empty byte-array series (`IN`, `SC`, `BB`) is read by no record, so its block is left out too.

r[cram_writer.version+2]
`CramWriterBuilder::version` MUST choose the CRAM version, 3.1 by default as htslib's, or 3.0 for readers without 3.1's codecs; the file definition's minor version MUST be the one chosen. CRAM 3.1 carries its codec profile (`r[cram_writer.profile]`), `CramVersion::V3_1(CramProfile::Normal)` by default, so a 3.1 codec cannot be asked of a 3.0 file.

r[cram_writer.profile]
`CramProfile::Normal` MUST try the codecs `r[cram_writer.compression+3]` lists for CRAM 3.1. `CramProfile::Archive` MUST also try, on every data block, the adaptive arithmetic coder (`r[cram_writer.arith]`) under htslib's eight combinations of order, STRIPE, RLE and PACK, and, on the `QS` block, fqzcomp (`r[cram_writer.fqzcomp]`) under each of htscodecs' four strategies, given each record's share of the block and its reverse (`0x10`) and second-read (`0x80`) flags. Codec learning (`r[cram_writer.compression.learning]`) picks among them as among the others, so the archive profile costs its trials and whatever its winners cost. htslib's archive profile also tries bzip2, lzma and the name tokeniser over the arithmetic coder, which the writer does not.

r[cram_writer.compression+3]
Each block's data MUST be stored by one of the version's and profile's codecs (`r[cram_writer.profile]`) or raw: gzip (at the builder's `compression_level`) for both; then rANS 4×8 order 0 and order 1 (`r[cram_writer.rans4x8]`) in CRAM 3.0, or rANS Nx16 under htslib's normal profile (`r[cram_writer.rans_nx16]`) and, for the `RN` block, the name tokeniser (`r[cram_writer.tok3]`) in CRAM 3.1. Which ones are tried for a block is `r[cram_writer.compression.learning]`'s; of those tried, the one giving the fewest bytes is stored, raw included, a tie going to the one listed first in this order: raw, rANS, gzip, the name tokeniser, the arithmetic coder, fqzcomp (roughly fastest to slowest to decode). Blocks shorter than 32 bytes are stored raw. The header block is gzip or raw, whichever is smaller. The choice is per block and changes no decoded byte.

r[cram_writer.compression.learning]
The writer MUST learn which codec suits each block, as htslib's `cram_metrics` does, keyed by the block's content id (a data series or a tag). A *trial* slice tries every codec on every block it compresses (not the raw ones under 32 bytes) and records each one's size; the first slice and every 70th after it (slice index `k` with `k % 70 == 0`, counting the file's data containers from 0) are trials. Each trial adds its sizes to a running score per content id and codec, halving the scores before it adds its own when it is not the first, so older trials weigh less; a codec that did not apply to a block (too short for order 1, say) scores the block's raw size. Every other slice tries, for each block, the codec with the lowest score (raw counting as a codec), and the runner-up too while its score is within 5% of the best's, after the trials up to the one a span before its own span's: slices 1–139 learn from trial 0, slices 141–209 from trials 0 and 70, and so on, so a compression thread waits for a trial only at the start of the file. A block whose content id no trial saw, whose length is far from the trials' average for it (`(len + 1000) / 4 > avg + 1000` or `len + 1000 < (avg + 1000) / 4`, htslib's test), or whose learned codecs all decline it, tries every codec, without changing the scores. `exhaustive_codec_search(true)` (doc-hidden, for tests and measurements) makes every slice a trial.

r[cram_writer.threads+2]
`CramWriterBuilder::compression_threads(n)` MUST assemble and compress containers on `n` worker threads while the calling thread encodes the next slice, and write them in order; `n = 0` (the default) does it on the calling thread. The file and its index MUST be byte for byte those of the single-threaded writer: a slice's content depends on its records only, and its codec choices on its records and the trial slices before it (`r[cram_writer.compression.learning]`), which a worker waits for.

r[cram_writer.rans4x8]
The rANS 4×8 encoder MUST write streams `r[cram.codec.rans4x8]` decodes to the input, order 0 and order 1, with each table's frequencies summing to 4095, as htscodecs normalizes them (its decoder takes 4095 or 4096), every present symbol at least 1. Order 1 needs four bytes at least; a block too small for it is not tried with it. The decoder is its oracle: every input MUST round-trip.

r[cram_writer.rans4x8.fast]
The rANS 4×8 encoder MUST count, build its tables and code without a division per byte: order-0 counts from histograms that do not serialize on repeated symbols; order-1 counts from every byte pair once, in eight sub-tables taking the pairs in turn and placed so that one pair's counts in two of them are not a multiple of 4 KiB apart, with the quarter starts then moved to context 0; order-1 tables for the block's contexts only, kept in per-thread scratch that is never a reason to panic (a fresh one while it is in use or being destroyed) and keeps no byte buffer past 32 MiB; and the renormalization bytes written backwards into a buffer of two a symbol. The reference encoder (each byte's context worked out on its own, counting and coding a byte at a time) stays as the oracle: the encoder MUST write its bytes for every input, order 0 and order 1.

r[cram_writer.rans_nx16]
The rANS Nx16 encoder MUST write streams `r[cram.codec.rans_nx16]` decodes to the input, laid out as htscodecs' `rans_compress_to_4x16` writes them, for every combination of order 0 or 1, 4 or 32 interleaved states, PACK, RLE and STRIPE it is asked for; a separate entry point writes a CAT stream (the bytes stored as they are). The transforms run as htscodecs runs them: PACK first, then RLE on what PACK left, then the entropy coder, and each one that does not apply is dropped from the flags rather than failing: PACK with more than 16 distinct symbols, RLE when no symbol repeats enough to pay or the runs and literals would not save 1%, order 1 on fewer than 8 bytes (order 0 instead), STRIPE on 20 bytes or fewer. An entropy-coded body no smaller than its input is stored as CAT under the same transforms. RLE's metadata is order-0 coded when that is smaller. STRIPE writes four substreams without sizes (`NOSZ`), each the smallest of the order the stream asks for, RLE and PACK if asked, and order 0. Every order-0 table's frequencies MUST sum to a power of two no larger than 4096, the decoder scaling them up to 4096; every order-1 row's MUST sum to a power of two no larger than `1 << shift`, with shift 10 or 12 (the only two htscodecs' decoders build tables for), and every symbol of the alphabet MUST have a row, which a symbol that never precedes another fills with symbol 0, so htscodecs' and seqair's decoders both take the row as a complete table. An order-1 table larger than 1000 bytes is order-0 coded when that saves more than 6 bytes. The writer's choice of combinations for a block MUST be htslib's normal profile (compression level 5): order 0, order 1, order 0 with RLE, order 1 with STRIPE, order 0 with PACK, and order 1 with RLE and PACK, each with 32 states from 50 000 bytes on except STRIPE's. The decoder is its oracle: every input MUST round-trip under every combination; so are htscodecs' decoders, scalar and SIMD, when a build of its `rANS_static4x16pr_test` is given (`HTSCODECS_RANS4X16PR`).

r[cram_writer.arith]
The adaptive arithmetic coder's encoder (CRAM 3.1 method 6) MUST write streams `r[cram.codec.arith]` decodes to the input, byte for byte as htscodecs' `arith_compress_to` writes them, for every combination of order 0 or 1, PACK, RLE and STRIPE it is asked for. The transforms run as htscodecs runs them, each one that does not apply dropped from the flags rather than failing: PACK first (dropped for more than 16 distinct symbols, or no input), then the range coder over what PACK left, with each literal's run length coded after it under RLE (dropped when nothing is left to code), and order 1 coded as order 0 below 8 bytes. A coded body no smaller than its input is stored as CAT, keeping the PACK and RLE bits but not the order. STRIPE (dropped at 20 bytes or fewer) writes four byte-interleaved substreams without sizes (`NOSZ`), each the smallest of htscodecs' candidates for its position, the first winning a tie: order 1, order 0 with RLE and order 0 for the first; order 1 and order 0 for the second; order 1 and order 0 with PACK for the others; order 1 only when the stream asks for it. EXT (bzip2) is not written. The combinations htslib tries above compression level 1, which the encoder offers as a profile, are order 0 plain, with RLE, with PACK and with both, and order 1 plain, with STRIPE, with PACK, and with RLE and PACK (`ARITH_PR0`, `PR64`, `PR128`, `PR192`, `PR1`, `PR9`, `PR129`, `PR193`). The decoder is its oracle: every input MUST round-trip under every combination; htscodecs, linked into the tests, MUST write the same bytes.

r[cram_writer.fqzcomp]
The fqzcomp encoder (CRAM 3.1 method 7) MUST write streams `r[cram.codec.fqzcomp]` decodes to its input: a `QS` block's qualities, given each record's count of them and whether the record is reverse-strand (`0x10`) and the second of a pair (`0x80`). It picks parameters as htscodecs' `fqz_pick_parameters` does for its four built-in strategies (htslib's `FQZ`, `FQZ_b`, `FQZ_c`, `FQZ_d`), for CRAM 3.x (`vers` 3): one parameter block; reverse-strand records' qualities coded reversed (`GFLAG_DO_REV`); dedup, a quality map for alphabets of at most 8 values, a fixed record length, and selectors from the records' average quality and READ1/READ2 split chosen by htscodecs' statistics, in the same floating-point operations. Those choices match htscodecs' up to floating-point rounding near a threshold: htscodecs may be compiled with fused multiply-adds (clang contracts them on aarch64) and the logarithm is the platform's, so a selector decision on the edge can differ; the stream still decodes. Given the same parameters it codes as `compress_block_fqz2f` does, byte for byte. An empty record cannot be coded: one before the last non-empty record, which htscodecs refuses, is left out; empty records after it are not coded but count in the statistics (record count, fixed length, selectors) as htscodecs counts them. Record counts that do not add up to the qualities MUST NOT be encoded. The decoder is its oracle: every input MUST round-trip under every strategy; htscodecs' `fqz_compress`, linked into the tests, MUST write the same bytes on the test data, on the platforms the tests run on.

r[cram_writer.rans_nx16.simd]
The rANS Nx16 encoder's entropy coding MUST run as SIMD kernels (`r[io.simd_portable]`) for order 0 and order 1 at 4 and 32 states, with htscodecs' exact reciprocal in place of the division (`q = ⌊x·m / 2^(31+s)⌋` with `s = ⌈log2 f⌉` and `m = ⌈2^(31+s) / f⌉`, which equals `⌊x / f⌋` for every state below `2^31`), each lane's renormalization words placed in decode order by a shuffle chosen by the lane mask, and the order-1 table looked up by the two bytes before and at each position. Its frequency counts MUST come from histograms that do not serialize on repeated symbols (several sub-tables for order 0; every byte pair once, then the stretch starts corrected, for order 1; for RLE's choice of run symbols, each symbol's count less twice its runs), and its order-1 tables are built over the alphabet's rows and columns only, in per-thread scratch that is never a reason to panic (a fresh one while it is in use or being destroyed) and keeps no word buffer past 32 MiB. The order-0 and RLE counts and PACK's packing have one implementation, which tests check against counting and packing a byte at a time; PACK's packing MUST run as a SIMD kernel too, each byte's index into the symbol map taken from a 16-entry table where the symbols span fewer than 16 values and otherwise counted as the map symbols at or below it, and pairs of indices merged in 16-bit lanes and narrowed until a byte holds two, four or eight; and PACK's symbols MUST be found without counting where no other combination has counted the block, each 32 bytes compared with the symbols seen so far and the scan given up at a 17th. The scalar encoder (counting each byte of order 1 with its context, dividing, one state at a time) stays as the reference for the rest: at every SIMD level, the scalar `Fallback` included, the encoder MUST write the scalar encoder's bytes for every input and every combination `r[cram_writer.rans_nx16]` takes.

r[cram_writer.tok3]
The tok3 encoder (CRAM 3.1 name tokeniser, method 8, [CRAMcodecs] §6) MUST write blocks that `r[cram.codec.tok3]` decodes, and htscodecs' `tok3_decode_names` too, to its input: the `RN` block as the writer stores it, every name followed by a NUL. It tokenises as htscodecs' `tokenise_name3.c` encoder does, so its token streams are the ones htscodecs would compress: a non-empty name identical to an earlier one in the block is a `DUP` of the most recent one; any other is a `DIFF` against the previous name or, when it fits one of htscodecs' platform naming schemes, against the most recent name sharing that scheme's prefix (none: no name to diff against). The schemes are htscodecs': PacBio movie names (a 60-byte prefix), Ion Torrent (6 bytes), ONT read UUIDs (37 bytes), and Illumina, four colons before the first space, whose prefix ends at the fourth colon from that end; but for PacBio's, the prefix is the name's first token. The rest splits into runs of letters and punctuation (`ALPHA`, or `CHAR` for one byte), digit runs of at most nine (`DIGITS0` with its length when it starts with `0` or the previous name's token there was a `DIGITS0` of the same length, else `DIGITS`) and single other bytes (`CHAR`); a token equal to the previous name's at its position is `MATCH`, and a number above that token's by less than 256 is `DELTA0` (same length) or `DELTA` (while deltas at that position have not fallen five behind plain values). Where htscodecs gives up, the encoder does not: a name of more than 126 tokens stores its rest in its 126th, one `ALPHA`, and bytes outside ASCII are tokens like any other. Empty input, which htscodecs' decoder refuses, and input whose last name has no NUL MUST NOT be encoded. As in htscodecs, a position's type stream that is its first type followed only by `MATCH` is left out (the decoder regenerates it), and a stream identical to an earlier one is written as a copy of it. Every stream is an rANS Nx16 stream (`use_arith` 0), compressed as htscodecs' `compress` does for its token type. The decoder is the oracle: every input MUST round-trip.

## Compression header

r[cram_writer.preservation_map+2]
The preservation map MUST give `RN` true, `AP` true (positions are deltas), `RR` true for a slice written against the FASTA (multi-reference slices included) and false otherwise — false for a slice that embeds its reference, as htslib writes it, since noodles reads the embedded block only when `RR` is false, `SM` the substitution matrix `r[cram_writer.substitution]` uses, and `TD` the slice's tag lines in the order of their first use.

r[cram_writer.data_series+2]
A data series of integers or bytes whose values in the slice are all one value MUST be encoded as a one-symbol `HUFFMAN` code of that value, which reads no bits and needs no block; one with no values, as a one-symbol code of 0. Every other data series MUST be encoded `EXTERNAL`, each in a block of its own content id. `RN`, `IN` and `SC` MUST be `BYTE_ARRAY_STOP` with stop byte 0; `BB` MUST be `BYTE_ARRAY_LEN` with external length and external bytes, in two blocks.

r[cram_writer.tags]
A tag is keyed by its name and BAM type, `name[0] << 16 | name[1] << 8 | type`, and stored in the external block of that key, as htslib does. `A`, `c`, `C`, `s`, `S`, `i`, `I` and `f` values MUST be `BYTE_ARRAY_LEN` with a constant length (a one-symbol `HUFFMAN` code of the type's size) and external bytes; `Z` and `H` values (with their NUL) `BYTE_ARRAY_STOP` with stop byte TAB; `B` arrays `BYTE_ARRAY_LEN` with external length and bytes in the same block. A `Z` or `H` value holding a TAB, which SAM does not allow and the stop byte could not frame, MUST be refused (`TabInStringTag`). A `cF` tag (htslib's private CRAM flags, which readers drop) MUST be refused (`ReservedTag`), and so MUST a tag name other than a letter followed by a letter or digit (`InvalidTagName`): the tag dictionary is NUL-separated.

## Records

r[cram_writer.record.order]
A record's data series MUST be written in the order `r[cram.record.decode_order]` reads them: `BF`, `CF`, `RL`, `AP`, `RG`, `RN`, the mate series, `TL`, the tags, then for a mapped record `FN`, the features, `MQ` and `QS`, and for an unmapped one `BA` and `QS`.

r[cram_writer.record.flags]
`BF` MUST be the record's FLAG. `CF` MUST have bit `0x1` (qualities stored as an array) set always, `0x8` when SEQ is `*`, and the mate bits of `r[cram_writer.mates]`.

r[cram_writer.record.position]
`RL` MUST be the read's length: SEQ's, or for SEQ `*` the CIGAR's query length (htslib's `fake_qual`), whose bases the decode rebuilds from the reference and then drops. `AP` MUST be the 1-based position minus the previous record's in the slice (the first record's minus the slice's alignment start); an unplaced record's position is 0.

r[cram_writer.record.read_group]
`RG` MUST be the index of the record's `RG:Z` value among the header's `@RG` lines, with the tag left out of the stored tags; `-1` if it has none. A value no `@RG` line has stays a stored tag, as does a value naming a last `@RG` line whose ID is `UNKNOWN`, which readers drop as Java CRAM's (`r[cram.record.rg_tag+2]`).

r[cram_writer.mates+2]
Two records of a slice MUST be stored attached when they are the only two records of the slice with their read name, both are mapped on the same reference (htslib's decode gives a pair across references TLEN 0, noodles computes one) (for a pair with an unmapped read noodles computes a TLEN where htslib and seqair give 0), and the decode rebuilds each one's mate fields as they are (`r[cram.record.mate_attached]`): its FLAG already has every bit the decode adds (`0x1`; `0x8` if the other is unmapped; `0x20` if the other is reverse), its RNEXT and PNEXT are the other's reference and position, and its TLEN is the one the decode computes for the pair. The first gets CRAM flag `0x4` and `NF`, the number of records between them; the second neither `0x2` nor `0x4`. Every other record is detached. CRAM 3 has no other way to store a TLEN: htslib's CRAM flag `0x10`, which htslib and seqair read, is CRAM 4's, and noodles and htsjdk would take the next record's `TS` for it, so a pair whose TLEN the decode would compute otherwise is detached. `CF` and the mate series are written once the slice is complete, when its pairs are known.

r[cram_writer.record.mate]
A detached record (CRAM flag `0x2`) MUST store `MF` (bit `0x1` for FLAG `0x20`, bit `0x2` for `0x8`), `NS` (the mate's reference id, `-1` for none), `NP` (the mate's 1-based position, 0 for none) and `TS` (the template length) carry the BAM mate fields as they are.

r[cram_writer.record.name]
`RN` MUST be the read name. An empty read name MUST be refused (`EmptyReadName`): BAM stores none as `*`, and a CRAM reader names a record with an empty name itself. A name holding a NUL, which would end it early in `RN`, MUST be refused (`NulInReadName`).

r[cram_writer.record.unmapped]
An unmapped record's bases MUST go to `BA`, one byte each, as the ASCII of its IUPAC code. CRAM stores no CIGAR and no mapping quality for an unmapped record: the writer MUST drop them, as htslib does (some aligners give unmapped reads a MAPQ), and MUST log a warning the first time it does.

r[cram_writer.record.quality]
`QS` MUST hold the record's qualities as they are, `0xFF` for each base of a record without, after the `QS` values of its `B` features.

r[cram_writer.record.mapq]
`MQ` MUST be the mapping quality of a mapped record.

## Read features

r[cram_writer.features]
A mapped record's CIGAR MUST become read features (`FC`, and `FP` the 1-based read position minus the previous feature's), walking the CIGAR against the reference: `I` of one base `i` (`BA`) and longer `I` (`IN`); `D` (`DL`); `N` (`RS`); `S` (`SC`); `H` (`HC`); `P` (`PD`). `M`, `=` and `X` compare each read base with the reference base: a base that is the same letter (the reference uppercased) is an implicit match and writes nothing, any other is a substitution (`r[cram_writer.substitution]`), as htslib compares them. A reference letter that is no IUPAC code, or `U`, matches no read base, though a decode copying it would read it as `N` or `T`. The CIGAR's query length MUST be SEQ's (`CigarSeqLengthMismatch`).

r[cram_writer.substitution+2]
A mismatching read base MUST be an `X` feature (`BS`) when the read base is `A`, `C`, `G` or `T`, or is `N` against an `A`, `C`, `G` or `T`, and otherwise a `B` feature (`BA` and its quality in `QS`), as htslib chooses. The substitution matrix MUST be the slice's own: for each reference base, codes 0–3 name the other bases from the one the slice substitutes for it most often to the least, ties in the order `ACGTN` lists them, so the commonest substitutions (a C→T conversion, say) all store code 0 and `BS` compresses to little. A reference base other than `ACGT` uses the `N` row, as the decoder does. htslib writes one fixed matrix.

r[cram_writer.past_reference_end]
A read base aligned past the end of the reference MUST be a `B` feature: there is no reference base to match.

r[cram_writer.embedded_reference]
A slice on a reference the FASTA lacks MUST embed a reference built from its records, as htslib's `embed_ref=2` does (`cram_generate_reference`), when its mapped records' `M`, `=` and `X` bases (SEQ `*` has none) number at least 1% of its extent, as htslib requires, and its extent is at most 2^24 positions (the counts below take memory per position; htslib gives up near 22 Mi). Above that density it is smaller than storing every base in `b` features: the reference costs each covered position once, and a run of `N` next to nothing; below it, the slice is written without a reference (`r[cram_writer.reference_free+2]`). The reference covers the slice's extent (`r[cram_writer.slice_extent+2]`). Where a record's `MD` fits its CIGAR (each match run and mismatch consumes aligned bases, each `^` run deleted ones, and none are left over), the bases it implies are the reference's there, a match's the read's base and a mismatch's or deletion's the letter, if it is one of `ACGTN` (the later record's where two disagree), as htslib's `cram_add_to_ref_MD` takes them; every other base is the one of `A`, `C`, `G`, `T` the most of those aligned bases at its position are (on a tie, the first in that order), or `N` where none is. It is stored in an external block of its own content id, which the slice header names (`r[cram_writer.slice_header+2]`), and its records are written against it as against the FASTA (`r[cram_writer.features]`), their `MD` and `NM` dropped where the decode regenerates them from it (`r[cram_writer.md_nm+2]`).

r[cram_writer.embedded_reference.md_nm]
A mapped record with bases in a slice that embeds its reference MUST flag the `MD` and `NM` tags it lacks in a `cF:C` tag, as htslib does under `embed_ref=2`: a decode would otherwise generate them against bases built from the reads rather than the reference they were aligned to.

r[cram_writer.reference_free+2]
A slice on a reference the FASTA lacks that does not embed one MUST be written without a reference (`RR` false): every `M`, `=` and `X` run becomes a `b` feature holding its bases (`BB`), as htslib does without a reference.

r[cram_writer.md_nm+2]
The writer MUST drop a record's `MD:Z` tag when it is byte for byte the `MD` a decode generates from the stored features against the slice's reference (the FASTA's, or the embedded one) (`r[cram.stream.md_nm]`), and its `NM` tag when its value is the generated `NM`; otherwise it keeps them. No tag is generated, so none is dropped, for an unmapped record, SEQ `*`, or a slice without a reference. Neither is dropped where the read's reference span holds a byte that is no IUPAC code: seqair's decode writes it into `MD` as `N`, htslib's as it is.

r[cram_writer.md_nm.past_reference_end]
A record with a deletion that runs past its reference's end MUST keep `MD` and `NM` as they are, and MUST flag the ones it lacks in a `cF:C` tag (bit 1 `MD`, bit 2 `NM`), htslib's private flags, which every reader drops (`r[cram.record.cf_tag]`). htslib's decoder fails such a record whenever it generates either tag ("CRAM CIGAR extends beyond slice reference extents"), even in files htslib wrote.

## Fidelity

r[cram_writer.fidelity]
A record read back MUST equal the record written, field by field, except for these normalizations, each one htslib's:

- `=` and `X` CIGAR operations read back as `M`, adjacent operations of one kind as one, and zero-length operations not at all.
- Tags read back in another order: a dropped `MD`/`NM` is regenerated after the stored tags, and `RG` comes last (`r[cram.record.rg_tag+2]`). A mapped record with bases and without `MD` or `NM` reads back with the one it lacks generated, as htslib's decoder generates it, except in a slice that embeds its reference (`r[cram_writer.embedded_reference.md_nm]`); an `NM` reads back as the smallest unsigned type holding it.
- The mate reference of a record without FLAG `0x1` reads back as none.
- An unmapped record reads back with MAPQ 0 and no CIGAR.

The oracles are seqair's CRAM reader and `samtools view`: the records of a CRAM the writer wrote MUST print, through both, as `samtools view` prints the input after the same normalizations.

## Index

r[cram_writer.crai+2]
With `write_index(true)`, the writer MUST build a `CramIndex` with one entry per slice: reference id (`-1` unplaced), alignment start and span as in the slice header; for a multi-reference slice, as `samtools index` does, one entry per run of its records on one reference id, with the run's first position and the span to the furthest end of its records (`max(pos, pos + reference length - 1)`, not cut at the reference's end; an unmapped record ends at its position); the container's byte offset in the file, the slice's offset from the container's data start, and the slice's size in bytes. On `to_writer` the offsets count from the first byte the writer wrote. `CramIndex::write_crai` MUST write it as `samtools index` does: gzip-compressed, one tab-separated line per entry.

## Testing

r[cram_writer.test_roundtrip]
Round-trip tests MUST write records through each input, read them back with seqair's reader and with `samtools view`, and compare against the input under `r[cram_writer.fidelity]`. Property tests MUST draw reads with substitutions, IUPAC bases, indels, clips, padding, skips, reads past the reference end, unmapped and unplaced reads, SEQ `*` and missing qualities, on references the FASTA has and lacks, sparse and deep, and on enough references with few reads that the writer switches to multi-reference slices and back.

r[cram_writer.test_size]
A test SHOULD compare the size of a written file with `samtools view -C --output-fmt-option version=3.0` of the same input.
