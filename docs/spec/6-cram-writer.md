# CRAM Writer

The CRAM writer serializes alignment records into CRAM 3.0 against a reference FASTA. It is the inverse of seqair's CRAM reader ([CRAM Reader](./2-cram-1-reader.md), [CRAM Stream](./2-cram-2-stream.md)): every structure it writes is one the reader parses, and every record it writes decodes, through seqair's reader and through htslib's, to the record it was given, except for the normalizations `r[cram_writer.fidelity]` lists. It sits next to the BAM writer ([BAM Writer](./6-bam-writer.md)) and takes the same inputs.

The purpose is smaller output for rastair/postalt than BAM at close to `BamWriter` speed. The writer is reference-based: it stores each aligned read as its differences from the reference.

> **Sources:** [CRAM3] §6 "File definition", §7 "Container structure", §8 "Block structure", §8.4 "Compression header block", §8.5 "Slice header block", §9 "End of file container", §10 "Record structure", §10.6 "Read feature records", §12 "Indexing"; [CRAMcodecs] §2 "rANS 4×8". htslib's `cram/cram_encode.c` (`process_one_read`, `cram_encode_aux`) for the choices the format leaves open. See [References](./99-references.md).

## Scope

r[cram_writer.scope]
The writer MUST write CRAM version 3.1, or 3.0 when asked (`r[cram_writer.version]`). Its input MUST be coordinate-sorted (`r[cram_writer.sort_order]`). Every slice holds the records of one reference, or only unplaced records; there are no multi-reference slices and no embedded references. Read names are always stored (`RN` true). Mates in one slice are attached where the decode rebuilds their mate fields exactly (`r[cram_writer.mates]`); every other record's are stored as they are (`r[cram_writer.record.mate]`).

## Writer

r[cram_writer.builder]
`CramWriterBuilder::to_path(path, header, reference)` and `CramWriterBuilder::to_writer(writer, header, reference)` MUST configure a writer of the records of `header` (a `BamHeader`) against `reference` (an `IndexedFastaReader` the writer owns). `build()` MUST write the file definition (`r[cram_writer.file_definition]`) and the header container (`r[cram_writer.header_container]`) before it returns, so a header the writer cannot write fails the build, not a later record. `to_path` creates (truncating) the file at `path`.

r[cram_writer.inputs]
`write(&OwnedBamRecord)`, `write_raw(&RawBamRecord)` and `write_store_record(&RecordStore, RecordIdx)` MUST encode the same record to the same bytes, whichever form it arrives in. A raw record whose CIGAR is held in a `CG` tag (`r[bam.record.long_cigar]`) is written with its real CIGAR and without the `CG` tag; CRAM stores CIGARs of any length.

r[cram_writer.validate_before_write]
A record the writer refuses MUST leave the writer as it was: it is not part of the output and the next record is written as if it had not been offered. Records are buffered until their slice is complete, so a refusal costs nothing. An error writing a container to the sink (I/O, or a codec failure) MUST poison the writer: every later call returns `CramWriteError::Poisoned`, as `r[bam_writer.error_poisoning]`.

r[cram_writer.finish]
`finish()` MUST write the pending slice, then the EOF container (`r[cram_writer.eof]`), flush, and return the inner writer and, when `write_index(true)` was set, the `CramIndex` of what it wrote (`r[cram_writer.crai]`). The writer does not write the index file; `CramIndex::write_crai` does, where the caller chooses. Dropping the writer without `finish()` leaves a file without its last slice and without the EOF container, which htslib reports as truncated.

## File structure

r[cram_writer.file_definition]
The file MUST begin with `CRAM`, major version 3, minor version 0, and a 20-byte file id: the file name of the path, cut to 20 bytes and padded with NUL, as htslib writes it; all NUL for `to_writer`.

r[cram_writer.header_container]
The second container MUST hold one block, of content type `FILE_HEADER`, whose data is the header text's length (i32, little-endian) followed by the text. The container's reference id, alignment start and span, record count, record counter and base count are 0, and it has one landmark, 0.

r[cram_writer.header_text]
CRAM has no binary reference list: the header text is the only place the reader learns the references. If the header's text has no `@SQ` lines, the writer MUST write one per target (`SN`, `LN`) after the `@HD` line (or first, if there is none). If it has `@SQ` lines, their names MUST be the targets', in order, or the build fails with `CramWriteError::HeaderTargetsMismatch`.

r[cram_writer.header_m5]
Every `@SQ` line of a reference the FASTA holds MUST carry `M5`, the lowercase hex MD5 of the reference's bases uppercased, which is what htslib checks a reference against. The writer MUST add `M5` where the line has none. An `LN` that is not the FASTA's length MUST fail the build (`ReferenceLengthMismatch`). An `M5` that is not the FASTA's MUST refuse the first record on that reference and every later one (`ReferenceMd5Mismatch`) without reading the reference again: the file would claim a reference it was not written against. It is checked when the writer first reads the reference rather than at build, so a header that carries every `M5` (one copied from a CRAM) costs no pass over the FASTA. A reference the FASTA lacks keeps its line as it is, and its records are written without the reference (`r[cram_writer.reference_free]`). Computing the missing MD5s reads those references once, at `build()`.

r[cram_writer.eof]
The file MUST end with the CRAM 3 EOF container, the 38 bytes htslib's `cram_check_EOF` matches.

## Containers and slices

r[cram_writer.sort_order]
Records MUST arrive ordered by reference id, then position, with every unplaced record (no reference) after all placed ones; unplaced records are in no order among themselves, though they keep the `POS` they have (htslib writes 1 for many). A record that breaks the order MUST be refused with `CramWriteError::Unsorted`. A record with a reference but no position, and a mapped record without a reference, MUST be refused (`InvalidPlacement`).

r[cram_writer.slicing]
Each container MUST hold one slice. A slice MUST end when the next record is on another reference (or is the first unplaced one), and when it holds 10,000 records or 5,000,000 bases, htslib's defaults. Unplaced records go in slices whose reference id is `-1`.

r[cram_writer.container_header]
A data container's header MUST give its slice's reference id, alignment start and span, its record count, the number of records before it in the file (the record counter), its total read length (`bases`), its block count, one landmark (the offset of the slice header block from the container's data start, which is the size of the compression header block), and the CRC32 of the header's bytes before it.

r[cram_writer.slice_extent]
A slice's alignment start MUST be the 1-based position of its first record, and its span MUST reach the end of every record's reference span (an unmapped record spans its own position), cut at the reference's end. A slice of unplaced records has start 0 and span 0.

r[cram_writer.slice_header]
The slice header MUST give the slice's reference id, alignment start and span, record count, record counter, block count (the core block and every external block), the content ids of its external blocks, embedded reference `-1`, and the MD5 of the uppercased reference bases over its extent — sixteen zero bytes for a slice of unplaced records or one written without the reference.

r[cram_writer.blocks]
Every block MUST be written as `r[cram.block.structure]` reads it: method, content type, content id, compressed and uncompressed size, data, CRC32 of all the bytes before it. Each slice MUST have a core data block (content id 0), empty, as every series is external. The compression header MUST list an encoding for every data series, and the slice MUST have an external block for every one of them, empty or not, and one for each tag it holds: htslib's decoder refuses a slice that lacks the encoding of a series it reads, and noodles one that lacks its block, even when the read is of nothing (`BA` and `QS` for a record without bases).

r[cram_writer.version]
`CramWriterBuilder::version` MUST choose the CRAM version, 3.1 by default as htslib's, or 3.0 for readers without 3.1's codecs; the file definition's minor version MUST be the one chosen.

r[cram_writer.compression]
Each block's data MUST be stored by whichever of the version's codecs gives the fewest bytes, raw included: gzip (at the builder's `compression_level`) for both; then rANS 4×8 order 0 and order 1 (`r[cram_writer.rans4x8]`) in CRAM 3.0, or rANS Nx16 under htslib's normal profile (`r[cram_writer.rans_nx16]`) and, for the `RN` block, the name tokeniser (`r[cram_writer.tok3]`) in CRAM 3.1. The header block is gzip or raw. The choice is per block and changes no decoded byte.

r[cram_writer.threads]
`CramWriterBuilder::compression_threads(n)` MUST assemble and compress containers on `n` worker threads while the calling thread encodes the next slice, and write them in order; `n = 0` (the default) does it on the calling thread. The file and its index MUST be byte for byte those of the single-threaded writer: a slice's content and its codec choices depend on its records only.

r[cram_writer.rans4x8]
The rANS 4×8 encoder MUST write streams `r[cram.codec.rans4x8]` decodes to the input, order 0 and order 1, with each table's frequencies summing to 4095, as htscodecs normalizes them (its decoder takes 4095 or 4096), every present symbol at least 1. Order 1 needs four bytes at least; a block too small for it is not tried with it. The decoder is its oracle: every input MUST round-trip.

r[cram_writer.rans_nx16]
The rANS Nx16 encoder MUST write streams `r[cram.codec.rans_nx16]` decodes to the input, laid out as htscodecs' `rans_compress_to_4x16` writes them, for every combination of order 0 or 1, 4 or 32 interleaved states, PACK, RLE and STRIPE it is asked for; a separate entry point writes a CAT stream (the bytes stored as they are). The transforms run as htscodecs runs them: PACK first, then RLE on what PACK left, then the entropy coder, and each one that does not apply is dropped from the flags rather than failing: PACK with more than 16 distinct symbols, RLE when no symbol repeats enough to pay or the runs and literals would not save 1%, order 1 on fewer than 8 bytes (order 0 instead), STRIPE on 20 bytes or fewer. An entropy-coded body no smaller than its input is stored as CAT under the same transforms. RLE's metadata is order-0 coded when that is smaller. STRIPE writes four substreams without sizes (`NOSZ`), each the smallest of the order the stream asks for, RLE and PACK if asked, and order 0. Every order-0 table's frequencies MUST sum to a power of two no larger than 4096, the decoder scaling them up to 4096; every order-1 row's MUST sum to a power of two no larger than `1 << shift`, with shift 10 or 12 (the only two htscodecs' decoders build tables for), and every symbol of the alphabet MUST have a row, which a symbol that never precedes another fills with symbol 0, so htscodecs' and seqair's decoders both take the row as a complete table. An order-1 table larger than 1000 bytes is order-0 coded when that saves more than 6 bytes. The writer's choice of combinations for a block MUST be htslib's normal profile (compression level 5): order 0, order 1, order 0 with RLE, order 1 with STRIPE, order 0 with PACK, and order 1 with RLE and PACK, each with 32 states from 50 000 bytes on except STRIPE's. The decoder is its oracle: every input MUST round-trip under every combination; so are htscodecs' decoders, scalar and SIMD, when a build of its `rANS_static4x16pr_test` is given (`HTSCODECS_RANS4X16PR`).

r[cram_writer.tok3]
The tok3 encoder (CRAM 3.1 name tokeniser, method 8, [CRAMcodecs] §6) MUST write blocks that `r[cram.codec.tok3]` decodes, and htscodecs' `tok3_decode_names` too, to its input: the `RN` block as the writer stores it, every name followed by a NUL. It tokenises as htscodecs' `tokenise_name3.c` encoder does, so its token streams are the ones htscodecs would compress: a non-empty name identical to an earlier one in the block is a `DUP` of the most recent one; any other is a `DIFF` against the previous name or, when it fits one of htscodecs' platform naming schemes, against the most recent name sharing that scheme's prefix (none: no name to diff against). The schemes are htscodecs': PacBio movie names (a 60-byte prefix), Ion Torrent (6 bytes), ONT read UUIDs (37 bytes), and Illumina, four colons before the first space, whose prefix ends at the fourth colon from that end; but for PacBio's, the prefix is the name's first token. The rest splits into runs of letters and punctuation (`ALPHA`, or `CHAR` for one byte), digit runs of at most nine (`DIGITS0` with its length when it starts with `0` or the previous name's token there was a `DIGITS0` of the same length, else `DIGITS`) and single other bytes (`CHAR`); a token equal to the previous name's at its position is `MATCH`, and a number above that token's by less than 256 is `DELTA0` (same length) or `DELTA` (while deltas at that position have not fallen five behind plain values). Where htscodecs gives up, the encoder does not: a name of more than 126 tokens stores its rest in its 126th, one `ALPHA`, and bytes outside ASCII are tokens like any other. Empty input, which htscodecs' decoder refuses, and input whose last name has no NUL MUST NOT be encoded. As in htscodecs, a position's type stream that is its first type followed only by `MATCH` is left out (the decoder regenerates it), and a stream identical to an earlier one is written as a copy of it. Every stream is an rANS Nx16 stream (`use_arith` 0), compressed as htscodecs' `compress` does for its token type. The decoder is the oracle: every input MUST round-trip.

## Compression header

r[cram_writer.preservation_map]
The preservation map MUST give `RN` true, `AP` true (positions are deltas), `RR` true for a slice written against the reference and false otherwise, `SM` the substitution matrix `r[cram_writer.substitution]` uses, and `TD` the slice's tag lines in the order of their first use.

r[cram_writer.data_series]
Every data series MUST be encoded `EXTERNAL`, each in a block of its own content id. `RN`, `IN` and `SC` MUST be `BYTE_ARRAY_STOP` with stop byte 0; `BB` MUST be `BYTE_ARRAY_LEN` with external length and external bytes, in two blocks.

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

r[cram_writer.mates]
Two records of a slice MUST be stored attached when they are the only two records of the slice with their read name, both are mapped (for a pair with an unmapped read noodles computes a TLEN where htslib and seqair give 0), and the decode rebuilds each one's mate fields as they are (`r[cram.record.mate_attached]`): its FLAG already has every bit the decode adds (`0x1`; `0x8` if the other is unmapped; `0x20` if the other is reverse), its RNEXT and PNEXT are the other's reference and position, and its TLEN is the one the decode computes for the pair. The first gets CRAM flag `0x4` and `NF`, the number of records between them; the second neither `0x2` nor `0x4`. Every other record is detached. CRAM 3 has no other way to store a TLEN: htslib's CRAM flag `0x10`, which htslib and seqair read, is CRAM 4's, and noodles and htsjdk would take the next record's `TS` for it, so a pair whose TLEN the decode would compute otherwise is detached. `CF` and the mate series are written once the slice is complete, when its pairs are known.

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

r[cram_writer.substitution]
A mismatching read base MUST be an `X` feature (`BS`) when the read base is `A`, `C`, `G` or `T`, or is `N` against an `A`, `C`, `G` or `T`, and otherwise a `B` feature (`BA` and its quality in `QS`), as htslib chooses. The substitution matrix is fixed: for each reference base the codes 0–3 name the other bases in the order `ACGTN` lists them (`SM` bytes `0x1B` each); a reference base other than `ACGT` uses the `N` row, as the decoder does.

r[cram_writer.past_reference_end]
A read base aligned past the end of the reference MUST be a `B` feature: there is no reference base to match.

r[cram_writer.reference_free]
A slice on a reference the FASTA lacks MUST be written without it (`RR` false): every `M`, `=` and `X` run becomes a `b` feature holding its bases (`BB`), as htslib does without a reference.

r[cram_writer.md_nm]
The writer MUST drop a record's `MD:Z` tag when it is byte for byte the `MD` a decode generates from the stored features (`r[cram.stream.md_nm]`), and its `NM` tag when its value is the generated `NM`; otherwise it keeps them. No tag is generated, so none is dropped, for an unmapped record, SEQ `*`, or a slice without the reference. Neither is dropped where the read's reference span holds a byte that is no IUPAC code: seqair's decode writes it into `MD` as `N`, htslib's as it is.

r[cram_writer.md_nm.past_reference_end]
A record with a deletion that runs past its reference's end MUST keep `MD` and `NM` as they are, and MUST flag the ones it lacks in a `cF:C` tag (bit 1 `MD`, bit 2 `NM`), htslib's private flags, which every reader drops (`r[cram.record.cf_tag]`). htslib's decoder fails such a record whenever it generates either tag ("CRAM CIGAR extends beyond slice reference extents"), even in files htslib wrote.

## Fidelity

r[cram_writer.fidelity]
A record read back MUST equal the record written, field by field, except for these normalizations, each one htslib's:

- `=` and `X` CIGAR operations read back as `M`, adjacent operations of one kind as one, and zero-length operations not at all.
- Tags read back in another order: a dropped `MD`/`NM` is regenerated after the stored tags, and `RG` comes last (`r[cram.record.rg_tag+2]`). A mapped record with bases and without `MD` or `NM` reads back with the one it lacks generated, as htslib's decoder generates it; an `NM` reads back as the smallest unsigned type holding it.
- The mate reference of a record without FLAG `0x1` reads back as none.
- An unmapped record reads back with MAPQ 0 and no CIGAR.

The oracles are seqair's CRAM reader and `samtools view`: the records of a CRAM the writer wrote MUST print, through both, as `samtools view` prints the input after the same normalizations.

## Index

r[cram_writer.crai]
With `write_index(true)`, the writer MUST build a `CramIndex` with one entry per slice: reference id (`-1` unplaced), alignment start and span as in the slice header, the container's byte offset in the file, the slice's offset from the container's data start, and the slice's size in bytes. On `to_writer` the offsets count from the first byte the writer wrote. `CramIndex::write_crai` MUST write it as `samtools index` does: gzip-compressed, one tab-separated line per entry.

## Testing

r[cram_writer.test_roundtrip]
Round-trip tests MUST write records through each input, read them back with seqair's reader and with `samtools view`, and compare against the input under `r[cram_writer.fidelity]`. Property tests MUST draw reads with substitutions, IUPAC bases, indels, clips, padding, skips, reads past the reference end, unmapped and unplaced reads, SEQ `*` and missing qualities.

r[cram_writer.test_size]
A test SHOULD compare the size of a written file with `samtools view -C --output-fmt-option version=3.0` of the same input.
