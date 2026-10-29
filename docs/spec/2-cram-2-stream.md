# Streaming CRAM

The indexed CRAM reader (`r[cram.index.query+2]`) answers region queries: it needs a `.crai`, a seekable file, and it keeps the records a query asks for. A tool that rewrites a whole alignment file — adding tags, rewriting SEQ, converting formats — instead reads it front to back, often from a pipe, and needs every record, unplaced and unmapped ones included, in file order, as `BamStreamReader` gives them for BAM (`r[bam.stream]`).

CRAM decoding is the expensive part, and it splits cleanly: reading a container is a handful of small reads and one copy, decoding it is all CPU. So the stream is two types. `CramStreamReader` reads containers off the input, serially and cheaply; `CramContainerDecoder` turns one container into records and can run on any thread, one decoder per thread.

> **Sources:** [CRAM3] §6 "File definition", §7 "Container header structure", §8 "Block structure", §8.4 "Slice header block", §9 "End of file container", §10 "Record structure". htslib's `cram_read_container` (`cram_io.c`), `cram_decode_slice` and `cram_to_bam` (`cram_decode.c`) are the behavioural reference. See [References](./99-references.md).

## Reading containers

r[cram.stream+2]
`CramStreamReader<R: Read>` MUST read a CRAM v3.0 or v3.1 file front to back from any `io::Read` without seeking: the 26-byte file definition (`r[cram.file.magic]`, `r[cram.scope.versions]`; a major version other than 3 is `CramError::UnsupportedVersion`), then the file header container (`r[cram.file.header_container]`), whose header text is parsed as the indexed reader's is (`BamHeader::from_sam_text`) but whose sort order is not checked — a stream reads every record in file order, so a `queryname`-sorted or unsorted file (`samtools sort -n`, `samtools collate`) streams like any other — then one container per `next_container` call, in file order. Every byte of a container's header and body is consumed, so whatever follows is the next container; padding a writer leaves in the header container (samtools reserves room for in-place header edits) is consumed with it.

r[cram.stream.container_header]
A container header MUST be read field by field — the fixed `length`, each ITF8/LTF8 field one prefix byte and then as many bytes as the prefix says, the landmark count and that many landmarks, the CRC32 — and only then parsed and its CRC32 verified (`r[cram.container.header]`). Reading a fixed-size chunk and parsing whatever arrived is not safe on a pipe, which may deliver any number of bytes per read.

r[cram.stream.framing]
`next_container` MUST return `None` after the EOF container (`r[cram.file.eof]`) and on every later call. Input that ends cleanly where a container would start, with no EOF container read, MUST also end the stream, with one `warn!` (`r[cram.file.eof]`: a missing EOF container is not an error). Input that ends inside a container header or body MUST fail with `CramStreamError::TruncatedContainer` carrying the container's 0-based number; a container whose CRC32 does not match, or whose length is negative or above the allocation limit (`r[io.fuzz.alloc_limits]`), MUST fail with the `CramError` that says so. Any error leaves the input misaligned, so every later call MUST fail with `CramStreamError::Poisoned`. Containers with zero records that are not the EOF container hold nothing to decode and are not returned.

r[cram.stream.container]
A `CramContainer` MUST own its header and body bytes and be `Send`, so a reading thread can hand containers to decoding threads. It carries its 0-based number in the file, for error messages.

## Decoding containers

r[cram.stream.decoder]
`CramContainerDecoder` MUST hold everything one thread needs to decode containers of one file: the file header, its `@RG` ids, a FASTA reader (`r[cram.scope.reference_required]`) and the scratch buffers and codec tables the indexed reader keeps per reader. `fork` MUST give a decoder for another thread that shares the parsed header and FASTA index. `decode(container, out)` MUST append every record of the container to `out` as a `RawBamRecord`, slice by slice in landmark order and record by record within a slice, so decoding containers in file order yields the file's records in file order whichever thread decodes which container.

r[cram.stream.all_records]
Every record of a slice MUST be decoded: in a multi-reference slice (`r[cram.slice.multi_ref]`) the records of every reference and the unplaced ones, in an unmapped slice (reference `-1`) every record, and placed-unmapped reads wherever they are. Nothing is filtered by reference or position. The indexed path keeps decoding one reference per slice and skipping the reconstruction of the others (`r[cram.slice.multi_ref_skip]`).

r[cram.stream.reference_window]
The reference a slice decodes against MUST come from the slice itself, as there is no index. A single-reference slice decodes against its header's range `[alignment_start, alignment_start + alignment_span)` of that reference, clamped to the reference's length; MD5 verification (`r[cram.edge.reference_mismatch]`), embedded references (`r[cram.slice.embedded_ref]`) and the missing-reference error with its `REF_PATH` hint (`r[cram.edge.missing_reference+2]`) apply as in the indexed reader. A slice of unmapped records needs no reference. A multi-reference slice carries no usable range in its header, so it MUST be decoded twice: a first pass reconstructs nothing (as `r[cram.slice.multi_ref_skip]` does for foreign records) and collects, per reference, the lowest start and highest end of the mapped records on it; the second pass decodes every record against those windows, fetched once per reference. Fetching whole references instead (htslib's choice, with a reference cache that outlives the slice) would read a whole chromosome for one read at its end; decoding twice costs a pass over a slice that is small by construction, since writers emit multi-reference slices only when the references are short of records.

r[cram.stream.raw_records]
Each decoded record MUST become a `RawBamRecord` as htslib's `cram_to_bam` builds the BAM record, validated against the file header with `RawBamRecord::from_payload` (`r[bam.raw_record.validation]`): reference and position as decoded, a position of `-1` where `AP` is below 1 on an unmapped read; FLAG and MAPQ as decoded (MAPQ 0 for unmapped reads, which store none); mate reference, position and TLEN as the slice resolves them (`r[cram.record.mate_detached]`, `r[cram.record.mate_tlen_reconstruction]`); the name, CIGAR, sequence, qualities (`0xFF` throughout when the record stores none) and aux tags, `RG` last (`r[cram.record.rg_tag]`); `bin` as `hts_reg2bin(pos, pos + rlen, 14, 5)` with `rlen` the CIGAR's reference length, or 1 when that is 0 or the read is unmapped (`4680` when unplaced); a CIGAR of more than 65535 operations as BAM's `CG` placeholder with the real CIGAR in a `CG:B:I` tag, as `bam_write1` writes it. The encoding is shared with the SAM encoder (`r[sam.encode.bin]`, `r[sam.encode.long_cigar]`), not a second implementation. A record that fails validation MUST fail the container with `CramStreamError::InvalidRecord`, carrying the record's 0-based number in the file (the slice's `record_counter` plus its index).

## Sequence

r[cram.stream.iupac]
The stream MUST rebuild SEQ as `IupacBase`s, keeping every code htslib keeps: a base an implicit match copies from the reference is the reference's code (`R` stays `R`, either case reads as upper case), and a base a feature stores as a byte (`B`, `b`, `i`, `I`, `S`, an unmapped read's `BA`) is that byte's code as `bam_set1` reads it (`=`, `Y`, ...), `N` for anything that is no IUPAC letter. An `X` substitution of a reference base other than A, C, G or T uses the substitution matrix's `N` row. The indexed reader keeps decoding to `Base`, which its pileup consumers expect (`r[base_decode.decode+2]`), so its path is not slowed by the wider type.

## MD and NM

r[cram.stream.md_nm]
`CramContainerDecoder` MUST generate `MD:Z` and `NM` as htslib does when it decodes (`decode_md`, on by default): for a mapped record that has bases (not SEQ `*`, CRAM flag `0x8`), in a slice that has a reference (an external one — the compression header's `RR` is true — or one embedded in the slice), `MD` unless the record stores an `MD` tag and `NM` unless it stores an `NM` tag. A `cF:C` tag's bit 1 counts as a stored `MD` and bit 2 as a stored `NM` (`r[cram.record.cf_tag]`). Both are computed from the read features as the record is rebuilt, not by comparing bases: a run of implicit matches adds to the match count, except that each reference `N` in it is a mismatch; `X`, `B` and each base of `b` are mismatches against the reference base there (even when a `B` or `b` base equals it); `D` writes `^` and the deleted reference bases; `I` and `i` count only towards `NM`; `S`, `H`, `P`, `Q`, `q` and the reference skip `N` leave the match count running, so `5M100N5M` continues one run across the skip. `NM` is the number of mismatches plus inserted and deleted bases, as the smallest unsigned integer type that holds it (`C`, `S`, else `I`). A generated `MD` follows every stored tag and `NM` follows `MD`; both precede the `RG` tag (`r[cram.record.rg_tag]`). Features past the end of the reference window stop the count as htslib's do (no more numbers in `MD`).

r[cram.stream.md_nm.option]
`CramContainerDecoder::decode_md(false)` MUST turn the generation off, as `samtools view --input-fmt-option decode_md=0` does; forks keep the setting. The indexed reader (`r[cram.index.query+2]`) does not generate either tag: its records feed the pileup, which reads neither.

## Fidelity

r[cram.stream.fidelity+2]
The oracle is samtools: each streamed record MUST be byte for byte the record `samtools view -b` writes for it, `format_sam_record` of it MUST equal the line `samtools view --no-PG` prints, and the stream MUST yield as many records as `samtools view -c` counts; with `decode_md(false)` the same holds against `--input-fmt-option decode_md=0`. The known differences, each pinned by name in the test (`KNOWN_MISMATCHES` in `tests/cram_stream.rs`) so that fixing one fails it until the entry goes:

- Read names when `RN` is false (`lossy_names`): htslib keeps a detached record's stored name and gives an attached record its mate's name or `<prefix>:<record number>`; seqair leaves both empty.
- TLEN of an attached pair with an unmapped mate: htslib sets it to 0, seqair computes the span. Only `lossy_names` files attach such pairs.

Differences no fixture shows, from htslib's source: an attached chain over several references gets TLEN 0 in htslib (seqair has no such rule) and an attached mate at `AP` 0 gets mate position `-1` (seqair 0); the qualities of `B`, `q` and `Q` features are not applied when a record stores no quality array (`r[cram.feature.quality_score]`, `r[cram.feature.quality_block]`; htslib fills 30 first); mate flags are taken from `BF` as stored, where htslib also ORs in `0x1`/`0x8`/`0x20` from an attached mate or `MF`; the last `@RG` named `UNKNOWN` is not dropped; CRAM flag `0x10` (explicit TLEN), `QO=0`, `MD*`/`NM*` placeholders and the `no_ref` `=` prefill are not handled. htslib writes none of these in CRAM 3.x except as noted.
