# Streaming SAM

The indexed SAM reader (`r[sam.reader.fetch_into]`) answers region queries into a `RecordStore`. A tool that consumes a whole SAM front to back — the output of an aligner on `stdin`, a gzipped SAM, a file without `@SQ` lines whose reads are all unmapped — instead needs every line as a BAM record, exactly as htslib would turn it into one, so that it can be modified and written as BAM (`BamWriter::write_raw`) or compared byte for byte with what `samtools view -b` produces.

> **Sources:** [SAM1] §1.3 "The header section", §1.4 "The alignment section: mandatory fields", §1.5 "The alignment section: optional fields"; §4.2 "The BAM format" (fixed fields, `bin`, packed SEQ, `0xFF` QUAL); §4.2.2 "N_CIGAR_OP field" (the `CG` placeholder); §5.3 (`reg2bin`). htslib's `sam_parse1`, `aux_parse`, `parse_sam_flag`, `bam_parse_cigar`, `bam_tag2cigar`, `bam_write1` and `sam_read1` (`sam.c`), `hts_str2uint`/`hts_str2int` (`textutils_internal.h`) and `seq_nt16_table` (`hts.c`) are the behavioural reference. See [References](./99-references.md).

## Encoding one line

r[sam.encode]
`SamRecordEncoder::encode(line, header, out)` MUST turn one SAM alignment line (without its line terminator) into a `RawBamRecord` (`r[bam.raw_record]`) validated against `header`. For every line htslib accepts and that this spec does not list as a deviation, the payload MUST be byte for byte what htslib's `sam_parse1` followed by `bam_write1` writes after the `block_size`. Every line htslib rejects MUST be rejected with a typed `SamRecordError`. The encoder keeps its scratch buffers and writes into `out`'s buffer, so encoding line after line allocates nothing once the buffers have grown; on error `out` is left as the default record (`r[bam.raw_record.validation]`).

r[sam.encode.deviations]
Where htslib accepts a line but silently loses or misreads part of it, the encoder MUST reject the line instead. Each case is named in the rule it belongs to, marked **Deviation**, and has a test showing htslib accepts and seqair rejects it. The encoder never accepts a line htslib rejects.

r[sam.encode.fields]
A line MUST have at least 11 TAB-separated fields (`TooFewFields`); an empty line MUST fail with `EmptyLine` (htslib fails it too: `sam_read1` hands it to `sam_parse1`). A NUL byte anywhere in the line MUST fail with `NulByte` (htslib works on C strings: a NUL ends a field early or the line is cut there). QNAME is kept verbatim and MUST be at most 254 bytes (`QnameTooLong`); an empty QNAME is accepted, as htslib accepts it (with a warning), and is stored as `l_read_name = 1`.

r[sam.encode.numbers]
FLAG, POS, MAPQ and PNEXT MUST match `[0-9]+` and TLEN `[-+]?[0-9]+` (the SAM spec's regular expressions); leading zeros are allowed except in FLAG. **Deviation:** htslib reads an empty field as 0 and accepts a leading `+` in POS, MAPQ and PNEXT; seqair rejects both (`InvalidFlag`, `InvalidPos`, `InvalidMapq`, `InvalidPnext`, `InvalidTlen`). Ranges: FLAG ≤ 65535, MAPQ ≤ 255, POS and PNEXT ≤ 2^31 (so the 0-based position fits BAM's `int32`, which `bam_write1` checks), TLEN within `int32`.

r[sam.encode.flag]
A FLAG with a leading zero other than `0` itself MUST fail with `InvalidFlag`. **Deviation:** htslib parses such a FLAG with `strtoul(…, 0)`, so `016` is octal 14 and `0x10` is 16, where the SAM spec reads `016` as 16.

r[sam.encode.references]
RNAME `*` MUST give no reference; any other name MUST be a target of the header (`UnknownReference`). RNEXT `*` MUST give no mate reference, `=` the record's reference (after `r[sam.encode.fixups]`), and any other name a target of the header (`UnknownMateReference`). **Deviation:** htslib warns and makes the read unmapped (RNAME) or drops the mate reference (RNEXT) for a name its header lacks, losing the name. Alternative names (`@SQ AN:`), which htslib also resolves, are not resolved: seqair's header knows only `SN`.

r[sam.encode.fixups]
The encoder MUST apply htslib's fix-ups, in htslib's order: a record with a reference but POS 0 loses its reference; a record without a reference gets FLAG `0x4`; a record whose CIGAR is `*` gets FLAG `0x4`; a mate reference with PNEXT 0 is dropped. These follow the SAM spec ("no assumptions can be made about RNAME and CIGAR" when POS is 0) rather than lose information.

r[sam.encode.bin]
`bin` MUST be htslib's `hts_reg2bin(pos, pos + rlen, 14, 5)` stored in 16 bits, where `pos` is the 0-based position (`-1` when unplaced) and `rlen` is the CIGAR's reference length, or 1 when the record is unmapped (FLAG `0x4`, after fix-ups) or the reference length is 0. An unplaced record therefore gets 4680. For a position past 2^29 the value does not fit 16 bits and is truncated as htslib truncates it.

## CIGAR

r[sam.encode.cigar]
`parse_cigar(text, out)` MUST clear `out` and parse a SAM CIGAR string: `*` gives no operations; otherwise one or more `[0-9]+[MIDNSHP=XB]`, each length below 2^28. An empty string (`EmptyCigar`), a missing or over-long length (`InvalidCigarLength`), an operator outside `MIDNSHP=XB` (`UnknownCigarOp`) and a length without an operator at the end (`MissingCigarOp`) MUST fail. `B` (htslib's `BAM_CBACK`, code 9) is accepted because htslib accepts it. The line encoder MUST use this parser. **Deviation:** `sam_parse1` takes a CIGAR field that starts with `*` as `*` whatever follows it; seqair rejects `*junk`.

r[sam.encode.cigar_seq]
When a record has a CIGAR and a sequence, the CIGAR's query length (M, I, S, =, X) MUST equal the sequence length, mapped or not, as `sam_parse1` checks (`CigarSeqLenMismatch`). The check runs after `r[sam.encode.cg_tag]`.

r[sam.encode.long_cigar]
A CIGAR of more than 65535 operations MUST be written as `bam_write1` writes it: the CIGAR field holds the two placeholder operations `<l_seq>S<rlen>N` (`n_cigar_op = 2`), where `rlen` is the reference length of the real CIGAR, and a `CG:B:I` tag holding the real CIGAR is appended after every other tag. A reference length of 2^28 or more cannot be written this way and MUST fail with `CigarRefLenTooLong`.

r[sam.encode.cg_tag]
As `bam_tag2cigar` does in `sam_parse1`, when a placed record (reference and POS set) has a CIGAR whose first operation is `<l_seq>S` and its first `CG` tag is a `B:I` or `B:i` array of at least as many elements as the CIGAR has operations and fewer than 2^29, the CG tag MUST become the CIGAR and be removed from the tags; `bin` is computed from the new CIGAR. The record is then written by `r[sam.encode.long_cigar]`, so a CG of more than 65535 operations goes back to the end of the tags.

## SEQ and QUAL

r[sam.encode.seq]
SEQ `*` MUST give an empty sequence. Otherwise every byte MUST be printable ASCII other than space (`!`..=`~`, `InvalidSeq`) and is packed with htslib's `seq_nt16_table`: IUPAC codes in either case, `=`, `U`/`u` as `T`, `0`–`3` as `A`/`C`/`G`/`T`, anything else `N`; an odd-length sequence leaves the last byte's low nibble 0. **Deviation:** htslib packs any byte; seqair rejects whitespace, control and non-ASCII bytes.

r[sam.encode.qual]
QUAL `*` MUST give `0xFF` for every base. Otherwise it MUST have exactly as many bytes as SEQ (`QualLenMismatch`, which includes a QUAL with SEQ `*`) and every byte MUST be in `33..=160` (`InvalidQual`), stored minus 33 — the range htslib's `COPY_MINUS_N` check accepts.

r[sam.encode.limits]
The encoded record MUST be at most the record size limit (`r[bam.record.max_size]`, `RecordTooLarge`), as for every BAM record seqair reads or writes. **Deviation:** htslib has no such limit.

## Optional fields

r[sam.aux.parse]
One strict parser MUST turn the optional fields into BAM aux bytes, for both the indexed reader (`r[sam.record.aux_tags+2]`) and the encoder. Each field MUST be `TG:T:VALUE`: at least 5 bytes, `:` at offsets 2 and 4, tag bytes printable ASCII other than space (`InvalidAuxField`). An empty field MUST fail, except that one trailing TAB at the end of the line is accepted, as htslib accepts it. Duplicate tags are kept, in order, as htslib keeps them. **Deviation:** htslib does not check the two `:` and reads `NMxiy5` as `NM:i:5`.

r[sam.aux.types]
The value MUST be non-empty except for `Z` and `H`, and is encoded by type (`UnknownAuxType` for any other type letter):

- `A`, and htslib's synonyms `a`, `c`, `C`: exactly one byte other than TAB, stored as type `A`. **Deviation:** htslib keeps the first byte of a longer value and drops the rest.
- `i`, and its synonym `I`: `[-+]?[0-9]+`. A value with a leading `-` is stored in the smallest signed type that holds it (`c`, `s`, `i`; so `-0` is `c`), any other value in the smallest unsigned type (`C`, `S`, `I`) — htslib's unsigned-first choice. Outside `-2^31..=2^32-1` MUST fail with `AuxIntOutOfRange`.
- `f` and `d`: a decimal floating-point number (`[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?`) or `[-+]?` followed by `inf`, `infinity` or `nan` in any case, which is what htslib's `%g` writes. It is read as a `double` and, for `f`, rounded to `float`, as `strtod` and the cast do; NaN is the default quiet NaN with the written sign.
- `Z`: any bytes (TAB and NUL cannot occur). `H`: an even number of hex digits (`0-9A-Fa-f`). **Deviation:** htslib accepts any bytes in `H`.
- `B`: a subtype in `cCsSiIf`, then nothing (an empty array) or `,` and comma-separated elements: `[-+]?[0-9]+` within the subtype's range for the integer subtypes, a number as for `f` for `f`. An unknown subtype MUST fail with `InvalidAuxSubtype`, an element out of range or malformed with `InvalidAuxValue`.

**Deviation** (all types): htslib skips anything after a number up to the next TAB (`NM:i:5abc` is 5), reads an empty array element as 0, accepts an empty `B:A` array, and retypes an integer array whose elements do not fit its subtype; seqair rejects each.

## Streaming

r[sam.stream]
`SamStreamReader<R: BufRead>` MUST read a SAM front to back from any `BufRead` — plain text, a `BgzfStream` over bgzipped SAM (`r[bgzf.stream]`), or a `flate2::read::MultiGzDecoder` over plain gzip — and yield every alignment line as a `RawBamRecord` through `read_record(&mut RawBamRecord) -> Result<bool, SamStreamError>` (`false` at the end) and an iterator.

r[sam.stream.header]
The header MUST be every line at the start of the input that begins with `@`, kept as text with each line ending in `\n` (a `\r\n` ending normalised). A header line that is not UTF-8 MUST fail with `HeaderNotUtf8`, naming its line. Unlike the indexed reader (`r[sam.header.sq_required]`), a header without `@SQ` lines — or no header at all — MUST be accepted: its records can only be unplaced, and a record naming a reference fails with `UnknownReference`.

r[sam.stream.lines]
Lines end at `\n`; one `\r` before it is removed, and a last line without `\n` is a line. A line may span any number of reads from the underlying reader. An empty line MUST fail (`r[sam.encode.fields]`), as it fails in htslib's `sam_read1` — unlike the indexed reader, which skips it (`r[sam.edge.empty_lines]`). A line longer than 16 MiB MUST fail with `LineTooLong` before it is buffered whole, so input without newlines cannot exhaust memory.

r[sam.stream.errors]
A line that fails to encode MUST fail with `SamStreamError::InvalidRecord` carrying its 1-based line number (counting header lines) and the `SamRecordError`; the stream stays usable and the next read returns the next line. An I/O error from the underlying reader (including a truncated gzip or BGZF stream) and `LineTooLong` MUST fail the stream: every later read MUST return `SamStreamError::Poisoned`.
