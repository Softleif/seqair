# Formatting SAM

Turning a `RawBamRecord` (`r[bam.raw_record]`) back into a SAM text line: what `samtools view` prints for a BAM record. The target is htslib's output byte for byte, so a pipeline can swap `samtools view` for seqair and diff clean.

> **Sources:** [SAM1] §1.4 "The alignment section: mandatory fields", §1.5 "The alignment section: optional fields", §4.2 "The BAM format", §4.2.2 "N_CIGAR_OP field". htslib is the behavioural reference: `sam_format1` (`sam.c`), `sam_format_aux1` (`htslib/sam.h`), `kputd` (`kstring.c`), `bam_tag2cigar` (`sam.c`) and `sam_hdr_write` (`sam.c`). C's `%g` is ISO C §7.21.6.1. See [References](./99-references.md).

## Records

r[sam.format]
`sam::format::format_sam_record(record, header, out)` MUST append to `out` the SAM line htslib's `sam_format1` produces for the same record read from a BAM file with the same header, byte for byte, without a trailing newline. On error, `out` MUST be left as it was.

r[sam.format.qname]
QNAME MUST be the stored read name without its NUL terminator: all `l_read_name - 1` bytes, including any NUL inside them, as htslib copies them.

r[sam.format.fixed_fields]
FLAG, MAPQ and TLEN MUST be printed as decimal integers (MAPQ 255 included). POS and PNEXT MUST be the 0-based position plus one, so an unplaced record (`-1`) prints `0`. RNAME MUST be the header's name for the reference id, or `*` for none. RNEXT MUST be `*` for no mate reference, `=` when the mate reference id equals the record's own, and the header's name otherwise.

r[sam.format.unknown_reference]
A reference or mate reference id the given header has no target for MUST fail with `SamFormatError::UnknownReference`, naming which field. (The record was validated against some header; it need not be this one.)

r[sam.format.cigar]
CIGAR MUST be `*` for a record without operations, and otherwise each operation's length in decimal followed by its letter from `MIDNSHP=XB` for codes 0–9 and `?` for the reserved codes 10–15, as htslib's `bam_cigar_opchr` prints them.

r[sam.format.long_cigar]
A record carrying a long CIGAR (`r[bam.raw_record.long_cigar]`) MUST print the real CIGAR from its `CG` tag in place of the stored placeholder, and that `CG` tag MUST NOT be printed; every other tag keeps its place. This is the record htslib's `bam_read1` hands `sam_format1` after `bam_tag2cigar` moved the CIGAR out of the tag.

r[sam.format.seq_qual]
A record without a sequence MUST print `*` for both SEQ and QUAL. Otherwise SEQ MUST be every base's 4-bit code as its letter from `=ACMGRSVTWYHKDBN`, and QUAL MUST be `*` when the first quality byte is `0xFF`, and otherwise every quality byte plus 33, wrapping modulo 256 as htslib's `add33` does (bytes above 222 do not give printable characters, and that is what htslib prints).

## Tags

r[sam.format.aux]
Every aux tag MUST be printed, in stored order, as a tab followed by `TAG:TYPE:VALUE`: every integer type (`c`, `C`, `s`, `S`, `i`, `I`) as `i` and its value in decimal; `A` as its byte; `Z` and `H` as their bytes up to the NUL terminator; `f` and `d` as C's `%g` (`r[sam.format.float_g]`) of the value; a `B` array as `B:` and its subtype letter, then `,` and each element, integers in decimal and `f` elements through htslib's `kputd` (`r[sam.format.float_kputd]`). An empty array prints as `B:` and its subtype with nothing after it.

r[sam.format.aux_walk]
The tags MUST be printed by one native-width SIMD kernel (`r[io.simd_portable]`) per record, checking each tag as `Aux::iter_checked` does (`r[bam.record.aux_checked]`) and finding each `Z`/`H` string's NUL once.

r[sam.format.aux_errors]
A malformed aux block MUST fail with `SamFormatError::Aux`, carrying the `AuxLayoutError` (`r[bam.record.aux_checked]`): an unknown type or array subtype (`B:A` included, which htslib also rejects), a truncated value, or a `Z`/`H` string without its NUL. Where htslib prints what it can and stops quietly — 1 to 3 bytes left over after the last whole tag — seqair MUST fail too, since dropping bytes without a word is the bug this formatter exists to avoid.

## Floats

r[sam.format.float_g]
`format_g` MUST write what C's `printf("%g", value)` writes for a `f64`: six significant digits, correctly rounded (ties to even on the exact binary value); the exponent `X` of that rounded value chooses fixed notation when `-4 <= X < 6`, else `d.ddddde±XX` with at least two exponent digits; trailing zeros of the fraction, and then a bare decimal point, removed. `-0.0` MUST print `-0`, infinities `inf`/`-inf`, and NaN `nan`, or `-nan` when its sign bit is set. That is glibc's and musl's output; macOS's libc differs in two places — it prints `nan` for every NaN, and it keeps the mantissa's trailing zeros when an exact tie rounds down in exponent form (`1000005` prints `1.00000e+06`, glibc `1e+06`) — so htslib's own output differs between platforms there, and seqair takes Linux's. It MUST NOT allocate or call libc.

r[sam.format.float_kputd]
`format_kputd` MUST be a port of htslib's `kputd` (what `B:f` elements go through): `0`/`-0` for zeros; for magnitudes in `[0.0001, 999999]` the value scaled to a six-digit integer with C's `rint` (round half to even) on the `f64` product, and the decimal point placed by the magnitude's band, with trailing zeros and a bare point removed; every other value — NaN and infinities included — as `%g` (`r[sam.format.float_g]`), its sign printed first. The scaling can round differently from `%g` (it rounds a rounded product: `1.000005` prints `1`, `%g` gives `1.00001`), and the port MUST keep that. (No `f32` in the scaled range is affected, so for the `B:f` elements htslib actually prints, `kputd` and `%g` agree.)

## Header

r[sam.format.header]
`sam::format::format_sam_header(header, out)` MUST append the header as `samtools view -h --no-PG` prints a BAM file's: the stored header text up to its first NUL (BAM writers pad `l_text` with NULs), a newline added when the text does not end with one, and, when no line starts with `@SQ\t`, one `@SQ\tSN:<name>\tLN:<length>` line per target, as htslib's `sam_hdr_write` adds them. A line that does not start with `@` — an empty one included — MUST fail with `SamFormatError::MalformedHeaderLine`, as htslib's `sam_hdr_sanitise` rejects the header. htslib writes the NUL padding itself and any text after it into its SAM output; seqair stops at the first NUL.
