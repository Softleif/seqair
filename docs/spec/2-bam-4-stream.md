# Streaming BAM

The indexed reader (`r[bam.reader.fetch_into+2]`) answers region queries: it needs an index, a seekable file, and it keeps placed records only. Many tools instead consume a whole alignment file front to back — postprocessing a BAM piped out of an aligner, converting formats, filtering — and need every record exactly as it was written, unplaced and unmapped ones included, from input that cannot seek (`stdin`, a pipe).

> **Sources:** [SAM1] §4.2 "The BAM format" — `block_size`, the fixed record fields and the variable-length fields that follow; §4.2.2 "N_CIGAR_OP field" — the `CG` placeholder for long CIGARs; §4.1 — BGZF. htslib's `bam_read1` (`sam.c`) is the behavioural reference for what a reader accepts. See [References](./99-references.md).

## Raw records

r[bam.raw_record]
`RawBamRecord` MUST hold one record in its wire form — the bytes after `block_size` — and MUST NOT decode, normalise or re-encode them: what a streaming reader yields, `BamWriter::write_raw` writes back byte for byte. Accessors decode on demand: reference (`Option<Tid>`), position (`Option<Pos0>`), `MateFields`, `BamFlags`, `MappingQuality`, the stored `bin`, read name, CIGAR operations, sequence (packed, as `=ACMGRSVTWYHKDBN` ASCII with every IUPAC code kept, or as lossy `Base`s), qualities (`None` when the record has none — no sequence, or a first quality byte of `0xFF`) and aux tags. FLAG and MAPQ MAY be set in place, since neither moves another field; setting FLAG MUST re-run the CIGAR/sequence check (`r[bam.raw_record.cigar_seq]`), which depends on `0x4`, and on failure leave FLAG unchanged — else a record valid when read would be written out as one htslib rejects.

r[bam.raw_record.validation]
A `RawBamRecord` MUST only be constructible through validation against a header, and every value MUST satisfy it, so that no accessor can fail or panic: the payload is at least 32 bytes and the variable-length fields it declares fit inside it (`r[bam.record.checked_offsets]`); `l_read_name` is at least 1 and the name's last byte is NUL (htslib repairs a missing terminator; seqair rejects the record rather than rewrite it); the reference and mate reference ids are `-1` or name a target of the header (`r[bam.record.ref_ids]`); the position and mate position are `-1` or non-negative (`r[bam.record.unplaced_pos]`); and the payload is no longer than the record size limit (`r[bam.record.max_size]`, `DecodeError::RecordTooLarge`), so every value `write_raw` can write. A reference id with an unset position (`-1`) is accepted, as htslib's `bam_read1` accepts it, though SAM parsing never produces one. A failed validation MUST leave the value it was refilling as the default record (unplaced, unmapped, named `*`), never half-updated.

r[bam.raw_record.long_cigar]
A record carries a long CIGAR (`[SAM1] §4.2.2`) exactly when htslib's `bam_tag2cigar` would move one out of its `CG` tag on read: it has a CIGAR whose first operation is exactly `<seq_len>S`, its reference id and position are both set, and its first `CG` tag is a `B:I` or `B:i` array with at least as many elements as the stored CIGAR has operations and fewer than 2^29. The array's elements are then the record's real CIGAR, in BAM's packed form. A malformed aux block before the first `CG` tag means there is none (htslib fails the read there; seqair does not validate aux on read, `r[bam.raw_record.aux_lazy]`, and every interpretation of the tags reports the malformed one). `RawBamRecord::has_cg_placeholder` answers whether a record carries one.

r[bam.raw_record.cigar_seq]
As htslib's `bam_read1` does, a record with a sequence MUST fail with `DecodeError::CigarQueryLenMismatch` when its CIGAR's query-consuming operations (M, I, S, =, X) do not sum to the sequence length, where the CIGAR checked is the long one if the record carries one (`r[bam.raw_record.long_cigar]`, checked whatever FLAG says), and otherwise the stored one, checked only when the record is mapped (FLAG `0x4` clear) and has CIGAR operations. (This is htslib 1.25's check; htslib before 1.25 did not check the `CG` CIGAR.)

r[bam.raw_record.cigar_access]
`cigar()` and `cigar_len()` MUST give the record's real CIGAR, as htslib's reader returns it: for a record with a long CIGAR (`r[bam.raw_record.long_cigar]`) the operations from its `CG` tag, otherwise the stored ones. `stored_cigar_bytes()` gives the CIGAR field exactly as stored (the placeholder, for a long CIGAR).

r[bam.raw_record.aux_lazy]
The aux block MUST NOT be validated when a record is read: htslib passes a malformed aux block through `samtools view -b` untouched and only fails when it has to interpret it, and a streaming copy must not reject files htslib copies. Interpreting the tags is where a malformed block MUST surface — `Aux::iter_checked` (`r[bam.record.aux_checked]`), and the SAM formatter, which uses it.

r[bam.raw_record.edit_seq]
`RawBamRecord::set_seq(&[IupacBase])` and `set_base(QPos, IupacBase)` MUST replace the sequence, or one base of it, in place in the packed SEQ field, storing each base's 4-bit code (so every IUPAC code survives). The sequence MUST keep its length: `set_seq` with a different number of bases, and `set_base` past the last base, MUST fail with a `RawRecordEditError` and leave the record unchanged. Nothing else MUST change — not QUAL, the CIGAR, the aux block, `bin`, nor the unused low nibble of an odd-length sequence's last byte. Because the length and the other fields are untouched, the record still satisfies `r[bam.raw_record.validation]` and `r[bam.raw_record.cigar_seq]` without re-validation.

r[bam.raw_record.iupac_bases]
`RawBamRecord::iupac_bases()` MUST yield the record's bases as `IupacBase`, one per base, every 4-bit code kept; it MUST agree with `decode_seq_into` (the ASCII spelling, used by the SAM formatter) letter for letter.

r[bam.raw_record.edit_aux]
`RawBamRecord::edit_aux(f)` MUST hand `f` an `AuxData` holding the record's aux block and replace the block with what `f` leaves. Tags `f` does not touch MUST keep their bytes, order and integer width (`AuxData::remove` and the setters only move the tag they name). The edit MUST fail, leaving the record unchanged, if `f` fails (`RawRecordEditError::Aux`), if the payload would exceed the record size limit (`r[bam.record.max_size]`), or if it would change the record's long CIGAR (`r[bam.raw_record.long_cigar]`): removing, replacing or adding a `CG` tag so that the real CIGAR differs from before (a record without a long CIGAR gaining one included) is `RawRecordEditError::LongCigarChanged`. Rewriting `CG` to the same bytes is allowed. The existing aux block MUST be validated first (`aux::validate_aux`): a malformed block fails with `RawRecordEditError::MalformedAux` and the record is unchanged, because an edit has to find tag boundaries and cannot do so in a block that does not parse. Reading and copying a record still does not validate it (`r[bam.raw_record.aux_lazy]`).

r[bam.raw_record.bin_verbatim]
The stored `bin` MUST be kept as the file has it. (htslib recomputes it on read, so `samtools view -b` may rewrite a wrong one; a raw copy does not.)

## Reading a stream

r[bam.stream]
`BamStreamReader<R: Read>` MUST read a BGZF-compressed BAM front to back from any `io::Read` through a `BgzfStream` (`r[bgzf.stream]`): first the header (the same parser as the indexed reader's, `r[bam.header.magic]` and following), then every record in file order, unplaced ones included, as `RawBamRecord`s. Threads and EOF-marker handling are the `BgzfStream`'s, configured by building one and passing it in.

r[bam.stream.framing]
Each record is a little-endian `i32` `block_size` and that many bytes. Zero bytes where a record would start MUST be the end of the stream. A `block_size` below 32 or above the record size limit (`r[bam.record.max_size]`) MUST fail with `BamStreamError::RecordSizeOutOfRange`; input that ends inside a `block_size` or a payload MUST fail with `BamStreamError::TruncatedRecord`. Both leave the stream misaligned, so every later read MUST fail with `BamStreamError::Poisoned`. A record whose bytes are framed correctly but fail validation (`r[bam.raw_record.validation]`) MUST fail with `BamStreamError::InvalidRecord`, carrying its 1-based record number and the `DecodeError`, and the stream MUST stay usable: the next read returns the next record.

r[bam.stream.payloads]
`read_payload_into` MUST frame the next record exactly as `read_record` does (`r[bam.stream.framing]`) but append its payload bytes unvalidated, so a caller can validate on other threads (`RawBamRecord::from_payload`). Any error MUST end the stream (`BamStreamError::Poisoned` afterwards) and leave the caller's buffer as it was.
