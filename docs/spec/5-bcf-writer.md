# BCF Binary Writer

> **Sources:** [BCF2] — magic bytes, header encoding (l_text + NUL-terminated header text), typed value encoding (type byte: count << 4 | type_code, type codes 0–7), missing/end-of-vector sentinel values, record binary layout (l_shared/l_indiv, 24-byte fixed fields, variable shared/individual sections), FILTER/INFO/FORMAT encoding, GT binary encoding ((allele+1) << 1 | phased), string dictionary. BGZF compression follows [SAM1] §4.1. See [references.md](99-references.md).

> *[BCF2] — magic `BCF\2\1`, l_text header length, NUL-terminated VCF header text*

r[bcf_writer.magic]
The BCF stream MUST begin with the 5-byte magic `BCF\x02\x02` (BCF version 2.2, the current version expected by htslib/bcftools), followed by `l_text: u32` (header text length including NUL), followed by the NUL-terminated VCF header text. The entire stream MUST be BGZF-compressed.

r[bcf_writer.coordinate_system]
BCF positions are 0-based. The writer MUST subtract 1 from the 1-based VcfRecord position when encoding. CHROM MUST be encoded as the contig's integer index from the header dictionary.

> *[BCF2] — record layout: l_shared, l_indiv, shared data (site), individual data (samples)*

r[bcf_writer.record_layout]
Each record MUST be encoded as `l_shared: u32`, `l_indiv: u32`, followed by `l_shared` bytes of shared data and `l_indiv` bytes of individual (per-sample) data. The conversion of buffer lengths to `u32` MUST use checked conversion (`u32::try_from`) and return a typed error on overflow rather than silently truncating.

r[bcf_writer.fixed_fields]
The first 24 bytes of shared data MUST contain: CHROM (i32), POS (i32, 0-based), rlen (i32), QUAL (f32), n_info|n_allele (u32 packed: n_allele << 16 | n_info), n_fmt|n_sample (u32 packed: n_fmt << 24 | n_sample). All multi-byte values are little-endian.

> *[BCF2] — typed value encoding: type byte (count << 4 | type_code), type codes 1=int8, 2=int16, 3=int32, 5=float, 7=char*

r[bcf_writer.typed_values]
All variable-length fields MUST use BCF typed value encoding. The type byte encodes `(count << 4) | type_code`. Type codes: 0=missing, 1=int8, 2=int16, 3=int32, 5=float, 7=char. For count >= 15, the type byte uses count=15 followed by a typed integer with the actual count.

r[bcf_writer.smallest_int_type]
Integer values MUST be encoded using the smallest type that fits all concrete (non-missing) values: int8 for [-120, 127], int16 for [-32760, 32767], int32 otherwise. The 8 most-negative values of each type are reserved for sentinels. Missing values within integer arrays MUST use the per-type sentinel matching the selected type (int8=0x80, int16=0x8000, int32=0x80000000), NOT a fixed i32::MIN.

> *[BCF2] — missing and end-of-vector sentinel values table. See also [hts-specs issue #145](https://github.com/samtools/hts-specs/issues/145) on signaling NaN handling*

r[bcf_writer.missing_sentinels]
Missing values MUST use the type-specific sentinel: int8=0x80, int16=0x8000 (LE), int32=0x80000000 (LE), float=0x7F800001 (LE, signaling NaN). Float sentinels MUST be written as raw bytes via bit manipulation, never through float arithmetic (signaling NaN would be silently converted to quiet NaN).

r[bcf_writer.end_of_vector]
When samples have variable-length values for a FORMAT field, shorter vectors MUST be padded with end-of-vector sentinels: int8=0x81, int16=0x8001, int32=0x80000001, float=0x7F800002, char=0x00 (NUL).

r[bcf_writer.shared_variable]
After the 24 fixed bytes, shared data MUST contain in order: ID (typed string, `.` if missing), alleles (n_allele typed strings, REF first), FILTER (typed integer vector of dictionary indices; PASS=[0], not-applied=type 0/count 0), INFO key-value pairs (typed int key = dictionary index, typed value).

> *[BCF2] — FORMAT/sample encoding: field-major layout, key + type + n_sample values per FORMAT field*

r[bcf_writer.indiv_field_major]
Individual data MUST be encoded in field-major order: for each FORMAT field, emit the key (typed int = dictionary index), the type descriptor, then all n_sample values contiguously. This layout is NOT sample-major.

> *[BCF2] — GT encoding: `(allele+1) << 1 | phased`, missing allele = 0x00*

r[bcf_writer.gt_encoding]
Genotype alleles MUST be encoded as `(allele_index + 1) << 1 | phased_bit`. Missing allele = 0x00. Alleles are stored as int8/int16/int32 depending on the maximum allele index. Haploid samples in a diploid context MUST be padded with end-of-vector sentinels.

r[bcf_writer.filter_pass]
PASS MUST always be assigned BCF dictionary index 0. A record that passes all filters MUST encode FILTER as a single-element typed int8 vector containing 0. A record with filters not applied MUST encode FILTER with type=0, count=0.

r[bcf_writer.qual_missing]
Missing QUAL MUST be encoded as the float missing sentinel (0x7F800001 as raw LE bytes). Present QUAL MUST be encoded as IEEE 754 single-precision float in little-endian.

r[bcf_writer.buffer_reuse]
The writer MUST reuse internal buffers for shared_data and indiv_data across records to avoid per-record allocation. Buffers are cleared (not deallocated) before each record.

r[bcf_writer.bgzf_blocks]
The writer MUST flush to a new BGZF block when the uncompressed buffer approaches 64 KB. A single BCF record SHOULD NOT span block boundaries when possible (for indexing compatibility).

r[bcf_writer.finish]
`finish()` MUST flush all pending data through the BGZF layer and write the 28-byte BGZF EOF marker block.

r[bcf_writer.string_encoding]
Strings in BCF MUST be encoded as typed char vectors (type code 7). They are NOT NUL-terminated within the typed value encoding. Empty strings MUST have count=0.

r[bcf_writer.flag_encoding]
Flag-type INFO fields MUST be encoded with type=0, count=0 (the key's presence indicates the flag is set). Absent flags are simply not included in the INFO key-value pairs.
