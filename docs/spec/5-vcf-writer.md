# VCF Text Writer

> **Sources:** [VCF43] §1.3 "Data lines" — tab-delimited format, field order, missing value `.` syntax; §1.3.1 "Fixed fields" — CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO serialization; §1.3.2 "Genotype fields" — FORMAT:sample colon-separated encoding, GT allele/phase syntax; §1.0.2 "Character encoding" — percent-encoding of special characters. BGZF compression for `.vcf.gz` follows [SAM1] §4.1. See [references.md](99-references.md).

r[vcf_writer.output_formats]
The text writer MUST support three output targets: uncompressed VCF (plain `io::Write`), BGZF-compressed VCF (`.vcf.gz` via `BgzfWriter`), and stdout.

r[vcf_writer.header_first]
`write_header()` MUST be called exactly once before any `write_record()` call. Calling `write_record` without a prior `write_header` MUST return an error.

> *[VCF43] §1.3 "Data lines" — tab-delimited, 8 fixed columns + FORMAT + samples*

r[vcf_writer.tab_delimited]
Data lines MUST be tab-delimited with exactly 8 fixed columns (CHROM through INFO) plus FORMAT and one column per sample when genotype data is present. Lines MUST end with `\n`.

r[vcf_writer.missing_dot]
Missing values MUST be serialized as `.`: missing QUAL, missing ID, missing INFO (entire column), missing ALT, and missing individual sample values.

r[vcf_writer.info_serialization]
INFO fields MUST be serialized as semicolon-separated `key=value` pairs. Flag-type fields emit the key only (no `=`). Multiple values within a field MUST be comma-separated. An empty INFO column MUST be written as `.`.

r[vcf_writer.format_serialization]
FORMAT keys MUST be colon-separated. Per-sample values MUST be colon-separated in the same order as FORMAT keys. Trailing missing values at the end of a sample MAY be omitted.

> *[VCF43] §1.3.2 "Genotype fields" — GT: `allele[sep allele]*`, sep is `/` or `|`, missing allele `.`*

r[vcf_writer.genotype_serialization]
GT values MUST be serialized as allele indices separated by `/` (unphased) or `|` (phased). Missing alleles MUST be `.`. Examples: `0/1`, `1|0`, `./.`, `0|0|1` (triploid).

r[vcf_writer.float_precision]
Float values MUST be written with sufficient precision to round-trip through parsing (at least 6 significant digits). Trailing zeros after the decimal point SHOULD be omitted for compactness.

r[vcf_writer.integer_format]
Integer values MUST be written as decimal without leading zeros. Negative values MUST use `-` prefix. The `itoa` crate or equivalent fast formatting SHOULD be used for performance.

> *[VCF43] §1.0.2 "Character encoding" — percent-encoding for `:;=%,` TAB/LF/CR*

r[vcf_writer.percent_encoding]
Special characters in field values MUST be percent-encoded per the VCF spec: `:` → `%3A`, `;` → `%3B`, `=` → `%3D`, `%` → `%25`, `,` → `%2C`, TAB → `%09`, LF → `%0A`, CR → `%0D`.

r[vcf_writer.buffer_reuse]
The writer MUST reuse an internal line buffer across records to avoid per-record allocation. The buffer is cleared (not deallocated) before each record.

r[vcf_writer.finish]
`finish()` MUST flush all buffered data and, for BGZF output, write the EOF marker block. The inner writer MUST be returned to allow the caller to perform additional operations (e.g., creating an index).

r[vcf_writer.validation]
The writer MUST validate that record contig names and INFO/FORMAT/FILTER IDs are declared in the header. Undeclared fields MUST produce a typed error, not silent output.
