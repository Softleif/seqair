# VCF Text Format Rules

> **Sources:** [VCF43] §1.3 "Data lines" — tab-delimited format, field order, missing value `.` syntax; §1.3.1 "Fixed fields" — CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO serialization; §1.3.2 "Genotype fields" — FORMAT:sample colon-separated encoding, GT allele/phase syntax; §1.0.2 "Character encoding" — percent-encoding of special characters. BGZF compression for `.vcf.gz` follows [SAM1] §4.1. See [References](./99-references.md).

These rules apply to the VCF text arm of the unified [`Writer`](./5-record-encoder.md). The `Writer` selects VCF text encoding when constructed with `OutputFormat::Vcf` or `OutputFormat::VcfGz`.

r[vcf_writer.output_formats]
The unified `Writer` MUST support three output modes: uncompressed VCF (plain `io::Write`), BGZF-compressed VCF (`.vcf.gz` via `BgzfWriter`), and BCF binary. `OutputFormat::from_path` MUST return an error for unrecognized file extensions rather than silently defaulting.

> _[VCF43] §1.3 "Data lines" — tab-delimited, 8 fixed columns + FORMAT + samples_

r[vcf_writer.tab_delimited]
Data lines MUST be tab-delimited with exactly 8 fixed columns (CHROM through INFO) plus FORMAT and one column per sample when genotype data is present. Lines MUST end with `\n`.

r[vcf_writer.missing_dot]
Missing values MUST be serialized as `.`: missing QUAL, missing ID, missing INFO (entire column), missing ALT, and missing individual sample values.

r[vcf_writer.info_serialization]
INFO fields MUST be serialized as semicolon-separated `key=value` pairs. Flag-type fields emit the key only (no `=`). Multiple values within a field MUST be comma-separated. An empty INFO column MUST be written as `.`.

r[vcf_writer.format_serialization]
FORMAT keys MUST be colon-separated. Per-sample values MUST be colon-separated in the same order as FORMAT keys. Trailing missing values at the end of a sample MAY be omitted.

> _[VCF43] §1.3.2 "Genotype fields" — GT: `allele[sep allele]_`, sep is `/`or`|`, missing allele `.`\*

r[vcf_writer.genotype_serialization]
GT values MUST be serialized as allele indices separated by `/` (unphased) or `|` (phased). Missing alleles MUST be `.`. Examples: `0/1`, `1|0`, `./.`, `0|0|1` (triploid).

r[vcf_writer.float_precision]
Float values MUST be written exactly as C's `%g` with 6 significant digits writes them, which is what htslib and bcftools emit, so seqair's VCF text can be diffed against theirs. Concretely: trailing zeros after the decimal point and a bare trailing decimal point MUST be omitted; the form is scientific when the value's decimal exponent — *after* rounding to 6 significant digits, so `0.0001` is `0.0001` and not `1e-04` — is below `-4` or at least `6`, and fixed otherwise; and a scientific exponent MUST carry a sign and at least two digits (`1e-05`, `1.23457e+06`). Negative zero MUST be written `-0`. Every finite value MUST be writable — p-values in INFO fields reach magnitudes whose fixed form would not fit six significant digits in any fixed-size buffer, and the scientific form always does.

r[vcf_writer.integer_format]
Integer values MUST be written as decimal without leading zeros. Negative values MUST use `-` prefix. The `itoa` crate or equivalent fast formatting SHOULD be used for performance.

> _[VCF43] §1.0.2 "Character encoding" — percent-encoding for `:;=%,` TAB/LF/CR_

r[vcf_writer.percent_encoding]
Special characters in field values MUST be percent-encoded per the VCF spec: `:` → `%3A`, `;` → `%3B`, `=` → `%3D`, `%` → `%25`, `,` → `%2C`, TAB → `%09`, LF → `%0A`, CR → `%0D`.

r[vcf_writer.buffer_reuse]
The writer MUST reuse an internal line buffer across records to avoid per-record allocation. The buffer is cleared (not deallocated) before each record.

r[vcf_writer.finish]
`finish()` MUST flush all buffered data and, for BGZF output, write the EOF marker block and finalize the TBI index.
