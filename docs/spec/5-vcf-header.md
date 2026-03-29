# VCF Header

> **Sources:** [VCF43] §1.1 "Meta-information lines" (fileformat, INFO/FILTER/FORMAT/contig field definitions, structured meta-information syntax), §1.2 "Header line syntax" (#CHROM line and sample columns). BCF dictionary mapping follows [BCF2] (string-to-index assignment for FILTER/INFO/FORMAT IDs, contig integer indices). See [references.md](99-references.md).

> *[VCF43] §1.1.1 "File format" — `##fileformat=VCFvX.Y` required first line*

r[vcf_header.file_format]
The header MUST begin with a `##fileformat=VCFvX.Y` line. The default version is VCFv4.3.

r[vcf_header.builder]
Headers MUST be constructed via a builder pattern. The builder validates all constraints at `build()` time and returns a typed error on violation.

> *[VCF43] §1.1.7 "Contig field format" — `##contig=<ID=name,length=N>`. [BCF2] — contig lines required for BCF, define integer mapping*

r[vcf_header.contig_required]
Every contig referenced by a record MUST be declared in the header via `##contig=<ID=name[,length=N]>`. For BCF output, contig declarations define the integer-to-name mapping (index = insertion order, 0-based).

> *[VCF43] §1.1.2 "Information field format" — ID, Number (A/R/G/./integer), Type (Integer/Float/Flag/Character/String), Description*

r[vcf_header.info_def]
Each INFO field MUST be declared with ID, Number, Type, and Description. Number is one of: a fixed count, `A` (one per ALT), `R` (one per allele), `G` (one per genotype), or `.` (unknown/variable). Type is one of: Integer, Float, Flag, Character, String. Flag type requires Number=0.

> *[VCF43] §1.1.4 "Individual format field format" — same Number/Type as INFO except no Flag; GT must be first*

r[vcf_header.format_def]
Each FORMAT field MUST be declared with ID, Number, Type, and Description. Same Number/Type rules as INFO except Flag type is not permitted. GT, if declared, MUST always appear first in the FORMAT column.

r[vcf_header.filter_def]
Each FILTER MUST be declared with ID and Description. PASS MUST always be present (implicitly or explicitly) and MUST map to BCF dictionary index 0.

> *[VCF43] §1.2 "Header line syntax" — sample names in #CHROM line must be unique*

r[vcf_header.sample_names]
Sample names MUST be unique. The order of samples in the header defines the order of sample columns in records.

r[vcf_header.no_duplicates]
Duplicate IDs within the same field category (INFO, FORMAT, FILTER, contig) MUST be rejected with a typed error at build time.

r[vcf_header.string_map]
For BCF output, the header MUST provide a string-to-index dictionary mapping all FILTER, INFO, and FORMAT IDs to integer indices. Contig names use a separate index namespace matching insertion order.

r[vcf_header.serialization]
`to_vcf_text()` MUST emit all meta-information lines (`##`) followed by the `#CHROM` header line. Lines MUST be ordered: fileformat first, then INFO, FORMAT, FILTER, contig, other lines, then the #CHROM line.

r[vcf_header.from_bam_header]
A `from_bam_header(&BamHeader)` constructor MUST copy contig names and lengths from the BAM header's @SQ lines, preserving order (and thus tid mapping).

r[vcf_header.smolstr]
All string fields (IDs, descriptions, contig names, sample names) MUST use `SmolStr` to avoid heap allocation for short strings.

r[vcf_header.ordered_maps]
INFO, FORMAT, FILTER, and contig maps MUST preserve insertion order. This is required for BCF dictionary index assignment and for deterministic VCF output.
