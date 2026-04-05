# VCF Record

> **Sources:** [VCF43] §1.3 "Data lines" — fixed fields CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO (§1.3.1), genotype fields and FORMAT column (§1.3.2), Number semantics (A/R/G/.), missing value `.` convention, genotype ordering. [BCF2] — 0-based POS, rlen computation, n_allele counting, sentinel values. See [references.md](99-references.md).

> *[VCF43] §1.3.1 "Fixed fields" — CHROM, POS (1-based), ID, REF, ALT, QUAL, FILTER, INFO*

r[vcf_record.fields]
A VCF record MUST contain: contig name, 1-based position, reference allele, and zero or more alternate alleles. Optional fields: ID, QUAL, FILTER, INFO, FORMAT + sample data.

> *[VCF43] §1.3.1 — POS is 1-based. [BCF2] — BCF POS field is 0-based*

r[vcf_record.pos_one_based]
Positions in VcfRecord are 1-based, matching VCF text format convention. Writers converting to BCF MUST subtract 1 to produce the 0-based BCF POS.

> *[VCF43] §1.3.1 — missing value `.`. [BCF2] — typed missing sentinels per type (int8=0x80, float=0x7F800001)*

r[vcf_record.missing_values]
Missing scalar values MUST be represented as `None` in the record model. Writers MUST serialize `None` as `.` (VCF text) or the appropriate sentinel value (BCF binary: int8=0x80, int16=0x8000, int32=0x80000000, float=0x7F800001).

r[vcf_record.filters]
Filters MUST be one of three states: Pass (all filters passed), Failed (list of filter IDs), or NotApplied (filters not run). Writers MUST serialize Pass as `PASS`, Failed as semicolon-separated IDs, and NotApplied as `.`.

r[vcf_record.info_fields]
INFO fields MUST be stored as key-value pairs where the key matches a declared INFO ID in the header. Flag-type INFO fields carry no value. Writers MUST validate that values match the declared Number and Type.

r[vcf_record.format_gt_first]
If GT is present in the FORMAT keys, it MUST be the first key. The writer MUST enforce this ordering.

> *[VCF43] §1.3.2 "Genotype fields" — GT allele indexing (0=REF, 1+=ALT), `/` unphased, `|` phased, `.` missing*

r[vcf_record.genotype_encoding]
Genotypes MUST store per-allele indices (0=REF, 1+=ALT, None=missing) and per-separator phasing (true=phased `|`, false=unphased `/`). The phase bit on the first allele separator is ignored by convention but MUST be written as unphased.

r[vcf_record.sample_count]
The number of sample value sets MUST equal the number of samples declared in the header. Writers MUST return an error on mismatch.

r[vcf_record.builder]
Records MUST support construction via a builder pattern with required fields (contig, pos, ref_allele) and optional fields. The builder MUST validate field consistency at build time.

r[vcf_record.rlen]
The reference length (rlen, used in BCF) MUST be computed as `ref_allele.len()` for concrete alleles. For symbolic alleles with an END INFO field, rlen MUST be `END - POS + 1`.

r[vcf_record.allele_count]
`n_allele` counts ALL alleles including REF. A record with REF=A and ALT=T has n_allele=2. A record with no ALT (monomorphic/reference site) has n_allele=1.

r[vcf_record.smallvec]
Alt alleles, INFO fields, FORMAT keys, and per-sample values MUST use `SmallVec` with inline capacity tuned for common cases (2 alt alleles, 8 INFO fields, 6 FORMAT keys) to avoid heap allocation for typical records.

## Type-safe alleles

> *[VCF43] §1.3.1 — REF: reference base(s). ALT: comma-separated list of alternate non-reference alleles. Rules for representing SNVs, insertions, deletions, and complex events via REF/ALT padding conventions.*

r[vcf_record.alleles_typed]
Alleles MUST be represented via a type-safe enum with variants for common VCF allele patterns: reference-only (REF with no ALT), SNV (single-base REF and ALT), insertion (anchor base + inserted bases), deletion (anchor base + deleted bases), and complex (arbitrary REF/ALT strings for MNVs, symbolic alleles, or mixed multi-allelic sites). Construction MUST validate structural invariants: SNV alts must differ from ref, insertion/deletion sequences must be non-empty.

r[vcf_record.alleles_rlen]
The reference length (rlen) MUST be derived from the alleles variant: 1 for Reference, SNV, and Insertion; `1 + deleted.len()` for Deletion; `ref_allele.len()` for Complex. This replaces manual rlen computation.

r[vcf_record.alleles_serialization]
VCF REF/ALT text serialization MUST be handled by the alleles type. Each variant produces the correct VCF encoding: Insertion emits `anchor` as REF and `anchor+inserted` as ALT; Deletion emits `anchor+deleted` as REF and `anchor` as ALT; Reference emits the ref base with `.` as ALT.
