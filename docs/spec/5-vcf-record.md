# VCF Record Types

> **Sources:** [VCF43] §1.3 "Data lines" — fixed fields, genotype fields. [BCF2] — 0-based POS, rlen computation, n_allele counting. See [References](./99-references.md).

## Genotype

> _[VCF43] §1.3.2 "Genotype fields" — GT allele indexing (0=REF, 1+=ALT), `/` unphased, `|` phased, `.` missing_

r[vcf_record.gt_first_phase]
The first allele of a `GT` carries a phase bit like any other, and from VCF 4.4
it is read rather than ignored (4.3 and earlier "should treat the phasing of the
first allele as missing"). The encoder MUST set it per the spec's rule for the
leading indicator: `|` when every separator in the genotype is phased, `/`
otherwise. Hardcoding it unphased makes htslib render a fully phased `0|1` as
`/0|1` once the header declares 4.4 or later.

A reader MAY omit a leading indicator it can infer, so a partially phased
genotype round-trips as `0|1/1` rather than `/0|1/1`; both spell the same
genotype and only the fully phased case distinguishes the bit.

r[vcf_record.genotype_encoding]
Genotypes MUST store per-allele indices (0=REF, 1+=ALT, None=missing) and per-separator phasing (true=phased `|`, false=unphased `/`). The phase bit on the first allele separator is ignored by convention but MUST be written as unphased.

## Type-safe alleles

> _[VCF43] §1.3.1 — REF: reference base(s). ALT: comma-separated list of alternate non-reference alleles._

r[vcf_record.alleles_typed]
Alleles MUST be represented via a type-safe enum with variants for common VCF allele patterns: reference-only (REF with no ALT), SNV (single-base REF and ALT), insertion (anchor base + inserted bases), deletion (anchor base + deleted bases), and complex (arbitrary REF/ALT strings for MNVs, symbolic alleles, or mixed multi-allelic sites). Construction MUST validate structural invariants: SNV alts must differ from ref, insertion/deletion sequences must be non-empty.

r[vcf_record.alleles_rlen]
The reference length (rlen) MUST be derived from the alleles variant: 1 for Reference, SNV, and Insertion; `1 + deleted.len()` for Deletion; `ref_allele.len()` for Complex.

r[vcf_record.alleles_serialization]
VCF REF/ALT text serialization MUST be handled by the alleles type. Each variant produces the correct VCF encoding: Insertion emits `anchor` as REF and `anchor+inserted` as ALT; Deletion emits `anchor+deleted` as REF and `anchor` as ALT; Reference emits the ref base with `.` as ALT.

r[vcf_record.allele_count]
`n_allele` counts ALL alleles including REF. A record with REF=A and ALT=T has n_allele=2. A record with no ALT (monomorphic/reference site) has n_allele=1.

r[vcf_record.pos_one_based]
Positions are 1-based, matching VCF text format convention. The BCF encoder MUST subtract 1 to produce the 0-based BCF POS.

## Validation

r[vcf_record.sample_count]
The encoder MUST validate that the number of samples provided for each FORMAT field matches the sample count declared in the VCF header. A mismatch MUST be reported as a typed error.

r[vcf_record.format_gt_first]
When a GT (genotype) FORMAT field is present, it MUST be the first FORMAT key. If GT appears at any other position, the encoder MUST reject the record with a typed error.
