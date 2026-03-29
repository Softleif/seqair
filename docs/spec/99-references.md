# References

The seqair spec is based on the following upstream specifications and reference implementations.

## Format specifications

- **[SAM1]** *Sequence Alignment/Map Format Specification*, The SAM/BAM Format Specification Working Group, last modified 2025-08-10. <https://samtools.github.io/hts-specs/SAMv1.pdf>
  Covers: SAM text format (§1), mandatory alignment fields and FLAG bits (§1.4), BAM binary encoding (§4 "The BAM Format Specification"), BGZF block compression (§4.1 "The BGZF compression format"), random access and virtual offsets (§4.1 "Random access"), EOF marker (§4.1 "End-of-file marker"), BAM record layout (§4.2 "The BAM format"), BIN field calculation (§4.2.1), SEQ/QUAL encoding (§4.2.4), auxiliary data encoding (§4.2.5), BAI index algorithm (§5 "Indexing BAM"), basic binning index (§5.1.1), BAI index format (§5.2 "The BAI index format for BAM files"), `reg2bin`/`reg2bins` C code (§5.3).

- **[CRAM3]** *CRAM format specification (version 3.1)*, samtools-devel@lists.sourceforge.net, last modified 2025-06-04. <https://samtools.github.io/hts-specs/CRAMv3.pdf>
  Covers: data types and ITF8/LTF8 variable-length integers (§2), bit stream reading (§2.2), encoding types (§3 "Encodings"), checksums (§4), file structure overview (§5), file definition and magic (§6), container header structure (§7), CRAM header container (§7.1), block structure (§8), block content types (§8.1), compression header block (§8.3) including preservation map, data series encodings, and substitution matrix, slice header block (§8.4), core data block (§8.5), external data blocks (§8.6), EOF container (§9), record structure and CRAM record decode order (§10), BAM and CRAM flags (§10.1), positional data and AP delta coding (§10.2), read names (§10.3), mate records (§10.4), auxiliary tags (§10.5), mapped reads and read features (§10.6), substitution codes (§10.6.1), unmapped reads (§10.7), reference sequences (§11), CRAM index / CRAI format (§12), encoding codec IDs (§13).

- **[CRAMcodecs]** *CRAM codec specification (version 3.1)*, samtools-devel@lists.sourceforge.net, last modified 2025-03-24. <https://samtools.github.io/hts-specs/CRAMcodecs.pdf>
  Covers: rANS 4×8 codec with order-0 and order-1 frequency tables and interleaving (§2), rANS Nx16 codec with transforms (RLE, PACK, STRIPE, CAT, N32) (§3), range/arithmetic coder (§4), name tokeniser tok3 (§5), fqzcomp quality codec (§6).

- **[VCF43]** *Variant Call Format Specification (version 4.3)*, The SAM/BAM Format Specification Working Group. <https://samtools.github.io/hts-specs/VCFv4.3.pdf>
  Covers: meta-information lines (§1.1 "Meta-information lines"): fileformat, INFO field format (§1.1.2), FILTER field format (§1.1.3), FORMAT field format (§1.1.4), ALT field format (§1.1.5), contig field format (§1.1.7), sample field format (§1.1.8); header line syntax (§1.2); data lines (§1.3): fixed fields CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO (§1.3.1), genotype fields and FORMAT column (§1.3.2), Number values (A/R/G/.), Type values (Integer/Float/Flag/Character/String), genotype ordering for GL/PL/GP/PP, reserved INFO keys, reserved FORMAT keys; character encoding and percent-escaping (§1.0.2 "Character encoding").

- **[BCF2]** *BCF2 Quick Reference (r198)*, The SAM/BAM Format Specification Working Group. <https://samtools.github.io/hts-specs/BCFv2_qref.pdf>
  Covers: BCF magic (`BCF\2\1`), header encoding (l_text + NUL-terminated VCF header text), typed value encoding (type byte layout: count << 4 | type_code, type codes 0–7), missing and end-of-vector sentinel values per type, record binary layout (l_shared/l_indiv, 24-byte fixed fields CHROM/POS/rlen/QUAL/n_info_allele/n_fmt_sample, variable shared and individual sections), allele encoding, FILTER encoding, INFO key-value encoding, FORMAT field-major sample encoding, genotype (GT) binary encoding ((allele+1) << 1 | phased), string dictionary mapping.

- **[TABIX]** *The Tabix index file format*, Heng Li. <https://samtools.github.io/hts-specs/tabix.pdf>
  Covers: TBI magic, header fields (n_ref, format, col_seq, col_beg, col_end, meta, skip, l_nm, names), bin/chunk/linear-index structure identical to BAI, `n_no_coor` trailing field, `reg2bin`/`reg2bins` C functions.

- **[CSI]** *Coordinate Sorted Index (CSI) format specification*, SAM/BAM Format Specification Working Group, last modified 2020-05-25. <https://samtools.github.io/hts-specs/CSIv1.pdf>
  Covers: CSI magic, header fields (min_shift, depth, l_aux), per-bin `loffset` field replacing BAI's linear index, parameterised `reg2bin`/`reg2bins`/`bin_limit` functions.

## Reference implementations

- **[htslib]** *htslib*, Wellcome Sanger Institute. <https://github.com/samtools/htslib>
  The canonical C implementation of SAM/BAM/CRAM/BGZF/BAI. Used to verify correctness of seqair's output: pileup positions, depths, qpos values, CRAM record fields, and flag handling are all cross-checked against htslib.

- **[noodles]** *noodles*, Michael Macias. <https://github.com/zaeleus/noodles>
  A pure-Rust implementation of genomics formats. Referenced for Rust-idiomatic design patterns in handling BAM/SAM/CRAM records.
