# VCF Header

> **Sources:** [VCF43] §1.1 "Meta-information lines" (fileformat, INFO/FILTER/FORMAT/contig field definitions, structured meta-information syntax), §1.2 "Header line syntax" (#CHROM line and sample columns). BCF dictionary mapping follows [BCF2] (string-to-index assignment for FILTER/INFO/FORMAT IDs, contig integer indices). See [References](./99-references.md).

> _[VCF43] §1.1.1 "File format" — `##fileformat=VCFvX.Y` required first line_

r[vcf_header.file_format+2]
Every header MUST carry its version as a typed, ordered `VcfVersion` (`major`, `minor`, with constants `V4_0` … `V4_5`), exposed by `VcfHeader::version()`, so that "is this header new enough" is a comparison and not a string parse. `to_vcf_text()` MUST emit `##fileformat=VCFv{major}.{minor}` as the first line. There is no string-returning accessor.

r[vcf_header.file_format.single_version+2]
A header built with the builder MUST declare `VcfVersion::V4_5` (`VcfHeader::WRITE_VERSION`), with no setter. A cardinality is versioned — `Number=M` is defined by VCF 4.5 and by no earlier version — so a built header that could declare an older version could promise a grammar its fields violate. Emitting one version removes that class of mismatch; backwards compatibility with readers older than 4.5 is not a goal.

A parsed header (`r[vcf_header.parse]`) MUST keep the version its file declared, and writing it back out MUST declare that version: its definitions are the file's own, and the file's version is the one that describes them. Rules whose meaning depends on the version — the first allele's phase bit (`r[vcf_record.gt_first_phase]`), the reference span (`r[vcf_values.rlen]`) — MUST take it from the header of the record being read.

r[vcf_header.builder]
Headers MUST be constructed via a typestate builder that enforces field registration order at compile time. The builder progresses through phases — Contigs → Filters → Infos → Formats → Samples — and each phase only exposes methods appropriate for that stage. Phase transitions consume the builder and return the next phase; phases may be skipped. `build()` is available from any phase. Duplicate-ID and type-constraint errors are returned as typed errors.

> _[VCF43] §1.1.7 "Contig field format" — `##contig=<ID=name,length=N>`. [BCF2] — contig lines required for BCF, define integer mapping_

r[vcf_header.contig_required]
Every contig referenced by a record MUST be declared in the header via `##contig=<ID=name[,length=N]>`. For BCF output, contig declarations define the integer-to-name mapping (index = insertion order, 0-based). Looking up a contig that was never declared MUST return `VcfHeaderError::MissingContig` naming the contig that was asked for, not a generic lookup failure — the caller's next move is usually to add that contig, so the name has to survive the error.

> _[VCF43] §1.1.2 "Information field format" — ID, Number (A/R/G/./integer), Type (Integer/Float/Flag/Character/String), Description_

r[vcf_header.info_def]
Each INFO field MUST be declared with ID, Number, Type, and Description. Number is one of: a fixed count, `A` (one per ALT), `R` (one per allele), `G` (one per genotype), or `.` (unknown/variable). Type is one of: Integer, Float, Flag, Character, String. Flag type requires Number=0.

> _[VCF43] §1.1.4 "Individual format field format" — same Number/Type as INFO except no Flag; GT must be first_

r[vcf_header.format_def]
Each FORMAT field MUST be declared with ID, Number, Type, and Description. Same Number/Type rules as INFO except Flag type is not permitted. GT, if declared, MUST always appear first in the FORMAT column.

r[vcf_header.filter_def]
Each FILTER MUST be declared with ID and Description. PASS MUST always be present (implicitly or explicitly) and MUST map to BCF dictionary index 0.

> _[VCF43] §1.2 "Header line syntax" — sample names in #CHROM line must be unique_

r[vcf_header.sample_names]
Sample names MUST be unique. The order of samples in the header defines the order of sample columns in records.

r[vcf_header.no_duplicates]
Duplicate IDs within the same field category (INFO, FORMAT, FILTER, contig) MUST be rejected with a typed error at build time.

r[vcf_header.string_map]
For BCF output, the header MUST provide a string-to-index dictionary mapping all FILTER, INFO, and FORMAT IDs to integer indices. The dictionary is built incrementally during builder construction; the typestate phase ordering guarantees entries appear in canonical order (PASS at index 0, then remaining FILTERs, then INFO, then FORMAT) without a runtime verification pass. Contig names use a separate index namespace matching insertion order.

r[vcf_header.dict_capacity+2]
The FILTER/INFO/FORMAT string dictionary MUST hold up to `MAX_DICT_ENTRIES` = 65 536 entries, so every index fits a `u16` `DictId`. Real headers carry far more than the 64 entries an earlier version allowed (gnomAD 4.1 declares 244 ids, the 1000 Genomes high-coverage panel 83). Registering (writer) or parsing (reader) a field past the limit MUST return `TooManyFields { max }`. The limit is the same for both paths (`r[io.writer_limits]`). The cap does not apply to the separate contig namespace.

r[vcf_header.field_set]
Per-record "already seen" tracking keyed by `DictId` (the writer's INFO and FORMAT duplicate trackers, the reader's warn-once sets) MUST use one crate-internal `FieldSet`: a bitset over `Box<[u64]>` with `insert(id) -> bool` (true when the id was not yet present), `contains`, `clear` (zeroes the words, keeping the allocation) and an iterator over the set ids. An id beyond the current capacity MUST grow the set rather than panic or fail, so a `FieldSet` needs no header to be constructed. No dependency: no crate measured faster for these operations, and inline-storage variants (`SmallVec<u64, N>`, an inline/heap enum) measured 3–5× slower at small sizes because the storage branch on every access keeps the words out of registers.

r[vcf_header.serialization]
`to_vcf_text()` MUST emit all meta-information lines (`##`) followed by the `#CHROM` header line. Lines MUST be ordered: fileformat first, then other lines (metadata), then FILTER (PASS first), INFO, FORMAT, contig, then the #CHROM line. This ordering ensures the BCF string dictionary (built by scanning header lines in order) assigns PASS to index 0 and matches the dictionary indices used during encoding.

r[vcf_header.from_bam_header]
A `from_bam_header(&BamHeader)` constructor MUST copy contig names and lengths from the BAM header's @SQ lines, preserving order (and thus tid mapping).

r[vcf_header.smolstr]
All string fields (IDs, descriptions, contig names, sample names) MUST use `SmolStr` to avoid heap allocation for short strings.

r[vcf_header.ordered_maps]
INFO, FORMAT, FILTER, and contig maps MUST preserve insertion order. This is required for BCF dictionary index assignment and for deterministic VCF output.

## Parsing

A header read from a file is the same `VcfHeader` type the writer builds, so a reader's header can be handed to a writer (`r[vcf_reader.passthrough_header]`).

> _[VCF45] §1.4 "Meta-information lines", §1.5 "Header line syntax"; [BCF2] "Dictionary of strings", "Dictionary of contigs"_

r[vcf_header.parse]
`VcfHeader::parse(text: &str)` MUST parse VCF header text — the `##` meta-information lines followed by the `#CHROM` line — as found at the start of a VCF file or in a BCF header's `l_text` block (trailing NUL bytes stripped). Structured lines (`##INFO=<…>`, `##FORMAT=<…>`, `##FILTER=<…>`, `##contig=<…>`, and any other `##key=<…>`) MUST be parsed as comma-separated `key=value` pairs where a value may be double-quoted and a quoted value may contain `\"` and `\\` escapes, commas and `=`. The `#CHROM` line MUST be tab-separated; the eight fixed column names MUST match exactly; a `FORMAT` column followed by sample names is optional.

r[vcf_header.parse.version]
The first line MUST be `##fileformat=VCFvMAJOR.MINOR`. It MUST be parsed into a `VcfVersion` (`r[vcf_header.file_format+2]`). A missing or malformed line, or a major version other than 4, MUST fail with a typed error. A 4.x minor version newer than 4.5 MUST be accepted with a `warn!`, read with the 4.5 rules.

r[vcf_header.parse.defs]
INFO and FORMAT lines MUST yield `InfoDef`/`FormatDef` with a parsed `Number` (integer, `A`, `R`, `G`, `M`, `.`) and `ValueType`. An unknown `Number` or `Type` token MUST fail with a typed error naming the field. Constraints the writer enforces but real files violate (a `Flag` with a `Number` other than 0, a `Number` that disagrees with the spec's reserved definition) MUST be accepted with a `warn!`, since rejecting them would refuse files htslib and bcftools read.

r[vcf_header.parse.dictionary]
The string dictionary MUST be assigned as htslib assigns it: `PASS` at index 0; every other FILTER, INFO and FORMAT id in order of first appearance across those lines, regardless of category; an id declared by both an INFO and a FORMAT line shares one index. An `IDX=n` attribute (htslib writes it into BCF headers) MUST be honoured and MAY leave gaps. Two ids claiming the same `IDX`, or one id given two different `IDX` values, MUST fail with a typed error. Contigs have their own dictionary: order of appearance, or `IDX` when present.

r[vcf_header.parse.limits]
Parsing MUST enforce practical limits before allocating (`r[io.fuzz.alloc_limits]`): header text ≤ 256 MiB, dictionary ≤ `MAX_DICT_ENTRIES`, contigs ≤ 1 048 576, samples ≤ 10 000 000. Duplicate sample names MUST fail with a typed error.

r[vcf_header.parse.preserve]
Lines the model does not interpret MUST be kept verbatim, in order, so `to_vcf_text()` of a parsed header reproduces them. When a parsed header's dictionary order differs from the canonical emission order of `r[vcf_header.serialization]`, or has gaps, serializing it MUST emit `IDX=` attributes so a re-read assigns the same indices.

r[vcf_header.parse.round_trip]
For every header the builder can produce, `VcfHeader::parse(&h.to_vcf_text())` MUST yield a header equal to `h` in its definitions, dictionary, contigs and samples. Tests MUST also parse headers bcftools wrote (with and without `IDX=`) and compare the dictionary with the one bcftools uses.

r[vcf_header.brand]
Every `VcfHeader` MUST carry a `HeaderId` distinct from every other header's in the process, drawn from an atomic counter at `build()` or `parse()`. Keys resolved from a header carry its `HeaderId` (`r[vcf_values.header_brand]`). Cloning a header keeps its id.
