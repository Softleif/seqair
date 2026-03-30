# VCF/BCF Writer Implementation Plan

## Goal

Add VCF (text) and BCF (binary) **writing** to seqair. This gives rastair (and other downstream consumers) a zero-dependency path from pileup columns to variant output files. Reading VCF/BCF is explicitly out of scope for now.

---

## Phase 0: BgzfWriter (prerequisite)

BCF files are BGZF-compressed, and bgzipped VCF (`.vcf.gz`) needs the same layer. seqair already has `BgzfReader` but no writer.

### Tasks

1. **`src/bam/bgzf_writer.rs`** — new module, re-exported from `bam.rs`
   - `BgzfWriter<W: Write>` wrapping a buffered inner writer
   - Accumulates uncompressed data in a 64 KB buffer
   - On flush/overflow: compress with `libdeflate::Compressor`, emit a valid BGZF block (gzip header with `BC` extra field, compressed payload, CRC32 + ISIZE trailer)
   - `finish()` writes the 28-byte EOF marker block and flushes
   - Returns `VirtualOffset` after each block write (needed for indexing later)
   - Error type: extend existing `BgzfError` with write variants (`CompressFailed`, `WriteFailed`)

2. **Tests**
   - Round-trip: write N blocks with `BgzfWriter`, read them back with `BgzfReader`, compare
   - EOF marker present and valid
   - Block boundary handling (data exactly 64 KB, data spanning two blocks)
   - Proptest: random byte payloads round-trip correctly

### Spec

- `docs/spec/1-2-bgzf.md` — add rules for `bgzf.writer.*`

---

## Phase 1: VCF Header Model

A shared header representation used by both VCF and BCF writers.

### Types

```
src/vcf.rs           — module declaration
src/vcf/header.rs    — VcfHeader, VcfHeaderBuilder, field definitions
src/vcf/error.rs     — VcfError hierarchy
```

### VcfHeader

```rust
pub struct VcfHeader {
    file_format: SmolStr,                          // "VCFv4.3"
    infos: IndexMap<SmolStr, InfoDef>,             // ordered by insertion
    formats: IndexMap<SmolStr, FormatDef>,
    filters: IndexMap<SmolStr, FilterDef>,
    contigs: IndexMap<SmolStr, ContigDef>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,                     // ##source, ##reference, etc.
    // BCF dictionary (built lazily or on demand)
    string_map: Option<StringMap>,
}
```

### Field definition types

```rust
pub enum Number {
    Count(u32),          // fixed count
    AlternateBases,      // A — one per ALT
    ReferenceAlternateBases, // R — one per allele (REF + ALTs)
    Genotypes,           // G — one per genotype combination
    Unknown,             // . — variable
}

pub enum ValueType { Integer, Float, Flag, Character, String }

pub struct InfoDef    { pub number: Number, pub typ: ValueType, pub description: SmolStr }
pub struct FormatDef  { pub number: Number, pub typ: ValueType, pub description: SmolStr }
pub struct FilterDef  { pub description: SmolStr }
pub struct ContigDef  { pub length: Option<u64> }
```

### Builder

```rust
VcfHeader::builder()
    .file_format("VCFv4.3")
    .add_contig("chr1", ContigDef { length: Some(248_956_422) })
    .add_info("DP", InfoDef { number: Number::Count(1), typ: ValueType::Integer, description: "Total Depth".into() })
    .add_format("GT", FormatDef { number: Number::Count(1), typ: ValueType::String, description: "Genotype".into() })
    .add_filter("PASS", FilterDef { description: "All filters passed".into() })
    .add_sample("sample0")
    .build()
```

### StringMap (for BCF)

Bidirectional mapping of header string IDs to integer indices. Built from the header's ordered maps. PASS is always index 0 in the filter namespace. Contig indices are a separate namespace (same as BAM tids).

### Error hierarchy

```rust
pub enum VcfError {
    Bgzf { #[from] source: BgzfError },
    Header { #[from] source: VcfHeaderError },
    Encode { #[from] source: VcfEncodeError },
    Io(#[from] std::io::Error),
}

pub enum VcfHeaderError {
    DuplicateContig { name: SmolStr },
    DuplicateInfo { id: SmolStr },
    DuplicateFormat { id: SmolStr },
    DuplicateFilter { id: SmolStr },
    DuplicateSample { name: SmolStr },
    MissingContig { name: SmolStr },       // record references undeclared contig
    MissingInfo { id: SmolStr },
    MissingFormat { id: SmolStr },
    MissingFilter { id: SmolStr },
}

pub enum VcfEncodeError {
    /// Value doesn't match declared Number/Type
    TypeMismatch { field: SmolStr, expected: ValueType, got: ValueType },
    /// Integer value out of BCF range
    IntegerOverflow { field: SmolStr, value: i64 },
}
```

### Spec

- `docs/spec/5-vcf-header.md` — rules for header construction, validation, field definitions

---

## Phase 2: VCF Record Model

A record type that can be serialized to both text VCF and binary BCF.

### VcfRecord (builder-style, owned)

```rust
pub struct VcfRecord {
    pub contig: SmolStr,
    pub pos: i64,                        // 1-based (VCF convention)
    pub id: Option<SmolStr>,             // "." if None
    pub ref_allele: SmolStr,
    pub alt_alleles: SmallVec<[SmolStr; 2]>,
    pub qual: Option<f32>,               // "." if None
    pub filters: Filters,
    pub info: InfoFields,
    pub samples: SampleFields,
}

pub enum Filters {
    Pass,
    Failed(SmallVec<[SmolStr; 2]>),
    NotApplied,                           // "." in VCF
}
```

### InfoFields

```rust
pub struct InfoFields {
    fields: SmallVec<[(SmolStr, InfoValue); 8]>,
}

pub enum InfoValue {
    Integer(i32),
    Float(f32),
    Flag,
    String(SmolStr),
    IntegerArray(SmallVec<[Option<i32>; 4]>),
    FloatArray(SmallVec<[Option<f32>; 4]>),
    StringArray(SmallVec<[Option<SmolStr>; 2]>),
}
```

### SampleFields

```rust
pub struct SampleFields {
    /// FORMAT keys in order (e.g., ["GT", "GQ", "DP", "AD", "PL"])
    format_keys: SmallVec<[SmolStr; 6]>,
    /// Per-sample values, outer = sample index, inner = format field index
    values: Vec<SmallVec<[SampleValue; 6]>>,
}

pub enum SampleValue {
    Missing,
    Integer(i32),
    Float(f32),
    String(SmolStr),
    Genotype(Genotype),
    IntegerArray(SmallVec<[Option<i32>; 4]>),
    FloatArray(SmallVec<[Option<f32>; 4]>),
}

pub struct Genotype {
    pub alleles: SmallVec<[Option<u16>; 2]>,  // None = missing allele
    pub phased: SmallVec<[bool; 2]>,          // per-separator phasing
}
```

### VcfRecordBuilder

Ergonomic builder for constructing records from pileup data:

```rust
VcfRecord::builder("chr1", 12345)
    .ref_allele("A")
    .alt_allele("T")
    .qual(30.0)
    .filter_pass()
    .info_integer("DP", 50)
    .info_float_array("AF", &[0.3])
    .info_flag("DB")
    .format_keys(&["GT", "GQ", "DP", "AD"])
    .sample(0, |s| s
        .genotype_unphased(&[0, 1])
        .integer("GQ", 40)
        .integer("DP", 30)
        .integer_array("AD", &[15, 15])
    )
    .build()
```

### Spec

- `docs/spec/5-vcf-record.md` — record structure, field validation, missing value semantics

---

## Phase 3: VCF Text Writer

### `src/vcf/writer.rs`

```rust
pub struct VcfWriter<W: io::Write> {
    inner: W,                    // could be BgzfWriter<File> for .vcf.gz
    header: Arc<VcfHeader>,
    buf: Vec<u8>,                // reusable line buffer
}
```

### API

```rust
impl<W: io::Write> VcfWriter<W> {
    pub fn new(inner: W, header: &VcfHeader) -> Result<Self, VcfError>;
    pub fn write_header(&mut self) -> Result<(), VcfError>;
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError>;
    pub fn finish(self) -> Result<W, VcfError>;
}
```

### Encoding rules

- Tab-delimited columns
- Missing values: `.`
- INFO: `key=value;key=value` (Flag: key only, no `=`)
- FORMAT/samples: colon-separated, trailing missing fields omittable
- GT: `allele[/|]allele`, missing allele = `.`
- Float formatting: sufficient precision, no trailing zeros (use `ryu` crate or manual)
- Integer formatting: use `itoa` crate for performance
- Percent-encoding for special characters in field values (`:`, `;`, `=`, `%`, `,`)

### Spec

- `docs/spec/5-vcf-writer.md` — text serialization rules, escaping, missing values

---

## Phase 4: BCF Binary Writer

### `src/vcf/bcf_writer.rs`

```rust
pub struct BcfWriter<W: io::Write> {
    bgzf: BgzfWriter<W>,
    header: Arc<VcfHeader>,
    string_map: StringMap,
    shared_buf: Vec<u8>,        // reusable buffer for shared section
    indiv_buf: Vec<u8>,         // reusable buffer for individual section
}
```

### API

Same as VcfWriter:

```rust
impl<W: io::Write> BcfWriter<W> {
    pub fn new(inner: W, header: &VcfHeader) -> Result<Self, VcfError>;
    pub fn write_header(&mut self) -> Result<(), VcfError>;
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError>;
    pub fn finish(self) -> Result<W, VcfError>;
}
```

### Binary encoding details

**Header:**
- Magic: `BCF\x02\x01` (5 bytes)
- `l_text: u32` — length of header text including NUL
- VCF header text, NUL-terminated

**Record layout:**
- `l_shared: u32` + `l_indiv: u32` (length prefixes)
- Fixed 24 bytes: CHROM (i32, 0-based tid), POS (i32, **0-based**), rlen (i32), QUAL (f32), n_info|n_allele (u32 packed), n_fmt|n_sample (u32 packed)
- Variable shared: ID (typed string), alleles (typed strings), FILTER (typed int vec), INFO key-value pairs (typed)
- Variable individual: FORMAT fields in field-major order (key, type, n_sample values)

**Typed value encoding:**
- Type byte: `(count << 4) | type_code`
- Type codes: 0=missing, 1=int8, 2=int16, 3=int32, 5=float, 7=char
- Count ≥ 15: type byte has count=15, followed by typed integer with actual count
- Missing sentinels: int8=0x80, int16=0x8000, int32=0x80000000, float=0x7F800001
- End-of-vector: int8=0x81, int16=0x8001, int32=0x80000001, float=0x7F800002

**GT encoding:**
- Per-allele: `(allele + 1) << 1 | phased`
- Missing: 0x00
- Stored as int8/int16/int32 depending on allele count

**Integer type selection:**
- Scan all values, pick smallest type that fits (int8 if all in [-120, 127], int16 if [-32760, 32767], else int32)

**Critical pitfalls:**
- VCF POS is 1-based, BCF POS is 0-based (subtract 1 when writing BCF)
- Float sentinels are signaling NaNs — write as raw bytes via `u32::to_le_bytes()`, never go through float arithmetic
- PASS filter = dictionary index 0, single-element vector `[0]`; not-applied filter = type=0/count=0

### Spec

- `docs/spec/5-bcf-writer.md` — binary encoding rules, type system, sentinel values

---

## Phase 5: Write Trait Unification

### `src/vcf/mod.rs`

```rust
pub trait VcfWrite {
    fn write_header(&mut self) -> Result<(), VcfError>;
    fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError>;
    fn finish(self) -> Result<(), VcfError>;
}
```

Both `VcfWriter` and `BcfWriter` implement this trait. Factory function:

```rust
pub fn open_writer(path: &Path, header: &VcfHeader) -> Result<Box<dyn VcfWrite>, VcfError> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("bcf") => Ok(Box::new(BcfWriter::from_path(path, header)?)),
        Some("gz") => Ok(Box::new(VcfWriter::bgzf(path, header)?)),
        _ => Ok(Box::new(VcfWriter::from_path(path, header)?)),
    }
}
```

---

## Module Layout

```
crates/seqair/src/
├── bam/
│   ├── bgzf.rs              (existing reader)
│   └── bgzf_writer.rs       (NEW — Phase 0)
├── vcf.rs                    (NEW — module declaration)
└── vcf/
    ├── header.rs             (NEW — Phase 1)
    ├── record.rs             (NEW — Phase 2)
    ├── writer.rs             (NEW — Phase 3, text VCF)
    ├── bcf_writer.rs         (NEW — Phase 4, binary BCF)
    └── error.rs              (NEW — Phase 1)

docs/spec/
├── 1-2-bgzf.md              (EXTEND — writer rules)
├── 5-vcf-header.md           (NEW)
├── 5-vcf-record.md           (NEW)
├── 5-vcf-writer.md           (NEW)
└── 5-bcf-writer.md           (NEW)
```

---

## Dependencies

| Crate | Purpose | Status |
|-------|---------|--------|
| `libdeflate` | BGZF compression (already used for decompression) | existing |
| `smallvec` | Inline allele/info/format arrays | existing |
| `smol_str` | Short string optimization | existing (via seqair-types) |
| `indexmap` | Ordered maps for header fields | **new** |
| `itoa` | Fast integer-to-string (VCF text) | **new, optional** |
| `ryu` | Fast float-to-string (VCF text) | **new, optional** |

---

## Implementation Order

| Phase | Deliverable | Depends on | Effort |
|-------|------------|------------|--------|
| 0 | BgzfWriter | — | Small |
| 1 | VcfHeader + errors | — | Medium |
| 2 | VcfRecord + builder | Phase 1 | Medium |
| 3 | VcfWriter (text) | Phase 1, 2 | Medium |
| 4 | BcfWriter (binary) | Phase 0, 1, 2 | Large |
| 5 | Write trait + factory | Phase 3, 4 | Small |

Phases 0 and 1 can be done in parallel. Phases 3 and 4 can be done in parallel once Phase 2 is complete.

---

## Testing Strategy

### Unit tests (per module)
- Header builder validation (duplicates, missing fields)
- Record builder (missing values, edge cases)
- BCF typed value encoding (all type codes, sentinels, overflow to wider types)
- BCF GT encoding (haploid, diploid, phased/unphased, missing)
- BgzfWriter round-trip

### Integration tests
- Write VCF, read back with `bcftools view` or noodles, compare
- Write BCF, read back with `bcftools view`, compare
- Write BCF, verify with `bcftools view -H` that all fields decode correctly

### Proptest
- Random VcfRecord → serialize to VCF text → parse back → compare
- Random VcfRecord → serialize to BCF → read back with noodles → compare
- BCF typed integer encoding: random i32 values → encode → decode → compare

### Comparison tests
- Generate VCF/BCF with rust-htslib and seqair from same input data, diff output

---

## Open Questions

1. **Indexing (TBI/CSI)**: out of scope for initial implementation? Can be added later since `bcftools index` can create indices externally.
2. **Multi-threaded writing**: not needed initially. BgzfWriter could support parallel compression later.
3. **VCF reading**: explicitly deferred. The header model should be designed to support future parsing, but no reader in this plan.
4. **From BamHeader**: should `VcfHeader` have a `from_bam_header(&BamHeader)` constructor that copies contig info? Seems useful.
5. **itoa/ryu**: add as dependencies or use `write!` formatting? Performance matters for large VCF files.
