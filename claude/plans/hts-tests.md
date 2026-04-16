# htslib Test Suite Analysis & Seqair Test Parity Plan

## Part 1: How htslib's Test Suite Works

### Architecture

htslib's test suite has three layers:

1. **Perl orchestrator** (`test/test.pl`, 1700 lines) — integration tests that run htslib CLI tools (`htsfile`, `bgzip`) and compare outputs against expected files. 25 test functions covering format conversions, indexing, and round-trips. Supports multi-threaded variants (0 and 4 threads), fail-fast, output regeneration, and selective test filtering.

2. **C unit test programs** (14 binaries) — test individual modules: `test_bgzf.c` (36KB), `sam.c` (84KB, the largest — covers all aux tag types), `hfile.c`, `test_faidx.c`, `test_expr.c`, `test_khash.c`, `test_kstring.c`, `test_nibbles.c`, `test_realn.c`, `test_mod.c`, plus threading stress tests (`thrash_threads1-7`).

3. **Shell script suites** (7 subdirectories) — each with a `.tst` file listing test cases and a `simple_test_driver.sh` harness:
   - `faidx/` — FASTA index creation, retrieval, line length, BGZIP support
   - `fastq/` — FASTQ parsing, interleaved pairs, CASAVA headers (~50 files)
   - `mpileup/` — pileup generation, indels, overlapping reads (~30 files)
   - `base_mods/` — MM/ML tag parsing, 27 test cases
   - `sam_filter/` — filter expression evaluation
   - `tabix/` — BED/GFF/VCF indexing, large chromosomes
   - `tlen/` — template length validation, 65 CRAM/SAM pairs

### Test Data Inventory

| Category | Count | Key Files |
|----------|-------|-----------|
| SAM files | 56 | `ce#*.sam` (C. elegans), `c1#*.sam`, `xx#*.sam` (synthetic) |
| CRAM files | 4 | `*_java.cram` (cross-implementation from Java htsjdk) |
| BAM files | 3 | `colons.bam`, `range.bam`, `no_hdr_sq_1.bam` |
| VCF files | 14 | `vcf44_1.vcf`, `formatcols.vcf`, `modhdr.vcf.gz` |
| FASTA refs | 9 | `ce.fa`, `c1.fa`, `c2.fa`, `xx.fa`, `auxf.fa`, `md.fa` |
| Index files | 10+ | `.bai`, `.csi`, `.tbi`, `.crai`, `.gzi` |
| Annotation TSV | 64+ | `annot-tsv/` |

### What htslib Tests (Comprehensive List)

#### A. BGZF Compression
- Round-trip compression/decompression (single + multi-threaded)
- Text mode vs binary mode
- Block boundary extraction with offset (`bgzip -b`)
- Named index creation (`-I` flag)
- Multi-file compression
- Re-indexing (`bgzip -g -r`) for alternate deflate implementations
- Virtual offset operations (tell, seek)
- Character-level I/O (`bgzf_getc`)

#### B. SAM/BAM Parsing & Conversion
- SAM -> BAM -> SAM round-trip (byte-identical via `compare_sam.pl`)
- Header parsing (SQ lines, NUL padding 0/1/678 bytes, DOS `\r\n` line endings)
- Missing @SQ headers in BAM (targets present but no SQ lines)
- Auxiliary tags: all types (A, c, C, s, S, i, I, f, d, Z, H, B arrays)
- Aux tag operations: get, delete, next iteration, update_int, update_array, update_str
- Integer range validation (int8 through uint32, int64)
- Float arrays, string updates with variable lengths
- CIGAR operations (M/I/D/N/S/H/P/X/=)
- Large sequences (2.1 MB in `ce#large_seq.sam`)
- Large aux tags (`xx#large_aux*.sam`)
- Records spanning BGZF block boundaries (>64KB records with ~32KB seq + 32KB CIGAR)
- Empty/blank reads, sequence-less reads, unknown bases
- Unmapped reads (multiple variants)
- Paired-end reads, read groups, supplementary alignments
- Unsorted records

#### C. CRAM
- Versions: 2.1, 3.0, 3.1 (experimental), 4.0 (experimental)
- SAM -> CRAM -> SAM round-trips (all versions)
- BAM -> CRAM -> SAM
- CRAM -> CRAM re-encoding
- Multi-ref containers (`ref_seq_id == -2`)
- Multi-slice options (`-1`, `0`, `1`)
- Embedded references (`embed_ref=2`)
- Compression profiles: fast/normal/small/archive (v3.1, v4.0)
- Java-generated CRAM interoperability (`*_java.cram`)
- Template length (TLEN) validation: 65 CRAM/SAM paired files
- Container skipping during indexed access
- Invalid UR field handling (forces reference reloading)

#### D. Indexing & Range Queries
- BAI index creation and querying
- CSI index (resolution 14) for BAM, BCF, gzip SAM
- CRAI index for CRAM
- TBI (Tabix) index for VCF/GFF/BED
- Multi-region iterator
- Large positions (>2GB) with CSI
- Region specification: `CHR`, `CHR:POS`, `CHR:POS-POS`, `CHR:POS,CHR:POS`, exact `{CHR:POS}`
- Custom index names
- Records spanning container boundaries

#### E. FASTA Indexing & Retrieval
- Index creation (SAM spec format)
- Sequence retrieval by name and numeric ID
- Range queries (32-bit and 64-bit positions)
- Quality string retrieval (FASTQ)
- Line length detection
- Trailing blank handling
- BGZIP-compressed FASTA with GZI

#### F. VCF/BCF
- BCF <-> VCF conversion
- VCF API (read, write, manipulate)
- Sweep iterator
- Header parsing edge cases (excess spaces, custom IDX fields)
- Format column validation
- Missing value handling
- VCF 4.4 spec compliance (implicit/explicit phasing)
- BCF sync reader: merge multiple files, allele pairing, randomized (10 seeds)
- Sync reader without indexes (error detection for out-of-order)
- Range queries with unusual chromosome names
- Variant type classification

#### G. Pileup / Mpileup
- Coverage depth calculation
- Quality filtering
- Indel handling in pileup output
- Overlapping reads
- Base modification display in pileup
- Padded alignment pileup

#### H. Base Modifications (MM/ML Tags)
- MM tag parsing and validation (27 test cases)
- Explicit vs implicit modification calls
- Position/probability fields
- Multi-modification per base
- ChEBI code handling
- Pileup integration

#### I. MD Tag Calculation
- Store MD/NM without auto-generation
- Skip auto-generation; check only if present
- Force auto-generation from reference
- Perfect alignment (no MD/NM added)

#### J. BAQ (Base Alignment Quality)
- Calculate and apply BAQ to reads
- Extended BAQ calculation
- Recalculation from SAM
- Apply from existing BQ tags
- Revert quality from ZQ tags

#### K. Filter Expressions
- Arithmetic, comparison, logical, bitwise operators
- String regex matching (`=~`, `!~`)
- Functions: `log`, `exp`, `pow`, `sqrt`, `default`, `exists`
- NULL/NaN propagation
- Operator precedence

#### L. Low-Level Infrastructure
- Hash tables (`khash`), string manipulation (`kstring`)
- Nibble encoding/decoding
- String-to-integer conversion
- Time/date parsing
- Endianness detection
- hfile I/O abstraction (local, memory, data: URLs, preload)
- Plugin loading

#### M. Thread Safety
- 7 `thrash_threads` stress tests for concurrent access patterns
- All major operations tested at 0 and 4 threads

---

## Part 2: What Seqair Already Has

### Test Infrastructure
- **37 integration test files**, ~277 test functions
- **49 source files** with `#[cfg(test)]` unit tests
- **23 fuzz targets** with 6,389 corpus + 3,959 seed + 100 artifact files
- **27+ proptest suites** for property-based testing
- **13 comparison test files** validating against noodles and rust-htslib
- **12 test data files** (~23MB): `test.bam`, `test.cram` (v2.1, v3.0, gzip), `test.fasta.gz`, indexes

### Seqair's Strengths (vs htslib)
- **Property-based testing** (proptest) — far more systematic than htslib's hand-crafted cases
- **Triple oracle comparison** — noodles, rust-htslib, and bcftools as independent oracles
- **VCF/BCF round-trip testing** — deep write->read validation with multi-allelic, indels, missing values, phasing
- **Fuzz testing** — 23 targets from atomic codec level to full-stack pileup pipeline
- **Structure-aware fuzzing** — `Arbitrary`-generated valid BAM records for pileup

### Seqair's Current Gaps

#### Closed by Tier 1 (hts-eqv-tests branch)

| Gap | Status |
|-----|--------|
| SAM/BAM writing round-trips | **Done** — `bam_write_roundtrip.rs`: 6 tests (simple, complex CIGARs, aux tags, multi-contig, paired-end, index co-production) |
| BGZF block boundary records | **Done** — `bgzf_block_boundary.rs`: 3 tests (large seqs, long CIGARs, mixed sizes) |
| DOS line endings | **Done** — `htslib_edge_cases.rs`: both SAM reader and BAM conversion paths |
| Sequence-less reads | **Done** — `htslib_edge_cases.rs`: c1#noseq.sam, secondary with SEQ=* |
| Unknown bases | **Done** — `htslib_edge_cases.rs`: c1#unknown.sam SEQ/QUAL combos |
| Read groups | **Done** — `htslib_sam_parity.rs`: xx#rg.sam |
| Colon characters in names | **Done** — `htslib_edge_cases.rs`: colons.bam with chr1:100, chr1,chr3 |
| CIGAR clipping edge cases | **Done** — `htslib_sam_parity.rs`: c1#clip.sam (S, H, N, I ops) |
| Padding CIGAR ops | **Done** — `htslib_edge_cases.rs`: c1#pad1/pad2/pad3.sam |
| Unmapped read variants | **Done** — `htslib_edge_cases.rs`: ce#unmap.sam, supplementary/secondary tests |
| BAM writer index co-production | **Done** — `index_roundtrip.rs`: 5 tests, samtools idxstats validation |
| htslib test data import | **Done** — `htslib_sam_parity.rs`: 8 htslib SAMs parsed + compared vs noodles |
| CRAM version matrix | **Done** — `cram_version_matrix.rs`: v3.0 (3 tests), v3.1 (4 tests) |
| Multi-ref CRAM containers | **Done** — `cram_version_matrix.rs`: multi_seq_per_slice=1, ce#5b + ce#1000 |
| Embedded CRAM references | **Done** — `cram_version_matrix.rs`: embed_ref=2, verified vs samtools |
| Template length (TLEN) | **Done** — `cram_tlen.rs`: all 30 htslib pairs via BAM + 6 CRAM field checks |
| BAM writer stress tests | **Done** — `bam_writer_stress.rs`: max fields, poisoning, unmapped dispatch, dense records |
| VCF/BCF bcftools validation | **Done** — `vcf_bcftools_roundtrip.rs`: integer boundaries, genotypes, indels, stats |
| BGZF writer validation | **Done** — `bgzf_writer_validation.rs`: bgzip decompresses small/large/chunked/empty output |

#### Still Open

| Gap | htslib Coverage | Seqair Status | Priority |
|-----|----------------|---------------|----------|
| Large aux tags | `xx#large_aux.sam` (>64KB) | Not tested with dedicated edge cases | Medium |
| CRAM Java interop | `*_java.cram` files from htsjdk | Not tested — files available in tests/htslib/cram/ | Medium |
| **CRAM TLEN for attached mates** | htslib reconstructs TLEN for non-detached mates | **seqair returns 0** — only detached mates get TLEN | Bug/Gap |
| Missing @SQ text header | `no_hdr_sq_1.bam` (binary ref list but @CO in text) | Not tested — needs BAI creation | Low |
| MD tag generation | 5 test modes (store/skip/force) | Not tested — seqair doesn't generate MD | Low |
| Large positions (>2GB) | `longrefs/` directory with CSI indexes | Not tested — seqair uses BAI (not CSI) | Low |
| Pileup vs samtools mpileup | htslib mpileup/ dir with expected outputs | Not compared — files in tests/htslib/mpileup/ | High (Tier 4) |
| **Base modifications (MM/ML)** | 27 test cases in base_mods/ | **Not implemented** — feature missing | Future |
| **BAQ calculation** | 5 modes (calculate, extended, recalc, apply, revert) | **Not implemented** — feature missing | Future |
| **SAM filter expressions** | Full expression engine (arithmetic, regex, functions) | **Not implemented** — feature missing | Future |
| **FASTQ format** | ~50 test files (interleaved, CASAVA, filtering) | **Not implemented** — no FASTQ reader | Future |
| **Multi-threaded I/O** | All ops tested at 0 and 4 threads | **Single-threaded only** — no threading | Future |
| **CSI index format** | CSI indexes for BAM, BCF, SAM.gz | **Not implemented** — BAI only | Future |
| VCF sync reader | Randomized multi-file merge | Not applicable (no sync reader) | N/A |
| Reference caching (MD5) | Multi-tier cache server | Not applicable (no cache server) | N/A |
| Annotation TSV | 64+ test files | Not applicable | N/A |
| BGZF re-indexing | `bgzip -g -r` for alternate deflate impls | Not applicable | N/A |
| Plugin loading | Shared library plugin tests | Not applicable | N/A |
| Thread stress tests | 7 `thrash_threads` concurrent access patterns | Not applicable (single-threaded) | N/A |

---

## Part 3: Test Parity Plan

### Priority Tiers

**Tier 1 — Core Format Correctness (High Value, Moderate Effort)**

These tests validate that seqair reads/writes the same data as htslib. Failures here mean real bugs.

#### 1.1 BAM Write Round-Trip
**Goal**: SAM -> BAM (seqair) -> read back (noodles/htslib) produces identical records.
- Write synthetic BAM with all aux tag types (A, c, C, s, S, i, I, f, d, Z, H, B:c, B:C, B:s, B:S, B:i, B:I, B:f)
- Write BAM with edge-case CIGARs (long soft clips, hard clips + match, pure insertions)
- Write BAM -> read with noodles -> compare field-by-field
- Write BAM -> `samtools view` -> compare SAM text output
- Proptest: random OwnedBamRecord -> write -> read -> compare

#### 1.2 htslib Test Data Import
**Goal**: Parse htslib's own test files and compare output against htslib.
- Copy/reference key SAM files from htslib: `ce#1.sam`, `ce#5.sam`, `ce#5b.sam`, `xx#minimal.sam`, `xx#pair.sam`, `xx#blank.sam`, `xx#large_aux.sam`, `c1#clip.sam`, `c1#noseq.sam`, `c1#unknown.sam`, `ce#unmap.sam`
- Copy FASTA refs: `ce.fa`, `c1.fa`, `xx.fa`
- Parse each -> compare record fields against noodles
- Run pileup on each -> compare against `samtools mpileup` output

#### 1.3 BGZF Block Boundary Records
**Goal**: Correctly handle BAM records that straddle BGZF 64KB block boundaries.
- Generate synthetic BAM with ~32KB sequences + ~32KB CIGARs (matching htslib's approach)
- Verify all records parse correctly
- Verify index queries work across boundaries
- Verify virtual offsets are correct

#### 1.4 Index Round-Trip
**Goal**: Indexes we write are readable by samtools/noodles, and vice versa.
- Write BAM + BAI -> `samtools view -c` region queries match seqair's query counts
- Write VCF + TBI -> `bcftools view -r` matches seqair
- Parse htslib-generated indexes -> query -> compare results

#### 1.5 SAM Edge Cases from htslib
**Goal**: Handle the specific edge cases htslib tests for.
- DOS line endings (`\r\n`) in SAM
- Missing @SQ headers in BAM (targets exist but no SQ lines)
- Sequence-less reads (`*` in SEQ field)
- Unknown bases in sequences
- Colon characters in reference names
- Large positions (>2GB, requires CSI)

---

**Tier 2 — CRAM Completeness (High Value, Higher Effort)**

#### 2.1 CRAM Version Matrix
- v3.0 and v3.1 round-trips (we already have v2.1 and v3.0)
- Generate test CRAMs with `samtools view -C --output-fmt-option version=3.1`
- Multi-ref containers: create SAM with reads mapping to many contigs in one slice
- Embedded references: `samtools view -C -o embed.cram --embed-ref`

#### 2.2 CRAM Java Interoperability
- Download/generate CRAM files from htsjdk (Java implementation)
- Verify seqair decodes them identically to samtools
- This catches codec interpretation differences between implementations

#### 2.3 CRAM TLEN Validation
- Port htslib's 65 TLEN test pairs: various read lengths, orientations, overlaps
- SAM -> CRAM (samtools) -> read with seqair -> verify TLEN matches original SAM

---

**Tier 3 — Write Path Hardening (Medium Value, Moderate Effort)**

#### 3.1 BAM Writer Stress Tests
- Records at maximum field sizes (2 MiB record, 256 KiB name, etc.)
- Writer error poisoning: verify writes fail cleanly after an error
- Index co-production: mapped, placed-unmapped, fully-unmapped dispatch
- `samtools quickcheck` on every written BAM

#### 3.2 VCF/BCF Writer Expanded Tests
- Port htslib's VCF edge cases: excess spaces in headers, format column variations, missing values
- VCF 4.4 features (if applicable)
- BCF magic version validation (`BCF\x02\x02` not `\x02\x01`)
- Integer type selection edge cases (int8/int16/int32 boundaries with missing sentinels)

#### 3.3 BGZF Writer Tests
- Virtual offset at buffer boundary (exactly 65536 bytes)
- Multi-block writes with index building
- Empty block (EOF marker) correctness

---

**Tier 4 — Pileup Engine (Medium Value, Low Effort — Much Already Exists)**

#### 4.1 Pileup vs samtools mpileup
- Run `samtools mpileup` on htslib test SAMs -> compare column-by-column with seqair pileup
- Focus on indel representation, depth counting, base quality
- Edge cases: overlapping mate pairs, supplementary alignments, unmapped in region

#### 4.2 Pileup Edge Cases
- Zero-coverage regions (no reads)
- All-deletion columns
- Reads with only soft-clips in region
- Maximum depth handling
- Region boundaries (reads partially overlapping region start/end)

---

**Tier 5 — Features Not Yet Implemented (Future Work)**

These are features htslib tests but seqair doesn't implement yet. Track for future:

| Feature | htslib Tests | Notes |
|---------|-------------|-------|
| Base modifications (MM/ML) | 27 cases | New feature needed |
| MD tag generation | 5 modes | Useful for CRAM |
| BAQ calculation | 5 modes | Quality recalibration |
| FASTQ format | ~50 files | New reader needed |
| Filter expressions | Full engine | Not planned? |
| Multi-threaded I/O | 0+4 threads | Performance feature |
| Reference caching | Cache server | CRAM optimization |

---

## Part 4: Implementation Order

### Phase 1: Quick Wins (1-2 days)
1. Copy htslib SAM test files + FASTA refs into `tests/data/htslib/`
2. Write `compare_htslib_sams.rs` — parse each with seqair, compare against noodles
3. Add DOS line endings test
4. Add sequence-less reads test
5. Add unknown bases test
6. Add colon-in-name test

### Phase 2: Write Path (2-3 days)
1. BAM write round-trip integration tests
2. BAM write -> samtools quickcheck validation
3. BAM write -> noodles read comparison
4. BGZF block boundary test (synthetic large records)
5. Index write -> samtools region query comparison

### Phase 3: CRAM Expansion (2-3 days)
1. Generate CRAM v3.1 test files with samtools
2. Multi-ref container test
3. Embedded reference test
4. TLEN validation subset (pick ~10 representative pairs from htslib's 65)

### Phase 4: Cross-Validation Suite (1-2 days)
1. Pileup vs `samtools mpileup` on htslib test files
2. VCF/BCF writer edge cases from htslib
3. Large position (>2GB) tests

### Phase 5: Continuous Enrichment
- As new features land (MM/ML tags, FASTQ, etc.), add corresponding htslib parity tests
- Periodically re-sync with htslib test data for new edge cases
- Add htsjdk-generated CRAM files when available

---

## Appendix: htslib Test File Reference

### Files to Import (Priority)

```
# SAM files (core edge cases)
htslib/test/ce#1.sam          — C. elegans basic
htslib/test/ce#5.sam          — C. elegans multi-record
htslib/test/ce#5b.sam         — C. elegans variant
htslib/test/ce#large_seq.sam  — 2.1 MB large sequence
htslib/test/ce#unmap.sam      — unmapped reads
htslib/test/ce#supp.sam       — supplementary alignments
htslib/test/xx#minimal.sam    — minimal valid record
htslib/test/xx#pair.sam       — paired-end
htslib/test/xx#blank.sam      — empty reads
htslib/test/xx#large_aux.sam  — large aux tags
htslib/test/xx#rg.sam         — read groups
htslib/test/xx#unsorted.sam   — unsorted records
htslib/test/c1#clip.sam       — CIGAR clipping
htslib/test/c1#noseq.sam      — sequence-less
htslib/test/c1#unknown.sam    — unknown bases
htslib/test/c1#bounds.sam     — boundary testing
htslib/test/c1#pad1.sam       — padding variation 1
htslib/test/c1#pad2.sam       — padding variation 2
htslib/test/c1#pad3.sam       — padding variation 3
htslib/test/index_dos.sam     — DOS line endings
htslib/test/auxf#values.sam   — all aux tag types

# FASTA references
htslib/test/ce.fa + ce.fa.fai
htslib/test/c1.fa + c1.fa.fai
htslib/test/xx.fa + xx.fa.fai
htslib/test/auxf.fa + auxf.fa.fai

# BAM edge cases
htslib/test/no_hdr_sq_1.bam   — missing SQ headers
htslib/test/colons.bam         — colons in names
htslib/test/range.bam          — range queries

# CRAM interop
htslib/test/auxf#values_java.cram   — Java-encoded
htslib/test/ce#5b_java.cram         — Java-encoded
htslib/test/xx#large_aux_java.cram  — Java large aux

# Indexes
htslib/test/range.bam.bai
htslib/test/colons.bam.bai
```

### htslib Test Execution Flow
```
make check
 |-- C unit tests (hts_endian, test_expr, test_kstring, ...)
 |-- test_bgzf, test_faidx, fieldarith, hfile
 |-- Shell suites: faidx/, sam_filter/, tabix/, mpileup/, fastq/, base_mods/, tlen/
 +-- test.pl (Perl orchestrator)
      |-- test_bgzip (0 + 4 threads)
      |-- test_index (0 + 4 threads)
      |-- test_view (0 + 4 threads): SAM/BAM/CRAM conversions
      |-- test_multi_ref (0 + 4 threads)
      |-- test_MD
      |-- test_vcf_* (6 VCF test functions)
      |-- test_bcf_sr_* (4 BCF sync reader functions)
      |-- test_realn, test_rebgzip, test_convert_padded_header
      |-- test_annot_tsv, test_ref_cache, test_plugin_loading
      +-- test_logging
```
