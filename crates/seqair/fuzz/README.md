# Fuzz Testing

Fuzz targets for all seqair file format readers and parsers, using [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz) (libFuzzer).

## Quick start

```sh
# Run all targets (30s each, default)
cd crates/seqair && ./fuzz/run_all.sh

# Shorter smoke test
./fuzz/run_all.sh 10

# Single target with seeds
cargo +nightly fuzz run fuzz_reader_indexed fuzz/corpus/fuzz_reader_indexed fuzz/seeds/fuzz_reader_indexed

# Multi-threaded
THREADS=8 ./fuzz/run_all.sh 60

# x86_64 via Docker (tests SSSE3 SIMD paths)
docker build --platform linux/amd64 -f fuzz/Dockerfile.x86_64 -t seqair-fuzz-x86 ../..
docker run --platform linux/amd64 --tmpfs /tmp --rm seqair-fuzz-x86 \
  "cd crates/seqair && ./fuzz/run_all.sh 30"
```

## Targets

### Full-stack (reader → fetch → pileup)

| Target | Description | Coverage | Seeds |
|--------|-------------|----------|-------|
| `fuzz_reader_indexed` | BAM/CRAM + index → `FuzzReaders` → `fetch_into` → pileup | 1273 | BAM+BAI, CRAM+CRAI |
| `fuzz_readers_pileup` | Full `FuzzReaders::pileup()` with FASTA reference | 981 | BAM+BAI+FASTA+FAI+GZI |
| `fuzz_reader_bam` | BGZF cursor → header → sequential records → pileup | 505 | test.bam |
| `fuzz_structured_pileup` | Arbitrary-generated valid BAM records → pileup | 918 | — |
| `fuzz_pileup_full` | Random bytes → `push_raw` → pileup | 141 | BAM records |

### BAM parsers

| Target | Description | Seeds |
|--------|-------------|-------|
| `fuzz_bam_record` | `BamRecord::decode` + `RecordStore::push_raw` | BAM records |
| `fuzz_bam_index` | BAI binary index parsing | test.bam.bai |
| `fuzz_bam_fields` | CIGAR ops, sequence decode (NEON/SSSE3), aux tags | — |
| `fuzz_bgzf` | BGZF block decompression | BGZF blocks |
| `fuzz_sam_header` | `BamHeader::from_sam_text` | SAM header lines |

### CRAM parsers

| Target | Description | Seeds |
|--------|-------------|-------|
| `fuzz_cram_block` | Block parsing + all compression codecs | 3 CRAM files |
| `fuzz_cram_container` | Container header parsing | 3 CRAM files |
| `fuzz_cram_codecs` | rANS 4x8, rANS Nx16, tok3 | noodles vector |
| `fuzz_cram_header` | Compression header + slice header + encodings | 3 CRAM files |
| `fuzz_cram_varint` | ITF8/LTF8 variable-length integers | 6 encodings |
| `fuzz_cram_bitstream` | BitReader operations | — |
| `fuzz_cram_encoding_decode` | HuffmanTable + ExternalCursor | — |
| `fuzz_cram_substitution` | SubstitutionMatrix | — |
| `fuzz_cram_decode_full` | CompressionHeader → encodings → DecodeContext | CRAM data |

### FASTA / index parsers

| Target | Description | Seeds |
|--------|-------------|-------|
| `fuzz_fasta_index` | FAI index parsing | 2 FAI files |
| `fuzz_gzi` | GZI index parsing | trimmed + full |

### Type-level

| Target | Description | Seeds |
|--------|-------------|-------|
| `fuzz_base_convert` | `Base::convert_ascii_in_place` / `from_ascii_vec` | — |
| `fuzz_region_string` | `RegionString::from_str` | 5 patterns |

## Seeds

Seeds in `fuzz/seeds/<target>/` are tracked in git. Auto-generated corpus goes to `fuzz/corpus/` (gitignored).

Regenerate seeds from test data:
```sh
cargo run --example gen_indexed_seeds    # for fuzz_reader_indexed
python3 fuzz/generate_seeds.py           # for all other targets
```

## Bugs found

| # | Bug | Severity | File | Found by |
|---|-----|----------|------|----------|
| 1 | CIGAR `matches += len` overflow panic | panic | `bam/cigar.rs` | `fuzz_bam_record` |
| 2 | SEGV in NEON/SSSE3 seq decode on short buffers | memory safety | `bam/seq.rs` | `fuzz_bam_fields` |
| 3 | OOM in CRAM block decompression (unbounded alloc) | DoS | `cram/block.rs` | `fuzz_cram_block` |
| 4 | OOM in container landmark array (unbounded alloc) | DoS | `cram/container.rs` | `fuzz_cram_container` |
| 5 | OOM in rANS/tok3 output buffers (unbounded alloc) | DoS | `cram/rans*.rs`, `cram/tok3.rs` | `fuzz_cram_codecs` |
| 6 | rANS state multiply overflow | panic | `cram/rans.rs`, `cram/rans_nx16.rs` | `fuzz_cram_codecs` |
| 7 | CRAM slice `matching_bases += 1` overflow | panic | `cram/slice.rs` | proactive audit |
| 8 | Huffman canonical code shift overflow | panic | `cram/encoding.rs` | `fuzz_cram_encoding_decode` |
| 9 | BAM header `l_text` / `n_ref` unbounded alloc | DoS | `bam/header.rs` | `fuzz_reader_bam` |
| 10 | RegionBuf OOM from corrupt BAI chunk ranges | DoS | `bam/region_buf.rs` | `fuzz_readers_pileup` |
| 11 | CRAM container negative length → capacity overflow | panic | `cram/reader.rs` | `fuzz_reader_indexed` |
| 12 | RegionBuf `.sum()` overflow on merged ranges | panic | `bam/region_buf.rs` | `fuzz_reader_indexed` |
| 13 | RegionBuf per-range alloc overflow | panic | `bam/region_buf.rs` | `fuzz_reader_indexed` |
| 14 | RegionBuf `debug_assert` panic on corrupt BGZF bsize | panic | `bam/region_buf.rs` | `fuzz_reader_indexed` |
| 15 | CRAI index query add overflow on large alignment_start + span | panic | `cram/index.rs` | `fuzz_reader_indexed` |
| 16 | OOM: validate BGZF ISIZE in RegionBuf::read_block | DoS | `bam/region_buf.rs` | `fuzz_readers_pileup` |

### Proactive fixes (applied based on patterns found)

- All `i32`/`i64` → `usize` casts in CRAM/BAM replaced with `usize::try_from().map_err()`
- All accumulation loops in CRAM slice reconstruction use `saturating_add`
- Allocation size cap (256 MiB) added to all parsed-size-driven allocations
- `debug_assert!`s added for all `#[allow(clippy::indexing_slicing)]` claimed invariants

## Architecture

The fuzz crate has a shared library (`src/lib.rs`) for types used by both fuzz targets and seed generators:

```
fuzz/
├── src/
│   ├── lib.rs                    # shared crate root
│   └── indexed_reader.rs         # Input/Format types for fuzz_reader_indexed
├── examples/
│   └── gen_indexed_seeds.rs      # seed generator using Input::encode
├── fuzz_targets/                 # 23 fuzz target binaries
├── seeds/                        # hand-crafted seeds (git-tracked)
├── corpus/                       # auto-generated corpus (gitignored)
├── artifacts/                    # crash/oom artifacts (gitignored)
├── generate_seeds.py             # seed generator for non-indexed targets
├── run_all.sh                    # CI runner
└── Dockerfile.x86_64             # x86 cross-fuzzing
```
