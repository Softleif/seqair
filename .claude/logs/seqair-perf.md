# Where rastair's time goes on the seqair backend

Companion to `rastair2/.claude/logs/seqair-perf.md`. That file measures the two
backends end to end; this one measures seqair on its own, which is the only way
to tell a seqair cost from a rastair one.

- Machine: macOS arm64 (Apple Silicon), 14 cores, 48 GB.
- Data: `NA12878_aa_chr12.bam` (1.6 GB compressed / 11.6 GB decompressed).
- Harness: `cargo run --release --example tiled_pileup -- <bam> --region
  chr12:1-20000000 --tile 10000`, single-threaded, which is rastair's access
  pattern (10 kb tiles, 200 bp overlap) with a trivial consumer.

## The shape of the seqair path

20 Mb of chr12, 2000 tiles, 4.4 M records, 20.3 M columns, 638 M
alignment-entries:

| phase | wall | share |
| --- | ---: | ---: |
| fetch (index query → inflate → decode → store) | 1.90 s | 33 % |
| prepare (sort + `link_mates`) | 0.04 s | **0.7 %** |
| columns (`PileupEngine::pileups`) | 3.80 s | 66 % |

Two things follow immediately.

**Mate linking is not a cost.** It was the most suspicious recent addition — a
hash table built per store, a qname hash per record, a byte comparison per
probe — and it is 0.7 % of the path. `qname_hash` itself is 0.3 % of total
samples. Nothing here is worth tuning; a faster hash function would buy
nothing measurable.

**The column loop is the cost**, and it is not the consumer's fault: running
the same harness with `--columns-only` (produce every column, never look at an
alignment) costs 3.61 s against 3.80 s. The engine spends 3.6 s producing what
the consumer reads in 0.19 s.

## Inside the column loop

`samply record --save-only --unstable-presymbolicate` + `tools/samply_hot.py`,
then leaf addresses joined against `llvm-objdump -d --line-numbers`. Shares are
of the whole run.

| share | what |
| ---: | --- |
| 56.2 % | `PileupEngine::pileups` (everything inlines into it) |
| 22.6 % | `libdeflate_deflate_decompress_ex` |
| 5.9 % | `CigarMapping::deletion_after_at` — **not inlined**, called per match column per read |
| 4.0 % | the harness's own consumer |
| 2.2 % | `build_decode_table` (libdeflate, per BGZF block) |
| 1.7 % | `RecordStore::push_raw` |
| 0.8 % | `BamIndex::query_split` |
| 0.7 % | `prepare_for_pileup` |
| 0.3 % | `qname_hash` |

And inside `pileups`, by source line:

| share | line | what |
| ---: | --- | --- |
| 10.7 % | `slice/iter/macros.rs:180` | iterating `&self.active` |
| 9.4 % | `vec/mod.rs:1873` | `Vec::push` (capacity check + len update) |
| 4.4 % | `ptr/mod.rs:1939` | writing the `PileupAlignment` |
| 4.1 % | `raw_vec/mod.rs:620` | the push's buffer pointer |
| 3.6 % + 3.2 % | `ptr/mod.rs:549`, `vec/mod.rs:1044` | the rest of the push |
| 6.0 % | `cigar.rs:466-469` | `pos_info_at`, `Linear` arm |
| 1.7 % | `pileup.rs:1092` | the alignment loop body |
| 1.3 % | `pileup.rs:931` | the eviction retain |

**~26 % of the entire run is `Vec::push` of a `PileupAlignment`.** The struct
is 56 bytes and is written once per read per column — 638 M times for 20 Mb of
chr12 — and most of its fields (`mapq`, `flags`, `seq_len`, `matching_bases`,
`indel_bases`, `record_idx`, `qname_hash`, `mate_idx`) are properties of the
*record*, not of the position: they are copied out of `ActiveRecord`, which the
engine is already holding, into every column entry.

## What that means end to end (measured in rastair)

Whole chr12, `--gpu`, `-@ 8`, hyperfine N=3, same machine state:

| backend | `--segment-max-length` | wall | user | sys |
| --- | ---: | ---: | ---: | ---: |
| seqair | 10 000 (default) | 30.34 s | 218.1 s | 20.1 s |
| seqair | 100 000 | **20.41 s** | 172.0 s | 5.7 s |
| htslib | 10 000 (default) | 33.55 s | 257.5 s | 15.0 s |
| htslib | 100 000 | 26.66 s | 228.5 s | 5.8 s |

seqair is **1.09× faster than htslib at the default tile size and 1.31× at
100 kb** — the tile size, not the backend, is what most of the current gap in
either direction is made of. Records are byte-identical across tile sizes
(`bcftools view -H | md5`, chr12:20–30 Mb).

The reason seqair gains more (1.49× vs htslib's 1.26×) is that it pays more per
query: a BAI `query_split`, a `RegionBuf` plan, and a window read, all of which
happen 13 k times per chromosome at 10 kb and 1.3 k times at 100 kb.

## Where the tile-size win actually comes from

Profiles of the same 30 Mb region at both tile sizes, worker threads only:

| | 10 kb | 100 kb |
| --- | ---: | ---: |
| `__semwait_signal` | **18.2 %** | 4.6 % |
| `libsystem_kernel` (total) | 22.0 % | 5.6 % |
| `libdeflate_deflate_decompress_ex` | 5.6 % | 3.5 % |
| `BamIndex::query_split` | 2.3 % | 0.7 % |

Most of it is workers *blocked*, not workers working: rastair dispatches its ML
batch per region and waits for it, so a 10 kb region buys one dispatch latency
for a tenth of the rows. Only the smaller part (`query_split`, inflate) is
seqair's per-query cost. That is a rastair-side finding, but it dominates
everything else measured here.

## Worker CPU at 100 kb, the configuration worth optimising

| share | function |
| ---: | --- |
| 22.6 % | `rastair …::get_pileups` (rastair's `from_seqair`, inlined) |
| 11.3 % | `PileupEngine::pileups` |
| 9.9 % | `rastair …::drops_overlapping_mate` |
| 8.1 % | `_platform_memmove` |
| 6.1 % | `add_fields` |
| 5.2 % | `rastair …::observed` |
| 4.3 % | `entropy_at` |
| 3.5 % | `libdeflate_deflate_decompress_ex` |
| 2.6 % | `build_indel_observation` |
| 1.9 % | `counted_mate` |
| 1.0 % | `deletion_after_at` |

seqair is ~18 % of worker CPU; rastair's per-column consumer is ~48 %.

## Where it ended up

Four changes, three of them here and one in rastair's segment default:

| build | chr12 wall | vs htslib |
| --- | ---: | ---: |
| seqair, session start | 30.34 s | 1.09× |
| **seqair, now** | **17.79 s** | **1.51×** |
| htslib, now | 26.91 s | — |

The seqair-side share of that: dropping the column entry's `qname_hash`,
`#[inline]` on `deletion_after_at`, and a lookup table for `Base::known_index`
took `examples/tiled_pileup` from 5.74 s to 5.10 s over 20 Mb, the column phase
from 3.80 s to 3.25 s.

`benches/bam.rs::pileup_tiled` reproduces the tiling cost in-repo, on
`tests/data/test.bam`, without needing a 1.6 GB BAM:

| tile | time | vs one query |
| ---: | ---: | ---: |
| 1 kb | 20.8 ms | 2.71× |
| 10 kb | 9.28 ms | 1.21× |
| 100 kb | 7.76 ms | 1.01× |
| one query | 7.68 ms | — |

## What is left, in order

1. **The column entry's remaining per-record fields.** `mapq`, `flags`,
   `seq_len`, `matching_bases`, `indel_bases`, `record_idx`, `mate_idx` are all
   in `ActiveRecord` already and get copied into every column entry. Behind
   `AlignmentView` accessors the entry goes from 48 bytes to ~16;
   `Vec::push` of it is ~26 % of the pure read+pileup path and the active-set
   walk another 11 %. Breaking change; rastair follows mechanically.
2. **Inflate is ~23 % of the path** and the residual over-read is the tiling's,
   shared with htslib. A small cross-query decompressed-block cache would take
   the tile-boundary blocks, which is most of what is left.
3. Nothing else here is above a few percent.
