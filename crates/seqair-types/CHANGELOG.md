# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0 (2026-05-08)

### Other

 - <csr-id-25d1791b0eceb3caa47cbba27690c6116fd4e04b/> post-BamRecord-removal cleanup
   Re-attach the spec rules orphaned by the BamRecord -> SlimRecord
   migration to their new code homes (BamFlags predicates, parse_header,
   decode_bases_into, RecordStore::{seq_at,aux}, push_raw SIMD path,
   matching/indel pre-compute fields). Rename remaining BamRecord prose
   references to OwnedBamRecord across the BAM-record-builder/writer/flags
   specs, point r[bam_writer.index_coproduction] at the new
   BamWriterBuilder::write_index API, and replace the BamRecord mention in
   r[perf.precompute_matches_indels] / base-quality field-usage with
   SlimRecord. Drop the contradictory "no universally correct value"
   paragraph from SegmentOptions::new now that Default exists, and note on
   BaseModState::from_record that we reject MM-with-deltas records lacking
   ML (stricter than SAM 1.6). Extract the duplicated push-raw test
   scaffold into bam::test_util.
 - <csr-id-ea65683bebe1edec40250801bc4042ca9e8845a0/> fix warnings in BAM/FASTA/types test modules
 - <csr-id-0643d430fe7174ebdb23e88ae56ac97902e5b1b4/> fix cast warnings in seqair-types (phred, pos)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 37 commits contributed to the release.
 - 3 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#1](https://github.com/Softleif/seqaire/issues/1)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#1](https://github.com/Softleif/seqaire/issues/1)**
    - Fixes after review ([`536263b`](https://github.com/Softleif/seqaire/commit/536263b1d779a5c491d29bbb036708723f65ee6e))
 * **Uncategorized**
    - Add readmes ([`983b64b`](https://github.com/Softleif/seqaire/commit/983b64bdcec1962268f4195aab1053e0c0ce5676))
    - Post-BamRecord-removal cleanup ([`25d1791`](https://github.com/Softleif/seqaire/commit/25d1791b0eceb3caa47cbba27690c6116fd4e04b))
    - Fallible Segment::new, buffer reuse, typed errors, conigs ([`e52ddb6`](https://github.com/Softleif/seqaire/commit/e52ddb6539f2acead01c1f8d87c6dc3889e7d52e))
    - Bump versions, use crates.io rust-htslib ([`4c21ac2`](https://github.com/Softleif/seqaire/commit/4c21ac29f46a95507e82ffb5c09c34082ba6c47c))
    - Fix doc refs ([`d35f27b`](https://github.com/Softleif/seqaire/commit/d35f27b039b70c447360e77e863e86669ee8efd7))
    - Three follow-ups from /critique: TLEN sentinel, set_filter removed, oracles ([`6f9a4a6`](https://github.com/Softleif/seqaire/commit/6f9a4a61ed4733045ef4590f96f78b4d6ce28abd))
    - Rework pileup + Readers API around lending iterator and extras provider ([`5138b3f`](https://github.com/Softleif/seqaire/commit/5138b3f42117472d582651c9d7a47a47238f9899))
    - Add #[non_exhaustive] to all public error enums, fix PathBuf in VcfError ([`cf1c926`](https://github.com/Softleif/seqaire/commit/cf1c9265fbbc59f5eb626b731e3b08232cd1e43f))
    - Add missing #[must_use] to Pos conversion and arithmetic methods ([`ebe6f6c`](https://github.com/Softleif/seqaire/commit/ebe6f6c54d9b870f757c7d98e163a6fb442bd4d7))
    - Tracey spec review: fix spec/code mismatches and annotation gaps ([`3422266`](https://github.com/Softleif/seqaire/commit/34222669970f9f1209a0172196d197a2392e0253))
    - Add BaseQuality newtype for 'quality unavailable' sentinel ([`d24dc9c`](https://github.com/Softleif/seqaire/commit/d24dc9c8009dc48fe7c0fc2650a883f87652650f))
    - Simplify pos ([`e6b7bc1`](https://github.com/Softleif/seqaire/commit/e6b7bc154eedaaf07d4372bd817494178bc273eb))
    - Fix second review round: region_string panic, TryFrom<u32>, checked_sub_offset ([`9e9e9bb`](https://github.com/Softleif/seqaire/commit/9e9e9bb8fc13d79754aef2fd09138f6e087e8565))
    - Fix review findings: use as_i32() in CIGAR hot path, silence clippy ([`9a4f19a`](https://github.com/Softleif/seqaire/commit/9a4f19ab7aa58bfa206ba181fbd34b5c642a3560))
    - Constrain Pos to i32::MAX, add standard TryFrom impls, use Pos0/Pos1 aliases ([`50b9c15`](https://github.com/Softleif/seqaire/commit/50b9c157b140815eb362bdf51c6fab4656da5dfd))
    - Some fixes ([`260dc47`](https://github.com/Softleif/seqaire/commit/260dc47970abcf00c90b3b834db4e60fe8bb07e4))
    - Add BamFlags newtype to replace raw u16 flag fields ([`ce84bc3`](https://github.com/Softleif/seqaire/commit/ce84bc3809e870878cbea958b5614bbe36122ece))
    - Fix unfulfilled expect(unreachable_code) on x86_64 and u8-as-i8 SIMD casts ([`e567ade`](https://github.com/Softleif/seqaire/commit/e567ade7ec25f75c8f5622b18bffc5c174fe71da))
    - Fix warnings in BAM/FASTA/types test modules ([`ea65683`](https://github.com/Softleif/seqaire/commit/ea65683bebe1edec40250801bc4042ca9e8845a0))
    - Fix cast warnings in seqair-types (phred, pos) ([`0643d43`](https://github.com/Softleif/seqaire/commit/0643d430fe7174ebdb23e88ae56ac97902e5b1b4))
    - Auto clippy fixes ([`8b00840`](https://github.com/Softleif/seqaire/commit/8b008403aedd8152b64d2097d4618702afeb6a60))
    - Update stale tracey references and annotate 28 uncovered spec rules ([`ffa2585`](https://github.com/Softleif/seqaire/commit/ffa2585dba60d7c5bf10f34cfff4986073c95697))
    - Use proper smallvec ([`2a7a8a9`](https://github.com/Softleif/seqaire/commit/2a7a8a90b7a5629468a73b172f5da09a09f33f71))
    - Add Pos0/Pos1 aliases, clean docs, TBI index comparison tests ([`a8ece81`](https://github.com/Softleif/seqaire/commit/a8ece8191ab85651c5a4fccc73ee401972caef68))
    - Replace all unsafe arithmetic with checked_*, saturating_*, wrapping_* ([`881d3f1`](https://github.com/Softleif/seqaire/commit/881d3f1d7a60d6f7ba364ab07d11d2ebb2aed922))
    - Fix/allow some arithmetic side effects ([`b663e65`](https://github.com/Softleif/seqaire/commit/b663e6539bc1b99d119107703b69c8058bdef1e0))
    - Enable clippy::arithmetic_side_effects warn ([`9ed26f7`](https://github.com/Softleif/seqaire/commit/9ed26f7dbf721de44e5178c9ed4d3ae889e9dc13))
    - Add Tracey r[impl/verify] annotations to pos.rs ([`e20c8d8`](https://github.com/Softleif/seqaire/commit/e20c8d8e5882340f7e00b71798673157c4f3751f))
    - Remove panicking Add/Sub operator impls from Pos and Offset ([`61073fb`](https://github.com/Softleif/seqaire/commit/61073fbff3e0cf4269be9045fc7c463e3717112c))
    - Add tests for invalid position error paths ([`11c5e66`](https://github.com/Softleif/seqaire/commit/11c5e667be87b98c73a22ef50b607d2a40dec338))
    - Harden Pos<S>: no panics on untrusted data, proper error propagation ([`eb8301d`](https://github.com/Softleif/seqaire/commit/eb8301d388be32ace94fe0a3d79307b64e0f8a8d))
    - Fix Pos<S> soundness: remove all unsafe, add proptests ([`6a0a63d`](https://github.com/Softleif/seqaire/commit/6a0a63dc9dbdd9ed46aeb36276b0cce922cfed5c))
    - Complete Pos<S> migration: Pos<One> conversions, API cleanup, nonmax niche ([`c7e41a1`](https://github.com/Softleif/seqaire/commit/c7e41a1c07a536b0722f0119282fc6fc95d3dade))
    - Add Pos<S> genomic position newtype with coordinate system safety ([`23099c4`](https://github.com/Softleif/seqaire/commit/23099c42a62a9e0fa05e654c2f0673e04bebdd8f))
    - Proptest Audit: Tautological vs. Valid Tests ([`49e446e`](https://github.com/Softleif/seqaire/commit/49e446ed28c951770c381d35dfc7dab061f70d43))
    - Initial import from local repo ([`1681005`](https://github.com/Softleif/seqaire/commit/1681005cf8e5b97089e50b224b6efabf7566787a))
</details>

