# Test data provenance

## `pairhmm-testdata.txt`

104 read/haplotype vectors with GATK's own expected log10 likelihoods, from Intel's
Genomics Kernel Library (<https://github.com/IntelLabs/GKL>, Apache-2.0), which in
turn took them from the Genome Analysis Toolkit
(<https://github.com/broadinstitute/gatk>, Apache-2.0/BSD).

## `10s.in`

The "10s" PairHMM benchmark: 7 region groups, 439 reads, 3,550 read-haplotype
pairs, 62,380,634 matrix cells. Read lengths 10-247 (mean 55), haplotype lengths
41-263 (mean 211). Base qualities are already floor-clamped at 6, the
gap-continuation penalty is 10 everywhere, and reads contain `N`.

Taken from <https://github.com/philipc/gkl-rs> (Apache-2.0), `tests/data/`, which
derives it from Intel GKL and GATK as above.

**`10s.out` is deliberately not vendored.** It ships next to `10s.in` in gkl-rs but
no test or bench in that repository reads it, and `10s.in` is shared there with the
Smith-Waterman bench. Checked against an independent `f64` DP, `10s.out` scores the
long reads under the *old* GATK mismatch prior `eps` rather than today's `eps/3`
(the gap is an exact multiple of log10(3) per mismatched cell), and its short-read
scores match no pairing of this input at all — a 19 bp read it scores at -21.50
scores -2.47 under every convention tried. It is not a conformance target for the
PairHMM. `10s.in` is used here for timing only; correctness is pinned by
`pairhmm-testdata.txt`, by `align_full`, and on x86_64 by gkl itself.

This is the same dataset the gpuPairHMM paper (Schmidt et al., arXiv:2411.11547,
Table 2) and the Endeavor paper (Graca & Ilic, HPDC '26, arXiv:2606.25738) report
CPU and GPU times on, which is why it is vendored here: it is the only published
PairHMM benchmark small enough to ship and widely enough used to compare against.
