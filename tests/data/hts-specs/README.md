# hts-specs CRAM test data

Copied verbatim from [samtools/hts-specs](https://github.com/samtools/hts-specs)
`test/cram/` at commit `da617203a9527537746e200abda2885bec3a822c` (2026-06-04).

- `cram/codecs/` — every CRAM 3.0/3.1 block codec on its own, next to the
  gzip-compressed original. See `cram/codecs/README.md` for the naming scheme
  and how to compare the decoded output. Apache 2.0 / BSD (`cram/codecs/LICENSE`).
- `cram/3.1/passed/level-{1,2,3,4}.cram` — the same records at successive
  compression levels; level 3 adds fqzcomp, level 4 the arithmetic coder.
  They reference `cram/ce.fa`. Apache 2.0 (`cram/3.1/LICENSE`).
