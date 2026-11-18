# VCF/BCF fixtures

Small slices of real files, used by the reading guide's doctests
(`crates/seqair/src/vcf/reading.md`) and the reader's tests.

| file | source | records | samples |
|---|---|---|---|
| `rastair.bcf` (+ `.csi`), `rastair.vcf.gz` (+ `.tbi`) | rastair calls (SV168), first 200 records of chr1 and of chr2; same records in both | 400 | 1 |
| `truth.vcf.gz` (+ `.tbi`) | GIAB HG002 GRCh38 v4.2.1 benchmark, `chr1:1-1500000,chr2:1-400000` | 1237 | 1 |
| `panel.bcf` (+ `.csi`) | 1000 Genomes high-coverage phased panel, chr22:20000000-20300000, samples HG00096, HG00097, HG00099, NA12878, first 150 sites with an ALT call | 150 | 4 |

Regenerate with `bcftools view` from the sources (indexes with `--write-index`).
The panel's header still declares `AC` as `Number=1` (bcftools warns); kept as is,
since real files carry such declarations.
