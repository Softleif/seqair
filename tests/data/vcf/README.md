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

`headers/` holds header text only, for `VcfHeader::parse` (`crates/seqair/tests/vcf_header_parse.rs`):

| file | source | trimmed |
|---|---|---|
| `kgp_chr22.bcf.header` | 1000 Genomes panel BCF `l_text` (with `IDX=`; FORMAT GT declared between INFO lines) | samples cut to 4 |
| `gnomad41_chr22.bcf.header` | gnomAD 4.1 chr22 BCF `l_text` (with `IDX=`) | 4 of every 5 INFO lines dropped, leaving IDX gaps |
| `rastair_sv168.bcf.header` | rastair BCF `l_text` (no `IDX=`, VCFv4.3, `Number=M`) | first 25 contigs |
| `giab_hg002.vcf.header` | GIAB HG002 v4.2.1 VCF text | first 25 contigs |
| `clinvar.vcf.header` | ClinVar VCF text | none |
