# CRAM stream fixtures

`fidelity.sam` (19 records) and `fidelity.fa` exercise the places where
htslib builds a BAM field the CRAM does not store verbatim: attached and
detached mates, TLEN ties, mapped SEQ `*`, QUAL `*`, IUPAC reference and read
bases, reads without MD/NM, a read group missing from the header, placed and
unplaced unmapped reads. They come from the CRAM decode fidelity notes
(`plans/cram-decode-fidelity-2026-10-29.md`), which say what each record is for.

Made with samtools 1.25 (htslib 1.25):

```sh
samtools faidx fidelity.fa
samtools view --no-PG -C -T fidelity.fa -o fidelity_v31.cram fidelity.sam
samtools view --no-PG -C -T fidelity.fa --output-fmt-option version=3.0 -o fidelity_v30.cram fidelity.sam
# several containers, several slices per container
samtools view --no-PG -C -T fidelity.fa --output-fmt-option seqs_per_slice=3 \
    --output-fmt-option slices_per_container=2 -o fidelity_slices.cram fidelity.sam
# multi-reference slices (samtools 1.25 needs multi_seq_per_slice to make them)
samtools view --no-PG -C -T fidelity.fa --output-fmt-option multi_seq_per_slice=1 \
    --output-fmt-option seqs_per_slice=6 -o fidelity_multiref.cram fidelity.sam
# RN=false: names of attached pairs are not stored
samtools view --no-PG -C -T fidelity.fa --output-fmt-option lossy_names=1 -o fidelity_lossy_names.cram fidelity.sam
# consensus embedded reference
samtools view --no-PG -C --output-fmt-option embed_ref=2 -o fidelity_embed_ref2.cram fidelity.sam
# no reference: every base stored as a feature
samtools view --no-PG -C --output-fmt-option no_ref=1 -o fidelity_no_ref.cram fidelity.sam
# a slice holding only an unplaced read with SEQ and QUAL `*` (RL 0)
samtools view --no-PG -C -T fidelity.fa -o onlyu.cram onlyu.sam
for f in *.cram; do samtools index "$f"; done
```

## Shared blocks (htsjdk)

`ce_5b_java.cram` and `auxf_values_java.cram` (with `auxf.fa`) are htslib's
test files `test/ce#5b_java.cram` and `test/auxf#values_java.cram`, written by
htsjdk (CRAM 3.0) and copied from htslib 1.25 (MIT/Expat). Unlike samtools,
htsjdk puts several series in one block: `BA`, `IN` and `SC` share content id
0, `NP` and `TS` share 3, and tags of one type share a block. They read
against `../hts-specs/cram/ce.fa` and `auxf.fa`; they test that a decode that
skips one of the series in a shared block still reads the others right.
