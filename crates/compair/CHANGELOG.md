# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- **`PcrIndelModel`**, GATK's PCR indel gap model: the gap-open quality falls inside tandem
  repeats (units of 1 to `max_unit` bases) as `gap_open - exp(repeats / (rate * ln 2)) + 1`,
  rounded and kept between `floor` and `gap_open`. `Default` is GATK's conservative setting (Q45,
  floor Q10, rate 3, units up to 8, 20 repetitions). `gap_open_qual(repeats)` and
  `fill_gap_open(bases, track)` compute it, and **`Read::with_pcr_indel_model(bases, base_quals,
  &model, gap_qual, strand)`** builds a read from it, so a caller passes no per-base gap tracks.
