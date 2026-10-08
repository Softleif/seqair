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

### Changed

- **A `Read` keeps its four quality tracks in one allocation**, and `Read::uniform` and
  `Read::with_pcr_indel_model` write the indel tracks straight into it: building a read allocates
  four times (bases, qualities, error probabilities, transitions) instead of seven, plus three
  temporary tracks for `uniform`. Scores are unchanged. When a read has more than one problem,
  `Read::new` now reports a wrong-length track before a missing quality in an earlier one.
