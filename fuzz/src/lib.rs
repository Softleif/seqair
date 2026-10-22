//! Shared types for fuzz targets and seed generation.

pub mod indexed_reader;

use seqair::bam::header::BamHeader;
use std::sync::LazyLock;

/// A header with 64 references, `c0`..`c63`, for targets that build records
/// without one: `push_raw` checks a record's reference ids against it, so a
/// seed's small ids decode and the fuzzer still reaches the rejection.
pub fn header() -> &'static BamHeader {
    static HEADER: LazyLock<BamHeader> = LazyLock::new(|| BamHeader::synthetic(64, 1_000_000_000));
    &HEADER
}
