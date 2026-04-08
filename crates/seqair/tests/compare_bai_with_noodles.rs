//! Compares BAI index query results between seqair and noodles-bam.
//!
//! TODO: These tests require noodles-csi for indexed BAM reading.
//! Add noodles-csi as a dev-dependency and implement:
//! - Query the same region with both seqair and noodles indexed readers
//! - Compare the set of records returned (should be identical)
//! - Compare chunk offsets from BAI queries
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

// TODO: Add BAI index query comparison tests
// - bai_query_returns_same_records_as_noodles
// - bai_chunk_offsets_match_noodles
