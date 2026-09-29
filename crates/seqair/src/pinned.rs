//! Case counts for property tests that pin their own `test_cases`.
//!
//! Since hegeltest 0.47 a count compiled into `#[hegel::test(test_cases = N)]`
//! wins over `HEGEL_TEST_CASES`, so that variable no longer reaches the pinned
//! tests — the expensive subprocess oracles, which are exactly the ones the
//! nightly `deep-oracles` job exists to run deeper. They are written
//! `test_cases = pinned::cases(N)` instead, and `SEQAIR_PINNED_CASES`
//! overrides them all.
//!
//! Shared with the integration tests via `#[path = "../src/pinned.rs"]`.

/// `n`, or the value of `SEQAIR_PINNED_CASES` when it is set.
pub fn cases(n: u64) -> u64 {
    match std::env::var("SEQAIR_PINNED_CASES") {
        Ok(v) if !v.is_empty() => {
            v.parse().expect("SEQAIR_PINNED_CASES must be a positive integer")
        }
        _ => n,
    }
}
