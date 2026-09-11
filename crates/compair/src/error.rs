/// Everything this crate rejects at construction time.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    #[error("read field `{field}` has {actual} entries but the read has {expected} bases")]
    ReadLengthMismatch { field: &'static str, expected: usize, actual: usize },

    #[error("read field `{field}` has no quality at index {index}")]
    MissingQuality { field: &'static str, index: usize },

    #[error("a read must report a known strand, got `Unknown`")]
    UnknownStrand,

    #[error("band width must be at least 2, got {0}")]
    BandTooNarrow(u32),

    #[error("band width must be at most {max}, got {got}")]
    BandTooWide { got: u32, max: u32 },
}
