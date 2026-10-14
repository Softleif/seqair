//! The reusable decoder state a CRAM reader keeps for its block codecs.

use super::{rans::Rans4x8Buf, rans_nx16::Nx16Order1Buf};

/// Everything a block decode reuses from one block to the next: the rANS
/// order-1 tables and the bzip2 decoder's working memory (several MiB each).
///
/// Every member allocates on first use, so a fresh scratch costs nothing and
/// one that never sees an order-1 or bzip2 block stays empty. The reader owns
/// one and passes it down to each codec as a required argument, so no decode
/// path falls back to allocating per call.
// r[impl cram.codec.scratch]
#[derive(Default)]
pub(crate) struct CodecScratch {
    pub(crate) rans_4x8: Rans4x8Buf,
    pub(crate) nx16_order1: Nx16Order1Buf,
    pub(crate) bz2: bz2ngrs::Decoder,
}
