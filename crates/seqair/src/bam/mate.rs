//! A record's mate fields: [`MateFields`].

use crate::reader::Tid;
use seqair_types::Pos0;
use std::num::NonZeroI32;

// r[impl bam.record.mate_fields+2]
/// A record's mate fields — SAM's RNEXT, PNEXT and TLEN, bytes 20–32 of a BAM
/// record: where the template's next segment is, and how long the template
/// is. Each is `None` where the format says "unavailable" (`*`/`-1`, `0`/`-1`,
/// `0`).
///
/// They are the values the file gives the record: BAM's stored ones, or for
/// CRAM those resolved from the slice's mate chain. Whether the mate is in a
/// store is answered by the store
/// ([`RecordStore::link_mates`](super::RecordStore::link_mates)), and
/// [`RecordRef::mate`](super::RecordRef::mate) is the linked record.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MateFields {
    /// The next segment's reference.
    pub tid: Option<Tid>,
    /// The next segment's 0-based position.
    pub pos: Option<Pos0>,
    /// Signed template length; `None` when unavailable or for a
    /// single-segment template.
    pub template_len: Option<NonZeroI32>,
}
