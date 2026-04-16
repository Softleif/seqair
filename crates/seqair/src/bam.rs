//! BAM reading and pileup. Open a BAM file with [`IndexedBamReader`], fetch records into a
//! [`RecordStore`], then drive [`PileupEngine`] to iterate [`PileupColumn`]s.

// # Slab storage vs zero-copy references into decompressed blocks
//
// We considered skipping the copy in `RecordStore::push_raw()` and instead
// holding references directly into decompressed BGZF blocks. Our *hypothesis*
// (not yet validated with profiling) is that the slab design wins for pileup
// performance. The reasoning:
//
// 1. **Cache locality** — `Vec<Base>` is likely the hottest structure (accessed
//    every position for every active read). Contiguous and sequential, it should
//    prefetch well. In decompressed BAM blocks the 4-bit packed seq is
//    interleaved with headers, cigar, qual, and aux of neighboring records,
//    which would pull in irrelevant cache lines.
//
// 2. **Dense qual access** — qual bytes are packed tightly in the data slab. In
//    raw BAM, qual sits at a variable offset within each record, scattered across
//    64 KB decompressed blocks.
//
// 3. **Discarded overhead** — each BAM record carries a 36-byte fixed header
//    (refID, pos, bin, mq, flags, next_*, tlen) parsed once and never touched
//    again. Slabs store only what the pileup engine needs.
//
// 4. **Cold copy vs hot iteration** — `push_raw()` runs once per record during
//    region loading. The pileup hot loop (thousands of positions × dozens of
//    active reads) would benefit from dense packing on every iteration.
//
// 5. **Block-spanning records** — ~1% of records straddle BGZF block boundaries.
//    Zero-copy would need a fallback allocation for these, adding branches to
//    every accessor.
//
// 6. **Seq still needs decoding** — 4-bit BAM nibbles must be converted to `Base`
//    regardless, so the bases slab cannot be eliminated.
//
// TODO: validate with samply profiling on a real-world BAM to confirm that
// push_raw() is not a significant bottleneck and that slab access patterns
// are indeed cache-friendly.

// r[impl io.minimal_public_api]
pub mod aux;
pub mod aux_data;
pub mod bgzf;
pub mod bgzf_writer;
pub mod cigar;
pub mod csi_index;
pub mod header;
pub mod index;
pub mod owned_record;
pub mod pileup;
pub mod reader;
pub mod record;
pub mod record_store;
pub mod region_buf;
pub mod seq;
pub mod writer;

pub use aux_data::{AuxData, AuxDataError};
pub use bgzf::BgzfError;
pub use cigar::CigarOp;
pub use csi_index::{CsiError, CsiIndex};
pub use header::{BamHeader, BamHeaderError, ContigInfo};
pub use index::{AlignmentIndex, BaiError, BamIndex};
pub use owned_record::{OwnedBamRecord, OwnedRecordError};
pub use pileup::{PileupColumn, PileupEngine, PileupOp, RefSeq};
pub use reader::{BamError, BamShared, IndexedBamReader};
pub use record::BamRecord;
pub use record_store::RecordStore;
pub use seqair_types::bam_flags as flags;
pub use seqair_types::{Offset, One, Pos, Pos0, Pos1, Zero};
pub use writer::{BamWriteError, BamWriter};
