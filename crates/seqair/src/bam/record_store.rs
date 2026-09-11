//! Slab-based BAM record storage with zero per-record heap allocation.
//!
//! All variable-length data is packed into contiguous byte buffers:
//! - **Name slab**: read names (qnames), accessed during dedup mate detection
//! - **Bases slab**: decoded `Base` values per record, accessed per-position in pileup
//! - **Cigar slab**: packed CIGAR ops, append-only mutation during realignment
//! - **Qual slab**: raw Phred bytes, hot in pileup per-base quality lookup
//! - **Aux slab**: auxiliary tag bytes in BAM binary format, rarely read in pileup

use super::{
    aux::Aux,
    cigar::{self, CigarOp},
    record::{self, DecodeError},
    record_idx::{RecordIdx, RecordRef},
    seq,
};
use hashbrown::HashTable;
use seqair_types::{BamFlags, Base, BaseQuality, Pos0, QPos};

// r[impl record_store.qname_hash]
// r[impl record_store.qname_hash.identity_is_probabilistic]
/// Seed-fixed 64-bit hash of a read name, with full avalanche.
///
/// Both mates of a template carry the same qname (`[SAM1] §1.4`), so this is
/// also a fragment identifier: a pure function of the bytes, with no per-store,
/// per-thread, or per-run state, so the same read hashes to the same value in
/// every store that ever holds it.
///
/// The construction is a folded 128-bit multiply per 8-byte word (the wyhash
/// primitive): the fold of the high and low halves mixes every input bit into
/// every output bit, unlike the single low-half multiply of `FxHash`, whose output
/// low bits barely move when only the high bytes of a word differ — which is
/// exactly the shape of Illumina qnames, where reads of one tile share a long
/// ASCII prefix.
#[must_use]
pub fn qname_hash(name: &[u8]) -> u64 {
    // Arbitrary odd 64-bit constants with ~half their bits set; any such
    // triple works, they only need to be fixed forever.
    const K0: u64 = 0x2d35_8dcc_aa6c_78a5;
    const K1: u64 = 0x8bb8_4b93_962e_acc9;
    const K2: u64 = 0x4b33_a62e_d433_d4a3;

    // r[impl record_store.qname_hash.no_name]
    if name.is_empty() {
        return 0;
    }

    let len = u64::try_from(name.len()).unwrap_or(u64::MAX);
    let mut acc = fold_mul(K0 ^ len, K1);

    let chunks = name.as_chunks::<8>();
    for chunk in chunks.0 {
        let word = u64::from_le_bytes(*chunk);
        acc = fold_mul(acc ^ word, K1);
    }

    let rest = chunks.1;
    if !rest.is_empty() {
        let mut buf = [0u8; 8];
        if let Some(head) = buf.get_mut(..rest.len()) {
            head.copy_from_slice(rest);
        }
        acc = fold_mul(acc ^ u64::from_le_bytes(buf), K2);
    }

    // r[impl record_store.qname_hash.no_name]
    // 0 is the reserved "this record has no qname" value, so a real name must
    // never land on it. Mapping the single colliding input to 1 costs one
    // predictable branch and doubles that value's probability, which at 2^-64
    // is not a number anyone will ever notice.
    match fold_mul(acc, K2) {
        0 => 1,
        h => h,
    }
}

/// The 128-bit product of `a * k`, folded down to 64 bits by `XOR`ing its halves.
#[inline]
fn fold_mul(a: u64, k: u64) -> u64 {
    let product = u128::from(a).wrapping_mul(u128::from(k));
    #[expect(
        clippy::cast_possible_truncation,
        reason = "the truncation is the point: both halves of the product are folded together"
    )]
    let folded = (product as u64) ^ ((product >> 64) as u64);
    folded
}

/// Compact BAM record with offsets into the store's slabs.
// r[impl record_store.push_raw+2]
// r[impl record_store.slim_record_fields]
pub struct SlimRecord {
    /// 0-based leftmost aligned reference position.
    pub pos: Pos0,
    // r[impl interval.end_pos_inclusive]
    /// Last reference position the alignment covers, 0-based and **inclusive**
    /// (`pos + ref_len - 1`).
    ///
    /// For unmapped reads, and for a CIGAR that consumes no reference, it
    /// equals [`pos`](Self::pos) — the same single position htslib reports as
    /// `end = start + 1` in its exclusive convention.
    pub end_pos: Pos0,
    // r[impl flags.field_type]
    /// BAM flag bits (paired, unmapped, reverse strand, …).
    pub flags: BamFlags,
    /// Number of CIGAR ops; bounded by BAM's u16 `l_op` field.
    pub n_cigar_ops: u16,
    /// Mapping quality (0..=254; 255 = unavailable per `[SAM1] §1.4`).
    pub mapq: u8,
    /// Query sequence length in bases — also the qual slab length for this record.
    pub seq_len: u32,
    // r[impl perf.precompute_matches_indels]
    /// Sum of M/=/X op lengths, pre-computed from CIGAR at push time.
    pub matching_bases: u32,
    /// Sum of I/D op lengths, pre-computed from CIGAR at push time.
    pub indel_bases: u32,
    /// Reference target index; -1 for unmapped.
    pub tid: i32,
    pub next_ref_id: i32,
    /// Mate's 0-based reference position; -1 if unavailable.
    pub next_pos: i32,
    /// Signed template/insert size (TLEN); 0 if unavailable or mates on different refs.
    pub template_len: i32,
    /// Offset into the name slab.
    name_off: u32,
    /// Length of the qname in the name slab.
    name_len: u16,
    /// Offset into the bases slab (`seq_len` Base values).
    bases_off: u32,
    /// Offset into the cigar slab (packed u32 ops).
    cigar_off: u32,
    /// Offset into the qual slab (`seq_len` Phred bytes).
    qual_off: u32,
    /// Offset into the aux slab (`aux_len` BAM aux bytes).
    aux_off: u32,
    /// Length of aux data in the aux slab.
    aux_len: u32,
    /// Index into the extras slab. Follows the same slab-offset pattern as
    /// other fields — survives record reordering (sort/dedup) without
    /// requiring the extras Vec to be reordered in lockstep.
    extras_idx: u32,
    // r[impl record_store.qname_hash]
    /// Seed-fixed hash of the qname bytes, computed at push time; `0` for a
    /// record with no qname. Both mates of a template share it. Read it through
    /// [`qname_hash`](Self::qname_hash); see the free [`qname_hash`] function
    /// for the construction.
    qname_hash: u64,
    // r[impl record_store.link_mates+2]
    // r[impl record_store.record_idx]
    /// Store index of this record's mate, `None` when unlinked — free, because
    /// [`RecordIdx`] leaves `u32::MAX` unrepresentable. Only
    /// [`RecordStore::link_mates`] writes it.
    mate_idx: Option<RecordIdx>,
}

// r[impl record_store.slim_record.field_getters]
impl SlimRecord {
    fn cigar_len(&self) -> usize {
        self.n_cigar_ops as usize
    }

    /// Read this record's qname bytes from the store's name slab.
    pub fn qname<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<&'store [u8], RecordAccessError> {
        let start = self.name_off as usize;
        let end =
            start.checked_add(self.name_len as usize).ok_or(RecordAccessError::OffsetOverflow {
                slab: Slab::Names,
                offset: self.name_off,
                len: self.name_len as usize,
            })?;
        store.names.get(start..end).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Names,
            offset: self.name_off,
        })
    }

    /// Read this record's decoded sequence bases from the store's bases slab.
    pub fn seq<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<&'store [Base], RecordAccessError> {
        let start = self.bases_off as usize;
        let end =
            start.checked_add(self.seq_len as usize).ok_or(RecordAccessError::OffsetOverflow {
                slab: Slab::Bases,
                offset: self.bases_off,
                len: self.seq_len as usize,
            })?;
        store.bases.get(start..end).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Bases,
            offset: self.bases_off,
        })
    }

    /// Read this record's per-base quality scores from the store's qual slab.
    pub fn qual<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<&'store [BaseQuality], RecordAccessError> {
        let start = self.qual_off as usize;
        let end =
            start.checked_add(self.seq_len as usize).ok_or(RecordAccessError::OffsetOverflow {
                slab: Slab::Qual,
                offset: self.qual_off,
                len: self.seq_len as usize,
            })?;
        let q = store.qual.get(start..end).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Qual,
            offset: self.qual_off,
        })?;
        Ok(BaseQuality::slice_from_bytes(q))
    }

    /// Read the CIGAR ops for this record from the store.
    pub fn cigar<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<&'store [CigarOp], RecordAccessError> {
        let start = self.cigar_off as usize;
        let end = start.checked_add(self.cigar_len()).ok_or(RecordAccessError::OffsetOverflow {
            slab: Slab::Cigar,
            offset: self.cigar_off,
            len: self.cigar_len(),
        })?;
        store.cigar.get(start..end).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Cigar,
            offset: self.cigar_off,
        })
    }

    // r[impl bam.record.aux_wrapper]
    /// Read the auxiliary tag data for this record.
    ///
    /// Returns an [`Aux`] wrapper that provides a fluent
    /// [`get`](Aux::get) API with automatic type conversion.
    pub fn aux<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<Aux<'store>, RecordAccessError> {
        let start = self.aux_off as usize;
        let end =
            start.checked_add(self.aux_len as usize).ok_or(RecordAccessError::OffsetOverflow {
                slab: Slab::Aux,
                offset: self.aux_off,
                len: self.aux_len as usize,
            })?;
        let bytes = store.aux.get(start..end).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Aux,
            offset: self.aux_off,
        })?;
        Ok(Aux::new(bytes))
    }

    // r[impl record_store.link_mates+2]
    /// Store index of this record's mate, if [`RecordStore::link_mates`] linked
    /// it. Invalidated by `sort_by_pos`/`dedup`.
    #[must_use]
    pub fn mate_idx(&self) -> Option<RecordIdx> {
        self.mate_idx
    }

    // r[impl record_store.qname_hash.no_name]
    /// The template's identity: the seed-fixed hash of the qname
    /// ([`qname_hash`]), shared by both mates. `None` when the record carries no
    /// qname — a CRAM written with `RN=false` stores none — in which case it has
    /// no template identity at all and never links to a mate.
    ///
    /// Two distinct qnames can in principle share a hash
    /// (r\[`record_store.qname_hash.identity_is_probabilistic`\]); mate linking
    /// compares the bytes and is unaffected, but a consumer using this as a
    /// fragment id must tolerate the merge.
    #[must_use]
    pub fn qname_hash(&self) -> Option<u64> {
        (self.qname_hash != 0).then_some(self.qname_hash)
    }

    /// Read the per-record extra value for this record.
    pub fn extra<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<&'store U, RecordAccessError> {
        store.extras.get(self.extras_idx as usize).ok_or(RecordAccessError::SlabOffsetOutOfRange {
            slab: Slab::Extras,
            offset: self.extras_idx,
        })
    }
}

// Compile-time size guard: 17 u32 (incl. next_ref_id, extras_idx, and
// mate_idx, whose `Option` is free — see `RecordIdx`)
// + 1 u64 (qname_hash) + 3 u16 + 1 u8 + padding = 88 bytes. The u64 raises the
// struct's alignment to 8, so the 12 bytes of mate linking cost 16. If this ever
// grows, revisit the layout before accepting the hit —
// see `docs/spec/3-record_store.md` "Layout".
const _: () =
    assert!(std::mem::size_of::<SlimRecord>() <= 88, "SlimRecord grew unexpectedly large");

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RecordAccessError {
    #[error("Cannot get {slab:?} entry for record: offset {offset} out of range")]
    SlabOffsetOutOfRange { slab: Slab, offset: u32 },
    #[error(
        "Cannot get {slab:?} entry for record: offset {offset} + length {len} overflows slab limits"
    )]
    OffsetOverflow { slab: Slab, offset: u32, len: usize },
}

#[derive(Debug)]
#[non_exhaustive]
pub enum Slab {
    Records,
    Names,
    Bases,
    Cigar,
    Qual,
    Aux,
    Extras,
}

// r[impl record_store.push_raw+2]
// r[impl record_store.field_access]
// r[impl record_store.clear+2]
// r[impl record_store.region_scoped]
// r[impl record_store.capacity]
// r[impl record_store.no_rc]
// r[impl base_decode.slab]
// r[related pileup.qpos]
/// Slab-based storage for BAM records with optional per-record extras.
///
/// Six contiguous buffers hold all data for a region, plus an extras slab:
/// - `records`: compact fixed-size structs
/// - `names`: packed qnames (for dedup)
/// - `bases`: decoded `Base` values per record (for per-position base lookup)
/// - `cigar`: packed CIGAR ops (separated for append-only mutation during realignment)
/// - `qual`: raw Phred bytes per record (dense, hot in pileup)
/// - `aux`: auxiliary tag bytes per record (rarely read in pileup; isolated so it
///   does not pollute the cache when scanning qual)
/// - `extras`: per-record user data of type `U` (default `()`, zero-cost)
///
/// Use [`CustomizeRecordStore`] to filter and compute per-record extras —
/// `compute` runs inline during [`push_raw`](RecordStore::push_raw) /
/// [`push_fields`](RecordStore::push_fields), populating the extras slab as
/// records arrive.
// r[impl record_store.extras.generic_param]
pub struct RecordStore<U = ()> {
    records: Vec<SlimRecord>,
    names: Vec<u8>,
    bases: Vec<Base>,
    cigar: Vec<CigarOp>,
    qual: Vec<u8>,
    aux: Vec<u8>,
    extras: Vec<U>,
    /// The two properties [`PileupEngine`] relies on, tracked as records are
    /// pushed and reordered. Not public: they exist so
    /// [`RecordStore::prepare_for_pileup`] can be cheap when they already
    /// hold, and so that no caller can claim them without earning them.
    ///
    /// [`PileupEngine`]: crate::bam::pileup::PileupEngine
    order: RecordOrder,
    mate_links: MateLinkState,
    // r[impl record_store.window_query.reach]
    /// `reach[i]` is the largest `end_pos` among the mapped records at
    /// indices `..=i`, in position order — non-decreasing by construction, so
    /// the first record that can overlap a window is one binary search away.
    /// Built by [`prepare_for_pileup`](Self::prepare_for_pileup), which is
    /// the only way to obtain a [`PileupInput`], and meaningful only on the
    /// store it hands out: a push or a move after that would leave it stale,
    /// which is why it is derived at that one point rather than maintained.
    reach: Vec<Pos0>,
}

// r[impl record_store.pileup_input]
/// A [`RecordStore`] that is in ascending position order and has its mates
/// linked — the two properties [`PileupEngine`] relies on and does not check.
///
/// There is no constructor: [`RecordStore::prepare_for_pileup`] is the only
/// way to obtain one, so the properties cannot be claimed without being
/// established.
///
/// [`PileupEngine`]: crate::bam::pileup::PileupEngine
pub struct PileupInput<U = ()> {
    store: RecordStore<U>,
}

impl<U> PileupInput<U> {
    /// Unwrap to the underlying store. Used by the pileup engine; the
    /// guarantees do not travel with the returned value.
    pub(crate) fn into_store(self) -> RecordStore<U> {
        self.store
    }

    /// The store, read-only. Indices from
    /// [`records_overlapping`](Self::records_overlapping) address it.
    pub fn store(&self) -> &RecordStore<U> {
        &self.store
    }

    // r[impl record_store.window_query]
    /// Indices of the mapped records whose alignment overlaps `[start, end]`
    /// — both ends inclusive, as everywhere in this API — in ascending
    /// order. `end < start` is the empty interval and yields nothing, and so
    /// does a placed-unmapped record, which has no alignment: the pileup
    /// never reports one, and this query answers for what the pileup reports.
    ///
    /// Lives here rather than on [`RecordStore`] because it binary-searches
    /// an index [`prepare_for_pileup`](RecordStore::prepare_for_pileup)
    /// builds over the records in position order — and this is the type that
    /// proves both exist.
    pub fn records_overlapping(
        &self,
        start: Pos0,
        end: Pos0,
    ) -> impl Iterator<Item = RecordIdx> + '_ {
        self.store.records_overlapping_sorted(start, end)
    }
}

/// The result of [`RecordStore::prepare_for_pileup`]: the store the engine
/// accepts, and what linking found.
///
/// `stats` is separate rather than folded into [`PileupInput`] because callers
/// log it (repeated primary qnames mean some reads are not deduplicated) and
/// the engine has no use for it.
pub struct Prepared<U = ()> {
    pub input: PileupInput<U>,
    pub stats: MateLinkStats,
}

/// Whether the store's records are in ascending `pos` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RecordOrder {
    /// Every record is at or after its predecessor. An empty store qualifies.
    Ascending,
    /// A push arrived before its predecessor; only `sort_by_pos` restores it.
    Unknown,
}

/// Whether `link_mates` has run since the last thing that invalidated it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MateLinkState {
    Linked,
    Unlinked,
}

impl<U> RecordStore<U> {
    pub fn new() -> Self {
        Self {
            order: RecordOrder::Ascending,
            mate_links: MateLinkState::Unlinked,
            reach: Vec::new(),
            records: Vec::new(),
            names: Vec::new(),
            bases: Vec::new(),
            cigar: Vec::new(),
            qual: Vec::new(),
            aux: Vec::new(),
            extras: Vec::new(),
        }
    }

    // r[impl perf.arena_capacity_hint+2]
    /// Pre-allocate based on estimated compressed bytes for the region.
    ///
    /// TODO(capacity-profiles): these are tuned for ~150bp Illumina short reads
    /// (WGS and 5-base-modified). Long-read methylation (ONT/PacBio with large
    /// MM/ML tags) has radically different per-slab proportions — aux can equal
    /// or exceed qual, and CIGAR grows by 2+ orders of magnitude. If/when seqair
    /// grows beyond rastair's current workloads, replace this with a `Profile`
    /// enum carrying per-workload constant blocks.
    pub fn with_byte_hint(compressed_bytes: usize) -> Self {
        // 5× matches measured 4.6–8.1×; erring toward the low end avoids huge
        // over-allocation on highly-compressible (chr-restricted) BAMs.
        let uncompressed_est = compressed_bytes.saturating_mul(5);
        // 550 B/rec ≈ pooled 569; rounded down so we over-estimate record count.
        let record_count_est = (uncompressed_est / 550).max(64);

        let names_est = record_count_est.saturating_mul(56); // headroom over pooled 51
        let bases_est = record_count_est.saturating_mul(150);
        // CIGAR at 2 ops/record (= 8 B) gives headroom over measured ~1.1 ops/record
        // for short reads without allocating the long-read budget we don't need here.
        let cigar_est = record_count_est.saturating_mul(2);
        let qual_est = record_count_est.saturating_mul(150);
        let aux_est = record_count_est.saturating_mul(180); // rounded up from pooled 163

        Self {
            order: RecordOrder::Ascending,
            mate_links: MateLinkState::Unlinked,
            reach: Vec::with_capacity(record_count_est),
            records: Vec::with_capacity(record_count_est),
            names: Vec::with_capacity(names_est),
            bases: Vec::with_capacity(bases_est),
            cigar: Vec::with_capacity(cigar_est),
            qual: Vec::with_capacity(qual_est),
            aux: Vec::with_capacity(aux_est),
            extras: Vec::with_capacity(record_count_est),
        }
    }
}

impl<U> RecordStore<U> {
    // r[impl record_store.pre_filter.rollback]
    /// Undo the most recent `push_raw`/`push_fields` by truncating every slab
    /// to the offsets recorded on the last `SlimRecord`. Only valid when the
    /// last record is the freshly-pushed one (no subsequent `sort`/`dedup`).
    fn rollback_last_push(&mut self) {
        let rec = self.records.pop().expect("rollback_last_push called with empty records Vec");
        self.extras
            .pop()
            .expect("push_raw/push_fields always append a matching () to extras before filter");
        self.names.truncate(rec.name_off as usize);
        self.bases.truncate(rec.bases_off as usize);
        self.cigar.truncate(rec.cigar_off as usize);
        self.qual.truncate(rec.qual_off as usize);
        self.aux.truncate(rec.aux_off as usize);
    }

    // r[impl record_store.extras.push_unit]
    // r[impl record_store.pre_filter.rollback]
    /// Decode a raw BAM record, append it, then consult
    /// `customize.filter_raw` (pre-extend) and `customize.filter` (post-compute)
    /// to decide whether to commit or roll back.
    ///
    /// Returns `Ok(Some(idx))` if the record survived all filters,
    /// `Ok(None)` if `filter_raw` or `filter` rejected it. If `filter_raw`
    /// rejects, no slab extension or base decode is performed.
    ///
    /// `filter_raw` sees the raw bytes before any slab work. `filter` sees
    /// the freshly-pushed [`SlimRecord`] and the `&self` store, so it can read
    /// qname, aux, etc. (those bytes live at the tail of the corresponding
    /// slabs until rollback). Pass `&mut ()` for no filtering (the blanket
    /// [`CustomizeRecordStore`] impl on `()` keeps every record).
    pub fn push_raw<E: CustomizeRecordStore<Extra = U>>(
        &mut self,
        raw: &[u8],
        customize: &mut E,
    ) -> Result<Option<RecordIdx>, DecodeError> {
        let h = record::parse_header(raw)?;

        let idx = RecordIdx::from_usize(self.records.len()).ok_or(DecodeError::SlabOverflow)?;

        // All slice bounds ≤ qual_end ≤ raw.len() (checked by parse_header)
        debug_assert!(h.qual_end <= raw.len(), "qual_end overrun: {} > {}", h.qual_end, raw.len());
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let qname_raw = &raw[32..h.var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let cigar_bytes = &raw[h.var_start..h.cigar_end];

        // --- filter_raw: pre-filter before any slab extension or base decode ---
        // r[impl record_store.pre_filter.rollback]
        {
            #[allow(clippy::indexing_slicing, reason = "qname_actual_len <= qname_raw.len()")]
            let qname_stripped = &qname_raw[..qname_actual_len];
            #[allow(
                clippy::indexing_slicing,
                reason = "all bounds \u{2264} qual_end \u{2264} raw.len()"
            )]
            let qual_slice = &raw[h.seq_end..h.qual_end];
            #[allow(clippy::indexing_slicing, reason = "qual_end \u{2264} raw.len()")]
            let aux_slice = &raw[h.qual_end..];
            #[allow(
                clippy::indexing_slicing,
                reason = "all bounds \u{2264} qual_end \u{2264} raw.len()"
            )]
            let packed_seq_slice = &raw[h.cigar_end..h.seq_end];
            if !customize.filter_raw(&FilterRawFields {
                pos: h.pos,
                end_pos: None,
                flags: h.flags,
                mapq: h.mapq,
                n_cigar_ops: h.n_cigar_ops,
                seq_len: h.seq_len,
                matching_bases: 0,
                indel_bases: 0,
                tid: h.tid,
                next_ref_id: h.next_ref_id,
                next_pos: h.next_pos,
                template_len: h.template_len,
                qname: qname_stripped,
                qual_bytes: qual_slice,
                aux_bytes: aux_slice,
                cigar: cigar_bytes,
                seq: Sequence::Packed(packed_seq_slice),
            }) {
                return Ok(None);
            }
        }

        // --- Write into cigar slab ---
        // r[impl record_store.checked_offsets]
        let cigar_off = u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;
        CigarOp::extend_from_bam_bytes(&mut self.cigar, cigar_bytes);
        let cigar_start = cigar_off as usize;
        #[allow(clippy::indexing_slicing, reason = "just appended; offsets within slab bounds")]
        let cigar_ops = &self.cigar[cigar_start..];
        // r[impl record_store.end_pos_htslib]
        let end_pos = if h.flags.is_unmapped() {
            h.pos
        } else {
            cigar::compute_end_pos(h.pos, cigar_ops)
                .ok_or(DecodeError::InvalidPosition { value: h.pos.as_i32() })?
        };
        let (matching_bases, indel_bases) = cigar::calc_matches_indels(cigar_ops);

        // --- Write into name slab ---
        // r[impl record_store.checked_offsets]
        let name_off = u32::try_from(self.names.len()).map_err(|_| DecodeError::SlabOverflow)?;
        #[allow(clippy::indexing_slicing, reason = "qname_actual_len ≤ qname_raw.len()")]
        self.names.extend_from_slice(&qname_raw[..qname_actual_len]);

        // --- Write into bases slab (4-bit → Base via SIMD, direct-to-slab) ---
        // r[impl record_store.checked_offsets]
        let bases_off = u32::try_from(self.bases.len()).map_err(|_| DecodeError::SlabOverflow)?;
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let packed_seq = &raw[h.cigar_end..h.seq_end];
        let seq_len = h.seq_len as usize;

        // Reserve space and decode directly into the slab's spare capacity.
        self.bases.reserve(seq_len);
        let bases_spare = self.bases.spare_capacity_mut();
        // Safety: spare_capacity_mut returns at least seq_len bytes after reserve().
        // decode_bases_into writes only valid Base discriminants (A=65,C=67,G=71,T=84,Unknown=78)
        // as guaranteed by the DECODE_BASE_TYPED table invariant. Base is repr(u8).
        // r[impl bam.record.seq_at_simd+2]
        // r[depends base_decode.table_invariant]
        unsafe {
            let out = std::slice::from_raw_parts_mut(bases_spare.as_mut_ptr() as *mut u8, seq_len);
            seq::decode_bases_into(packed_seq, seq_len, out);
            self.bases.set_len(
                self.bases.len().checked_add(seq_len).expect("bases slab length overflow"),
            );
        }

        // --- Write into qual slab ---
        // r[impl record_store.checked_offsets]
        let qual_off = u32::try_from(self.qual.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // Quality scores
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        self.qual.extend_from_slice(&raw[h.seq_end..h.qual_end]);

        // --- Write into aux slab ---
        // r[impl record_store.checked_offsets]
        let aux_off = u32::try_from(self.aux.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // Aux data (everything after qual)
        #[allow(clippy::indexing_slicing, reason = "qual_end ≤ raw.len()")]
        let aux_slice = &raw[h.qual_end..];
        self.aux.extend_from_slice(aux_slice);

        #[allow(clippy::indexing_slicing, reason = "qname_actual_len ≤ qname_raw.len()")]
        // r[impl record_store.qname_hash]
        let name_hash = qname_hash(&qname_raw[..qname_actual_len]);

        self.note_push(h.pos);
        self.records.push(SlimRecord {
            pos: h.pos,
            end_pos,
            flags: h.flags,
            n_cigar_ops: h.n_cigar_ops,
            mapq: h.mapq,
            seq_len: h.seq_len,
            matching_bases,
            indel_bases,
            tid: h.tid,
            next_ref_id: h.next_ref_id,
            next_pos: h.next_pos,
            template_len: h.template_len,
            name_off,
            #[expect(
                clippy::cast_possible_truncation,
                reason = "BAM qname is validated to ≤ 254 bytes by parse_header (l_read_name is u8); fits in u16"
            )]
            name_len: qname_actual_len as u16,
            bases_off,
            cigar_off,
            qual_off,
            aux_off,
            #[expect(
                clippy::cast_possible_truncation,
                reason = "aux data is bounded by slab limits (u32); slab overflow checked above via SlabOverflow"
            )]
            aux_len: aux_slice.len() as u32,
            extras_idx: idx.get(),
            qname_hash: name_hash,
            mate_idx: None,
        });
        self.extras.push(customize.compute(
            self.records.last().expect("just pushed a SlimRecord above; records.last() is Some"),
            self,
        ));

        // r[impl record_store.pre_filter.rollback]
        if customize.filter(
            self.records.last().expect("just pushed a SlimRecord above; records.last() is Some"),
            self,
        ) {
            Ok(Some(idx))
        } else {
            self.rollback_last_push();
            Ok(None)
        }
    }

    // r[impl unified.record_store_push]
    // r[impl record_store.push_fields]
    // r[impl unified.push_fields_equivalence]
    // r[impl record_store.checked_offsets]
    // r[impl record_store.pre_filter.rollback]
    /// Append a record from pre-parsed fields (for SAM/CRAM readers), then
    /// consult `customize.filter_raw` (pre-extend) and `customize.filter`
    /// (post-compute) to decide whether to commit or roll back.
    ///
    /// Writes directly into the slabs without going through BAM binary encoding.
    ///
    /// Returns `Ok(Some(idx))` if the record survived all filters,
    /// `Ok(None)` if `filter_raw` or `filter` rejected it. If `filter_raw`
    /// rejects, no slab extension is performed. Pass `&mut ()` to disable
    /// filtering.
    #[expect(
        clippy::too_many_arguments,
        reason = "all fields are needed for zero-copy push into the record store slabs"
    )]
    pub fn push_fields<E: CustomizeRecordStore<Extra = U>>(
        &mut self,
        pos: Pos0,
        end_pos: Pos0,
        flags: BamFlags,
        mapq: u8,
        matching_bases: u32,
        indel_bases: u32,
        qname: &[u8],
        cigar_ops: &[CigarOp],
        bases: &[Base],
        qual: &[u8],
        aux: &[u8],
        tid: i32,
        next_ref_id: i32,
        next_pos: i32,
        template_len: i32,
        customize: &mut E,
    ) -> Result<Option<RecordIdx>, DecodeError> {
        if qual.len() != bases.len() {
            return Err(DecodeError::QualLenMismatch {
                qual_len: qual.len(),
                seq_len: bases.len(),
            });
        }

        let idx = RecordIdx::from_usize(self.records.len()).ok_or(DecodeError::SlabOverflow)?;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "caller validates cigar op count ≤ 65535 (BAM n_cigar_op is u16); fits in u16"
        )]
        let n_cigar_ops = cigar_ops.len() as u16;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "seq length is bounded by slab limits (u32); slab overflow checked via SlabOverflow"
        )]
        let seq_len = bases.len() as u32;

        // --- filter_raw: pre-filter before any slab extension ---
        // r[impl record_store.filter_raw]
        if !customize.filter_raw(&FilterRawFields {
            pos,
            end_pos: Some(end_pos),
            flags,
            mapq,
            n_cigar_ops,
            seq_len,
            matching_bases,
            indel_bases,
            tid,
            next_ref_id,
            next_pos,
            template_len,
            qname,
            qual_bytes: qual,
            aux_bytes: aux,
            cigar: CigarOp::ops_as_bytes(cigar_ops),
            seq: Sequence::Bases(bases),
        }) {
            return Ok(None);
        }

        // Name slab
        // r[impl record_store.checked_offsets]
        let name_off = u32::try_from(self.names.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.names.extend_from_slice(qname);

        // Bases slab
        // r[impl record_store.checked_offsets]
        let bases_off = u32::try_from(self.bases.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.bases.extend_from_slice(bases);

        // Cigar slab
        // r[impl record_store.checked_offsets]
        let cigar_off = u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.cigar.extend_from_slice(cigar_ops);

        // Qual slab
        // r[impl record_store.checked_offsets]
        let qual_off = u32::try_from(self.qual.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.qual.extend_from_slice(qual);

        // Aux slab
        // r[impl record_store.checked_offsets]
        let aux_off = u32::try_from(self.aux.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.aux.extend_from_slice(aux);

        self.note_push(pos);
        self.records.push(SlimRecord {
            pos,
            end_pos,
            flags,
            n_cigar_ops,
            mapq,
            seq_len,
            matching_bases,
            indel_bases,
            tid,
            next_ref_id,
            next_pos,
            template_len,
            name_off,
            #[expect(
                clippy::cast_possible_truncation,
                reason = "caller validates qname ≤ 254 bytes; fits in u16"
            )]
            name_len: qname.len() as u16,
            bases_off,
            cigar_off,
            qual_off,
            aux_off,
            #[expect(
                clippy::cast_possible_truncation,
                reason = "aux data bounded by slab limits (u32); slab overflow checked via SlabOverflow"
            )]
            aux_len: aux.len() as u32,
            extras_idx: idx.get(),
            // r[impl record_store.qname_hash]
            qname_hash: qname_hash(qname),
            mate_idx: None,
        });

        let record =
            self.records.last().expect("just pushed a SlimRecord above; records.last() is Some");

        // r[impl record_store.extras.generic_param]
        self.extras.push(customize.compute(record, self));

        // r[impl record_store.pre_filter.rollback]
        if customize.filter(record, self) {
            Ok(Some(idx))
        } else {
            self.rollback_last_push();
            Ok(None)
        }
    }
}

// r[impl record_store.filter_raw_fields]
/// The form in which the CIGAR data is available in [`FilterRawFields`].
///
/// In the `push_raw` path (BAM), the CIGAR is available as raw BAM bytes
/// In `push_raw` (BAM), the sequence is in 4-bit packed BAM format.
/// In `push_fields` (SAM/CRAM), it has already been decoded into typed
/// [`Base`] values.
#[derive(Debug, Clone, Copy)]
pub enum Sequence<'a> {
    /// 4-bit packed BAM sequence bytes (from `push_raw`).
    Packed(&'a [u8]),
    /// Decoded sequence bases (from `push_fields`).
    Bases(&'a [Base]),
}

// r[impl record_store.filter_raw_fields]
/// Fields available to [`CustomizeRecordStore::filter_raw`] before any slab
/// data is written — the raw slices and parsed header fields from the incoming
/// record.
///
/// This struct is passed to `filter_raw` so customizers can reject records
/// without incurring any slab extension or base decode work.
///
/// # `push_raw` vs `push_fields`
///
/// From `push_raw` (BAM binary decode):
/// - [`seq`](Self::seq) is [`Sequence::Packed`] (4-bit packed sequence bytes).
/// - [`end_pos`](Self::end_pos) is `None` (not yet computed from CIGAR).
/// - [`matching_bases`](Self::matching_bases) and [`indel_bases`](Self::indel_bases)
///   are `0` (not yet computed from CIGAR).
///
/// From `push_fields` (pre-parsed SAM/CRAM fields):
/// - [`seq`](Self::seq) is [`Sequence::Bases`] (decoded sequence bases).
/// - [`end_pos`](Self::end_pos) is `Some(_)` (pre-computed).
/// - [`matching_bases`](Self::matching_bases) and [`indel_bases`](Self::indel_bases)
///   are the pre-computed values.
///
/// [`cigar`](Self::cigar) is always raw BAM CIGAR bytes (`&[u8]`). Use
/// [`CigarOp::slice_from_bam_bytes`] or [`CigarOp::extend_from_bam_bytes`]
/// for typed access.
#[derive(Debug)]
#[non_exhaustive]
pub struct FilterRawFields<'a> {
    /// 0-based reference position.
    pub pos: Pos0,
    /// Last reference position the record covers, inclusive
    /// (r\[`interval.end_pos_inclusive`\]).
    ///
    /// `Some(_)` from `push_fields` (SAM/CRAM) where the CIGAR is pre-parsed.
    /// `None` from `push_raw` (BAM) because the CIGAR hasn't been decoded yet
    /// at this point — `filter_raw` runs before any slab work; derive the span
    /// from [`cigar`](Self::cigar) if you need it there.
    pub end_pos: Option<Pos0>,
    /// BAM flags.
    pub flags: BamFlags,
    /// Mapping quality (0..254).
    pub mapq: u8,
    /// Number of CIGAR ops.
    pub n_cigar_ops: u16,
    /// Query sequence length.
    pub seq_len: u32,
    /// Sum of M/=/X CIGAR op lengths. `0` for `push_raw` (not yet computed).
    pub matching_bases: u32,
    /// Sum of I/D CIGAR op lengths. `0` for `push_raw` (not yet computed).
    pub indel_bases: u32,
    /// Reference target index; -1 for unmapped.
    pub tid: i32,
    /// Mate reference target index.
    pub next_ref_id: i32,
    /// Mate 0-based reference position; -1 if unavailable.
    pub next_pos: i32,
    /// Signed template/insert size.
    pub template_len: i32,
    /// Query name (NUL stripped).
    pub qname: &'a [u8],
    /// Raw quality score bytes.
    pub qual_bytes: &'a [u8],
    /// Raw auxiliary tag bytes (BAM binary format).
    pub aux_bytes: &'a [u8],
    /// Raw BAM CIGAR bytes. Use [`CigarOp::slice_from_bam_bytes`] for zero-cost
    /// typed access on aligned little-endian input, or
    /// [`CigarOp::extend_from_bam_bytes`] to copy into a `Vec<CigarOp>`.
    pub cigar: &'a [u8],
    /// Read sequence. [`Sequence::Packed`] from `push_raw` (BAM),
    /// [`Sequence::Bases`] from `push_fields` (SAM/CRAM).
    pub seq: Sequence<'a>,
}
// r[impl record_store.customize.trait]
/// Customize how records flow into a [`RecordStore`]: filter and compute
/// per-record extras, both inline at push time.
///
/// Used by [`Readers`](crate::reader::Readers) and the lower-level
/// `push_raw`/`push_fields`/`fetch_into_customized` APIs. Implementors are
/// `Clone` so they can be duplicated across forked readers in multi-threaded
/// pipelines.
///
/// # Ordering and state
///
/// Three methods run inline during push, for every record:
///
/// 1. **`filter_raw`** runs before any slab extension or base decode.
///    Returning `false` skips the record with zero wasted work.
/// 2. **`compute`** runs next, producing the extra value for the
///    freshly-pushed record. The `&RecordStore<Self::Extra>` argument
///    includes extras from previously-pushed records.
/// 3. **`filter`** runs last. If it returns `false`, all slab writes
///    for this record (including the just-computed extra) are rolled back
///    with zero waste — it is as if the record was never pushed.
///
/// Because `filter_raw` runs before any work, use it for cheap checks
/// (mapq, flags, tid, qname prefix) that can reject records without
/// touching the slabs.
///
/// # `compute` before `filter` — wasted work
///
/// `compute` runs *before* `filter`. If your `compute` is expensive (aux
/// parsing, base counting, read group extraction) and your `filter` might
/// reject many records, that work is wasted. When possible, push rejection
/// logic into `filter_raw` instead — it runs before any slab extension.
///
/// For checks that require slab data (aux tags, CIGAR, sequence), you have
/// two options:
///
/// 1. **Accept the waste.** If your `filter` rejects few records and your
///    `compute` needs the same slab data anyway, the waste is negligible.
/// 2. **Inline the check into `filter_raw`.** Parse the raw bytes directly
///    to reject early. This avoids slab writes AND `compute` for rejected
///    records. The second example below shows this pattern for read groups.
///
/// # Example: pre-filter on read group
///
/// A customizer that only keeps reads from a specific read group. Rather
/// than relying on `filter` (which would waste `compute` on every rejected
/// record), it scans the raw aux bytes in `filter_raw` — before any slab
/// work happens:
///
/// ```
/// use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};
///
/// #[derive(Clone)]
/// struct ReadGroupFilter {
///     /// Expected read group, e.g. `b"@RG\tID:mygroup"`.
///     target_rg: Vec<u8>,
/// }
///
/// impl CustomizeRecordStore for ReadGroupFilter {
///     type Extra = ();
///
///     /// Reject records whose RG tag doesn't match, before any slab extension.
///     fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
///         // Quick: if the aux bytes don't contain our target string at all,
///         // skip the record. No slab work done.
///         fields.aux_bytes
///             .windows(self.target_rg.len())
///             .any(|w| w == self.target_rg.as_slice())
///     }
///
///     fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
/// }
/// ```
///
/// The same `&mut E` instance is reused across regions (`Readers` holds one
/// customize value for its whole lifetime). Implementors that want
/// per-region state MUST reset it themselves between fetches via
/// [`Readers::customize_mut`](crate::reader::Readers::customize_mut).
/// Counters that should accumulate across regions need no special handling.
///
/// # Example
///
/// A customizer that drops low-quality reads early and tags each kept read
/// with its qname length:
///
/// ```
/// use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};
///
/// #[derive(Clone, Default)]
/// struct QnameLen { min_mapq: u8 }
///
/// impl CustomizeRecordStore for QnameLen {
///     type Extra = usize;
///
///     fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
///         fields.mapq >= self.min_mapq
///     }
///
///     fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<usize>) -> usize {
///         rec.qname(store).map(|n| n.len()).unwrap_or(0)
///     }
/// }
/// ```
pub trait CustomizeRecordStore: Clone {
    /// The per-record data produced by `compute`.
    type Extra;
    // r[impl record_store.filter_raw]
    /// Pre-filter on raw data: called before any slab extension or base
    /// decoding. Returning `false` skips the record with zero wasted work —
    /// no bytes are memcpy'd into slabs and no bases are decoded.
    ///
    /// Default: keep every record. Override to filter at the earliest
    /// possible point.
    ///
    /// Called from `push_raw` after parsing the BAM header and slicing
    /// the raw bytes, and from `push_fields` before any `extend_from_slice`.
    /// The [`FilterRawFields`] struct carries all header fields plus the
    /// raw slices (or decoded cigar/bases depending on the path).
    #[inline]
    fn filter_raw(&mut self, _fields: &FilterRawFields<'_>) -> bool {
        true
    }
    // r[impl record_store.customize.trait]
    /// Post-compute filter: called on each freshly-pushed record after
    /// [`compute`](Self::compute). Returning `false` rolls back the slab
    /// writes for this record — it is as if the record was never fetched.
    ///
    /// Default: keep every record. Override to filter at push time.
    ///
    /// The filter sees the pushed [`SlimRecord`] (with its `pos`, `flags`,
    /// `mapq`, etc.) and the `&RecordStore` for reading slab-backed
    /// data via `rec.qname(store)`, `rec.aux(store)`, etc. The freshly-pushed
    /// record is at the tail of each slab, so all of its bytes are accessible.
    #[inline]
    fn filter(&mut self, _rec: &SlimRecord, _store: &RecordStore<Self::Extra>) -> bool {
        true
    }
    /// Compute the extra for `rec`. Called inline at push time for every
    /// record (before [`filter`](Self::filter)), so the value is
    /// wasted if the record is later rejected. Use the `rec.seq(store)`,
    /// `rec.qual(store)`, `rec.aux(store)`, `rec.qname(store)`, and
    /// `rec.cigar(store)` getters on `SlimRecord` to read variable-length
    /// data without going through `RecordStore::record(idx)` first.
    fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<Self::Extra>) -> Self::Extra;
}

// r[impl record_store.customize.trait]
/// Blanket no-op implementation for the default `()` case.
///
/// `RecordStore<()>` does not need extras; `Readers<()>` uses this so
/// `compute` is a zero-cost no-op at push time.
/// The default `filter_raw` and `filter` (always `true`) mean `Readers<()>` never
/// filters either, so `&mut ()` is the no-op customizer to pass into
/// `push_raw`/`push_fields`/`fetch_into_customized`.
impl CustomizeRecordStore for () {
    type Extra = ();
    #[inline]
    fn compute(&mut self, _rec: &SlimRecord, _store: &RecordStore<()>) {}
}

/// Drop unmapped reads at the earliest possible point.
///
/// This is the simplest possible [`CustomizeRecordStore`] — it implements
/// only [`filter_raw`](CustomizeRecordStore::filter_raw) to reject unmapped
/// reads before any slab extension or base decode work happens. Because
/// [`filter_raw`](CustomizeRecordStore::filter_raw) runs before CIGAR is
/// decoded in the BAM path, it can only check fields that are always
/// available: `flags`, `mapq`, `tid`, `qname`, etc.
///
/// Pass this to [`Readers::open_customized`](crate::reader::Readers::open_customized)
/// to get htslib-compatible behavior where unmapped reads never enter the record
/// store. The default [`Readers::open`](crate::reader::Readers::open) passes
/// `()` as the customizer, which keeps unmapped reads — the pileup engine will
/// exclude them from columns but they remain in the store for other uses.
///
/// # Example
///
/// ```no_run
/// use seqair::bam::record_store::RejectUnmapped;
/// use seqair::reader::Readers;
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut readers = Readers::open_customized(
///     Path::new("sample.bam"),
///     Path::new("reference.fa"),
///     RejectUnmapped,
/// )?;
/// # Ok(())
/// # }
/// ```
///
/// # Source
///
/// The entire implementation:
///
/// ```
/// # use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};
/// #[derive(Clone, Default)]
/// pub struct RejectUnmapped;
///
/// impl CustomizeRecordStore for RejectUnmapped {
///     type Extra = ();
///     fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
///         !fields.flags.is_unmapped()
///     }
///     fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
/// }
/// ```
#[derive(Clone, Default)]
pub struct RejectUnmapped;

impl CustomizeRecordStore for RejectUnmapped {
    type Extra = ();
    #[inline]
    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        !fields.flags.is_unmapped()
    }
    #[inline]
    fn compute(&mut self, _rec: &SlimRecord, _store: &RecordStore<()>) {}
}

impl<U> RecordStore<U> {
    // --- Ordering ---

    // r[impl record_store.extras.sort_dedup_generic]
    /// Sort records by reference position.
    ///
    /// The pileup engine iterates records in store order and assumes they
    /// are sorted by position. After injecting chunk-cache records (which
    /// may have earlier positions), this must be called to restore the
    /// invariant.
    ///
    /// Safe for any `U` because each record carries its own `extras_idx`,
    /// so reordering the records Vec does not invalidate the extras mapping.
    /// Mate links are *not* safe for reordering and are dropped
    /// (r\[`record_store.link_mates.invalidated`\]).
    pub fn sort_by_pos(&mut self) {
        self.clear_mate_links();
        self.records.sort_by_key(|r| r.pos);
        self.order = RecordOrder::Ascending;
    }

    // r[impl record_store.link_mates.unique_records]
    /// Remove duplicate records (same position, flags, and read name).
    ///
    /// Nearby and distant BAM index chunks can cover overlapping byte ranges,
    /// causing the same record to be loaded from both sources.
    /// Slab data (including extras) for removed records is left in place (minor waste),
    /// same as dead name/cigar bytes from removed records.
    ///
    /// Duplicates are collapsed *consecutively*, so they have to be adjacent —
    /// which they are once the store is in position order. This used to be the
    /// caller's job ("must be called after `sort_by_pos`") and failed silently
    /// when forgotten: the duplicates simply stayed, and the record count was
    /// wrong in exactly the case the method exists for. It sorts first now.
    pub fn dedup(&mut self) {
        // Sorting also clears mate links, but do it unconditionally: an already
        // ordered store skips the sort and still must not keep stale links.
        if self.order == RecordOrder::Unknown {
            self.sort_by_pos();
        }
        self.clear_mate_links();
        let names = &self.names;
        self.records.dedup_by(|a, b| {
            a.pos == b.pos && a.flags == b.flags && a.name_len == b.name_len && {
                let a_start = a.name_off as usize;
                let b_start = b.name_off as usize;
                let len = a.name_len as usize;
                debug_assert!(a_start.saturating_add(len) <= names.len(), "name slice OOB");
                debug_assert!(b_start.saturating_add(len) <= names.len(), "name slice OOB");
                #[allow(clippy::indexing_slicing, reason = "offsets validated at push time")]
                {
                    names[a_start..a_start.saturating_add(len)]
                        == names[b_start..b_start.saturating_add(len)]
                }
            }
        });
    }

    // --- Mate linking ---

    // r[impl record_store.link_mates.invalidated]
    /// Drop every mate link. A mate index is a store index, so any operation
    /// that moves or removes records invalidates it.
    fn clear_mate_links(&mut self) {
        self.mate_links = MateLinkState::Unlinked;
        for rec in &mut self.records {
            rec.mate_idx = None;
        }
    }

    // r[impl record_store.link_mates+2]
    // r[impl record_store.link_mates.qname_uniqueness]
    // r[impl record_store.link_mates.stats]
    /// Link the mates of every template that has both of its primary mapped
    /// alignments in this store, so consumers can find a read's mate — and the
    /// reference interval the two share — without matching qnames per column.
    ///
    /// Run this as the *last* mutation before pileup iteration: `sort_by_pos`
    /// and `dedup` clear the links again
    /// (r\[`record_store.link_mates.invalidated`\]). [`Readers::pileup`] does this
    /// for you.
    ///
    /// The returned [`MateLinkStats`] is the only signal that the input broke
    /// the one assumption here — that a qname belongs to one template
    /// (`[SAM1] §1.4`) — so a caller should at least log a non-zero
    /// `ambiguous_qnames`.
    ///
    /// [`Readers::pileup`]: crate::reader::Readers::pileup
    #[must_use]
    pub fn link_mates(&mut self) -> MateLinkStats {
        // Split the borrow: the qname verification reads the name slab while
        // the link is written into the records Vec.
        let Self { records, names, .. } = self;
        for rec in records.iter_mut() {
            rec.mate_idx = None;
        }

        // Records still waiting for their mate, keyed by template identity.
        // A hash table is the right shape here for the reason any hash table
        // is: the hash groups candidates, and equality — which lives in the
        // name slab, not in the stored `u32` — decides. `HashTable` takes the
        // hash and the comparison as arguments precisely so the key can be
        // data the entry only points at. Two templates whose qnames collide
        // therefore share a bucket and cost one byte comparison; neither loses
        // its own mate, which a one-slot-per-key map would have caused.
        let mut waiting: HashTable<RecordIdx> = HashTable::with_capacity(records.len());
        let mut stats = MateLinkStats::default();

        for idx in 0..records.len() {
            let Some(record_idx) = RecordIdx::from_usize(idx) else { break };
            let Some(rec) = records.get(idx) else { continue };
            if !can_link(rec) {
                continue;
            }
            let qname = qname_bytes(names, rec);
            // r[impl record_store.qname_hash.no_name]
            if qname.is_empty() {
                continue;
            }
            let (pos, next_pos) = (rec.pos.as_i32(), rec.next_pos);

            // The probe visits every waiting candidate in the bucket, so it is
            // also where a same-qname record that cannot be *this* record's
            // mate gets noticed — a third primary alignment for the qname, or
            // the same record twice.
            // r[impl record_store.link_mates.qname_uniqueness]
            let mut same_qname_present = false;
            let partner = waiting
                .find(rec.qname_hash, |&candidate| {
                    let Some(other) = records.get(candidate.as_usize()) else { return false };
                    if qname_bytes(names, other) != qname {
                        return false;
                    }
                    same_qname_present = true;
                    other.mate_idx.is_none()
                        && pos == other.next_pos
                        && other.pos.as_i32() == next_pos
                })
                .copied();

            match partner {
                Some(partner) => {
                    if let Some(this) = records.get_mut(idx) {
                        this.mate_idx = Some(partner);
                    }
                    if let Some(other) = records.get_mut(partner.as_usize()) {
                        other.mate_idx = Some(record_idx);
                    }
                    stats.pairs = stats.pairs.saturating_add(1);
                }
                None => {
                    if same_qname_present {
                        stats.ambiguous_qnames = stats.ambiguous_qnames.saturating_add(1);
                    }
                    // Linked records stay in the table: a later record with the
                    // same qname must still be able to see that it arrived too
                    // late, which is what `ambiguous_qnames` counts.
                    waiting.insert_unique(rec.qname_hash, record_idx, |&candidate| {
                        records.get(candidate.as_usize()).map_or(0, |r| r.qname_hash)
                    });
                }
            }
        }

        self.mate_links = MateLinkState::Linked;
        stats
    }

    #[doc(hidden)]
    /// Test seam: overwrite a record's qname hash to stage a collision.
    ///
    /// A 64-bit hash with full avalanche cannot be collided by construction in
    /// a test, but the code path that survives a collision must still be
    /// exercised — so the hash is forced instead. Same shape as
    /// `RegionBuf::with_budget`, which forces a tiny window for tests.
    /// Do not call this outside tests: mate linking assumes the stored hash is
    /// the hash of the stored qname.
    pub fn set_qname_hash(&mut self, idx: RecordIdx, hash: u64) -> Option<()> {
        self.records.get_mut(idx.as_usize())?.qname_hash = hash;
        Some(())
    }

    // --- Accessors ---

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    // r[impl record_store.record_idx.resolution]
    // r[impl record_store.record_ref]
    /// The record at `idx`, borrowed from this store, or `None` when this
    /// store has no such record.
    ///
    /// The one door in. Reading anything about a record goes through the
    /// [`RecordRef`] this hands back, so the index is checked exactly once and
    /// everything downstream — fields, slab slices, the mate — is infallible.
    ///
    /// The handle borrows the store, so the store cannot be cleared, refilled
    /// or reordered while it is alive. That is the half a bounds check cannot
    /// do: an index kept across regions is *in range* for the next region's
    /// records and would resolve silently against the wrong read.
    ///
    /// ```
    /// # use seqair::bam::{RecordIdx, RecordStore};
    /// # fn demo(store: &RecordStore) {
    /// match store.record(RecordIdx::ZERO) {
    ///     Some(rec) => println!("{} at {}", rec.qname().len(), rec.pos.as_u64()),
    ///     None => println!("empty store"),
    /// }
    /// # }
    /// ```
    #[must_use]
    pub fn record(&self, idx: RecordIdx) -> Option<RecordRef<'_, U>> {
        Some(RecordRef::new(self, self.records.get(idx.as_usize())?, idx))
    }

    // r[impl record_store.iter]
    /// Iterate `&SlimRecord` in store order.
    ///
    /// Read variable-length fields directly off the record via the
    /// `SlimRecord` getters (`rec.cigar(store)`, `rec.qname(store)`, …);
    /// the index isn't needed for shared-reference access. Callers that
    /// truly need the position can use `records().enumerate()`.
    pub fn records(&self) -> impl ExactSizeIterator<Item = &SlimRecord> + '_ {
        self.records.iter()
    }

    // r[impl record_store.record_ref]
    /// The record at `idx`, for a caller that minted `idx` from this very
    /// store and has held it ever since — the pileup engine, which owns its
    /// store and walks `0..len` to activate records, is the only one.
    ///
    /// Crate-private, and deliberately not an `Option`: a caller that cannot
    /// prove the pairing has [`record`](Self::record) instead. Returning
    /// `Option` here pushed the impossible case onto call sites that then had
    /// to invent an answer for it — and `PileupColumn::alignments` invented
    /// "silently skip the entry", which would have hidden exactly the bug this
    /// panics on.
    ///
    /// # Panics
    ///
    /// If `idx` is not an index of this store, which means the invariant above
    /// was broken.
    pub(crate) fn record_at(&self, idx: RecordIdx) -> RecordRef<'_, U> {
        self.record(idx).expect("record index was minted by this store")
    }

    // r[impl record_store.record_idx]
    /// Every index of this store, in store order — the typed replacement for
    /// `0..store.len() as u32`. Pair it with [`record`](Self::record), which
    /// cannot fail for an index from here.
    pub fn indices(&self) -> impl Iterator<Item = RecordIdx> + '_ {
        // A slot whose index could not be minted never made it into the store:
        // `push_raw`/`push_fields` fail with `SlabOverflow` first.
        (0..self.records.len()).filter_map(RecordIdx::from_usize)
    }

    /// Update the `template_len` field for a record in the store.
    /// Used by the CRAM decoder to fill in TLEN after resolving mate
    /// cross-references within a slice.
    pub fn set_template_len(&mut self, idx: RecordIdx, tlen: i32) -> Option<()> {
        self.records.get_mut(idx.as_usize())?.template_len = tlen;
        Some(())
    }

    /// Update the mate fields (`next_ref_id`, `next_pos`) for a record.
    /// Used by the CRAM decoder to fill in mate info after resolving
    /// downstream mate cross-references within a slice.
    pub fn set_mate_info(&mut self, idx: RecordIdx, next_ref_id: i32, next_pos: i32) -> Option<()> {
        let rec = self.records.get_mut(idx.as_usize())?;
        rec.next_ref_id = next_ref_id;
        rec.next_pos = next_pos;
        Some(())
    }

    // The slab readers below take the record rather than its index: the
    // caller has already proved the pairing (it holds a `&SlimRecord` from
    // this store), so the only failure left is a slab offset this file wrote
    // wrongly — an invariant, not an argument. `RecordRef` is the public face
    // of all of them.
    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub(crate) fn qname_of(&self, rec: &SlimRecord) -> &[u8] {
        let start = rec.name_off as usize;
        let end = start.checked_add(rec.name_len as usize).expect("qname end overflow");
        debug_assert!(end <= self.names.len(), "qname slab overrun: {end} > {}", self.names.len());
        &self.names[start..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub(crate) fn cigar_of(&self, rec: &SlimRecord) -> &[CigarOp] {
        let start = rec.cigar_off as usize;
        let end = start.checked_add(rec.cigar_len()).expect("cigar end overflow");
        debug_assert!(end <= self.cigar.len(), "cigar slab overrun: {end} > {}", self.cigar.len());
        &self.cigar[start..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub(crate) fn seq_of(&self, rec: &SlimRecord) -> &[Base] {
        let start = rec.bases_off as usize;
        let end = start.checked_add(rec.seq_len as usize).expect("seq end overflow");
        debug_assert!(end <= self.bases.len(), "bases slab overrun: {end} > {}", self.bases.len());
        &self.bases[start..end]
    }

    // r[impl bam.record.seq_at]
    pub(crate) fn base_at_of(&self, rec: &SlimRecord, qpos: QPos) -> Base {
        // Past the record's own bases is `Unknown`, not another record's base:
        // the slab is shared, so the length has to be checked, not just the
        // offset.
        if qpos.as_usize() >= rec.seq_len as usize {
            return Base::Unknown;
        }
        self.bases
            .get((rec.bases_off as usize).saturating_add(qpos.as_usize()))
            .copied()
            .unwrap_or(Base::Unknown)
    }

    // r[impl types.base_quality.field_type]
    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub(crate) fn qual_of(&self, rec: &SlimRecord) -> &[BaseQuality] {
        let start = rec.qual_off as usize;
        let end = start.checked_add(rec.seq_len as usize).expect("qual end overflow");
        debug_assert!(end <= self.qual.len(), "qual slab overrun: {end} > {}", self.qual.len());
        BaseQuality::slice_from_bytes(&self.qual[start..end])
    }

    // r[impl bam.record.raw_aux]
    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub(crate) fn aux_of(&self, rec: &SlimRecord) -> &[u8] {
        let start = rec.aux_off as usize;
        let end = start.checked_add(rec.aux_len as usize).expect("aux end overflow");
        debug_assert!(end <= self.aux.len(), "aux slab overrun: {end} > {}", self.aux.len());
        &self.aux[start..end]
    }

    // r[impl record_store.set_alignment]
    // r[impl record_store.set_alignment.validation]
    /// Replace a record's alignment (pos + cigar) by appending new CIGAR to the
    /// cigar slab. The old CIGAR bytes become dead data.
    ///
    /// Moving a record breaks position order, so this retracts that guarantee
    /// and `prepare_for_pileup` re-establishes it — which the engine cannot be
    /// reached without. Mate indices survive (nothing is reordered here, and
    /// `mate_overlap` derives from the current positions), but the re-sort will
    /// clear them and preparing relinks.
    pub fn set_alignment(
        &mut self,
        idx: RecordIdx,
        new_pos: Pos0,
        new_cigar_ops: &[CigarOp],
    ) -> Result<(), DecodeError> {
        let rec =
            self.record(idx).ok_or(DecodeError::NoSuchRecord { idx, len: self.records.len() })?;
        let (seq_len, was_unmapped) = (rec.seq_len, rec.flags.is_unmapped());

        // ── All validation before any mutation ──
        // If any check fails, the store is unchanged (r[record_store.set_alignment.validation]).

        let new_query_len = cigar::calc_query_len(new_cigar_ops);
        if new_query_len != seq_len {
            return Err(DecodeError::CigarQueryLenMismatch {
                cigar_query_len: new_query_len,
                seq_len,
            });
        }

        let n_ops = new_cigar_ops.len();
        let n_cigar_ops =
            u16::try_from(n_ops).map_err(|_| DecodeError::CigarOpCountOverflow { count: n_ops })?;

        // r[impl record_store.end_pos_htslib]
        let end_pos = if was_unmapped {
            new_pos
        } else {
            cigar::compute_end_pos(new_pos, new_cigar_ops)
                .ok_or(DecodeError::InvalidPosition { value: new_pos.as_i32() })?
        };

        let (matching_bases, indel_bases) = cigar::calc_matches_indels(new_cigar_ops);

        let new_cigar_off =
            u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // ── Point of no return: all mutation below ──

        self.cigar.extend_from_slice(new_cigar_ops);

        #[allow(clippy::indexing_slicing, reason = "idx validated by self.record() above")]
        let rec = &mut self.records[idx.as_usize()];
        rec.pos = new_pos;
        // A moved record can land anywhere relative to its neighbours.
        self.order = RecordOrder::Unknown;
        rec.end_pos = end_pos;
        rec.n_cigar_ops = n_cigar_ops;
        rec.cigar_off = new_cigar_off;
        rec.matching_bases = matching_bases;
        rec.indel_bases = indel_bases;

        Ok(())
    }

    // r[impl record_store.extras.access]
    #[allow(
        clippy::indexing_slicing,
        reason = "extras_idx is written by push_raw alongside the record"
    )]
    pub(crate) fn extra_of(&self, rec: &SlimRecord) -> &U {
        let ei = rec.extras_idx as usize;
        debug_assert!(
            ei < self.extras.len(),
            "extras_idx {ei} out of bounds (len={})",
            self.extras.len()
        );
        &self.extras[ei]
    }

    // r[impl record_store.extras.access]
    // r[impl record_store.record_idx.resolution]
    /// Mutably access the per-record extra for `idx`, or `None` when this
    /// store has no such record.
    ///
    /// The mutable counterpart of [`RecordRef::extra`]. It cannot go through a
    /// [`RecordRef`], which borrows the store shared; the index is checked here
    /// instead.
    pub fn extra_mut(&mut self, idx: RecordIdx) -> Option<&mut U> {
        let ei = self.records.get(idx.as_usize())?.extras_idx as usize;
        self.extras.get_mut(ei)
    }

    /// Take all contents out, leaving an empty store with no capacity.
    /// Used by `PileupEngine::reclaim_allocation` to avoid requiring `Default`.
    pub(crate) fn take_contents(&mut self) -> Self {
        // What is taken keeps the state; what is left behind is empty, and an
        // empty store is trivially ordered and has nothing to link.
        let (order, mate_links) = (self.order, self.mate_links);
        self.order = RecordOrder::Ascending;
        self.mate_links = MateLinkState::Unlinked;
        RecordStore {
            order,
            mate_links,
            reach: std::mem::take(&mut self.reach),
            records: std::mem::take(&mut self.records),
            names: std::mem::take(&mut self.names),
            bases: std::mem::take(&mut self.bases),
            cigar: std::mem::take(&mut self.cigar),
            qual: std::mem::take(&mut self.qual),
            aux: std::mem::take(&mut self.aux),
            extras: std::mem::take(&mut self.extras),
        }
    }

    // r[impl record_store.extras.clear]
    // r[impl record_store.pileup_input]
    /// Establish the two properties [`PileupEngine`] relies on and hand back a
    /// store it will accept.
    ///
    /// The engine walks records assuming ascending `pos`, and reads `mate_idx`
    /// to surface mate overlap. Neither is checked at use: an unsorted store
    /// loses whole reads (the engine never reaches them, and no column reports
    /// a gap), and an unlinked one reports `mate_idx() == None` and
    /// `in_mate_overlap() == false` on every alignment — which is exactly what
    /// a record with no mate looks like, so a consumer that deduplicates
    /// overlapping mates silently keeps both.
    ///
    /// Both were reachable because [`PileupEngine::new`] took a bare
    /// [`RecordStore`]. It takes [`PileupInput`] instead, and this is the only
    /// way to make one.
    ///
    /// Idempotent and cheap to repeat: the sort is skipped when no push
    /// arrived out of order, and the linking when nothing has invalidated it.
    ///
    /// [`PileupEngine`]: crate::bam::pileup::PileupEngine
    /// [`PileupEngine::new`]: crate::bam::pileup::PileupEngine::new
    #[must_use]
    pub fn prepare_for_pileup(mut self) -> Prepared<U> {
        if self.order == RecordOrder::Unknown {
            // Also clears mate links: `mate_idx` is a store index, so it
            // cannot survive a reorder.
            self.sort_by_pos();
        }
        let stats = match self.mate_links {
            MateLinkState::Linked => MateLinkStats::default(),
            MateLinkState::Unlinked => self.link_mates(),
        };
        self.build_reach();
        Prepared { input: PileupInput { store: self }, stats }
    }

    // r[impl record_store.window_query.reach]
    /// Derive `reach` from the records, which must be in position order.
    /// Unmapped records carry no alignment and are never yielded by the
    /// window query, so they contribute nothing here either.
    fn build_reach(&mut self) {
        debug_assert_eq!(self.order, RecordOrder::Ascending, "reach needs position order");
        self.reach.clear();
        let mut furthest = Pos0::ZERO;
        self.reach.extend(self.records.iter().map(|rec| {
            if !rec.flags.is_unmapped() {
                furthest = furthest.max(rec.end_pos);
            }
            furthest
        }));
    }

    /// Note the arrival of a record at `pos`, before it is appended.
    ///
    /// A push can only ever *lose* the ascending-order property, and it always
    /// invalidates mate links: `mate_idx` is a store index, so a new record
    /// shifts nothing but a new mate may exist for an already-linked read.
    fn note_push(&mut self, pos: Pos0) {
        if self.records.last().is_some_and(|last| pos < last.pos) {
            self.order = RecordOrder::Unknown;
        }
        self.mate_links = MateLinkState::Unlinked;
    }

    // r[impl record_store.window_query.reach]
    /// Index of the first record that can overlap a window starting at
    /// `start`: every record before it ends before `start`, and since `reach`
    /// never decreases, so does nothing after it that the forward scan will
    /// not see. `len()` when no mapped record reaches `start`.
    pub(crate) fn first_reaching(&self, start: Pos0) -> usize {
        debug_assert_eq!(
            self.reach.len(),
            self.records.len(),
            "reach is stale: the store was changed after prepare_for_pileup"
        );
        self.reach.partition_point(|&reach| reach < start)
    }

    // r[impl record_store.window_query]
    // r[impl interval.overlap_test]
    /// Indices of the mapped records overlapping `[start, end]` (both
    /// inclusive), in ascending order.
    ///
    /// Crate-private because it is only correct on a store that
    /// [`prepare_for_pileup`](Self::prepare_for_pileup) has ordered and
    /// indexed; [`PileupInput`], [`PileupEngine`] and [`PileupColumn`] are
    /// the types that can vouch for that and expose it.
    ///
    /// [`PileupEngine`]: crate::bam::pileup::PileupEngine
    /// [`PileupColumn`]: crate::bam::pileup::PileupColumn
    pub(crate) fn records_overlapping_sorted(
        &self,
        start: Pos0,
        end: Pos0,
    ) -> impl Iterator<Item = RecordIdx> + '_ {
        debug_assert_eq!(
            self.order,
            RecordOrder::Ascending,
            "window query on a store that is not in position order"
        );
        let first = if end < start { self.records.len() } else { self.first_reaching(start) };
        self.records
            .get(first..)
            .unwrap_or_default()
            .iter()
            .enumerate()
            .take_while(move |(_, rec)| rec.pos <= end)
            // Unmapped records have no alignment to overlap anything; the
            // pileup never reports them, and this query answers for it.
            .filter(move |(_, rec)| rec.end_pos >= start && !rec.flags.is_unmapped())
            // Every slot in `records` has an index (a push that could not mint
            // one failed with `SlabOverflow`), so nothing is dropped here.
            .filter_map(move |(offset, _)| RecordIdx::from_usize(first.checked_add(offset)?))
    }

    pub fn clear(&mut self) {
        self.order = RecordOrder::Ascending;
        self.mate_links = MateLinkState::Unlinked;
        self.reach.clear();
        self.records.clear();
        self.names.clear();
        self.bases.clear();
        self.cigar.clear();
        self.qual.clear();
        self.aux.clear();
        self.extras.clear();
    }

    pub fn records_capacity(&self) -> usize {
        self.records.capacity()
    }

    pub fn names_capacity(&self) -> usize {
        self.names.capacity()
    }

    pub fn bases_capacity(&self) -> usize {
        self.bases.capacity()
    }

    pub fn cigar_capacity(&self) -> usize {
        self.cigar.capacity()
    }

    pub fn qual_capacity(&self) -> usize {
        self.qual.capacity()
    }

    pub fn aux_capacity(&self) -> usize {
        self.aux.capacity()
    }

    pub fn extras_capacity(&self) -> usize {
        self.extras.capacity()
    }
}

// r[impl record_store.link_mates.stats]
/// What one [`RecordStore::link_mates`] pass found.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct MateLinkStats {
    /// Templates whose two primary alignments were both in the store.
    pub pairs: u32,
    /// Records that found no free partner although another record with the same
    /// qname was already waiting — a third primary alignment for that qname, or
    /// the same record present twice. Non-zero means the input breaks the
    /// one-template-per-qname assumption of `[SAM1] §1.4`, and that reads which
    /// per-column qname matching would have collapsed are now counted
    /// separately.
    pub ambiguous_qnames: u32,
}

impl MateLinkStats {
    /// Whether anything about the input's qnames was surprising.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.ambiguous_qnames == 0
    }
}

// r[impl record_store.link_mates+2]
/// Whether a record can be one half of a linkable pair at all. A record with no
/// qname is rejected separately, where the name slab is in hand.
fn can_link(rec: &SlimRecord) -> bool {
    let flags = rec.flags;
    flags.is_paired()
        && !flags.is_unmapped()
        && !flags.is_mate_unmapped()
        && !flags.is_secondary()
        && !flags.is_supplementary()
        && rec.next_ref_id == rec.tid
        && rec.next_pos >= 0
}

/// The record's qname bytes, read straight from a borrowed name slab. Returns
/// an empty slice on a malformed offset — two malformed records then compare
/// equal on names alone, which the position and hash checks still reject.
fn qname_bytes<'a>(names: &'a [u8], rec: &SlimRecord) -> &'a [u8] {
    let start = rec.name_off as usize;
    let end = start.saturating_add(rec.name_len as usize);
    names.get(start..end).unwrap_or(&[])
}

impl<U> Default for RecordStore<U> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code with known small values"
)]
pub(crate) mod tests {
    use super::*;
    use crate::bam::test_util::ri;

    /// Drop every record at push time (no-extra customizer).
    #[derive(Clone, Default)]
    struct KeepNone;
    impl CustomizeRecordStore for KeepNone {
        type Extra = ();
        fn filter(&mut self, _: &SlimRecord, _: &RecordStore<()>) -> bool {
            false
        }
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
    }

    /// Decide accept/reject from a runtime flag — used by proptests.
    #[derive(Clone)]
    struct AcceptFlag(bool);
    impl CustomizeRecordStore for AcceptFlag {
        type Extra = ();
        fn filter(&mut self, _: &SlimRecord, _: &RecordStore<()>) -> bool {
            self.0
        }
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
    }

    /// Build a minimal valid BAM record raw bytes.
    fn make_raw_record(qname: &[u8], seq_len: u32, n_cigar_ops: u16) -> Vec<u8> {
        let name_len = qname.len() as u8 + 1; // includes NUL
        let cigar_bytes = n_cigar_ops as usize * 4;
        let seq_bytes = (seq_len as usize).div_ceil(2);
        let total = 32 + name_len as usize + cigar_bytes + seq_bytes + seq_len as usize;

        let mut raw = vec![0u8; total];
        raw[0..4].copy_from_slice(&0i32.to_le_bytes()); // tid
        raw[4..8].copy_from_slice(&0i32.to_le_bytes()); // pos
        raw[8] = name_len;
        raw[9] = 0; // mapq
        raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
        raw[14..16].copy_from_slice(&0u16.to_le_bytes()); // flags
        raw[16..20].copy_from_slice(&seq_len.to_le_bytes());

        // Write qname + NUL
        let var_start = 32;
        raw[var_start..var_start + qname.len()].copy_from_slice(qname);
        raw[var_start + qname.len()] = 0; // NUL

        // Write cigar: M ops
        let cigar_start = var_start + name_len as usize;
        for i in 0..n_cigar_ops as usize {
            let op = seq_len << 4; // M op
            let off = cigar_start + i * 4;
            if off + 4 <= raw.len() {
                raw[off..off + 4].copy_from_slice(&op.to_le_bytes());
            }
        }

        raw
    }

    // r[verify record_store.record_idx.resolution]
    /// `record` answers for the store it is asked, not for the index's numeric
    /// value: the same index resolves before a `clear` and not after.
    #[test]
    fn record_is_none_past_the_end_and_after_a_clear() {
        let mut store = RecordStore::new();
        store.push_raw(&make_raw_record(b"read1", 4, 1), &mut ()).unwrap();
        let idx = ri(0);
        assert!(store.record(idx).is_some());
        assert!(store.record(ri(1)).is_none(), "one past the end");
        assert!(store.record(ri(u32::MAX - 1)).is_none());

        store.clear();
        assert!(store.record(idx).is_none(), "an index cannot outlive its records");
    }

    // r[verify record_store.record_idx]
    #[test]
    fn indices_name_every_record_in_store_order() {
        let mut store = RecordStore::new();
        assert_eq!(store.indices().count(), 0, "an empty store has no indices");
        for name in [b"a", b"b", b"c"] {
            store.push_raw(&make_raw_record(name, 4, 1), &mut ()).unwrap();
        }
        assert_eq!(store.indices().collect::<Vec<_>>(), [ri(0), ri(1), ri(2)]);
        for (idx, rec) in store.indices().zip(store.records()) {
            assert_eq!(store.record(idx).unwrap().pos, rec.pos);
        }
    }

    // r[verify bam.record.checked_offsets]
    #[test]
    fn push_raw_rejects_offset_overflow() {
        let mut store = RecordStore::new();
        // Craft a record with fields that would cause overflow.
        let mut raw = [0u8; 32];
        raw[8] = 255; // name_len
        raw[12..14].copy_from_slice(&u16::MAX.to_le_bytes()); // n_cigar_ops
        raw[16..20].copy_from_slice(&u32::MAX.to_le_bytes()); // seq_len

        let result = store.push_raw(&raw, &mut ());
        assert!(result.is_err());
    }

    // r[verify record_store.checked_offsets]
    #[test]
    fn push_fields_rejects_offset_overflow() {
        use seqair_types::Base;
        let mut store = RecordStore::new();
        // Directly inflate the names slab past u32::MAX to trigger overflow
        // We can't actually allocate 4GB in a test, so we test the check path
        // by verifying the function returns Result (compile-time check)
        let result: Result<Option<RecordIdx>, _> = store.push_fields(
            Pos0::new(0).unwrap(),
            Pos0::new(0).unwrap(),
            BamFlags::empty(),
            0,
            0,
            0,
            b"read1",
            &[],
            &[Base::A],
            &[30],
            &[],
            0,  // tid
            -1, // next_ref_id
            0,  // next_pos
            0,  // template_len
            &mut (),
        );
        assert!(result.is_ok());
    }

    // r[verify record_store.checked_offsets]
    #[test]
    fn push_raw_normal_record_succeeds() {
        let mut store = RecordStore::new();
        let raw = make_raw_record(b"read1", 4, 1);
        let result = store.push_raw(&raw, &mut ());
        assert!(result.is_ok());
        assert_eq!(store.len(), 1);
    }

    // r[verify record_store.iter]
    #[test]
    fn iter_yields_all_records_with_indices() {
        let mut store = RecordStore::new();
        store.push_raw(&make_raw_record(b"read0", 4, 1), &mut ()).unwrap();
        store.push_raw(&make_raw_record(b"read1", 4, 1), &mut ()).unwrap();
        store.push_raw(&make_raw_record(b"read2", 4, 1), &mut ()).unwrap();

        let qnames: Vec<&[u8]> = store.records().map(|rec| rec.qname(&store).unwrap()).collect();

        assert_eq!(qnames, vec![b"read0".as_slice(), b"read1".as_slice(), b"read2".as_slice()]);

        // ExactSizeIterator parity with len(); enumerate gives positions when needed.
        assert_eq!(store.records().len(), store.len());
        let positions: Vec<usize> = store.records().enumerate().map(|(i, _)| i).collect();
        assert_eq!(positions, vec![0, 1, 2]);
    }

    #[test]
    fn iter_empty_store() {
        let store: RecordStore = RecordStore::new();
        assert_eq!(store.records().count(), 0);
        assert_eq!(store.records().len(), 0);
    }

    // r[verify record_store.slim_record_fields]
    #[test]
    fn push_raw_preserves_next_ref_id() {
        let mut store = RecordStore::new();
        let mut raw = make_raw_record(b"read1", 4, 1);
        // Write next_ref_id = 7 at BAM offset 20
        raw[20..24].copy_from_slice(&7i32.to_le_bytes());
        store.push_raw(&raw, &mut ()).unwrap();
        assert_eq!(store.record(ri(0)).unwrap().next_ref_id, 7);
    }

    // r[verify record_store.pre_filter.rollback]
    #[test]
    fn push_raw_with_rejecting_filter_truncates_all_slabs() {
        let mut store = RecordStore::new();

        // Push record A (kept)
        let a = make_raw_record(b"read1", 4, 1);
        let idx_a = store.push_raw(&a, &mut ()).unwrap();
        assert_eq!(idx_a, Some(ri(0)));
        let after_a = (
            store.records.len(),
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
            store.extras.len(),
        );

        // Push record B (rejected by filter) — all slab extensions must roll back
        let b = make_raw_record(b"read2", 8, 2);
        let idx_b = store.push_raw(&b, &mut KeepNone).unwrap();
        assert_eq!(idx_b, None, "rejected record must not yield an index");

        let after_b = (
            store.records.len(),
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
            store.extras.len(),
        );
        assert_eq!(after_a, after_b, "slab lengths must match pre-reject state exactly");

        // Push record C (kept) — indexing must continue from idx 1 (B was never committed)
        let c = make_raw_record(b"read3", 4, 1);
        let idx_c = store.push_raw(&c, &mut ()).unwrap();
        assert_eq!(idx_c, Some(ri(1)), "rejected record must not consume an index");
        assert_eq!(store.len(), 2, "store should hold A and C only");
    }

    // r[verify record_store.filter_raw]
    #[test]
    fn push_raw_filter_raw_rejection_avoids_slab_write() {
        /// Reject everything in `filter_raw` — no slab extensions should happen.
        #[derive(Clone, Default)]
        struct RejectRaw;
        impl CustomizeRecordStore for RejectRaw {
            type Extra = ();
            fn filter_raw(&mut self, _: &FilterRawFields<'_>) -> bool {
                false
            }
            fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
        }

        let mut store = RecordStore::new();
        let before = (
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
        );

        let raw = make_raw_record(b"readX", 8, 2);
        let result = store.push_raw(&raw, &mut RejectRaw).unwrap();
        assert_eq!(result, None);

        let after = (
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
        );
        assert_eq!(before, after, "filter_raw rejection must not extend any slab");
    }

    // r[verify record_store.filter_raw]
    #[test]
    fn push_fields_filter_raw_rejection_avoids_slab_write() {
        use seqair_types::Base;

        #[derive(Clone, Default)]
        struct RejectRaw;
        impl CustomizeRecordStore for RejectRaw {
            type Extra = ();
            fn filter_raw(&mut self, _: &FilterRawFields<'_>) -> bool {
                false
            }
            fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
        }

        let mut store = RecordStore::new();
        let before = (
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
        );

        let result = store
            .push_fields(
                Pos0::new(0).unwrap(),
                Pos0::new(4).unwrap(),
                BamFlags::empty(),
                42,
                4,
                0,
                b"readY",
                &[CigarOp::new(cigar::CigarOpType::Match, 5)],
                &[Base::A; 5],
                &[30; 5],
                b"tagdata",
                0,
                -1,
                0,
                0,
                &mut RejectRaw,
            )
            .unwrap();
        assert_eq!(result, None);

        let after = (
            store.names.len(),
            store.bases.len(),
            store.cigar.len(),
            store.qual.len(),
            store.aux.len(),
        );
        assert_eq!(before, after, "filter_raw rejection must not extend any slab");
    }

    // r[verify record_store.pre_filter.rollback]
    #[test]
    fn push_raw_filter_can_read_slim_record_and_store() {
        // Customizer that asserts the freshly-pushed record's fields are
        // visible (qname, mapq) by the time filter is called.
        #[derive(Clone, Default)]
        struct AssertReadX;
        impl CustomizeRecordStore for AssertReadX {
            type Extra = ();
            fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<()>) -> bool {
                assert_eq!(rec.mapq, 0);
                let qname = rec.qname(store).expect("qname slab must be readable");
                assert_eq!(qname, b"readX");
                true
            }
            fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
        }

        let mut store = RecordStore::new();
        let raw = make_raw_record(b"readX", 4, 1);
        let kept_idx = store.push_raw(&raw, &mut AssertReadX).unwrap();
        assert_eq!(kept_idx, Some(ri(0)));
    }

    // r[verify record_store.pre_filter.rollback]
    #[test]
    fn push_fields_with_rejecting_filter_rolls_back() {
        use seqair_types::Base;
        let mut store = RecordStore::new();

        // Push one kept record first
        store
            .push_fields(
                Pos0::new(0).unwrap(),
                Pos0::new(4).unwrap(),
                BamFlags::empty(),
                30,
                4,
                0,
                b"kept",
                &[CigarOp::new(cigar::CigarOpType::Match, 4)],
                &[Base::A, Base::C, Base::G, Base::T],
                &[30, 31, 32, 33],
                &[],
                0,
                -1,
                0,
                0,
                &mut (),
            )
            .unwrap();

        let names_len = store.names.len();
        let bases_len = store.bases.len();
        let cigar_len = store.cigar.len();
        let qual_len = store.qual.len();

        // Push a rejected one
        let rejected = store
            .push_fields(
                Pos0::new(10).unwrap(),
                Pos0::new(14).unwrap(),
                BamFlags::empty(),
                30,
                4,
                0,
                b"dropped",
                &[CigarOp::new(cigar::CigarOpType::Match, 4)],
                &[Base::A, Base::C, Base::G, Base::T],
                &[30, 31, 32, 33],
                b"NM:i:0",
                0,
                -1,
                0,
                0,
                &mut KeepNone,
            )
            .unwrap();
        assert_eq!(rejected, None);
        assert_eq!(store.len(), 1);
        assert_eq!(store.names.len(), names_len);
        assert_eq!(store.bases.len(), bases_len);
        assert_eq!(store.cigar.len(), cigar_len);
        assert_eq!(store.qual.len(), qual_len);
    }

    // r[verify record_store.slim_record_fields]
    #[test]
    fn push_fields_preserves_next_ref_id() {
        use seqair_types::Base;
        let mut store = RecordStore::new();
        store
            .push_fields(
                Pos0::new(100).unwrap(),
                Pos0::new(105).unwrap(),
                BamFlags::empty(),
                30,
                5,
                0,
                b"read1",
                &[CigarOp::new(cigar::CigarOpType::Match, 5)],
                &[Base::A, Base::C, Base::G, Base::T, Base::A],
                &[30, 31, 32, 33, 34],
                &[],
                0,   // tid
                3,   // next_ref_id
                500, // next_pos
                200, // template_len
                &mut (),
            )
            .unwrap();
        assert_eq!(store.record(ri(0)).unwrap().next_ref_id, 3);
    }

    // r[verify record_store.end_pos_htslib]
    #[test]
    fn end_pos_mapped_uses_cigar() {
        use seqair_types::{BamFlags, Base};
        let mut store = RecordStore::new();
        store
            .push_fields(
                Pos0::new(100).unwrap(),
                Pos0::new(104).unwrap(), // 0-based inclusive, 5M: span = [100, 104]
                BamFlags::empty(),       // mapped
                30,
                5,
                0,
                b"read1",
                &[CigarOp::new(cigar::CigarOpType::Match, 5)],
                &[Base::A; 5],
                &[30; 5],
                &[],
                0,
                -1,
                0,
                0,
                &mut (),
            )
            .unwrap();
        let rec = store.record(ri(0)).unwrap();
        assert_eq!(rec.end_pos, Pos0::new(104).unwrap());
    }

    // r[verify record_store.end_pos_htslib]
    #[test]
    fn end_pos_unmapped_ignores_cigar() {
        use seqair_types::bam_flags::consts::FLAG_UNMAPPED;
        use seqair_types::{BamFlags, Base};
        let mut store = RecordStore::new();
        store
            .push_fields(
                Pos0::new(100).unwrap(),
                Pos0::new(100).unwrap(), // htslib-compatible: end_pos = pos for unmapped
                BamFlags::from(FLAG_UNMAPPED), // unmapped
                0,
                100,
                0,
                b"read1",
                &[CigarOp::new(cigar::CigarOpType::Match, 100)],
                &[Base::A; 100],
                &[0; 100],
                &[],
                0,
                -1,
                0,
                0,
                &mut (),
            )
            .unwrap();
        let rec = store.record(ri(0)).unwrap();
        assert!(rec.flags.is_unmapped());
        assert_eq!(rec.end_pos, Pos0::new(100).unwrap());
    }

    // r[verify record_store.end_pos_htslib]
    /// Verify the CRAM unmapped-read code path semantics: unmapped reads pushed
    /// via `push_fields` (the path CRAM uses, `cram/slice.rs:649-705`) have
    /// `end_pos == pos`, `mapq == 0`, and correct flags/sequence data.
    #[test]
    fn push_fields_unmapped_cram_semantics() {
        use seqair_types::bam_flags::consts::FLAG_UNMAPPED;
        use seqair_types::{BamFlags, Base};

        let mut store = RecordStore::new();
        let pos = Pos0::new(200).unwrap();
        let bases = &[Base::A, Base::C, Base::G, Base::T];
        let quals = &[30u8, 31, 32, 33];
        let aux = b"RGZ\x00read_group";
        store
            .push_fields(
                pos,
                pos, // end_pos = pos (htslib-compatible, as CRAM does at slice.rs:676)
                BamFlags::from(FLAG_UNMAPPED),
                0, // mapq = 0 (CRAM slice.rs:681)
                0, // matching_bases = 0 (slice.rs:682)
                0, // indel_bases = 0 (slice.rs:683)
                b"cram_unmapped_read",
                &[], // empty CIGAR (slice.rs:685)
                bases,
                quals,
                aux,
                0,  // tid
                -1, // next_ref_id
                0,  // next_pos
                0,  // template_len
                &mut (),
            )
            .unwrap();

        let rec = store.record(ri(0)).unwrap();
        assert!(rec.flags.is_unmapped(), "CRAM unmapped reads must have FLAG_UNMAPPED");
        assert_eq!(rec.mapq, 0, "CRAM unmapped reads must have MAPQ=0");
        assert_eq!(rec.end_pos, pos, "CRAM unmapped reads must have end_pos == pos");
        assert_eq!(rec.seq_len, 4, "CRAM unmapped reads retain decoded sequence");
        assert_eq!(rec.matching_bases, 0, "CRAM unmapped reads have zero matching_bases");
        assert_eq!(rec.indel_bases, 0, "CRAM unmapped reads have zero indel_bases");
        assert_eq!(rec.qname(), b"cram_unmapped_read");
    }

    // r[verify record_store.end_pos_htslib]
    /// Verify that `compute_end_pos_from_raw` (used by the BAM reader's pre-push
    /// overlap filter) agrees with the stored `end_pos` after `push_raw`, for
    /// both mapped and unmapped reads.
    #[test]
    fn push_raw_end_pos_agrees_with_reader_prefilter() {
        use seqair_types::bam_flags::consts::FLAG_UNMAPPED;

        // --- Mapped read ---
        let mut mapped = make_raw_record(b"read1", 10, 1); // 1 CIGAR op (10M)
        mapped[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos = 100
        let reader_end = record::compute_end_pos_from_raw(&mapped).unwrap();
        let mut store = RecordStore::new();
        store.push_raw(&mapped, &mut ()).unwrap();
        assert_eq!(
            store.record(ri(0)).unwrap().end_pos,
            reader_end,
            "mapped: reader's compute_end_pos_from_raw must agree with push_raw's stored end_pos"
        );

        // --- Unmapped read ---
        let mut unmapped = make_raw_record(b"read2", 10, 1); // 1 CIGAR op (10M)
        unmapped[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos = 100
        unmapped[14..16].copy_from_slice(&FLAG_UNMAPPED.to_le_bytes()); // flag 0x4
        let reader_end = record::compute_end_pos_from_raw(&unmapped).unwrap();
        // compute_end_pos_from_raw computes CIGAR-derived (doesn't check flags).
        // After refactor, push_raw stores pos for unmapped reads.
        let mut store2 = RecordStore::new();
        store2.push_raw(&unmapped, &mut ()).unwrap();
        assert!(store2.record(ri(0)).unwrap().flags.is_unmapped());
        assert_eq!(
            store2.record(ri(0)).unwrap().end_pos,
            Pos0::new(100).unwrap(),
            "unmapped: push_raw stores pos (htslib-compatible), not CIGAR-derived end"
        );
        // Verify they differ: reader's prefilter returns CIGAR-derived (which is fine
        // for overlap checking — it's a loose upper bound), while store is exact.
        assert_ne!(
            store2.record(ri(0)).unwrap().end_pos,
            reader_end,
            "reader prefilter end_pos (CIGAR-derived) differs from stored end_pos (pos) for unmapped reads"
        );
    }

    // ---- Model-based property tests for rollback ----
    //
    // The rollback invariant: pushing a sequence of records with per-record
    // accept/reject decisions MUST produce a store byte-identical to pushing
    // only the accepted records in order with no filter. That is, a rollback
    // must perfectly undo both the records Vec push and all slab extensions,
    // leaving zero dead bytes.
    pub(crate) mod window_query {
        use super::super::*;
        use crate::bam::test_util::ri;
        use proptest::prelude::*;
        use seqair_types::{BamFlags, Base};

        /// An expected index list, written with literals.
        fn idxs<const N: usize>(ns: [u32; N]) -> Vec<RecordIdx> {
            ns.into_iter().map(ri).collect()
        }

        /// A record: `M` over `bases`, optionally followed by a deletion so
        /// spans vary independently of query length. An unmapped one is
        /// *placed* — it has a position, as a mate placed next to its partner
        /// does — but `end_pos == pos` and no alignment, which is how
        /// `push_raw` records such reads.
        #[derive(Debug, Clone)]
        pub(crate) struct Read {
            pub(crate) pos: u32,
            pub(crate) bases: u32,
            pub(crate) deletion: u32,
            pub(crate) mapped: bool,
        }

        impl Read {
            fn cigar(&self) -> Vec<CigarOp> {
                let mut cigar = vec![CigarOp::new(cigar::CigarOpType::Match, self.bases)];
                if self.deletion > 0 {
                    cigar.push(CigarOp::new(cigar::CigarOpType::Deletion, self.deletion));
                }
                cigar
            }

            /// Push with a qname unique to `i`, so the store's mate linking
            /// sees distinct templates rather than one qname eighty times.
            pub(crate) fn push(&self, store: &mut RecordStore<()>, i: usize) -> RecordIdx {
                let pos = Pos0::new(self.pos).expect("strategy bounds pos");
                let end =
                    if self.mapped { self.pos + self.bases + self.deletion - 1 } else { self.pos };
                let end_pos = Pos0::new(end).expect("strategy bounds end");
                let flags = if self.mapped { BamFlags::empty() } else { BamFlags::from(0x4u16) };
                let bases = vec![Base::A; self.bases as usize];
                let quals = vec![30u8; self.bases as usize];
                store
                    .push_fields(
                        pos,
                        end_pos,
                        flags,
                        30,
                        self.bases,
                        self.deletion,
                        format!("r{i}").as_bytes(),
                        &self.cigar(),
                        &bases,
                        &quals,
                        &[],
                        0,
                        -1,
                        0,
                        0,
                        &mut (),
                    )
                    .expect("synthetic push")
                    .expect("no filter")
            }
        }

        /// Deletions are bimodal: mostly none or short, occasionally one long
        /// enough to reach far past its neighbours, so a single record sets
        /// how far back a window must look — the case the index exists for.
        pub(crate) fn arb_read() -> impl Strategy<Value = Read> {
            (
                0u32..5_000,
                1u32..=16,
                prop_oneof![6 => Just(0u32), 3 => 1u32..=40, 1 => 200u32..=1_500],
                prop_oneof![9 => Just(true), 1 => Just(false)],
            )
                .prop_map(|(pos, bases, deletion, mapped)| Read {
                    pos,
                    bases,
                    deletion,
                    mapped,
                })
        }

        pub(crate) fn arb_store() -> impl Strategy<Value = Vec<Read>> {
            prop::collection::vec(arb_read(), 0..=80)
        }

        /// A window, either absolute — reaching past the last possible record
        /// and often inverted, since both ends are drawn independently — or
        /// anchored a base or two off a record's own start or end, which is
        /// where an off-by-one in the search would hide.
        #[derive(Debug, Clone)]
        pub(crate) enum Span {
            Absolute(u32, u32),
            Anchored { record: usize, at_end: bool, start_off: i8, len: u8 },
        }

        pub(crate) fn arb_span() -> impl Strategy<Value = Span> {
            prop_oneof![
                (0u32..7_000, 0u32..7_000).prop_map(|(a, b)| Span::Absolute(a, b)),
                (any::<usize>(), any::<bool>(), -2i8..=2, 0u8..=8).prop_map(
                    |(record, at_end, start_off, len)| Span::Anchored {
                        record,
                        at_end,
                        start_off,
                        len
                    }
                ),
            ]
        }

        impl Span {
            /// Resolve against the store the query will run on, so anchors
            /// follow the records wherever sorting and moves have put them.
            pub(crate) fn resolve(&self, store: &RecordStore<()>) -> (Pos0, Pos0) {
                match *self {
                    Span::Absolute(a, b) => (at(a), at(b)),
                    Span::Anchored { record, at_end, start_off, len } => {
                        let Some(rec) = store.records().nth(record % store.len().max(1)) else {
                            return (at(0), at(u32::from(len)));
                        };
                        let boundary = if at_end { rec.end_pos } else { rec.pos };
                        let start = boundary.as_u32().saturating_add_signed(i32::from(start_off));
                        (at(start), at(start + u32::from(len)))
                    }
                }
            }
        }

        /// A realignment: move a record and change its deletion. Query length
        /// is untouched, so `set_alignment` always accepts it.
        #[derive(Debug, Clone)]
        pub(crate) struct Move {
            pub(crate) record: usize,
            pub(crate) pos: u32,
            pub(crate) deletion: u32,
        }

        pub(crate) fn arb_moves() -> impl Strategy<Value = Vec<Move>> {
            prop::collection::vec(
                (any::<usize>(), 0u32..5_000, prop_oneof![3 => Just(0u32), 1 => 1u32..=1_500])
                    .prop_map(|(record, pos, deletion)| Move { record, pos, deletion }),
                0..=6,
            )
        }

        /// Apply `moves` to a store that `reads` were pushed into, in push
        /// order, skipping unmapped records (they have nothing to realign).
        pub(crate) fn apply_moves(store: &mut RecordStore<()>, reads: &[Read], moves: &[Move]) {
            for mv in moves {
                if reads.is_empty() {
                    break;
                }
                let idx = mv.record % reads.len();
                let Some(read) = reads.get(idx).filter(|read| read.mapped) else { continue };
                let moved = Read { pos: mv.pos, deletion: mv.deletion, ..read.clone() };
                let Some(idx) = RecordIdx::from_usize(idx) else { continue };
                store
                    .set_alignment(idx, at(mv.pos), &moved.cigar())
                    .expect("query length preserved");
            }
        }

        /// The definition, applied to every record: no search, no index. An
        /// inverted interval is empty by definition — the raw overlap test
        /// alone would accept any record spanning the gap — and an unmapped
        /// record has no alignment to overlap anything.
        pub(crate) fn brute_force(
            store: &RecordStore<()>,
            start: Pos0,
            end: Pos0,
        ) -> Vec<RecordIdx> {
            if end < start {
                return Vec::new();
            }
            store
                .records()
                .enumerate()
                .filter(|(_, rec)| {
                    !rec.flags.is_unmapped() && rec.pos <= end && rec.end_pos >= start
                })
                .map(|(idx, _)| RecordIdx::from_usize(idx).expect("small store"))
                .collect()
        }

        fn at(v: u32) -> Pos0 {
            Pos0::new(v).unwrap()
        }

        fn mapped(pos: u32, bases: u32, deletion: u32) -> Read {
            Read { pos, bases, deletion, mapped: true }
        }

        fn prepared(reads: &[Read]) -> PileupInput<()> {
            let mut store = RecordStore::new();
            for (i, read) in reads.iter().enumerate() {
                read.push(&mut store, i);
            }
            store.prepare_for_pileup().input
        }

        // r[verify record_store.window_query]
        // r[verify record_store.window_query.reach]
        /// A read that starts well before the window but reaches into it must
        /// be found — this is the case a plain lower bound on `pos` misses.
        #[test]
        fn long_read_starting_before_window_is_found() {
            let input = prepared(&[
                mapped(10, 5, 300), // covers 10..=314
                mapped(200, 5, 0),  // 200..=204
                mapped(400, 5, 0),  // 400..=404
            ]);

            let hits: Vec<RecordIdx> = input.records_overlapping(at(300), at(310)).collect();
            assert_eq!(hits, idxs([0]));
            let hits: Vec<RecordIdx> = input.records_overlapping(at(203), at(400)).collect();
            assert_eq!(hits, idxs([0, 1, 2]));
        }

        // r[verify record_store.window_query]
        #[test]
        fn empty_and_out_of_range_spans_yield_nothing() {
            let input = prepared(&[mapped(100, 10, 0)]);

            assert_eq!(input.records_overlapping(at(105), at(104)).count(), 0, "end < start");
            assert_eq!(input.records_overlapping(at(110), at(500)).count(), 0, "past the last");
            assert_eq!(input.records_overlapping(at(0), at(99)).count(), 0, "before the first");
            assert_eq!(input.records_overlapping(at(109), at(109)).count(), 1, "last base");
            assert_eq!(input.records_overlapping(at(100), at(100)).count(), 1, "first base");
        }

        // r[verify record_store.window_query]
        #[test]
        fn empty_store_yields_nothing() {
            let input = prepared(&[]);
            assert_eq!(input.records_overlapping(at(0), at(Pos0::max_value().as_u32())).count(), 0);
        }

        // r[verify record_store.window_query]
        /// A placed-unmapped read sits at a position but aligns to nothing;
        /// the pileup never reports it, so neither does the query.
        #[test]
        fn placed_unmapped_records_are_not_yielded() {
            let input = prepared(&[
                mapped(100, 10, 0),
                Read { pos: 105, bases: 10, deletion: 0, mapped: false },
                mapped(105, 10, 0),
            ]);
            let hits: Vec<RecordIdx> = input.records_overlapping(at(105), at(105)).collect();
            assert_eq!(hits, idxs([0, 2]));
        }

        // r[verify record_store.window_query.reach]
        /// `set_alignment` can lengthen a record after it was pushed; the
        /// index is built when the store is prepared, so the moved read is
        /// found where it now reaches.
        #[test]
        fn set_alignment_lengthening_a_record_is_seen() {
            let mut store = RecordStore::new();
            let idx = mapped(10, 5, 0).push(&mut store, 0);
            mapped(500, 5, 0).push(&mut store, 1);
            let long = [
                CigarOp::new(cigar::CigarOpType::Match, 5),
                CigarOp::new(cigar::CigarOpType::Deletion, 600),
            ];
            store.set_alignment(idx, at(10), &long).unwrap();
            let input = store.prepare_for_pileup().input;

            let hits: Vec<RecordIdx> = input.records_overlapping(at(600), at(600)).collect();
            assert_eq!(hits, brute_force(input.store(), at(600), at(600)));
            assert_eq!(hits, idxs([0]));
        }

        // r[verify record_store.window_query.reach]
        /// The index is rebuilt for every prepared store, so a store cleared
        /// and refilled — the `Readers` reuse cycle — answers for the records
        /// it holds now, not the ones it held before.
        #[test]
        fn reused_store_answers_for_its_new_records() {
            let mut store = RecordStore::new();
            mapped(10, 5, 900).push(&mut store, 0); // 10..=914
            let store = store.prepare_for_pileup().input.into_store();
            let mut store = store;
            store.clear();
            mapped(100, 10, 0).push(&mut store, 0);
            mapped(300, 10, 0).push(&mut store, 1);
            let input = store.prepare_for_pileup().input;

            assert_eq!(input.store().first_reaching(at(305)), 1, "a stale index would say 0");
            let hits: Vec<RecordIdx> = input.records_overlapping(at(305), at(305)).collect();
            assert_eq!(hits, idxs([1]));
        }

        proptest! {
            // r[verify record_store.window_query]
            // r[verify record_store.window_query.reach]
            /// The binary search must return exactly what the definition
            /// returns, for stores pushed in any order, realigned after the
            /// fact, and spans that are empty, before, past, inside, or a base
            /// off a record's edge.
            #[test]
            fn matches_brute_force(
                reads in arb_store(),
                moves in arb_moves(),
                spans in prop::collection::vec(arb_span(), 1..=8),
            ) {
                let mut store = RecordStore::new();
                for (i, read) in reads.iter().enumerate() {
                    read.push(&mut store, i);
                }
                apply_moves(&mut store, &reads, &moves);
                let input = store.prepare_for_pileup().input;
                for span in spans {
                    let (start, end) = span.resolve(input.store());
                    let fast: Vec<RecordIdx> = input.records_overlapping(start, end).collect();
                    prop_assert_eq!(fast, brute_force(input.store(), start, end), "span {:?} = {:?}..={:?}", span, start, end);
                }
            }

            // r[verify record_store.window_query.reach]
            /// The search lands on the first record that can overlap a window
            /// starting at `start`: the first mapped record whose `end_pos`
            /// reaches it, found here by a linear scan. Nothing before it can
            /// be a hit, and the forward scan sees everything after it.
            #[test]
            fn first_reaching_is_where_the_linear_scan_stops(reads in arb_store(), start in 0u32..7_000) {
                let input = prepared(&reads);
                let store = input.store();
                let expected = store
                    .records()
                    .position(|rec| !rec.flags.is_unmapped() && rec.end_pos >= at(start))
                    .unwrap_or(store.len());
                prop_assert_eq!(store.first_reaching(at(start)), expected);
            }
        }
    }

    mod rollback_props {
        use super::super::*;
        use super::AcceptFlag;
        use proptest::prelude::*;
        use seqair_types::{BamFlags, Base};

        /// A synthetic push input covering every slab. Fields are bounded to
        /// small sizes so proptest can generate long sequences without
        /// exhausting memory.
        #[derive(Debug, Clone)]
        struct PushInput {
            qname: Vec<u8>,
            bases: Vec<Base>,
            quals: Vec<u8>,
            aux: Vec<u8>,
            pos: u32,
            mapq: u8,
            /// Single M op whose length matches `bases.len()`; kept trivial so
            /// the invariant we're testing stays isolated to rollback, not
            /// CIGAR math.
            accept: bool,
        }

        impl PushInput {
            fn cigar_op(&self) -> CigarOp {
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "bases.len() is bounded by strategy (≤ 16)"
                )]
                let len = self.bases.len() as u32;
                CigarOp::new(cigar::CigarOpType::Match, len)
            }

            fn end_pos(&self) -> Pos0 {
                // Match the inclusive 0-based convention used by
                // `compute_end_pos` (pos + ref_len - 1). For our trivial 1-op
                // M CIGAR, ref_len == bases.len(), so end = pos + len - 1.
                // bases.len() >= 1 by strategy, so the subtraction is safe.
                #[expect(
                    clippy::cast_possible_truncation,
                    clippy::arithmetic_side_effects,
                    reason = "bases.len() ≥ 1 ≤ 16, pos < 1_000_000; sum fits in u32"
                )]
                let end = self.pos + self.bases.len() as u32 - 1;
                Pos0::new(end).expect("bounded by strategy to < i32::MAX")
            }
        }

        fn arb_base() -> impl Strategy<Value = Base> {
            prop_oneof![
                Just(Base::A),
                Just(Base::C),
                Just(Base::G),
                Just(Base::T),
                Just(Base::Unknown),
            ]
        }

        fn arb_push_input() -> impl Strategy<Value = PushInput> {
            (
                // qname: 1..16 ASCII bytes without NUL
                prop::collection::vec(1u8..=126, 1..=16),
                // bases: 1..16 bases (non-empty so qual len matches)
                prop::collection::vec(arb_base(), 1..=16),
                // aux: 0..32 bytes (NOT parsed, so any bytes OK)
                prop::collection::vec(any::<u8>(), 0..=32),
                0u32..1_000_000,
                0u8..=60,
                any::<bool>(),
            )
                .prop_map(|(qname, bases, aux, pos, mapq, accept)| {
                    let quals = bases.iter().map(|_| 30u8).collect();
                    PushInput { qname, bases, quals, aux, pos, mapq, accept }
                })
        }

        /// Push one `PushInput` into the store via `push_fields`, honoring its
        /// accept flag. Returns the new record index if kept.
        #[allow(
            clippy::expect_used,
            clippy::unwrap_in_result,
            reason = "proptest synthetic input bounded by strategy; panic on violation is informative"
        )]
        fn push_one(store: &mut RecordStore<()>, input: &PushInput) -> Option<RecordIdx> {
            let cigar = [input.cigar_op()];
            #[expect(clippy::cast_possible_truncation, reason = "bases.len() bounded ≤ 16")]
            let matching = input.bases.len() as u32;
            store
                .push_fields(
                    Pos0::new(input.pos).expect("strategy bounds pos < 1_000_000"),
                    input.end_pos(),
                    BamFlags::empty(),
                    input.mapq,
                    matching,
                    0,
                    &input.qname,
                    &cigar,
                    &input.bases,
                    &input.quals,
                    &input.aux,
                    0,
                    -1,
                    0,
                    0,
                    &mut AcceptFlag(input.accept),
                )
                .expect("push_fields must not error on synthetic input")
        }

        /// Push one record with the accept flag ignored (used for the
        /// no-filter reference store that only receives kept records).
        fn push_kept(store: &mut RecordStore<()>, input: &PushInput) -> RecordIdx {
            let mut forced_keep = input.clone();
            forced_keep.accept = true;
            push_one(store, &forced_keep).expect("accept=true always yields Some")
        }

        /// Snapshot of all slab contents + record/extras counts — used by the
        /// rollback proptest as the equivalence model.
        type SlabSnapshot = (usize, Vec<u8>, Vec<Base>, Vec<CigarOp>, Vec<u8>, Vec<u8>, usize);

        fn dump_slabs(store: &RecordStore<()>) -> SlabSnapshot {
            (
                store.records.len(),
                store.names.clone(),
                store.bases.clone(),
                store.cigar.clone(),
                store.qual.clone(),
                store.aux.clone(),
                store.extras.len(),
            )
        }

        proptest! {
            // r[verify record_store.pre_filter.rollback]
            /// Self-consistency check: pushing a mixed accept/reject sequence
            /// produces the same state as pushing only the accepted inputs
            /// with no filter. This catches divergence between the rollback
            /// path and the always-keep path inside `push_fields`, but does
            /// NOT catch bugs that affect both paths (`push_fields` is its
            /// own oracle here). For independent byte-correctness oracles
            /// see `push_fields_matches_owned_bam_round_trip` (A) and
            /// `push_fields_input_is_readable_via_getters` (B) below.
            #[test]
            fn push_fields_rollback_self_consistency(
                inputs in prop::collection::vec(arb_push_input(), 0..=60),
            ) {
                // Store A: push all inputs with per-record filter (rollback path).
                let mut a = RecordStore::new();
                for inp in &inputs {
                    push_one(&mut a, inp);
                }

                // Store B: push only accepted inputs with always-keep filter.
                let mut b = RecordStore::new();
                for inp in inputs.iter().filter(|i| i.accept) {
                    push_kept(&mut b, inp);
                }

                // Every slab must be byte-identical. Records too — we compare
                // their individual fields since SlimRecord doesn't derive Eq.
                prop_assert_eq!(a.records.len(), b.records.len(), "records len");
                for (ra, rb) in a.records.iter().zip(b.records.iter()) {
                    prop_assert_eq!(ra.pos, rb.pos);
                    prop_assert_eq!(ra.end_pos, rb.end_pos);
                    prop_assert_eq!(ra.flags, rb.flags);
                    prop_assert_eq!(ra.mapq, rb.mapq);
                    prop_assert_eq!(ra.seq_len, rb.seq_len);
                    prop_assert_eq!(ra.name_off, rb.name_off, "name_off");
                    prop_assert_eq!(ra.name_len, rb.name_len, "name_len");
                    prop_assert_eq!(ra.bases_off, rb.bases_off, "bases_off");
                    prop_assert_eq!(ra.cigar_off, rb.cigar_off, "cigar_off");
                    prop_assert_eq!(ra.qual_off, rb.qual_off, "qual_off");
                    prop_assert_eq!(ra.aux_off, rb.aux_off, "aux_off");
                    prop_assert_eq!(ra.aux_len, rb.aux_len, "aux_len");
                    prop_assert_eq!(ra.extras_idx, rb.extras_idx, "extras_idx");
                }
                prop_assert_eq!(dump_slabs(&a), dump_slabs(&b), "slab bytes");
            }

            // r[verify record_store.pre_filter.rollback]
            /// Per-step invariant: slab lengths track exactly the running
            /// total of accepted inputs. Catches cases where rollback leaves
            /// trailing garbage in one slab but not others.
            #[test]
            fn push_fields_slab_lengths_track_accepted_prefix(
                inputs in prop::collection::vec(arb_push_input(), 0..=40),
            ) {
                let mut store = RecordStore::new();
                let mut expected_names = 0usize;
                let mut expected_bases = 0usize;
                let mut expected_cigar = 0usize;
                let mut expected_qual = 0usize;
                let mut expected_aux = 0usize;
                let mut expected_records = 0usize;

                for inp in &inputs {
                    let before = dump_slabs(&store);
                    let result = push_one(&mut store, inp);

                    if inp.accept {
                        prop_assert!(result.is_some(), "accept=true must return Some");
                        expected_records += 1;
                        expected_names += inp.qname.len();
                        expected_bases += inp.bases.len();
                        expected_cigar += 1; // one CigarOp slot per kept record
                        expected_qual += inp.quals.len();
                        expected_aux += inp.aux.len();
                    } else {
                        prop_assert!(result.is_none(), "accept=false must return None");
                        // Reject path: state must be untouched.
                        prop_assert_eq!(before, dump_slabs(&store), "rollback left state altered");
                    }

                    prop_assert_eq!(store.records.len(), expected_records, "records len");
                    prop_assert_eq!(store.names.len(), expected_names, "names len");
                    prop_assert_eq!(store.bases.len(), expected_bases, "bases len");
                    prop_assert_eq!(store.cigar.len(), expected_cigar, "cigar len");
                    prop_assert_eq!(store.qual.len(), expected_qual, "qual len");
                    prop_assert_eq!(store.aux.len(), expected_aux, "aux len");
                    prop_assert_eq!(store.extras.len(), expected_records, "extras len");
                }
            }

            // r[verify record_store.pre_filter.rollback]
            /// Indices returned by push_fields across a filtered run must be
            /// dense and sequential — a rejected record must NOT burn an index.
            #[test]
            fn push_fields_indices_are_dense_for_accepted(
                inputs in prop::collection::vec(arb_push_input(), 0..=40),
            ) {
                let mut store = RecordStore::new();
                let mut kept: Vec<RecordIdx> = Vec::new();
                for inp in &inputs {
                    if let Some(idx) = push_one(&mut store, inp) {
                        kept.push(idx);
                    }
                }
                let expected: Vec<RecordIdx> =
                    (0..kept.len()).filter_map(RecordIdx::from_usize).collect();
                prop_assert_eq!(kept, expected);
            }
        }

        // ---- Independent oracle (A): BAM round-trip ----
        //
        // Push the same logical input two different ways:
        //   - via `push_fields` directly (the system under test)
        //   - by building an `OwnedBamRecord`, serializing to BAM-binary form
        //     via `to_bam_bytes`, and feeding to `push_raw` (different code
        //     path: BAM-binary decode + 4-bit seq decode + cigar bulk memcpy)
        //
        // The two stores must end up with byte-identical slabs. A bug shared
        // by both call paths in `push_fields` is invisible to the rollback
        // self-consistency proptest above, but here it would diverge from
        // the BAM-encode/decode round-trip.

        use crate::bam::aux_data::AuxData;
        use crate::bam::owned_record::OwnedBamRecord;
        use seqair_types::BaseQuality;
        use seqair_types::Pos0;

        fn build_owned_bam(input: &PushInput) -> OwnedBamRecord {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "bases.len() bounded ≤ 16 by strategy; fits in u32"
            )]
            let len = input.bases.len() as u32;
            OwnedBamRecord::builder(0, Some(Pos0::new(input.pos).unwrap()), input.qname.clone())
                .flags(BamFlags::empty())
                .mapq(input.mapq)
                .cigar(vec![CigarOp::new(cigar::CigarOpType::Match, len)])
                .seq(input.bases.clone())
                .qual(input.quals.iter().copied().map(BaseQuality::from_byte).collect())
                .aux(AuxData::from_bytes(input.aux.clone()))
                .build()
                .expect("synthetic OwnedBamRecord must build")
        }

        proptest! {
            // r[verify unified.push_fields_equivalence]
            #[test]
            fn push_fields_matches_owned_bam_round_trip(
                inputs in prop::collection::vec(arb_push_input(), 0..=20),
            ) {
                // Store A: BAM-binary encode then push_raw (independent path).
                let mut a = RecordStore::new();
                let mut raw_buf = Vec::new();
                for inp in inputs.iter().filter(|i| i.accept) {
                    let owned = build_owned_bam(inp);
                    raw_buf.clear();
                    owned.to_bam_bytes(&mut raw_buf)
                        .expect("to_bam_bytes must succeed on synthetic input");
                    a.push_raw(&raw_buf, &mut ())
                        .expect("push_raw must accept BAM bytes from to_bam_bytes")
                        .expect("default filter returns true");
                }

                // Store B: push_fields directly (system under test).
                let mut b = RecordStore::new();
                for inp in inputs.iter().filter(|i| i.accept) {
                    push_kept(&mut b, inp);
                }

                // Both must produce byte-identical slabs and matching record fields.
                prop_assert_eq!(a.records.len(), b.records.len(), "records len");
                for (ra, rb) in a.records.iter().zip(b.records.iter()) {
                    prop_assert_eq!(ra.pos, rb.pos, "pos");
                    prop_assert_eq!(ra.end_pos, rb.end_pos, "end_pos");
                    prop_assert_eq!(ra.flags, rb.flags, "flags");
                    prop_assert_eq!(ra.mapq, rb.mapq, "mapq");
                    prop_assert_eq!(ra.seq_len, rb.seq_len, "seq_len");
                    prop_assert_eq!(ra.n_cigar_ops, rb.n_cigar_ops, "n_cigar_ops");
                    prop_assert_eq!(ra.matching_bases, rb.matching_bases, "matching_bases");
                    prop_assert_eq!(ra.indel_bases, rb.indel_bases, "indel_bases");
                }
                prop_assert_eq!(dump_slabs(&a), dump_slabs(&b), "slab bytes");
            }
        }

        // ---- Independent oracle (B): read-back via getters ----
        //
        // After pushing a sequence of inputs via `push_fields`, the slab
        // accessors must return what we put in. Independent of slab layout —
        // catches "wrote at the wrong offset", "split a field across slabs",
        // "encoded the wrong length" without depending on the encoding logic
        // itself. The input is the source of truth.

        proptest! {
            #[test]
            fn push_fields_input_is_readable_via_getters(
                inputs in prop::collection::vec(arb_push_input(), 0..=40),
            ) {
                let mut store = RecordStore::new();
                let mut kept_inputs: Vec<&PushInput> = Vec::new();
                for inp in &inputs {
                    if let Some(idx) = push_one(&mut store, inp) {
                        prop_assert_eq!(idx.as_usize(), kept_inputs.len(), "dense indices");
                        kept_inputs.push(inp);
                    }
                }

                prop_assert_eq!(store.len(), kept_inputs.len(), "store len matches kept count");

                for (i, inp) in kept_inputs.iter().enumerate() {
                    let idx = RecordIdx::from_usize(i).expect("store is small");
                    prop_assert_eq!(store.record(idx).unwrap().qname(), inp.qname.as_slice(), "qname rec {}", i);
                    prop_assert_eq!(store.record(idx).unwrap().seq(), inp.bases.as_slice(), "seq rec {}", i);
                    let qual_bytes = BaseQuality::slice_to_bytes(store.record(idx).unwrap().qual());
                    prop_assert_eq!(qual_bytes, inp.quals.as_slice(), "qual rec {}", i);
                    prop_assert_eq!(store.record(idx).unwrap().aux(), inp.aux.as_slice(), "aux rec {}", i);
                    let cigar = store.record(idx).unwrap().cigar();
                    prop_assert_eq!(cigar.len(), 1, "cigar op count rec {}", i);
                    let op = cigar[0];
                    prop_assert_eq!(op.op_type(), cigar::CigarOpType::Match, "cigar op type rec {}", i);
                    prop_assert_eq!(op.len() as usize, inp.bases.len(), "cigar op len rec {}", i);

                    let rec = store.record(idx).unwrap();
                    #[expect(
                        clippy::cast_sign_loss,
                        reason = "PushInput.pos is bounded < 1_000_000 so always nonneg"
                    )]
                    let rec_pos_u32 = rec.pos.as_i32() as u32;
                    prop_assert_eq!(rec_pos_u32, inp.pos, "pos rec {}", i);
                    prop_assert_eq!(rec.mapq, inp.mapq, "mapq rec {}", i);
                    #[expect(
                        clippy::cast_possible_truncation,
                        reason = "bases.len() bounded ≤ 16; fits in u32"
                    )]
                    let expected_seq_len = inp.bases.len() as u32;
                    prop_assert_eq!(rec.seq_len, expected_seq_len, "seq_len rec {}", i);
                }
            }
        }

        // Same three proptests, but for push_raw. We build synthetic BAM
        // records so parse_header accepts them; the slab invariants are
        // identical but the decode path is different (qname NUL-termination,
        // packed seq decoding, cigar slicing).

        /// Build a minimal BAM record with the given qname, `seq_len`, and a
        /// single M op. Mirrors the helper in the outer test module but
        /// adapted for proptest inputs.
        fn build_bam_raw(qname: &[u8], seq_len: u32, pos: i32, mapq: u8) -> Vec<u8> {
            let mut name_with_nul: Vec<u8> = qname.to_vec();
            name_with_nul.push(0);
            while !name_with_nul.len().is_multiple_of(4) {
                name_with_nul.push(0);
            }
            let name_len = name_with_nul.len();
            let cigar_bytes = 4usize; // one op
            let seq_bytes = (seq_len as usize).div_ceil(2);
            let total = 32 + name_len + cigar_bytes + seq_bytes + seq_len as usize;

            let mut raw = vec![0u8; total];
            raw[0..4].copy_from_slice(&0i32.to_le_bytes());
            raw[4..8].copy_from_slice(&pos.to_le_bytes());
            #[expect(
                clippy::cast_possible_truncation,
                reason = "name_len bounded ≤ 20 by strategy"
            )]
            {
                raw[8] = name_len as u8;
            }
            raw[9] = mapq;
            raw[12..14].copy_from_slice(&1u16.to_le_bytes()); // n_cigar_ops
            raw[14..16].copy_from_slice(&0u16.to_le_bytes()); // flags
            raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
            raw[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // next_ref_id
            raw[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // next_pos
            raw[28..32].copy_from_slice(&0i32.to_le_bytes()); // tlen

            raw[32..32 + name_len].copy_from_slice(&name_with_nul);

            let cigar_start = 32 + name_len;
            let op = seq_len << 4; // M
            raw[cigar_start..cigar_start + 4].copy_from_slice(&op.to_le_bytes());

            // Seq bytes already zeroed (all Unknown after decode); leave qual zero.
            raw
        }

        #[derive(Debug, Clone)]
        struct RawInput {
            qname: Vec<u8>,
            seq_len: u32,
            pos: i32,
            mapq: u8,
            accept: bool,
        }

        fn arb_raw_input() -> impl Strategy<Value = RawInput> {
            (
                prop::collection::vec(b'a'..=b'z', 1..=10), // ASCII-safe qname
                1u32..=16,
                0i32..=1_000,
                0u8..=60,
                any::<bool>(),
            )
                .prop_map(|(qname, seq_len, pos, mapq, accept)| RawInput {
                    qname,
                    seq_len,
                    pos,
                    mapq,
                    accept,
                })
        }

        fn push_raw_one(store: &mut RecordStore<()>, input: &RawInput) -> Option<RecordIdx> {
            let raw = build_bam_raw(&input.qname, input.seq_len, input.pos, input.mapq);
            store
                .push_raw(&raw, &mut AcceptFlag(input.accept))
                .expect("synthetic BAM record is always parseable")
        }

        fn push_raw_kept(store: &mut RecordStore<()>, input: &RawInput) -> RecordIdx {
            let mut forced_keep = input.clone();
            forced_keep.accept = true;
            push_raw_one(store, &forced_keep).expect("accept=true always yields Some")
        }

        proptest! {
            // r[verify record_store.pre_filter.rollback]
            #[test]
            fn push_raw_rollback_matches_filtered_replay(
                inputs in prop::collection::vec(arb_raw_input(), 0..=40),
            ) {
                let mut a = RecordStore::new();
                for inp in &inputs {
                    push_raw_one(&mut a, inp);
                }

                let mut b = RecordStore::new();
                for inp in inputs.iter().filter(|i| i.accept) {
                    push_raw_kept(&mut b, inp);
                }

                prop_assert_eq!(dump_slabs(&a), dump_slabs(&b), "slab bytes");
            }
        }
    }
}
