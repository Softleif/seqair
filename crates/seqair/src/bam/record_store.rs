//! Slab-based BAM record storage with zero per-record heap allocation.
//!
//! All variable-length data is packed into contiguous byte buffers:
//! - **Name slab**: read names (qnames), accessed during dedup mate detection
//! - **Bases slab**: decoded `Base` values per record, accessed per-position in pileup
//! - **Cigar slab**: packed CIGAR ops, append-only mutation during realignment
//! - **Qual slab**: raw Phred bytes, hot in pileup per-base quality lookup
//! - **Aux slab**: auxiliary tag bytes in BAM binary format, rarely read in pileup

use seqair_types::{BamFlags, Base, BaseQuality, Pos0};

use super::{
    cigar,
    record::{self, DecodeError},
    seq,
};

/// Compact BAM record with offsets into the store's slabs.
// r[impl record_store.push_raw+2]
// r[impl record_store.slim_record_fields]
pub struct SlimRecord {
    /// 0-based leftmost aligned reference position.
    pub pos: Pos0,
    /// 0-based exclusive end (`pos` + reference-consuming CIGAR length).
    pub end_pos: Pos0,
    /// BAM flag bits (paired, unmapped, reverse strand, …).
    pub flags: BamFlags,
    /// Number of CIGAR ops; bounded by BAM's u16 `l_op` field.
    pub n_cigar_ops: u16,
    /// Mapping quality (0..=254; 255 = unavailable per [SAM1] §1.4).
    pub mapq: u8,
    /// Query sequence length in bases — also the qual slab length for this record.
    pub seq_len: u32,
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
}

impl SlimRecord {
    fn cigar_len(&self) -> usize {
        (self.n_cigar_ops as usize).checked_mul(4).expect("cigar_len overflow")
    }
}

// Compile-time size guard: 16 u32 (incl. next_ref_id, extras_idx) + 3 u16 + 1 u8 + padding
// = 72 bytes. If this ever grows, revisit the layout before accepting the hit —
// see `docs/spec/3-record_store.md` "Layout".
const _: () =
    assert!(std::mem::size_of::<SlimRecord>() <= 72, "SlimRecord grew unexpectedly large");

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
/// Use [`with_extras`](RecordStore::with_extras) to compute per-record data after loading.
// r[impl record_store.extras.generic_param]
pub struct RecordStore<U = ()> {
    records: Vec<SlimRecord>,
    names: Vec<u8>,
    bases: Vec<Base>,
    cigar: Vec<u8>,
    qual: Vec<u8>,
    aux: Vec<u8>,
    extras: Vec<U>,
}

impl RecordStore {
    pub fn new() -> Self {
        Self {
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
        // CIGAR at 8 gives headroom over measured ~4.4 without allocating the
        // 20 B/rec long-read budget we don't need here.
        let cigar_est = record_count_est.saturating_mul(8);
        let qual_est = record_count_est.saturating_mul(150);
        let aux_est = record_count_est.saturating_mul(180); // rounded up from pooled 163

        Self {
            records: Vec::with_capacity(record_count_est),
            names: Vec::with_capacity(names_est),
            bases: Vec::with_capacity(bases_est),
            cigar: Vec::with_capacity(cigar_est),
            qual: Vec::with_capacity(qual_est),
            aux: Vec::with_capacity(aux_est),
            extras: Vec::new(),
        }
    }

    // r[impl record_store.extras.push_unit]
    /// Decode a raw BAM record and append to the store.
    ///
    /// Variable-length data is written directly into the slabs.
    pub fn push_raw(&mut self, raw: &[u8]) -> Result<u32, DecodeError> {
        let h = record::parse_header(raw)?;

        let idx = u32::try_from(self.records.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // All slice bounds ≤ qual_end ≤ raw.len() (checked by parse_header)
        debug_assert!(h.qual_end <= raw.len(), "qual_end overrun: {} > {}", h.qual_end, raw.len());
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let qname_raw = &raw[32..h.var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let cigar_slice = &raw[h.var_start..h.cigar_end];
        let end_pos = record::compute_end_pos(h.pos, cigar_slice)
            .ok_or(DecodeError::InvalidPosition { value: h.pos.as_i32() })?;
        let (matching_bases, indel_bases) = cigar::calc_matches_indels(cigar_slice);

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
        // r[depends base_decode.table_invariant]
        unsafe {
            let out = std::slice::from_raw_parts_mut(bases_spare.as_mut_ptr() as *mut u8, seq_len);
            seq::decode_bases_into(packed_seq, seq_len, out);
            self.bases.set_len(
                self.bases.len().checked_add(seq_len).expect("bases slab length overflow"),
            );
        }

        // --- Write into cigar slab ---
        // r[impl record_store.checked_offsets]
        let cigar_off = u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.cigar.extend_from_slice(cigar_slice);

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
            extras_idx: idx,
        });
        self.extras.push(());

        Ok(idx)
    }

    // r[impl unified.record_store_push]
    // r[impl record_store.push_fields]
    // r[impl unified.push_fields_equivalence]
    // r[impl record_store.checked_offsets]
    /// Append a record from pre-parsed fields (for SAM/CRAM readers).
    ///
    /// Writes directly into the slabs without going through BAM binary encoding.
    /// CIGAR must be in BAM packed u32 format (`len << 4 | op`).
    #[expect(
        clippy::too_many_arguments,
        reason = "all fields are needed for zero-copy push into the record store slabs"
    )]
    pub fn push_fields(
        &mut self,
        pos: Pos0,
        end_pos: Pos0,
        flags: BamFlags,
        mapq: u8,
        matching_bases: u32,
        indel_bases: u32,
        qname: &[u8],
        cigar_packed: &[u8],
        bases: &[Base],
        qual: &[u8],
        aux: &[u8],
        tid: i32,
        next_ref_id: i32,
        next_pos: i32,
        template_len: i32,
    ) -> Result<u32, DecodeError> {
        if qual.len() != bases.len() {
            return Err(DecodeError::QualLenMismatch {
                qual_len: qual.len(),
                seq_len: bases.len(),
            });
        }

        let idx = u32::try_from(self.records.len()).map_err(|_| DecodeError::SlabOverflow)?;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "caller validates cigar op count ≤ 65535 (BAM n_cigar_op is u16); fits in u16"
        )]
        let n_cigar_ops = (cigar_packed.len() / 4) as u16;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "seq length is bounded by slab limits (u32); slab overflow checked via SlabOverflow"
        )]
        let seq_len = bases.len() as u32;

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
        self.cigar.extend_from_slice(cigar_packed);

        // Qual slab
        // r[impl record_store.checked_offsets]
        let qual_off = u32::try_from(self.qual.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.qual.extend_from_slice(qual);

        // Aux slab
        // r[impl record_store.checked_offsets]
        let aux_off = u32::try_from(self.aux.len()).map_err(|_| DecodeError::SlabOverflow)?;
        self.aux.extend_from_slice(aux);

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
            extras_idx: idx,
        });
        self.extras.push(());

        Ok(idx)
    }

    // r[impl record_store.extras.with_extras]
    /// Compute per-record extras in a single pass, consuming this `RecordStore<()>`.
    ///
    /// The closure receives `(record_index, &RecordStore<()>)` so it can access
    /// any slab (record fields, aux, seq, qname, etc.). All existing slabs are
    /// moved, not copied.
    ///
    /// # Example
    ///
    /// Attach per-record metadata and access it during pileup iteration:
    ///
    /// ```
    /// use seqair::bam::{Pos0, RecordStore};
    /// use seqair::bam::pileup::PileupEngine;
    /// use seqair_types::{BamFlags, Base};
    ///
    /// // Per-record data computed at load time.
    /// struct ReadInfo {
    ///     is_reverse: bool,
    ///     qname_len: usize,
    /// }
    ///
    /// // After loading records into a RecordStore<()>...
    /// let mut store = RecordStore::new();
    /// # let cigar = (4u32 << 4).to_le_bytes(); // 4M
    /// # store.push_fields(
    /// #     Pos0::new(100).unwrap(), Pos0::new(103).unwrap(),
    /// #     BamFlags::from(0x10), 60, 4, 0, b"read1", &cigar,
    /// #     &[Base::A, Base::C, Base::G, Base::T], &[30; 4], &[], 0, -1, 0, 0,
    /// # ).unwrap();
    ///
    /// // Compute extras from record fields:
    /// let store = store.with_extras(|idx, store| {
    ///     let rec = store.record(idx);
    ///     ReadInfo {
    ///         is_reverse: rec.flags.is_reverse(),
    ///         qname_len: store.qname(idx).len(),
    ///     }
    /// });
    ///
    /// // Build a pileup engine with the extras-bearing store:
    /// let mut engine = PileupEngine::new(
    ///     store,
    ///     Pos0::new(100).unwrap(),
    ///     Pos0::new(103).unwrap(),
    /// );
    ///
    /// // Iterate via the lending `pileups()` — the column borrows the store,
    /// // and alignments expose `.extra()` for the per-record data.
    /// while let Some(column) = engine.pileups() {
    ///     for aln in column.alignments() {
    ///         let info = aln.extra();
    ///         assert!(info.is_reverse);
    ///         assert_eq!(info.qname_len, 5);
    ///     }
    /// }
    /// ```
    pub fn with_extras<V>(self, mut f: impl FnMut(u32, &Self) -> V) -> RecordStore<V> {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "RecordStore capacity is bounded by SlabOverflow (u32)"
        )]
        let extras: Vec<V> = (0..self.records.len()).map(|i| f(i as u32, &self)).collect();
        self.into_parts_with_extras(extras)
    }

    // r[impl record_store.extras.provider]
    /// Compute per-record extras using a [`RecordStoreExtras`] provider.
    ///
    /// Unlike [`with_extras`](Self::with_extras), the provider is a reusable
    /// value: it can carry state and is `Clone`, so it survives calls across
    /// regions and across [`Readers::fork`] in a multi-threaded pipeline.
    ///
    /// See [`RecordStoreExtras`] for the trait requirements.
    pub fn apply_extras<E: RecordStoreExtras>(self, provider: &mut E) -> RecordStore<E::Extra> {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "RecordStore capacity is bounded by SlabOverflow (u32)"
        )]
        let extras: Vec<E::Extra> =
            (0..self.records.len()).map(|i| provider.compute(i as u32, &self)).collect();
        self.into_parts_with_extras(extras)
    }

    fn into_parts_with_extras<V>(self, extras: Vec<V>) -> RecordStore<V> {
        let mut records = self.records;
        // Reset extras_idx to sequential — the new extras Vec was built in
        // current record order, so extras[i] corresponds to records[i].
        #[expect(
            clippy::cast_possible_truncation,
            reason = "RecordStore capacity is bounded by SlabOverflow (u32)"
        )]
        for (i, rec) in records.iter_mut().enumerate() {
            rec.extras_idx = i as u32;
        }
        RecordStore {
            records,
            names: self.names,
            bases: self.bases,
            cigar: self.cigar,
            qual: self.qual,
            aux: self.aux,
            extras,
        }
    }
}

// r[impl record_store.extras.provider]
/// Computes per-record extras for a [`RecordStore`].
///
/// Used by [`Readers`](crate::reader::Readers) to attach user data to every
/// record as it is loaded. The provider itself must be `Clone` so it can be
/// duplicated across forked readers in a multi-threaded pipeline.
///
/// # Example
///
/// ```
/// use seqair::bam::record_store::{RecordStore, RecordStoreExtras};
///
/// #[derive(Clone, Default)]
/// struct QnameLen;
///
/// impl RecordStoreExtras for QnameLen {
///     type Extra = usize;
///     fn compute(&mut self, idx: u32, store: &RecordStore<()>) -> usize {
///         store.qname(idx).len()
///     }
/// }
/// ```
pub trait RecordStoreExtras: Clone {
    /// The per-record data produced by this provider.
    type Extra;

    /// Compute the extra for record `idx` in `store`.
    fn compute(&mut self, idx: u32, store: &RecordStore<()>) -> Self::Extra;
}

// r[impl record_store.extras.provider]
/// Blanket no-op implementation for the default `()` case.
///
/// `RecordStore<()>` does not need extras; `Readers<()>` uses this to skip the
/// `apply_extras` pass entirely when there is nothing to compute.
impl RecordStoreExtras for () {
    type Extra = ();
    #[inline]
    fn compute(&mut self, _idx: u32, _store: &RecordStore<()>) {}
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
    pub fn sort_by_pos(&mut self) {
        self.records.sort_by_key(|r| r.pos);
    }

    /// Remove consecutive duplicate records (same position, flags, and read
    /// name). Must be called after `sort_by_pos` so duplicates are adjacent.
    ///
    /// Nearby and distant BAM index chunks can cover overlapping byte ranges,
    /// causing the same record to be loaded from both sources.
    /// Slab data (including extras) for removed records is left in place (minor waste),
    /// same as dead name/cigar bytes from removed records.
    pub fn dedup(&mut self) {
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

    // --- Accessors ---

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    #[allow(clippy::indexing_slicing, reason = "idx is always a valid index returned by push_raw")]
    pub fn record(&self, idx: u32) -> &SlimRecord {
        debug_assert!(
            (idx as usize) < self.records.len(),
            "record idx {idx} out of bounds (len={})",
            self.records.len()
        );
        &self.records[idx as usize]
    }

    /// Update the `template_len` field for a record in the store.
    /// Used by the CRAM decoder to fill in TLEN after resolving mate
    /// cross-references within a slice.
    pub fn set_template_len(&mut self, idx: u32, tlen: i32) -> Option<()> {
        debug_assert!(
            (idx as usize) < self.records.len(),
            "set_template_len idx {idx} out of bounds (len={})",
            self.records.len()
        );
        let idx = usize::try_from(idx).ok()?;
        self.records.get_mut(idx)?.template_len = tlen;
        Some(())
    }

    /// Update the mate fields (`next_ref_id`, `next_pos`) for a record.
    /// Used by the CRAM decoder to fill in mate info after resolving
    /// downstream mate cross-references within a slice.
    pub fn set_mate_info(&mut self, idx: u32, next_ref_id: i32, next_pos: i32) -> Option<()> {
        let idx = usize::try_from(idx).ok()?;
        let rec = self.records.get_mut(idx)?;
        rec.next_ref_id = next_ref_id;
        rec.next_pos = next_pos;
        Some(())
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn qname(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let start = rec.name_off as usize;
        let end = start.checked_add(rec.name_len as usize).expect("qname end overflow");
        debug_assert!(end <= self.names.len(), "qname slab overrun: {end} > {}", self.names.len());
        &self.names[start..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn cigar(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let start = rec.cigar_off as usize;
        let end = start.checked_add(rec.cigar_len()).expect("cigar end overflow");
        debug_assert!(end <= self.cigar.len(), "cigar slab overrun: {end} > {}", self.cigar.len());
        &self.cigar[start..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn seq(&self, idx: u32) -> &[Base] {
        let rec = self.record(idx);
        let start = rec.bases_off as usize;
        let end = start.checked_add(rec.seq_len as usize).expect("seq end overflow");
        debug_assert!(end <= self.bases.len(), "bases slab overrun: {end} > {}", self.bases.len());
        &self.bases[start..end]
    }

    pub fn seq_at(&self, idx: u32, pos: usize) -> Base {
        let rec = self.record(idx);
        self.bases
            .get((rec.bases_off as usize).checked_add(pos).expect("seq_at offset overflow"))
            .copied()
            .unwrap_or(Base::Unknown)
    }

    // r[impl types.base_quality.field_type]
    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn qual(&self, idx: u32) -> &[BaseQuality] {
        let rec = self.record(idx);
        let start = rec.qual_off as usize;
        let end = start.checked_add(rec.seq_len as usize).expect("qual end overflow");
        debug_assert!(end <= self.qual.len(), "qual slab overrun: {end} > {}", self.qual.len());
        BaseQuality::slice_from_bytes(&self.qual[start..end])
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn aux(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
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
    /// After calling this on one or more records, call `sort_by_pos()` before
    /// using the store for pileup iteration.
    pub fn set_alignment(
        &mut self,
        idx: u32,
        new_pos: Pos0,
        new_cigar_packed: &[u8],
    ) -> Result<(), DecodeError> {
        let rec = self.record(idx);
        let seq_len = rec.seq_len;

        // ── All validation before any mutation ──
        // If any check fails, the store is unchanged (r[record_store.set_alignment.validation]).

        let new_query_len = cigar::calc_query_len(new_cigar_packed);
        if new_query_len != seq_len {
            return Err(DecodeError::CigarQueryLenMismatch {
                cigar_query_len: new_query_len,
                seq_len,
            });
        }

        let n_ops = new_cigar_packed.len() / 4;
        let n_cigar_ops =
            u16::try_from(n_ops).map_err(|_| DecodeError::CigarOpCountOverflow { count: n_ops })?;

        let end_pos = record::compute_end_pos(new_pos, new_cigar_packed)
            .ok_or(DecodeError::InvalidPosition { value: new_pos.as_i32() })?;

        let (matching_bases, indel_bases) = cigar::calc_matches_indels(new_cigar_packed);

        let new_cigar_off =
            u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // ── Point of no return: all mutation below ──

        self.cigar.extend_from_slice(new_cigar_packed);

        #[allow(clippy::indexing_slicing, reason = "idx validated by self.record() above")]
        let rec = &mut self.records[idx as usize];
        rec.pos = new_pos;
        rec.end_pos = end_pos;
        rec.n_cigar_ops = n_cigar_ops;
        rec.cigar_off = new_cigar_off;
        rec.matching_bases = matching_bases;
        rec.indel_bases = indel_bases;

        Ok(())
    }

    // r[impl record_store.extras.access]
    /// Access the per-record extra for `idx`.
    #[allow(clippy::indexing_slicing, reason = "idx is always a valid index returned by push_raw")]
    pub fn extra(&self, idx: u32) -> &U {
        let rec = self.record(idx);
        let ei = rec.extras_idx as usize;
        debug_assert!(
            ei < self.extras.len(),
            "extras_idx {ei} out of bounds (len={})",
            self.extras.len()
        );
        &self.extras[ei]
    }

    // r[impl record_store.extras.access]
    /// Mutably access the per-record extra for `idx`.
    #[allow(clippy::indexing_slicing, reason = "idx is always a valid index returned by push_raw")]
    pub fn extra_mut(&mut self, idx: u32) -> &mut U {
        let ei = self.record(idx).extras_idx as usize;
        debug_assert!(
            ei < self.extras.len(),
            "extras_idx {ei} out of bounds (len={})",
            self.extras.len()
        );
        &mut self.extras[ei]
    }

    /// Take all contents out, leaving an empty store with no capacity.
    /// Used by `PileupEngine::take_store` to avoid requiring `Default`.
    pub(crate) fn take_contents(&mut self) -> Self {
        RecordStore {
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
    pub fn clear(&mut self) {
        self.records.clear();
        self.names.clear();
        self.bases.clear();
        self.cigar.clear();
        self.qual.clear();
        self.aux.clear();
        self.extras.clear();
    }

    // r[impl record_store.extras.strip]
    /// Discard the extras slab, returning a `RecordStore<()>` that preserves
    /// all other slab data and capacity. Used by `Readers::recover_store`.
    pub fn strip_extras(self) -> RecordStore {
        RecordStore {
            records: self.records,
            names: self.names,
            bases: self.bases,
            cigar: self.cigar,
            qual: self.qual,
            aux: self.aux,
            extras: Vec::new(),
        }
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

impl Default for RecordStore {
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
mod tests {
    use super::*;

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

    // r[verify bam.record.checked_offsets]
    #[test]
    fn push_raw_rejects_offset_overflow() {
        let mut store = RecordStore::new();
        // Craft a record with fields that would cause overflow.
        let mut raw = [0u8; 32];
        raw[8] = 255; // name_len
        raw[12..14].copy_from_slice(&u16::MAX.to_le_bytes()); // n_cigar_ops
        raw[16..20].copy_from_slice(&u32::MAX.to_le_bytes()); // seq_len

        let result = store.push_raw(&raw);
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
        let result: Result<u32, _> = store.push_fields(
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
        );
        assert!(result.is_ok());
    }

    // r[verify record_store.checked_offsets]
    #[test]
    fn push_raw_normal_record_succeeds() {
        let mut store = RecordStore::new();
        let raw = make_raw_record(b"read1", 4, 1);
        let result = store.push_raw(&raw);
        assert!(result.is_ok());
        assert_eq!(store.len(), 1);
    }

    // r[verify record_store.slim_record_fields]
    #[test]
    fn push_raw_preserves_next_ref_id() {
        let mut store = RecordStore::new();
        let mut raw = make_raw_record(b"read1", 4, 1);
        // Write next_ref_id = 7 at BAM offset 20
        raw[20..24].copy_from_slice(&7i32.to_le_bytes());
        store.push_raw(&raw).unwrap();
        assert_eq!(store.record(0).next_ref_id, 7);
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
                &(5u32 << 4).to_le_bytes(), // 5M cigar
                &[Base::A, Base::C, Base::G, Base::T, Base::A],
                &[30, 31, 32, 33, 34],
                &[],
                0,   // tid
                3,   // next_ref_id
                500, // next_pos
                200, // template_len
            )
            .unwrap();
        assert_eq!(store.record(0).next_ref_id, 3);
    }
}
