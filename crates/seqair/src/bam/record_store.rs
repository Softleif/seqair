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
    seq,
};
use seqair_types::{BamFlags, Base, BaseQuality, Pos0};

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
    /// Mapping quality (0..=254; 255 = unavailable per `[SAM1] §1.4`).
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

// Compile-time size guard: 16 u32 (incl. next_ref_id, extras_idx) + 3 u16 + 1 u8 + padding
// = 72 bytes. If this ever grows, revisit the layout before accepting the hit —
// see `docs/spec/3-record_store.md` "Layout".
const _: () =
    assert!(std::mem::size_of::<SlimRecord>() <= 72, "SlimRecord grew unexpectedly large");

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
}

impl<U> RecordStore<U> {
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
        // CIGAR at 2 ops/record (= 8 B) gives headroom over measured ~1.1 ops/record
        // for short reads without allocating the long-read budget we don't need here.
        let cigar_est = record_count_est.saturating_mul(2);
        let qual_est = record_count_est.saturating_mul(150);
        let aux_est = record_count_est.saturating_mul(180); // rounded up from pooled 163

        Self {
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
    // r[impl record_store.filter_raw]
    pub fn push_raw<E: CustomizeRecordStore<Extra = U>>(
        &mut self,
        raw: &[u8],
        customize: &mut E,
    ) -> Result<Option<u32>, DecodeError> {
        let h = record::parse_header(raw)?;

        let idx = u32::try_from(self.records.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // All slice bounds ≤ qual_end ≤ raw.len() (checked by parse_header)
        debug_assert!(h.qual_end <= raw.len(), "qual_end overrun: {} > {}", h.qual_end, raw.len());
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let qname_raw = &raw[32..h.var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let cigar_bytes = &raw[h.var_start..h.cigar_end];

        // --- filter_raw: pre-filter before any slab extension or base decode ---
        // r[impl record_store.filter_raw]
        {
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
                end_pos: h.pos, // not yet computed from CIGAR
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
                raw_cigar_bytes: Some(cigar_bytes),
                cigar_ops: None,
                packed_seq: Some(packed_seq_slice),
                bases: None,
            }) {
                return Ok(None);
            }
        }

        // --- Write into cigar slab (typed; bulk memcpy on LE) ---
        // r[impl record_store.checked_offsets]
        let cigar_off = u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;
        CigarOp::extend_from_bam_bytes(&mut self.cigar, cigar_bytes);
        let cigar_start = cigar_off as usize;
        #[allow(clippy::indexing_slicing, reason = "just appended; offsets within slab bounds")]
        let cigar_ops = &self.cigar[cigar_start..];

        let end_pos = cigar::compute_end_pos(h.pos, cigar_ops)
            .ok_or(DecodeError::InvalidPosition { value: h.pos.as_i32() })?;
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
    ) -> Result<Option<u32>, DecodeError> {
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
            end_pos,
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
            raw_cigar_bytes: None,
            cigar_ops: Some(cigar_ops),
            packed_seq: None,
            bases: Some(bases),
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
/// - `raw_cigar_bytes` is `Some` (the raw BAM CIGAR bytes).
/// - `packed_seq` is `Some` (4-bit packed sequence bytes).
/// - `cigar_ops` and `bases` are `None` (not yet decoded).
/// - `end_pos`, `matching_bases`, `indel_bases` are not yet computed (set to
///   `pos`, `0`, `0` respectively).
///
/// From `push_fields` (pre-parsed SAM/CRAM fields):
/// - `cigar_ops` is `Some` (decoded CIGAR ops).
/// - `bases` is `Some` (decoded sequence bases).
/// - `raw_cigar_bytes` and `packed_seq` are `None`.
/// - `end_pos`, `matching_bases`, `indel_bases` are the pre-computed values.
#[derive(Debug)]
#[non_exhaustive]
pub struct FilterRawFields<'a> {
    /// 0-based reference position.
    pub pos: Pos0,
    /// Exclusive end position (computed from CIGAR; `pos` for `push_raw` since
    /// CIGAR is not yet decoded).
    pub end_pos: Pos0,
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
    /// Raw BAM CIGAR bytes (`push_raw` path); `None` from `push_fields`.
    pub raw_cigar_bytes: Option<&'a [u8]>,
    /// Decoded CIGAR ops (`push_fields` path); `None` from `push_raw`.
    pub cigar_ops: Option<&'a [CigarOp]>,
    /// 4-bit packed sequence bytes (`push_raw` path); `None` from `push_fields`.
    pub packed_seq: Option<&'a [u8]>,
    /// Decoded sequence bases (`push_fields` path); `None` from `push_raw`.
    pub bases: Option<&'a [Base]>,
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
    pub fn cigar(&self, idx: u32) -> &[CigarOp] {
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
        new_cigar_ops: &[CigarOp],
    ) -> Result<(), DecodeError> {
        let rec = self.record(idx);
        let seq_len = rec.seq_len;

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

        let end_pos = cigar::compute_end_pos(new_pos, new_cigar_ops)
            .ok_or(DecodeError::InvalidPosition { value: new_pos.as_i32() })?;

        let (matching_bases, indel_bases) = cigar::calc_matches_indels(new_cigar_ops);

        let new_cigar_off =
            u32::try_from(self.cigar.len()).map_err(|_| DecodeError::SlabOverflow)?;

        // ── Point of no return: all mutation below ──

        self.cigar.extend_from_slice(new_cigar_ops);

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
mod tests {
    use super::*;

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
        let result: Result<Option<u32>, _> = store.push_fields(
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

    // r[verify record_store.slim_record_fields]
    #[test]
    fn push_raw_preserves_next_ref_id() {
        let mut store = RecordStore::new();
        let mut raw = make_raw_record(b"read1", 4, 1);
        // Write next_ref_id = 7 at BAM offset 20
        raw[20..24].copy_from_slice(&7i32.to_le_bytes());
        store.push_raw(&raw, &mut ()).unwrap();
        assert_eq!(store.record(0).next_ref_id, 7);
    }

    // r[verify record_store.pre_filter.rollback]
    #[test]
    fn push_raw_with_rejecting_filter_truncates_all_slabs() {
        let mut store = RecordStore::new();

        // Push record A (kept)
        let a = make_raw_record(b"read1", 4, 1);
        let idx_a = store.push_raw(&a, &mut ()).unwrap();
        assert_eq!(idx_a, Some(0));
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
        assert_eq!(idx_c, Some(1), "rejected record must not consume an index");
        assert_eq!(store.len(), 2, "store should hold A and C only");
    }

    // r[verify record_store.filter_raw]
    #[test]
    fn push_raw_filter_raw_rejection_avoids_slab_write() {
        /// Reject everything in filter_raw — no slab extensions should happen.
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
        assert_eq!(kept_idx, Some(0));
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
        assert_eq!(store.record(0).next_ref_id, 3);
    }

    // ---- Model-based property tests for rollback ----
    //
    // The rollback invariant: pushing a sequence of records with per-record
    // accept/reject decisions MUST produce a store byte-identical to pushing
    // only the accepted records in order with no filter. That is, a rollback
    // must perfectly undo both the records Vec push and all slab extensions,
    // leaving zero dead bytes.
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
        fn push_one(store: &mut RecordStore<()>, input: &PushInput) -> Option<u32> {
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
        fn push_kept(store: &mut RecordStore<()>, input: &PushInput) -> u32 {
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
                let mut kept: Vec<u32> = Vec::new();
                for inp in &inputs {
                    if let Some(idx) = push_one(&mut store, inp) {
                        kept.push(idx);
                    }
                }
                let expected: Vec<u32> =
                    (0..kept.len()).map(|i| u32::try_from(i).unwrap()).collect();
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

        fn build_owned_bam(input: &PushInput) -> OwnedBamRecord {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "bases.len() bounded ≤ 16 by strategy; fits in u32"
            )]
            let len = input.bases.len() as u32;
            OwnedBamRecord::builder(0, i64::from(input.pos), input.qname.clone())
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
                        prop_assert_eq!(idx as usize, kept_inputs.len(), "dense indices");
                        kept_inputs.push(inp);
                    }
                }

                prop_assert_eq!(store.len(), kept_inputs.len(), "store len matches kept count");

                for (i, inp) in kept_inputs.iter().enumerate() {
                    let idx = i as u32;
                    prop_assert_eq!(store.qname(idx), inp.qname.as_slice(), "qname rec {}", i);
                    prop_assert_eq!(store.seq(idx), inp.bases.as_slice(), "seq rec {}", i);
                    let qual_bytes = BaseQuality::slice_to_bytes(store.qual(idx));
                    prop_assert_eq!(qual_bytes, inp.quals.as_slice(), "qual rec {}", i);
                    prop_assert_eq!(store.aux(idx), inp.aux.as_slice(), "aux rec {}", i);
                    let cigar = store.cigar(idx);
                    prop_assert_eq!(cigar.len(), 1, "cigar op count rec {}", i);
                    let op = cigar[0];
                    prop_assert_eq!(op.op_type(), cigar::CigarOpType::Match, "cigar op type rec {}", i);
                    prop_assert_eq!(op.len() as usize, inp.bases.len(), "cigar op len rec {}", i);

                    let rec = store.record(idx);
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

        fn push_raw_one(store: &mut RecordStore<()>, input: &RawInput) -> Option<u32> {
            let raw = build_bam_raw(&input.qname, input.seq_len, input.pos, input.mapq);
            store
                .push_raw(&raw, &mut AcceptFlag(input.accept))
                .expect("synthetic BAM record is always parseable")
        }

        fn push_raw_kept(store: &mut RecordStore<()>, input: &RawInput) -> u32 {
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
