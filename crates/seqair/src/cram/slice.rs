//! Decode CRAM slices into records. Reads data series blocks, applies the encodings from
//! [`CompressionHeader`], and pushes decoded records into a [`RecordStore`].

use std::ops::Neg;

use super::{
    block::{self, Block, ContentType},
    compression_header::CompressionHeader,
    encoding::{DecodeContext, ExternalCursor},
    reader::CramError,
    varint,
};
use crate::bam::{
    BamHeader,
    cigar::CigarOp,
    record_store::{CustomizeRecordStore, RecordStore},
};
use seqair_types::{BamFlags, Base, Pos0, Pos1, SmallVec, SmolStr};
use tracing::warn;

/// Mate cross-reference info collected during CRAM record decoding.
/// Used to reconstruct `template_len` for attached/downstream mates.
struct SliceMateInfo {
    /// 0-based alignment start position.
    pos: i64,
    /// 0-based exclusive alignment end position.
    end_pos: i64,
    /// Absolute index of the mate record within this slice, or -1 if detached/no mate.
    mate_line: i32,
    /// Index in the `RecordStore` if this record was pushed, or `None` if filtered out.
    store_idx: Option<u32>,
    /// BAM flags (needed for READ1/READ2 tie-breaking).
    bam_flags: u16,
    /// Reference sequence ID for this record (needed for `next_ref_id` resolution).
    ref_id: i32,
}

/// Parsed CRAM slice header.
#[derive(Debug)]
pub struct SliceHeader {
    pub ref_seq_id: i32,
    pub alignment_start: i32,
    pub alignment_span: i32,
    pub num_records: i32,
    pub record_counter: i64,
    pub num_blocks: i32,
    pub block_content_ids: Vec<i32>,
    pub embedded_reference: i32,
    pub reference_md5: [u8; 16],
}

impl SliceHeader {
    // r[impl cram.slice.header]
    pub fn parse(data: &[u8]) -> Result<Self, CramError> {
        let mut cursor: &[u8] = data;

        let ref_seq_id = read_itf8(&mut cursor)?.cast_signed();
        let alignment_start = read_itf8(&mut cursor)?.cast_signed();
        let alignment_span = read_itf8(&mut cursor)?.cast_signed();
        let num_records = read_itf8(&mut cursor)?.cast_signed();
        let record_counter = read_ltf8(&mut cursor)?.cast_signed();
        let num_blocks = read_itf8(&mut cursor)?.cast_signed();

        // r[impl cram.slice.validated_lengths]
        let num_content_ids_raw = read_itf8(&mut cursor)?;
        let num_content_ids = usize::try_from(num_content_ids_raw.cast_signed())
            .map_err(|_| CramError::InvalidLength { value: num_content_ids_raw.cast_signed() })?;
        super::reader::check_alloc_size(num_content_ids.saturating_mul(4), "slice content IDs")?;
        let mut block_content_ids = Vec::with_capacity(num_content_ids);
        for _ in 0..num_content_ids {
            block_content_ids.push(read_itf8(&mut cursor)?.cast_signed());
        }

        let embedded_reference = read_itf8(&mut cursor)?.cast_signed();

        let md5_bytes =
            cursor.get(..16).ok_or(CramError::Truncated { context: "slice reference MD5" })?;
        let mut reference_md5 = [0u8; 16];
        reference_md5.copy_from_slice(md5_bytes);

        Ok(SliceHeader {
            ref_seq_id,
            alignment_start,
            alignment_span,
            num_records,
            record_counter,
            num_blocks,
            block_content_ids,
            embedded_reference,
            reference_md5,
        })
    }
}

/// Decode all records from a slice's blocks and push them into a `RecordStore`.
///
/// `container_data` starts at the first block after the container header.
/// `slice_offset` is the byte offset from container data start to this slice's header block.
///
/// `customize` is a push-time customizer: each record that passes the
/// reader's built-in overlap/tid/unmapped checks is pushed and
/// `customize.filter` is consulted. When it returns `false`, the push
/// is rolled back with zero slab waste (same as BAM/SAM). Returns
/// `(fetched, kept)` where `fetched` counts records that reached the push
/// step and `kept` counts those that survived the filter.
// r[impl cram.fetch_into_customized.push_time]
#[expect(
    clippy::too_many_arguments,
    reason = "CRAM slice decoding requires all compression header, data, and offset parameters"
)]
pub(crate) fn decode_slice<E: CustomizeRecordStore>(
    ch: &CompressionHeader,
    container_data: &[u8],
    slice_offset: usize,
    reference_seq: &[u8],
    ref_start_0based: i64,
    header: &BamHeader,
    read_group_ids: &[SmolStr],
    tid: u32,
    query_start: Pos0,
    query_end: Pos0,
    store: &mut RecordStore<E::Extra>,
    cigar_buf: &mut Vec<CigarOp>,
    bases_buf: &mut Vec<Base>,
    qual_buf: &mut Vec<u8>,
    aux_buf: &mut Vec<u8>,
    name_buf: &mut Vec<u8>,
    feature_byte_buf: &mut Vec<u8>,
    cigar_ops_buf: &mut Vec<(u32, u8)>,
    customize: &mut E,
    rans_4x8_buf: &mut Option<super::rans::Rans4x8Buf>,
    nx16_order1_buf: &mut Option<super::rans_nx16::Nx16Order1Buf>,
) -> Result<(usize, usize), CramError> {
    let slice_data = container_data
        .get(slice_offset..)
        .ok_or(CramError::Truncated { context: "slice offset" })?;

    // Parse slice header block
    let (slice_header_block, mut pos) =
        block::parse_block_with_buf(slice_data, rans_4x8_buf.as_mut(), nx16_order1_buf.as_mut())?;
    if slice_header_block.content_type != ContentType::SliceHeader {
        return Err(CramError::ExpectedSliceHeader { found: slice_header_block.content_type });
    }
    let sh = SliceHeader::parse(&slice_header_block.data)?;

    let is_multi_ref = sh.ref_seq_id == -2;

    // r[impl cram.edge.unknown_read_names]
    if !ch.preservation.read_names_included {
        warn!(
            "CRAM file has RN=false (read names not stored). Mate-pair overlap dedup will be unreliable."
        );
    }

    // r[impl cram.edge.empty_slice]
    if sh.num_records == 0 {
        return Ok((0, 0));
    }

    // r[impl cram.edge.reference_mismatch]
    if !is_multi_ref && sh.reference_md5 != [0u8; 16] && !reference_seq.is_empty() {
        let slice_start_0based = Pos1::try_from(sh.alignment_start.max(1))
            .map_err(|_| CramError::InvalidPosition { value: i64::from(sh.alignment_start) })?
            .to_zero_based()
            .as_i64();
        let slice_ref_start_i64 = slice_start_0based.wrapping_sub(ref_start_0based);
        let slice_ref_start = usize::try_from(slice_ref_start_i64)
            .map_err(|_| CramError::InvalidPosition { value: slice_ref_start_i64 })?;
        // r[impl cram.slice.validated_lengths]
        let span = usize::try_from(sh.alignment_span)
            .map_err(|_| CramError::InvalidLength { value: sh.alignment_span })?;
        let slice_ref_end = slice_ref_start.wrapping_add(span);
        if let Some(slice_ref) = reference_seq.get(slice_ref_start..slice_ref_end) {
            let computed =
                md5::compute(slice_ref.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>());
            if *computed != sh.reference_md5 {
                return Err(CramError::ReferenceMd5Mismatch {
                    contig: header.target_name(tid).unwrap_or("?").into(),
                    start: sh.alignment_start as u64,
                    end: sh.alignment_start.wrapping_add(sh.alignment_span) as u64,
                });
            }
        }
    }

    // Parse remaining blocks: core data + external data
    let mut core_block: Option<Block> = None;
    let mut external_blocks: SmallVec<(i32, ExternalCursor), 4> = SmallVec::new();

    for _ in 0..sh.num_blocks {
        let remaining =
            slice_data.get(pos..).ok_or(CramError::Truncated { context: "slice block" })?;
        let (blk, consumed) = block::parse_block_with_buf(
            remaining,
            rans_4x8_buf.as_mut(),
            nx16_order1_buf.as_mut(),
        )?;
        pos = pos.wrapping_add(consumed);

        match blk.content_type {
            ContentType::CoreData => {
                core_block = Some(blk);
            }
            ContentType::ExternalData => {
                external_blocks.push((blk.content_id, ExternalCursor::new(blk.data)));
            }
            _ => {
                // Ignore unexpected block types
            }
        }
    }

    // r[impl cram.slice.embedded_ref]
    let embedded_ref_data = if sh.embedded_reference >= 0 {
        let idx = external_blocks.iter().position(|(id, _)| *id == sh.embedded_reference);
        match idx {
            Some(i) => external_blocks.swap_remove(i).1.into_data(),
            None => Vec::new(),
        }
    } else {
        Vec::new()
    };
    let effective_ref: &[u8] =
        if sh.embedded_reference >= 0 { &embedded_ref_data } else { reference_seq };

    let core_data = core_block.ok_or(CramError::MissingCoreDataBlock)?;

    let mut ctx = DecodeContext::new(&core_data.data, external_blocks);

    // Decode records, collecting mate cross-reference info for TLEN reconstruction.
    let mut alignment_pos = i64::from(sh.alignment_start);
    let mut fetched_count = 0usize;
    let mut kept_count = 0usize;
    let num_records = usize::try_from(sh.num_records)
        .map_err(|_| CramError::InvalidLength { value: sh.num_records })?;
    let mut mate_infos: Vec<SliceMateInfo> = Vec::with_capacity(num_records);

    for record_index in 0..num_records {
        // For embedded reference, ref_start is the slice's alignment_start
        let effective_ref_start = if sh.embedded_reference >= 0 {
            Pos1::try_from(sh.alignment_start.max(1))
                .map_err(|_| CramError::InvalidPosition { value: i64::from(sh.alignment_start) })?
                .to_zero_based()
                .as_i64()
        } else {
            ref_start_0based
        };

        let (fetched, mate_info) = decode_record(
            ch,
            &sh,
            is_multi_ref,
            &mut ctx,
            &mut alignment_pos,
            effective_ref,
            effective_ref_start,
            read_group_ids,
            tid,
            query_start,
            query_end,
            store,
            cigar_buf,
            bases_buf,
            qual_buf,
            aux_buf,
            name_buf,
            feature_byte_buf,
            cigar_ops_buf,
            record_index,
            customize,
        )?;
        if fetched {
            fetched_count = fetched_count.wrapping_add(1);
            if mate_info.store_idx.is_some() {
                kept_count = kept_count.wrapping_add(1);
            }
        }
        mate_infos.push(mate_info);
    }

    // Post-process: resolve template_len for attached/downstream mates.
    // Follow mate_line chains to compute alignment span, then assign signed TLEN.
    resolve_mate_tlen(&mate_infos, store);

    Ok((fetched_count, kept_count))
}

/// Decode a single CRAM record from the decode context.
///
/// Returns `(fetched, mate_info)` where `fetched` is `true` if the record
/// reached the push step (i.e., was not rejected by the reader's built-in
/// overlap/tid/unmapped checks) and `false` otherwise. `mate_info.store_idx`
/// is `Some(idx)` if the record survived both the reader check and the
/// user's `customize.filter`, and `None` if either dropped it — the
/// mate-resolution pass already handles the `None` case.
// r[impl cram.record.decode_order]
// r[impl cram.fetch_into_customized.push_time]
#[expect(
    clippy::too_many_arguments,
    reason = "CRAM record decoding requires compression header, slice header, context, and customize parameters"
)]
fn decode_record<E: CustomizeRecordStore>(
    ch: &CompressionHeader,
    sh: &SliceHeader,
    is_multi_ref: bool,
    ctx: &mut DecodeContext<'_>,
    prev_alignment_pos: &mut i64,
    reference_seq: &[u8],
    ref_start_0based: i64,
    read_group_ids: &[SmolStr],
    tid: u32,
    query_start: Pos0,
    query_end: Pos0,
    store: &mut RecordStore<E::Extra>,
    cigar_buf: &mut Vec<CigarOp>,
    bases_buf: &mut Vec<Base>,
    qual_buf: &mut Vec<u8>,
    aux_buf: &mut Vec<u8>,
    name_buf: &mut Vec<u8>,
    feature_byte_buf: &mut Vec<u8>,
    cigar_ops_buf: &mut Vec<(u32, u8)>,
    record_index: usize,
    customize: &mut E,
) -> Result<(bool, SliceMateInfo), CramError> {
    let ds = &ch.data_series;

    // r[impl cram.record.flags]
    // 1. BF (BAM flags)
    let raw_flags = ds.bam_flags.decode(ctx)?;
    let raw_flags: u16 = u16::try_from(raw_flags)
        .map_err(|source| CramError::InvalidBamFlags { value: raw_flags, source })?;
    let bam_flags = BamFlags::from(raw_flags);

    // 2. CF (CRAM flags)
    #[expect(
        clippy::cast_sign_loss,
        reason = "CRAM flags are small non-negative integers encoded as ITF8 i32"
    )]
    let cram_flags = ds.cram_flags.decode(ctx)? as u32;
    let quality_as_array = cram_flags & 0x1 != 0;
    let detached = cram_flags & 0x2 != 0;
    let mate_downstream = cram_flags & 0x4 != 0;
    let seq_unknown = cram_flags & 0x8 != 0;

    // r[impl cram.slice.multi_ref]
    // 3. RI (ref ID — only for multi-ref slices)
    let record_ref_id = if is_multi_ref { ds.ref_id.decode(ctx)? } else { sh.ref_seq_id };

    // r[impl cram.record.read_length]
    // 4. RL (read length)
    let read_length_i32 = ds.read_length.decode(ctx)?;
    let read_length = usize::try_from(read_length_i32)
        .map_err(|_| CramError::InvalidLength { value: read_length_i32 })?;

    // r[impl cram.record.position]
    // 5. AP (alignment position)
    let ap = i64::from(ds.alignment_pos.decode(ctx)?);
    let alignment_pos = if ch.preservation.ap_delta {
        *prev_alignment_pos = prev_alignment_pos.wrapping_add(ap);
        *prev_alignment_pos
    } else {
        *prev_alignment_pos = ap;
        ap
    };
    // r[impl cram.edge.position_overflow]
    // Convert from 1-based to 0-based
    let pos_0based = Pos1::try_from(alignment_pos)
        .map(|p| p.to_zero_based())
        .map_err(|_| super::reader::CramError::InvalidPosition { value: alignment_pos })?;

    // r[impl cram.record.read_group]
    // 6. RG (read group)
    let read_group = ds.read_group.decode(ctx)?;

    // r[impl cram.record.read_name]
    // 7. RN (read name) — decode straight into the caller's reusable
    // `name_buf` to avoid a per-record `Vec<u8>` allocation.
    name_buf.clear();
    if ch.preservation.read_names_included {
        ds.read_name.decode_into(ctx, name_buf)?;
    }

    // r[impl cram.record.mate_detached]
    // r[impl cram.record.mate_attached]
    // 8. Mate data
    let mut next_pos_val: i32 = -1;
    let mut next_ref_id_val: i32 = -1;
    let mut template_len_val: i32 = 0;
    let mut mate_line: i32 = -1;
    if detached {
        let _mate_flags = ds.mate_flags.decode(ctx)?;
        if !ch.preservation.read_names_included {
            // Detached + RN=false: the codec stream carries a read-name
            // (the mate's, per the CRAM spec) that we must consume but
            // also must discard — when read_names_included is false the
            // record's qname is unknown. Decode into `name_buf` (which is
            // already empty for this branch) and immediately clear so
            // qname stays empty downstream. No allocation in steady state.
            ds.read_name.decode_into(ctx, name_buf)?;
            name_buf.clear();
        }
        next_ref_id_val = ds.next_segment_ref.decode(ctx)?;
        next_pos_val = ds.next_mate_pos.decode(ctx)?;
        template_len_val = ds.template_size.decode(ctx)?;
    } else if mate_downstream {
        let nf = ds.next_fragment.decode(ctx)?;
        // NF is relative: absolute mate index = current + 1 + nf
        #[expect(
            clippy::cast_possible_wrap,
            clippy::cast_possible_truncation,
            reason = "record_index is bounded by num_records (i32), fits in i32"
        )]
        {
            mate_line = (record_index as i32).wrapping_add(1).wrapping_add(nf);
        }
    }

    // 9. TL (tag line index)
    let tag_line_idx_i32 = ds.tag_line.decode(ctx)?;
    let tag_line_idx = usize::try_from(tag_line_idx_i32)
        .map_err(|_| CramError::InvalidLength { value: tag_line_idx_i32 })?;

    // r[impl cram.record.aux_tags]
    // r[impl base_mod.passthrough.cram] — MM (Z) and ML (B:C) flow through the generic
    //   tag-encoding path without truncation; arbitrary lengths are emitted into the
    //   reconstructed BAM aux block as the encoded `tag_value` bytes.
    // 10. Decode tag values
    aux_buf.clear();
    if let Some(tag_set) = ch.preservation.tag_dictionary.get(tag_line_idx) {
        for entry in tag_set {
            let tag_key = (i32::from(entry.tag[0]) << 16)
                | (i32::from(entry.tag[1]) << 8)
                | i32::from(entry.bam_type);
            if let Some(enc) = ch.tag_encodings.get(&tag_key) {
                // Serialize to BAM binary aux format: tag[0] tag[1] type value_bytes.
                // Decoding the value directly into `aux_buf` skips the
                // per-tag `Vec<u8>` allocation that the old `decode` API
                // forced.
                aux_buf.push(entry.tag[0]);
                aux_buf.push(entry.tag[1]);
                aux_buf.push(entry.bam_type);
                enc.decode_into(ctx, aux_buf)?;
            }
        }
    }

    // r[impl cram.record.rg_tag]
    // Add RG tag if read group >= 0. RG IDs are pre-parsed from the header
    // once at open time (cached on `CramShared::read_group_ids`); used to
    // be a per-record header text scan + String allocation, which dominated
    // `decode_record` profiles.
    if let Ok(rg_idx) = usize::try_from(read_group)
        && let Some(rg_id) = read_group_ids.get(rg_idx)
    {
        aux_buf.extend_from_slice(b"RGZ");
        aux_buf.extend_from_slice(rg_id.as_bytes());
        aux_buf.push(0); // null terminator for Z-type
    }

    let is_unmapped = bam_flags.is_unmapped();

    // 11. Mapped read: features + MQ + quality
    cigar_buf.clear();
    bases_buf.clear();
    qual_buf.clear();
    // r[impl cram.edge.long_reads]
    bases_buf.reserve(read_length);
    qual_buf.reserve(read_length);

    if !is_unmapped {
        let feature_count_i32 = ds.feature_count.decode(ctx)?;
        let feature_count = usize::try_from(feature_count_i32)
            .map_err(|_| CramError::InvalidLength { value: feature_count_i32 })?;

        // Decode features and reconstruct sequence + CIGAR
        let result = decode_features_and_reconstruct(
            ch,
            ctx,
            feature_count,
            read_length,
            pos_0based.as_i64(),
            ref_start_0based,
            reference_seq,
            cigar_buf,
            bases_buf,
            feature_byte_buf,
            cigar_ops_buf,
        )?;

        // MQ
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "mapping quality is 0..=255 encoded as ITF8 i32; only lower 8 bits are valid"
        )]
        let mapq = ds.mapping_quality.decode(ctx)? as u8;

        // r[impl cram.record.quality]
        // Quality scores. Bulk decode for External (single external-block
        // lookup + one memcpy of `read_length` bytes) instead of the
        // per-byte `decode(ctx)` loop that dominated profiles.
        if quality_as_array {
            ds.quality_score.decode_n_into(ctx, read_length, qual_buf)?;
        } else {
            qual_buf.resize(read_length, 0xFF);
        }

        // `end_pos_raw` is the half-open exclusive reference end
        // (`pos + ref_consumed`). It is used for TLEN math
        // (`aright - aleft`), where half-open arithmetic is the natural
        // form. The value stored in `SlimRecord::end_pos` is the 0-based
        // **inclusive** last reference position covered, matching the
        // BAM/SAM path's `compute_end_pos` and the pileup engine's
        // eviction check (`end_pos < pos`). For zero-refspan reads we
        // store `pos` itself, the same fallback BAM uses.
        // r[impl cram.record.end_pos]
        let end_pos_raw = pos_0based.as_i64().wrapping_add(i64::from(result.ref_consumed));
        let end_pos_inclusive_raw = if result.ref_consumed == 0 {
            pos_0based.as_i64()
        } else {
            end_pos_raw.wrapping_sub(1)
        };
        let end_pos = Pos0::try_from(end_pos_inclusive_raw).map_err(|_| {
            super::reader::CramError::InvalidPosition { value: end_pos_inclusive_raw }
        })?;

        // r[impl cram.index.multi_ref_slices]
        // Skip records from different references in multi-ref slices
        #[expect(
            clippy::cast_possible_wrap,
            reason = "tid comes from BAM header, capped at MAX_REFERENCES (1M), well within i32"
        )]
        if record_ref_id != tid as i32 {
            return Ok((
                false,
                SliceMateInfo {
                    pos: pos_0based.as_i64(),
                    end_pos: end_pos_raw,
                    mate_line,
                    store_idx: None,
                    bam_flags: raw_flags,
                    ref_id: record_ref_id,
                },
            ));
        }

        // Check overlap with query region. Both `query_start` and
        // `query_end` are inclusive (matches the BAM IndexedReader's
        // convention: a record is kept iff `rec_pos <= query_end &&
        // rec_end >= query_start`). `end_pos` is the inclusive last
        // covered position. Mismatching this with the half-open
        // semantics used elsewhere drops boundary records (reads
        // starting exactly at the requested end position).
        if pos_0based > query_end || end_pos < query_start {
            return Ok((
                false,
                SliceMateInfo {
                    pos: pos_0based.as_i64(),
                    end_pos: end_pos_raw,
                    mate_line,
                    store_idx: None,
                    bam_flags: raw_flags,
                    ref_id: record_ref_id,
                },
            ));
        }

        let qname: &[u8] = name_buf;

        // r[impl cram.fetch_into_customized.push_time]
        // Push the record, then consult the user's filter. `push_fields`
        // returns `Ok(None)` when the filter rejects and rolls back slab
        // writes with zero waste (same as BAM/SAM).
        let store_idx = store.push_fields(
            pos_0based,
            end_pos,
            bam_flags,
            mapq,
            result.matching_bases,
            result.indel_bases,
            qname,
            cigar_buf,
            bases_buf,
            qual_buf,
            aux_buf,
            record_ref_id,
            next_ref_id_val,
            next_pos_val,
            template_len_val,
            customize,
        )?;

        return Ok((
            true,
            SliceMateInfo {
                pos: pos_0based.as_i64(),
                end_pos: end_pos_raw,
                mate_line,
                store_idx,
                bam_flags: raw_flags,
                ref_id: record_ref_id,
            },
        ));
    }

    // Unmapped read
    // r[impl cram.edge.seq_unknown]
    if !seq_unknown {
        for _ in 0..read_length {
            let base_byte = ds.base.decode(ctx)?;
            bases_buf.push(Base::from(base_byte));
        }
    }

    // Quality. See the mapped-read branch above for the bulk-decode
    // rationale.
    if quality_as_array {
        ds.quality_score.decode_n_into(ctx, read_length, qual_buf)?;
    } else {
        qual_buf.resize(read_length, 0xFF);
    }

    // r[impl cram.edge.unmapped_reads]
    // Unmapped reads — skip for fetch_into (same as BAM/SAM)
    Ok((
        false,
        SliceMateInfo {
            pos: -1,
            end_pos: -1,
            mate_line: -1,
            store_idx: None,
            bam_flags: raw_flags,
            ref_id: record_ref_id,
        },
    ))
}

// r[impl cram.record.mate_tlen_reconstruction]
/// Post-process mate cross-references to reconstruct `template_len` for
/// attached/downstream mates. Mirrors htslib's `cram_decode_slice_xref`.
///
/// For each record with a downstream mate link, follows the chain to find
/// the alignment span (min start to max end) across all fragments, then
/// assigns signed TLEN: positive for the leftmost fragment, negative for
/// the rightmost, with ties broken by the READ1 flag.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "everything is based on infos.len()"
)]
fn resolve_mate_tlen<U>(infos: &[SliceMateInfo], store: &mut RecordStore<U>) {
    let n = infos.len();
    // Track which records we've already resolved to avoid reprocessing.
    let mut resolved = vec![false; n];

    for i in 0..n {
        if resolved[i] || infos[i].mate_line < 0 {
            continue;
        }

        // Collect all records in this mate chain.
        let mut chain: Vec<usize> = Vec::with_capacity(4);
        let mut idx = i;
        loop {
            if idx >= n || resolved[idx] {
                break;
            }
            chain.push(idx);
            resolved[idx] = true;
            let ml = infos[idx].mate_line;
            if ml < 0 || ml as usize >= n {
                break;
            }
            idx = ml as usize;
            // Safety: break cycles
            if chain.len() > n {
                break;
            }
        }

        if chain.len() < 2 {
            continue;
        }

        // Compute alignment span and position counts across all fragments.
        let mut aleft = i64::MAX;
        let mut aright = i64::MIN;
        let mut left_cnt = 0u32;
        let mut right_cnt = 0u32;

        for &j in &chain {
            let info = &infos[j];
            if info.pos < 0 {
                continue; // unmapped
            }
            if info.pos < aleft {
                aleft = info.pos;
                left_cnt = 1;
            } else if info.pos == aleft {
                left_cnt = left_cnt.wrapping_add(1);
            }
            if info.end_pos > aright {
                aright = info.end_pos;
                right_cnt = 1;
            } else if info.end_pos == aright {
                right_cnt = right_cnt.wrapping_add(1);
            }
        }

        if aright <= aleft {
            continue; // unmapped or degenerate
        }

        // In 0-based half-open coordinates: tlen = max(end) - min(start).
        // Equivalent to htslib's (aright - aleft + 1) in 1-based inclusive.
        #[expect(
            clippy::cast_possible_truncation,
            reason = "template length fits in i32 for genomic data"
        )]
        let tlen_magnitude = (aright.wrapping_sub(aleft)) as i32;

        // Determine sign for the FIRST record in the chain, then assign
        // the opposite sign to all remaining records. This mirrors htslib's
        // cram_decode_slice_xref sign assignment logic.
        let first = &infos[chain[0]];
        let first_tlen = if first.pos == aleft && (first.end_pos < aright || left_cnt <= 1) {
            // Leftmost, and either not rightmost or unique leftmost
            tlen_magnitude
        } else if first.pos == aleft && first.end_pos == aright && left_cnt > 1 && right_cnt > 1 {
            // Both leftmost and rightmost with ties — break by READ1 flag
            if first.bam_flags & 0x40 != 0 { tlen_magnitude } else { tlen_magnitude.neg() }
        } else {
            // Rightmost or internal
            tlen_magnitude.neg()
        };

        // First record gets first_tlen; all remaining get the opposite.
        if let Some(idx) = first.store_idx {
            store.set_template_len(idx, first_tlen);
        }
        let rest_tlen = first_tlen.neg();
        for &j in &chain[1..] {
            if let Some(idx) = infos[j].store_idx {
                store.set_template_len(idx, rest_tlen);
            }
        }

        // Resolve next_ref_id and next_pos for each record in the chain.
        // Each record's mate is the next entry in the chain (last wraps to first).
        // r[impl cram.fetch_into_customized.filtered_mate_sentinel]
        // If the mate was rejected by the user's filter (store_idx == None),
        // null out the kept record's mate fields so the BAM "mate unavailable"
        // sentinel (-1, -1) reflects what's actually in the store. TLEN is left
        // as the full-chain span — it's a per-template property and matches what
        // an unfiltered fetch would compute.
        for ci in 0..chain.len() {
            let mate_ci = if ci + 1 < chain.len() { ci + 1 } else { 0 };
            let record_idx = infos[chain[ci]].store_idx;
            let mate = &infos[chain[mate_ci]];
            let Some(idx) = record_idx else { continue };
            if mate.store_idx.is_none() {
                // Mate was filtered out — set the BAM "mate unavailable" sentinel.
                store.set_mate_info(idx, -1, -1);
                continue;
            }
            #[expect(
                clippy::cast_possible_truncation,
                reason = "mate positions fit in i32 for genomic data"
            )]
            let mate_pos = if mate.pos < 0 { -1 } else { mate.pos as i32 };
            store.set_mate_info(idx, mate.ref_id, mate_pos);
        }
    }
}

struct ReconstructResult {
    ref_consumed: u32,
    matching_bases: u32,
    indel_bases: u32,
}

// r[impl cram.record.features]
// r[impl cram.record.sequence]
// r[impl cram.record.cigar_reconstruction]
// r[impl cram.slice.ref_bounds_warning+2]
/// Look up a reference base by index, logging a warning on the first out-of-bounds access
/// per slice and falling back to `b'N'`.
fn ref_base_at(reference_seq: &[u8], index: usize, warned: &mut bool) -> u8 {
    match reference_seq.get(index) {
        Some(&b) => b,
        None => {
            if !*warned {
                warn!(
                    index,
                    ref_len = reference_seq.len(),
                    "reference shorter than expected during CRAM sequence reconstruction; \
                     substituting N (may indicate reference/CRAM mismatch)"
                );
                *warned = true;
            }
            b'N'
        }
    }
}

/// Decode features and reconstruct sequence + CIGAR in a single fused
/// pass.
///
/// Earlier this was two passes: collect every feature into a local
/// `Vec<Feature>` (with `Vec<u8>` payloads for I/S/B/q), then walk the
/// read filling reference + applying features. That allocated three
/// times per record (`features` + each `FeatureData::Insertion(Vec)` /
/// `SoftClip(Vec)` / `Bases(Vec)`) and a fourth `cigar_ops: Vec::new()`.
///
/// CRAM features are emitted in *increasing* read-position order
/// (`feature_pos` is delta-encoded), so we can decode each feature just
/// before applying it: fill reference bases up to the feature's anchor,
/// decode the feature's payload into the caller's reusable
/// `feature_byte_buf`, apply, repeat. `cigar_ops_buf` is also threaded
/// in so the run-length encoder works in-place. Net effect: zero
/// allocations per record after the first record warms the buffers.
#[expect(
    clippy::too_many_arguments,
    reason = "sequence reconstruction requires all CRAM encoding handles, reference slice, and output buffers"
)]
fn decode_features_and_reconstruct(
    ch: &CompressionHeader,
    ctx: &mut DecodeContext<'_>,
    feature_count: usize,
    read_length: usize,
    pos_0based: i64,
    slice_start_0based: i64,
    reference_seq: &[u8],
    cigar_buf: &mut Vec<CigarOp>,
    bases_buf: &mut Vec<Base>,
    feature_byte_buf: &mut Vec<u8>,
    cigar_ops_buf: &mut Vec<(u32, u8)>,
) -> Result<ReconstructResult, CramError> {
    let ds = &ch.data_series;

    // pos_0based >= slice_start_0based is guaranteed by CRAM structure; a negative
    // difference would indicate a malformed file, so clamp to 0 and let ref_base_at
    // handle out-of-bounds via its warning path.
    let ref_offset = usize::try_from(pos_0based.wrapping_sub(slice_start_0based)).unwrap_or(0);
    let mut read_pos = 0usize; // 0-based position in the read
    let mut ref_pos = 0usize; // 0-based position relative to alignment start on reference
    let mut ref_warned = false; // log at most one warning per reconstruction
    let mut matching_bases = 0u32;
    let mut indel_bases = 0u32;

    cigar_ops_buf.clear();

    fn push_cigar_op(ops: &mut Vec<(u32, u8)>, len: u32, op: u8) {
        if len == 0 {
            return;
        }
        if let Some(last) = ops.last_mut()
            && last.1 == op
        {
            last.0 = last.0.saturating_add(len);
            return;
        }
        ops.push((len, op));
    }

    /// Emit one matching base copied from the reference, advancing both
    /// `read_pos` and `ref_pos`. Used both for fill-up-to-next-feature
    /// and for quality-only features (which don't modify the sequence).
    #[expect(clippy::too_many_arguments, reason = "shared helper for ref-match emit")]
    fn emit_ref_match(
        reference_seq: &[u8],
        ref_offset: usize,
        ref_pos: &mut usize,
        read_pos: &mut usize,
        ref_warned: &mut bool,
        bases_buf: &mut Vec<Base>,
        cigar_ops: &mut Vec<(u32, u8)>,
        matching_bases: &mut u32,
    ) {
        let ref_base = ref_base_at(reference_seq, ref_offset.wrapping_add(*ref_pos), ref_warned);
        bases_buf.push(Base::from(ref_base));
        push_cigar_op(cigar_ops, 1, 0); // M
        *matching_bases = matching_bases.saturating_add(1);
        *read_pos = read_pos.wrapping_add(1);
        *ref_pos = ref_pos.wrapping_add(1);
    }

    let mut feature_read_pos = 0u32; // 1-based, accumulator for delta-encoded positions

    for _ in 0..feature_count {
        let fc = ds.feature_code.decode(ctx)?;
        let fp = ds.feature_pos.decode(ctx)? as u32;
        feature_read_pos = feature_read_pos.wrapping_add(fp);
        // Feature positions are 1-based in CRAM; the 0-based read_pos
        // catches up to `feature_read_pos - 1` before applying.
        let feat_target = (feature_read_pos as usize).saturating_sub(1);

        // Fill reference matches until read_pos reaches the feature's anchor.
        while read_pos < feat_target && read_pos < read_length {
            emit_ref_match(
                reference_seq,
                ref_offset,
                &mut ref_pos,
                &mut read_pos,
                &mut ref_warned,
                bases_buf,
                cigar_ops_buf,
                &mut matching_bases,
            );
        }

        // Apply the feature inline. Payload-bearing variants decode straight
        // into `feature_byte_buf` (cleared, then filled) — no per-feature alloc.
        // r[related cram.slice.ref_bounds_warning+2]
        match fc {
            // X — substitution: 1 read base, 1 ref base, M cigar.
            b'X' => {
                let bs = ds.base_sub.decode(ctx)?;
                let ref_base =
                    ref_base_at(reference_seq, ref_offset.wrapping_add(ref_pos), &mut ref_warned);
                let read_base = ch.preservation.substitution_matrix.substitute(ref_base, bs);
                bases_buf.push(Base::from(read_base));
                push_cigar_op(cigar_ops_buf, 1, 0); // M
                matching_bases = matching_bases.saturating_add(1);
                read_pos = read_pos.wrapping_add(1);
                ref_pos = ref_pos.wrapping_add(1);
            }
            // I — insertion: N read bases, 0 ref, I cigar.
            b'I' => {
                feature_byte_buf.clear();
                ds.insertion.decode_into(ctx, feature_byte_buf)?;
                let len = feature_byte_buf.len();
                bases_buf.reserve(len);
                for &b in feature_byte_buf.iter() {
                    bases_buf.push(Base::from(b));
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "insertion length is bounded by read_length (validated); fits in u32"
                )]
                let len_u32 = len as u32;
                push_cigar_op(cigar_ops_buf, len_u32, 1); // I
                indel_bases = indel_bases.saturating_add(len_u32);
                read_pos = read_pos.wrapping_add(len);
            }
            // i — single-base insertion.
            b'i' => {
                let base = ds.base.decode(ctx)?;
                bases_buf.push(Base::from(base));
                push_cigar_op(cigar_ops_buf, 1, 1); // I
                indel_bases = indel_bases.saturating_add(1);
                read_pos = read_pos.wrapping_add(1);
            }
            // D — deletion: 0 read, N ref, D cigar.
            b'D' => {
                let len = ds.deletion_length.decode(ctx)? as u32;
                push_cigar_op(cigar_ops_buf, len, 2); // D
                indel_bases = indel_bases.saturating_add(len);
                ref_pos = ref_pos.wrapping_add(len as usize);
            }
            // S — soft clip: N read bases, 0 ref, S cigar.
            b'S' => {
                feature_byte_buf.clear();
                ds.soft_clip.decode_into(ctx, feature_byte_buf)?;
                let len = feature_byte_buf.len();
                bases_buf.reserve(len);
                for &b in feature_byte_buf.iter() {
                    bases_buf.push(Base::from(b));
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "soft-clip length is bounded by read_length (validated); fits in u32"
                )]
                let len_u32 = len as u32;
                push_cigar_op(cigar_ops_buf, len_u32, 4); // S
                read_pos = read_pos.wrapping_add(len);
            }
            // H — hard clip: emit cigar but consume neither read nor ref.
            b'H' => {
                let len = ds.hard_clip.decode(ctx)? as u32;
                push_cigar_op(cigar_ops_buf, len, 5); // H
            }
            // N — ref skip: 0 read, N ref, N cigar.
            b'N' => {
                let len = ds.ref_skip.decode(ctx)? as u32;
                push_cigar_op(cigar_ops_buf, len, 3); // N
                ref_pos = ref_pos.wrapping_add(len as usize);
            }
            // P — padding: emit cigar but consume neither read nor ref.
            b'P' => {
                let len = ds.padding.decode(ctx)? as u32;
                push_cigar_op(cigar_ops_buf, len, 6); // P
            }
            // B — base + quality (BA + QS): treated as a 1-base match.
            // The qual byte is decoded but currently ignored (matches previous
            // behaviour; see r[unimpl cram.feature.quality_score]).
            b'B' => {
                let base = ds.base.decode(ctx)?;
                let _qual = ds.quality_score.decode(ctx)?;
                bases_buf.push(Base::from(base));
                push_cigar_op(cigar_ops_buf, 1, 0); // M
                matching_bases = matching_bases.saturating_add(1);
                read_pos = read_pos.wrapping_add(1);
                ref_pos = ref_pos.wrapping_add(1);
            }
            // b — bases block: N read bases, N ref, M cigar.
            b'b' => {
                feature_byte_buf.clear();
                ds.bases_block.decode_into(ctx, feature_byte_buf)?;
                let len = feature_byte_buf.len();
                bases_buf.reserve(len);
                for &b in feature_byte_buf.iter() {
                    bases_buf.push(Base::from(b));
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "bases length is bounded by read_length (validated); fits in u32"
                )]
                let len_u32 = len as u32;
                push_cigar_op(cigar_ops_buf, len_u32, 0); // M
                matching_bases = matching_bases.saturating_add(len_u32);
                read_pos = read_pos.wrapping_add(len);
                ref_pos = ref_pos.wrapping_add(len);
            }
            // q — quality block: data is consumed (so the codec stream stays
            // aligned), value discarded. Anchors a 1-base reference match.
            // r[unimpl cram.feature.quality_block]
            b'q' => {
                feature_byte_buf.clear();
                ds.quality_block.decode_into(ctx, feature_byte_buf)?;
                emit_ref_match(
                    reference_seq,
                    ref_offset,
                    &mut ref_pos,
                    &mut read_pos,
                    &mut ref_warned,
                    bases_buf,
                    cigar_ops_buf,
                    &mut matching_bases,
                );
            }
            // Q — single quality score: same shape as q but no payload.
            // r[unimpl cram.feature.quality_score]
            b'Q' => {
                let _qual = ds.quality_score.decode(ctx)?;
                emit_ref_match(
                    reference_seq,
                    ref_offset,
                    &mut ref_pos,
                    &mut read_pos,
                    &mut ref_warned,
                    bases_buf,
                    cigar_ops_buf,
                    &mut matching_bases,
                );
            }
            _ => return Err(CramError::UnknownFeatureCode { feature_code: fc }),
        }
    }

    // Tail: fill reference matches until the read length is reached.
    while read_pos < read_length {
        emit_ref_match(
            reference_seq,
            ref_offset,
            &mut ref_pos,
            &mut read_pos,
            &mut ref_warned,
            bases_buf,
            cigar_ops_buf,
            &mut matching_bases,
        );
    }

    // Pack CIGAR ops into typed BAM-layout CigarOps
    cigar_buf.clear();
    cigar_buf.reserve(cigar_ops_buf.len());
    for &(len, op) in cigar_ops_buf.iter() {
        cigar_buf.push(CigarOp::from_bam_u32((len << 4) | u32::from(op)));
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "ref_pos is bounded by reference sequence length which fits in u32; BAM positions are capped at 2^29"
    )]
    Ok(ReconstructResult { ref_consumed: ref_pos as u32, matching_bases, indel_bases })
}

fn read_itf8(cursor: &mut &[u8]) -> Result<u32, CramError> {
    varint::read_itf8_from(cursor).ok_or(CramError::Truncated { context: "slice header itf8" })
}

fn read_ltf8(cursor: &mut &[u8]) -> Result<u64, CramError> {
    varint::read_ltf8_from(cursor).ok_or(CramError::Truncated { context: "slice header ltf8" })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_test_cram() -> Vec<u8> {
        std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram")).unwrap()
    }

    #[test]
    fn parse_slice_header_from_real_cram() {
        use super::super::container::ContainerHeader;

        let data = load_test_cram();
        let after_file_def = &data[26..];

        // Skip header container
        let hdr_container = ContainerHeader::parse(after_file_def).unwrap();
        let data_start = hdr_container.header_size + hdr_container.length as usize;

        // Parse first data container
        let dc_bytes = &after_file_def[data_start..];
        let dc = ContainerHeader::parse(dc_bytes).unwrap();

        // Parse compression header block
        let container_data = &dc_bytes[dc.header_size..];
        let (comp_block, _) = block::parse_block(container_data).unwrap();
        assert_eq!(comp_block.content_type, ContentType::CompressionHeader);

        // First slice header should be at the first landmark offset
        let first_landmark = *dc.landmarks.first().unwrap();
        let slice_data = &container_data[first_landmark as usize..];
        let (slice_header_block, _) = block::parse_block(slice_data).unwrap();
        assert_eq!(slice_header_block.content_type, ContentType::SliceHeader);

        let sh = SliceHeader::parse(&slice_header_block.data).unwrap();
        assert!(sh.num_records > 0, "slice should have records");
        assert!(sh.alignment_start > 0, "slice should have valid alignment start");
        assert!(sh.num_blocks > 0, "slice should have blocks");
    }

    // r[verify cram.edge.reference_mismatch]
    #[test]
    fn md5_mismatch_detected_with_wrong_reference() {
        // Verify that decode_slice catches a reference MD5 mismatch.
        // We do this by reading a real CRAM slice and passing a wrong reference.
        use super::super::{compression_header::CompressionHeader, container::ContainerHeader};

        let data = load_test_cram();
        let after_file_def = &data[26..];
        let hdr = ContainerHeader::parse(after_file_def).unwrap();
        let data_start = hdr.header_size + hdr.length as usize;
        let dc_bytes = &after_file_def[data_start..];
        let dc = ContainerHeader::parse(dc_bytes).unwrap();
        let container_data = &dc_bytes[dc.header_size..];
        let (comp_block, _) = block::parse_block(container_data).unwrap();
        let ch = CompressionHeader::parse(&comp_block.data).unwrap();

        let bam_header =
            crate::bam::BamHeader::from_sam_text("@HD\tVN:1.6\n@SQ\tSN:chr19\tLN:61431566\n")
                .unwrap();

        // Create a fake reference full of N's — MD5 will not match
        let fake_ref = vec![b'N'; 100_000];
        let ref_start = Pos1::try_from(dc.alignment_start.max(1))
            .map(|p| p.to_zero_based().as_i64())
            .unwrap_or(0);

        let first_landmark = *dc.landmarks.first().unwrap();
        let result = decode_slice(
            &ch,
            container_data,
            first_landmark as usize,
            &fake_ref,
            ref_start,
            &bam_header,
            &[],
            0,
            Pos0::new(0).unwrap(),
            Pos0::max_value(),
            &mut crate::bam::record_store::RecordStore::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut Vec::new(),
            &mut (),
            &mut None,
            &mut None,
        );

        assert!(
            matches!(result, Err(CramError::ReferenceMd5Mismatch { .. })),
            "should detect MD5 mismatch, got: {result:?}"
        );
    }

    #[test]
    fn expected_slice_header_wrong_block_type() {
        // Build a block with ExternalData instead of SliceHeader and check the error variant.
        let data = b"\x01\x02\x03";
        let block_bytes = block::build_test_block(4, 0, data); // 4 = ExternalData
        let (blk, _) = block::parse_block(&block_bytes).unwrap();
        assert_eq!(blk.content_type, ContentType::ExternalData);
        // The error is constructed in decode_slice when the first block isn't a SliceHeader.
        let err = CramError::ExpectedSliceHeader { found: blk.content_type };
        assert!(matches!(err, CramError::ExpectedSliceHeader { found: ContentType::ExternalData }));
    }

    #[test]
    fn missing_core_data_block_error_variant() {
        // Verify the MissingCoreDataBlock error variant exists and formats correctly.
        let err = CramError::MissingCoreDataBlock;
        let msg = format!("{err}");
        assert!(
            msg.contains("core data block"),
            "error message should mention core data block: {msg}"
        );
    }

    #[test]
    fn unknown_feature_code_error_variant() {
        // UnknownFeatureCode is returned from decode_features_and_reconstruct when
        // encountering a feature code byte that doesn't match any known code.
        // The known codes are: X, I, i, D, S, H, N, P, B, b, q, Q
        // Value 0xFF is not a valid feature code.
        let err = CramError::UnknownFeatureCode { feature_code: 0xFF };
        assert!(matches!(err, CramError::UnknownFeatureCode { feature_code: 0xFF }));
        let msg = format!("{err}");
        assert!(
            msg.contains("ff") || msg.contains("0xff") || msg.contains("255"),
            "error message should contain the code: {msg}"
        );
    }

    // r[verify cram.slice.ref_bounds_warning+2]
    #[test]
    fn ref_base_at_warns_on_out_of_bounds() {
        // When requesting a base beyond the reference length, ref_base_at should
        // return b'N' and log a warning.
        let reference = b"ACGT";
        let mut warned = false;

        // In-bounds access should return the actual base
        let base = ref_base_at(reference, 0, &mut warned);
        assert_eq!(base, b'A');
        assert!(!warned, "should not warn for in-bounds access");

        // Out-of-bounds access should return b'N' and set warned flag
        let base = ref_base_at(reference, 10, &mut warned);
        assert_eq!(base, b'N');
        assert!(warned, "should warn for out-of-bounds access");

        // Second out-of-bounds access should still return b'N' but not re-warn
        // (warned is already true, so the warn! is skipped)
        let base = ref_base_at(reference, 20, &mut warned);
        assert_eq!(base, b'N');
    }

    // r[verify cram.slice.validated_lengths]
    #[test]
    fn negative_itf8_num_content_ids_rejected() {
        // Build a slice header where num_content_ids encodes -1 (0xFFFFFFFF as ITF8).
        // Fields before num_content_ids: ref_seq_id, alignment_start, alignment_span,
        // num_records, record_counter (ltf8), num_blocks — all as ITF8/LTF8.
        let mut data = vec![
            0x00, // ref_seq_id = 0 (1 byte ITF8)
            0x01, // alignment_start = 1 (1 byte ITF8)
            0x01, // alignment_span = 1 (1 byte ITF8)
            0x00, // num_records = 0 (1 byte ITF8)
            0x00, // record_counter = 0 (1 byte LTF8)
            0x00, // num_blocks = 0 (1 byte ITF8)
        ];
        // num_content_ids = -1 (5 byte ITF8: 0xFF 0xFF 0xFF 0xFF 0xFF)
        data.extend_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF]);

        let result = SliceHeader::parse(&data);
        assert!(
            matches!(result, Err(CramError::InvalidLength { value: -1 })),
            "expected InvalidLength error for negative num_content_ids, got: {result:?}"
        );
    }
}
