//! Decode CRAM slices into records. Reads data series blocks, applies the encodings from
//! [`CompressionHeader`], and decodes a slice's records into a `DecodedSlice`, from
//! which a fetch copies the ones it wants into a [`RecordStore`].

// See rans.rs: lazy `ok_or_else(|| CramError::...)` keeps error construction and its
// `drop_in_place<CramError>` off the per-record path.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use core::range::Range;
use std::ops::Neg;

use super::{
    block::{self, Block, ContentType},
    codec_scratch::CodecScratch,
    compression_header::CompressionHeader,
    container::RefSeqId,
    encoding::{ByteArrayEncoding, ByteEncoding, DecodeContext, ExternalCursor, IntEncoding},
    reader::CramError,
    varint,
};
use crate::bam::{
    BamHeader, MateFields,
    cigar::CigarOp,
    record::DecodeError,
    record_store::{CustomizeRecordStore, RecordStore},
};
use crate::reader::{FetchedRecord, Tid};
use crate::sam::encode::{AlignmentFields, Qualities, Seq};
use seqair_types::{BamFlags, Base, IupacBase, MappingQuality, Pos0, Pos1, SmallVec, SmolStr};
use std::num::NonZeroI32;
use tracing::warn;

// r[impl cram.stream.iupac]
/// A base as a decode rebuilds SEQ: [`Base`] for the indexed reader, whose
/// records feed the pileup, [`IupacBase`] for the stream, which keeps every
/// code htslib keeps.
pub(crate) trait SeqBase: Copy + Default + PartialEq + std::fmt::Debug + 'static {
    /// A reference base as the decode keeps the slice's reference.
    type Ref: Copy + std::fmt::Debug + 'static;
    /// `N`, for a read base past the reference's end.
    const N: Self;
    /// `N`, for a reference base past its end.
    const REF_N: Self::Ref;
    /// Whether decodes of this kind may generate `MD`/`NM`; `false` lets
    /// the indexed reader's decode compile without the bookkeeping.
    const MD_NM: bool;
    /// A read base the CRAM stores as a byte.
    fn from_ascii(byte: u8) -> Self;
    /// A stretch of the reference, converted once per slice.
    fn reference(ascii: &[u8]) -> Vec<Self::Ref>;
    /// Implicit matches: the reference bases, as read bases.
    fn extend_from_ref(reference: &[Self::Ref], bases: &mut Vec<Self>);
    /// The upper-case letter, for substitutions and `MD`.
    fn ref_ascii(base: Self::Ref) -> u8;
    /// Whether a reference base is `N`, a mismatch in a generated `MD`.
    fn ref_is_n(base: Self::Ref) -> bool;
}

impl SeqBase for Base {
    type Ref = Base;
    const N: Self = Base::Unknown;
    const REF_N: Base = Base::Unknown;
    const MD_NM: bool = false;

    #[inline]
    fn from_ascii(byte: u8) -> Self {
        Base::from(byte)
    }

    fn reference(ascii: &[u8]) -> Vec<Base> {
        Base::from_ascii_vec(ascii.to_vec())
    }

    #[inline]
    fn extend_from_ref(reference: &[Base], bases: &mut Vec<Self>) {
        bases.extend_from_slice(reference);
    }

    #[inline]
    fn ref_ascii(base: Base) -> u8 {
        base as u8
    }

    #[inline]
    fn ref_is_n(base: Base) -> bool {
        base == Base::Unknown
    }
}

/// The stream keeps the reference as upper-case ASCII, as htslib does: its
/// `MD` counts a literal `N` only, not the bytes `bam_set1` reads as `N` (an
/// `embed_ref=2` consensus can hold a NUL).
impl SeqBase for IupacBase {
    type Ref = u8;
    const N: Self = IupacBase::N;
    const REF_N: u8 = b'N';
    const MD_NM: bool = true;

    #[inline]
    fn from_ascii(byte: u8) -> Self {
        IupacBase::from_ascii(byte)
    }

    fn reference(ascii: &[u8]) -> Vec<u8> {
        ascii.to_ascii_uppercase()
    }

    fn extend_from_ref(reference: &[u8], bases: &mut Vec<Self>) {
        bases.extend(reference.iter().map(|&b| IupacBase::from_ascii(b)));
    }

    /// A byte that is no IUPAC letter is written as `N`, where htslib would
    /// copy it into `MD` as it is.
    #[inline]
    fn ref_ascii(base: u8) -> u8 {
        IupacBase::from_ascii(base).to_ascii()
    }

    #[inline]
    fn ref_is_n(base: u8) -> bool {
        base == b'N'
    }
}

/// Mate cross-reference info collected during CRAM record decoding.
/// Used to reconstruct `template_len` for attached/downstream mates.
struct SliceMateInfo {
    /// 0-based alignment start position.
    pos: i64,
    /// 0-based exclusive alignment end position.
    end_pos: i64,
    /// Absolute index of the mate record within this slice, or -1 if detached/no mate.
    mate_line: i32,
    /// Index in the decoded slice's records, or `None` for a record of
    /// another reference in a multi-reference slice.
    slot: Option<u32>,
    /// BAM flags (needed for READ1/READ2 tie-breaking).
    bam_flags: u16,
    /// The record's reference (needed for `next_ref_id` resolution).
    tid: Option<Tid>,
}

// r[impl cram.slice_cache]
/// A slice's records for one reference, decoded once and independent of any
/// query: every record the slice holds for that reference, with TLEN and
/// attached-mate fields resolved over the whole slice. A fetch copies the
/// records it wants out of it with [`copy_into`](Self::copy_into).
#[derive(Debug, Default)]
pub(crate) struct DecodedSlice<B: SeqBase = Base> {
    records: Vec<SliceRecord>,
    names: Vec<u8>,
    cigar: Vec<CigarOp>,
    bases: Vec<B>,
    qual: Vec<u8>,
    aux: Vec<u8>,
}

/// One record of a [`DecodedSlice`]; its variable-length fields are ranges
/// into the slice's slabs.
#[derive(Debug, Clone, Copy)]
struct SliceRecord {
    pos: Pos0,
    end_pos: Pos0,
    /// An unmapped read whose `AP` is below 1: BAM position `-1`. Its `pos`
    /// is 0 so region filters see a position.
    unplaced: bool,
    flags: BamFlags,
    mapq: MappingQuality,
    matching_bases: u32,
    indel_bases: u32,
    tid: Option<Tid>,
    mate_fields: MateFields,
    name: Range<usize>,
    cigar: Range<usize>,
    seq: Range<usize>,
    qual: Range<usize>,
    aux: Range<usize>,
}

impl<B: SeqBase> DecodedSlice<B> {
    /// Heap bytes the slice holds, for the cache's budget.
    pub(crate) fn heap_bytes(&self) -> usize {
        [
            self.records.capacity().saturating_mul(size_of::<SliceRecord>()),
            self.names.capacity(),
            self.cigar.capacity().saturating_mul(size_of::<CigarOp>()),
            self.bases.capacity().saturating_mul(size_of::<B>()),
            self.qual.capacity(),
            self.aux.capacity(),
        ]
        .into_iter()
        .fold(0, usize::saturating_add)
    }

    /// Append a record whose variable-length fields are the given slices.
    #[expect(clippy::too_many_arguments, reason = "one argument per record field")]
    fn push(
        &mut self,
        pos: Pos0,
        end_pos: Pos0,
        unplaced: bool,
        flags: BamFlags,
        mapq: MappingQuality,
        matching_bases: u32,
        indel_bases: u32,
        name: &[u8],
        cigar: &[CigarOp],
        bases: &[B],
        qual: &[u8],
        aux: &[u8],
        tid: Option<Tid>,
        mate_fields: MateFields,
    ) -> Result<u32, CramError> {
        // r[impl cram.record.seq_qual_lengths]
        if qual.len() != bases.len() {
            return Err(CramError::QualityLengthMismatch {
                bases: bases.len(),
                qualities: qual.len(),
            });
        }
        let slot = u32::try_from(self.records.len())
            .map_err(|_| CramError::InvalidLength { value: i32::MAX })?;
        let append = |slab: &mut Vec<u8>, bytes: &[u8]| {
            let start = slab.len();
            slab.extend_from_slice(bytes);
            Range { start, end: slab.len() }
        };
        let name = append(&mut self.names, name);
        let qual = append(&mut self.qual, qual);
        let aux = append(&mut self.aux, aux);
        let cigar_start = self.cigar.len();
        self.cigar.extend_from_slice(cigar);
        let seq_start = self.bases.len();
        self.bases.extend_from_slice(bases);
        self.records.push(SliceRecord {
            pos,
            end_pos,
            unplaced,
            flags,
            mapq,
            matching_bases,
            indel_bases,
            tid,
            mate_fields,
            name,
            cigar: Range { start: cigar_start, end: self.cigar.len() },
            seq: Range { start: seq_start, end: self.bases.len() },
            qual,
            aux,
        });
        Ok(slot)
    }

    /// Give back the slabs' spare capacity: a cached slice lives long.
    fn shrink_to_fit(&mut self) {
        self.records.shrink_to_fit();
        self.names.shrink_to_fit();
        self.cigar.shrink_to_fit();
        self.bases.shrink_to_fit();
        self.qual.shrink_to_fit();
        self.aux.shrink_to_fit();
    }
}

impl DecodedSlice {
    // r[impl cram.fetch_into_customized.push_time+2]
    // r[impl cram.fetch_into_customized.filtered_mate_fields]
    // r[impl interval.overlap_test]
    /// Push the records overlapping `query_start..=query_end` into `store`,
    /// consulting `customize` for each and adding each fetched record to
    /// `counts`.
    ///
    /// Mate fields are copied as the slice resolved them, also for a record
    /// whose mate the fetch leaves out: they describe the template, not the
    /// fetch, as in htslib, noodles and our BAM reader.
    pub(crate) fn copy_into<E: CustomizeRecordStore>(
        &self,
        query_start: Pos0,
        query_end: Pos0,
        store: &mut RecordStore<E::Extra>,
        customize: &mut E,
        counts: &mut crate::reader::FetchCounts,
    ) -> Result<(), CramError> {
        for r in &self.records {
            // Both bounds are inclusive, as are `pos` and `end_pos`: a record
            // is kept iff `pos <= query_end && end_pos >= query_start`. An
            // unmapped record's `end_pos` is its `pos`, as in BAM.
            let wanted = r.pos <= query_end && r.end_pos >= query_start;
            if !wanted {
                continue;
            }
            let idx = store.push_fields(
                r.pos,
                r.end_pos,
                r.flags,
                r.mapq,
                r.matching_bases,
                r.indel_bases,
                field(&self.names, r.name),
                field(&self.cigar, r.cigar),
                field(&self.bases, r.seq),
                field(&self.qual, r.qual),
                field(&self.aux, r.aux),
                r.tid,
                r.mate_fields,
                customize,
            )?;
            counts.add(FetchedRecord {
                flags: r.flags,
                started: r.pos >= query_start,
                kept: idx.is_some(),
            });
        }
        Ok(())
    }
}

impl DecodedSlice<IupacBase> {
    /// Every record, in slice order, as the fields of its BAM record
    /// (`r[cram.stream.raw_records]`).
    pub(crate) fn alignments(&self) -> impl ExactSizeIterator<Item = AlignmentFields<'_>> {
        self.records.iter().map(|r| AlignmentFields {
            tid: r.tid,
            pos: (!r.unplaced).then_some(r.pos),
            flags: r.flags.raw(),
            mapq: r.mapq,
            mate_fields: r.mate_fields,
            qname: field(&self.names, r.name),
            cigar: field(&self.cigar, r.cigar),
            seq: Seq::Iupac(field(&self.bases, r.seq)),
            qual: Qualities::Raw(field(&self.qual, r.qual)),
            aux: field(&self.aux, r.aux),
        })
    }
}

/// A record's field in a [`DecodedSlice`] slab; its ranges were cut from
/// the slab as it grew, so they always fit.
fn field<T>(slab: &[T], range: Range<usize>) -> &[T] {
    slab.get(range).unwrap_or_default()
}

/// Parsed CRAM slice header.
#[derive(Debug)]
pub struct SliceHeader {
    /// The reference(s) the slice's records are on, as numbered in the file;
    /// not checked against the file header until the slice is decoded.
    pub ref_seq_id: RefSeqId,
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

        let ref_seq_id = RefSeqId::from_wire(read_itf8(&mut cursor)?.cast_signed())?;
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

        let md5_bytes = cursor
            .get(..16)
            .ok_or_else(|| CramError::Truncated { context: "slice reference MD5" })?;
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

/// Which of a slice's records a decode keeps; the others are read past
/// without reconstructing their sequence (`r[cram.slice.multi_ref_skip]`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Keep {
    /// The records on one reference: a query's.
    Tid(Tid),
    /// Every record: a stream's (`r[cram.stream.all_records]`).
    All,
    /// None, for collecting where the records lie.
    Nothing,
}

impl Keep {
    fn keeps(self, tid: Option<Tid>) -> bool {
        match self {
            Self::Tid(t) => tid == Some(t),
            Self::All => true,
            Self::Nothing => false,
        }
    }
}

/// A stretch of one reference, its bases converted, for the records of a
/// multi-reference slice on it.
#[derive(Debug)]
pub(crate) struct TidWindow<B> {
    pub(crate) tid: Tid,
    /// 0-based position of the first base.
    pub(crate) start: i64,
    pub(crate) bases: Vec<B>,
}

/// The reference bases a slice's records are reconstructed against.
#[derive(Debug, Clone, Copy)]
pub(crate) enum SliceReference<'a, B> {
    /// One window of the slice's reference starting at 0-based `start`;
    /// `None` when the FASTA lacks the contig, which only a slice that does
    /// not embed its reference minds.
    Window { bases: Option<&'a [u8]>, start: i64 },
    /// A window per reference, sorted by tid, for every reference a
    /// multi-reference slice keeps records of (`r[cram.stream.reference_window]`).
    PerTid(&'a [TidWindow<B>]),
}

/// The reference a record is reconstructed against.
enum RecordReference<'a, B> {
    One { bases: &'a [B], start: i64 },
    PerTid(&'a [TidWindow<B>]),
}

impl<B> RecordReference<'_, B> {
    fn for_tid(&self, tid: Option<Tid>) -> (&[B], i64) {
        match self {
            Self::One { bases, start } => (bases, *start),
            Self::PerTid(windows) => tid
                .and_then(|tid| windows.binary_search_by_key(&tid, |w| w.tid).ok())
                .and_then(|i| windows.get(i))
                .map_or((&[][..], 0), |w| (w.bases.as_slice(), w.start)),
        }
    }
}

/// Buffers a reader reuses across slices and records, so decoding does not
/// allocate per record once they have grown.
#[derive(Default)]
pub(crate) struct SliceScratch<B = Base> {
    pub(crate) cigar: Vec<CigarOp>,
    pub(crate) bases: Vec<B>,
    pub(crate) qual: Vec<u8>,
    pub(crate) aux: Vec<u8>,
    /// Read name (qname) bytes.
    pub(crate) name: Vec<u8>,
    /// Per-feature byte payloads (insertion bases, soft-clip bases,
    /// bases-block, quality-block).
    pub(crate) feature_bytes: Vec<u8>,
    /// `(length, op_code)` CIGAR ops before they are packed into
    /// BAM-layout `CigarOp`s.
    pub(crate) cigar_ops: Vec<(u32, u8)>,
    /// rANS order-1 tables and the bzip2 decoder, each allocated on first use.
    pub(crate) codec: CodecScratch,
    pub(crate) md: MdNm,
}

/// Decode a slice's records — those `keep` names — into a [`DecodedSlice`].
///
/// `container_data` starts at the first block after the container header.
/// `slice_offset` is the byte offset from container data start to this slice's header block.
#[expect(
    clippy::too_many_arguments,
    reason = "CRAM slice decoding requires all compression header, data, and offset parameters"
)]
pub(crate) fn decode_slice<B: SeqBase>(
    ch: &CompressionHeader,
    container_data: &[u8],
    slice_offset: usize,
    reference: SliceReference<'_, B::Ref>,
    header: &BamHeader,
    read_group_ids: &[SmolStr],
    keep: Keep,
    decode_md: bool,
    scratch: &mut SliceScratch<B>,
) -> Result<DecodedSlice<B>, CramError> {
    decode_slice_records(
        ch,
        container_data,
        slice_offset,
        reference,
        header,
        read_group_ids,
        keep,
        decode_md,
        scratch,
    )
    .map(|(slice, _)| slice)
}

// r[impl cram.stream.reference_window]
/// Where the mapped records of a multi-reference slice lie: per reference,
/// sorted by tid, the 0-based half-open span from the lowest start to the
/// highest end. Reconstructs nothing, so it needs no reference.
pub(crate) fn multi_ref_spans<B: SeqBase>(
    ch: &CompressionHeader,
    container_data: &[u8],
    slice_offset: usize,
    header: &BamHeader,
    read_group_ids: &[SmolStr],
    scratch: &mut SliceScratch<B>,
) -> Result<Vec<(Tid, Range<i64>)>, CramError> {
    let (_, infos) = decode_slice_records(
        ch,
        container_data,
        slice_offset,
        SliceReference::PerTid(&[]),
        header,
        read_group_ids,
        Keep::Nothing,
        false,
        scratch,
    )?;
    let mut spans: Vec<(Tid, Range<i64>)> = Vec::new();
    for info in infos.iter().filter(|i| !BamFlags::from(i.bam_flags).is_unmapped()) {
        let Some(tid) = info.tid else { continue };
        match spans.binary_search_by_key(&tid, |(t, _)| *t) {
            Ok(i) => {
                if let Some((_, span)) = spans.get_mut(i) {
                    span.start = span.start.min(info.pos);
                    span.end = span.end.max(info.end_pos);
                }
            }
            Err(i) => spans.insert(i, (tid, Range { start: info.pos, end: info.end_pos })),
        }
    }
    Ok(spans)
}

#[expect(
    clippy::too_many_arguments,
    reason = "CRAM slice decoding requires all compression header, data, and offset parameters"
)]
fn decode_slice_records<B: SeqBase>(
    ch: &CompressionHeader,
    container_data: &[u8],
    slice_offset: usize,
    reference: SliceReference<'_, B::Ref>,
    header: &BamHeader,
    read_group_ids: &[SmolStr],
    keep: Keep,
    decode_md: bool,
    scratch: &mut SliceScratch<B>,
) -> Result<(DecodedSlice<B>, Vec<SliceMateInfo>), CramError> {
    let slice_data = container_data
        .get(slice_offset..)
        .ok_or_else(|| CramError::Truncated { context: "slice offset" })?;

    // Parse slice header block
    // The scratch lives on the reader so every block of every slice reuses
    // one set of codec tables instead of allocating its own.
    let (slice_header_block, mut pos) = block::parse_block_with(slice_data, &mut scratch.codec)?;
    if slice_header_block.content_type != ContentType::SliceHeader {
        return Err(CramError::ExpectedSliceHeader { found: slice_header_block.content_type });
    }
    let sh = SliceHeader::parse(&slice_header_block.data)?;
    let slice_ref = sh.ref_seq_id.resolve(header)?;

    // r[impl cram.edge.missing_reference+2]
    let (reference_seq, ref_start_0based): (&[u8], i64) = match reference {
        SliceReference::Window { bases: Some(seq), start } => (seq, start),
        SliceReference::Window { bases: None, start } if sh.embedded_reference >= 0 => (&[], start),
        SliceReference::Window { bases: None, .. } => {
            let tid = match (slice_ref, keep) {
                (RefSeqId::Ref(t), _) | (_, Keep::Tid(t)) => Some(t),
                _ => None,
            };
            return Err(CramError::MissingReference {
                contig: tid.and_then(|t| header.target_name(t)).unwrap_or("?").into(),
            });
        }
        SliceReference::PerTid(_) => (&[], 0),
    };

    // r[impl cram.edge.unknown_read_names]
    if !ch.preservation.read_names_included {
        warn!(
            "CRAM file has RN=false (read names not stored). Mate-pair overlap dedup will be unreliable."
        );
    }

    // r[impl cram.edge.empty_slice]
    if sh.num_records == 0 {
        return Ok((DecodedSlice::default(), Vec::new()));
    }

    // r[impl cram.edge.reference_mismatch]
    // A slice carrying its own reference is not making a claim about the FASTA.
    // Under `embed_ref=2` htslib embeds a *consensus* computed from the reads
    // and MD5s that, so the digest matches the external reference only where
    // the two happen to agree — which is to say, not at low coverage and not
    // wherever the reads carry a real difference. Checking it there rejects
    // perfectly good files.
    let claims_the_fasta = sh.embedded_reference < 0;
    if slice_ref != RefSeqId::Multi
        && claims_the_fasta
        && sh.reference_md5 != [0u8; 16]
        && !reference_seq.is_empty()
    {
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
        if let Some(slice_bases) = reference_seq.get(slice_ref_start..slice_ref_end) {
            let computed = md5_uppercase_streaming(slice_bases);
            if computed != sh.reference_md5 {
                let contig = match slice_ref {
                    RefSeqId::Ref(t) => header.target_name(t),
                    _ => None,
                };
                return Err(CramError::ReferenceMd5Mismatch {
                    contig: contig.unwrap_or("?").into(),
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
            slice_data.get(pos..).ok_or_else(|| CramError::Truncated { context: "slice block" })?;
        let (blk, consumed) = block::parse_block_with(remaining, &mut scratch.codec)?;
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
    // Converted once, so reconstruction copies runs of reference bases
    // instead of converting base by base.
    let effective_ref = B::reference(effective_ref);
    // For embedded reference, ref_start is the slice's alignment_start
    let effective_ref_start = if sh.embedded_reference >= 0 {
        Pos1::try_from(sh.alignment_start.max(1))
            .map_err(|_| CramError::InvalidPosition { value: i64::from(sh.alignment_start) })?
            .to_zero_based()
            .as_i64()
    } else {
        ref_start_0based
    };
    let record_reference = match reference {
        SliceReference::PerTid(windows) => RecordReference::PerTid(windows),
        SliceReference::Window { .. } => {
            RecordReference::One { bases: &effective_ref, start: effective_ref_start }
        }
    };

    let core_data = core_block.ok_or_else(|| CramError::MissingCoreDataBlock)?;
    // r[impl cram.stream.md_nm]
    // htslib has no reference for a slice of a file written without one
    // (`RR` false) unless the slice embeds it.
    let decode_md = B::MD_NM
        && decode_md
        && (ch.preservation.reference_required || sh.embedded_reference >= 0);

    order_external_blocks(ch, &mut external_blocks);
    let mut ctx = DecodeContext::new(&core_data.data, external_blocks);
    let tag_lines = resolve_tag_lines(ch);

    // Decode records, collecting mate cross-reference info for TLEN reconstruction.
    let mut alignment_pos = i64::from(sh.alignment_start);
    let num_records = usize::try_from(sh.num_records)
        .map_err(|_| CramError::InvalidLength { value: sh.num_records })?;
    let mut mate_infos: Vec<SliceMateInfo> = Vec::with_capacity(num_records);
    let mut out = DecodedSlice::default();

    for record_index in 0..num_records {
        let mate_info = decode_record(
            ch,
            &tag_lines,
            &mut ctx,
            &mut alignment_pos,
            &record_reference,
            header,
            slice_ref,
            read_group_ids,
            keep,
            decode_md,
            &mut out,
            scratch,
            record_index,
        )?;
        mate_infos.push(mate_info);
    }

    // Post-process: resolve template_len for attached/downstream mates.
    // Follow mate_line chains to compute alignment span, then assign signed TLEN.
    resolve_mate_tlen(&mate_infos, &mut out);
    out.shrink_to_fit();
    Ok((out, mate_infos))
}

// r[impl cram.perf.external_order]
/// Orders a slice's external blocks by when a record first reads them —
/// the data series in decode order, then the tags — so that
/// [`DecodeContext`]'s linear content-id scan finds the per-record series
/// in a step or two instead of wading through every tag's block.
fn order_external_blocks(ch: &CompressionHeader, blocks: &mut [(i32, ExternalCursor)]) {
    let ds = &ch.data_series;
    let mut order: SmallVec<i32, 64> = SmallVec::new();
    for e in [
        &ds.bam_flags,
        &ds.cram_flags,
        &ds.ref_id,
        &ds.read_length,
        &ds.alignment_pos,
        &ds.read_group,
    ] {
        int_ids(e, &mut order);
    }
    byte_array_ids(&ds.read_name, &mut order);
    for e in [
        &ds.mate_flags,
        &ds.next_segment_ref,
        &ds.next_mate_pos,
        &ds.template_size,
        &ds.next_fragment,
        &ds.tag_line,
    ] {
        int_ids(e, &mut order);
    }
    for enc in ch.preservation.tag_dictionary.iter().flatten().filter_map(|entry| {
        let key = (i32::from(entry.tag[0]) << 16)
            | (i32::from(entry.tag[1]) << 8)
            | i32::from(entry.bam_type);
        ch.tag_encodings.get(&key)
    }) {
        byte_array_ids(enc, &mut order);
    }
    int_ids(&ds.feature_count, &mut order);
    byte_ids(&ds.feature_code, &mut order);
    int_ids(&ds.feature_pos, &mut order);
    for e in [&ds.mapping_quality, &ds.deletion_length, &ds.ref_skip, &ds.padding, &ds.hard_clip] {
        int_ids(e, &mut order);
    }
    for e in [&ds.base, &ds.quality_score, &ds.base_sub] {
        byte_ids(e, &mut order);
    }
    for e in [&ds.insertion, &ds.soft_clip, &ds.bases_block, &ds.quality_block] {
        byte_array_ids(e, &mut order);
    }
    blocks.sort_by_key(|(id, _)| order.iter().position(|o| o == id).unwrap_or(usize::MAX));
}

fn int_ids(e: &IntEncoding, order: &mut SmallVec<i32, 64>) {
    if let IntEncoding::External { content_id } = e {
        order.push(*content_id);
    }
}

fn byte_ids(e: &ByteEncoding, order: &mut SmallVec<i32, 64>) {
    if let ByteEncoding::External { content_id } = e {
        order.push(*content_id);
    }
}

fn byte_array_ids(e: &ByteArrayEncoding, order: &mut SmallVec<i32, 64>) {
    match e {
        ByteArrayEncoding::External { content_id }
        | ByteArrayEncoding::ByteArrayStop { content_id, .. } => order.push(*content_id),
        ByteArrayEncoding::ByteArrayLen { len_encoding, val_encoding } => {
            int_ids(len_encoding, order);
            byte_ids(val_encoding, order);
        }
        ByteArrayEncoding::Null => {}
    }
}

/// A tag dictionary line with each tag's BAM head (`tag[0] tag[1] type`)
/// and encoding, looked up once per slice; a tag without an encoding is
/// left out, as the per-record lookup skipped it.
type ResolvedTagLine<'a> = SmallVec<([u8; 3], &'a ByteArrayEncoding), 8>;

// r[impl cram.perf.tag_lines_resolved]
/// The compression header's tag lines with their encodings resolved, so a
/// record's tags cost no hash lookups.
fn resolve_tag_lines(ch: &CompressionHeader) -> Vec<ResolvedTagLine<'_>> {
    ch.preservation
        .tag_dictionary
        .iter()
        .map(|line| {
            line.iter()
                .filter_map(|entry| {
                    let key = (i32::from(entry.tag[0]) << 16)
                        | (i32::from(entry.tag[1]) << 8)
                        | i32::from(entry.bam_type);
                    let enc = ch.tag_encodings.get(&key)?;
                    Some(([entry.tag[0], entry.tag[1], entry.bam_type], enc))
                })
                .collect()
        })
        .collect()
}

/// Decode a single CRAM record from the decode context into `out`, unless
/// it belongs to another reference of a multi-reference slice.
///
/// Returns the record's mate cross-reference info; its `slot` is `None`
/// for a record left out.
// r[impl cram.record.decode_order]
#[expect(
    clippy::too_many_arguments,
    reason = "CRAM record decoding requires compression header, slice header, context, and scratch parameters"
)]
fn decode_record<B: SeqBase>(
    ch: &CompressionHeader,
    tag_lines: &[ResolvedTagLine<'_>],
    ctx: &mut DecodeContext<'_>,
    prev_alignment_pos: &mut i64,
    reference: &RecordReference<'_, B::Ref>,
    header: &BamHeader,
    slice_ref: RefSeqId<Tid>,
    read_group_ids: &[SmolStr],
    keep: Keep,
    decode_md: bool,
    out: &mut DecodedSlice<B>,
    scratch: &mut SliceScratch<B>,
    record_index: usize,
) -> Result<SliceMateInfo, CramError> {
    let ds = &ch.data_series;
    let SliceScratch {
        cigar: cigar_buf,
        bases: bases_buf,
        qual: qual_buf,
        aux: aux_buf,
        name: name_buf,
        feature_bytes: feature_byte_buf,
        cigar_ops: cigar_ops_buf,
        codec: _,
        md: md_buf,
    } = scratch;

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
    // r[impl cram.record.ref_ids]
    // 3. RI (ref ID — only for multi-ref slices)
    let record_ref_id = match slice_ref {
        RefSeqId::Multi => Tid::from_wire(ds.ref_id.decode(ctx)?, header)
            .map_err(|source| DecodeError::InvalidTid { source })?,
        RefSeqId::Ref(t) => Some(t),
        RefSeqId::Unmapped => None,
    };
    // A record the decode does not want: for a query, in a multi-reference
    // slice, one of another reference or none. A query decodes only slices
    // that list its reference, so in any other slice every record is on it.
    let foreign_tid = !keep.keeps(record_ref_id);

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
    // Convert from 1-based to 0-based. Unmapped/unplaced reads carry AP <= 0
    // (htslib uses 0 for "no position"). They have no real coordinate, so map
    // them to 0-based 0 rather than erroring — the unmapped path below still
    // pushes them to the store and lets the region filter decide. A *mapped*
    // read must always have a valid (>= 1) AP, so it still errors.
    let unplaced = bam_flags.is_unmapped() && alignment_pos < 1;
    let pos_0based = if unplaced {
        Pos0::ZERO
    } else {
        Pos1::try_from(alignment_pos)
            .map(|p| p.to_zero_based())
            .map_err(|_| super::reader::CramError::InvalidPosition { value: alignment_pos })?
    };

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
    let mut mate_fields = MateFields::default();
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
        // r[impl cram.record.ref_ids]
        mate_fields.tid = Tid::from_wire(ds.next_segment_ref.decode(ctx)?, header)
            .map_err(|source| DecodeError::InvalidMateTid { source })?;
        // CRAM `NP` is 1-based per the spec, NP=0 is "no mate position";
        // matching htslib's `cr->mate_pos - 1` at the BAM emit
        // (cram_decode.c:3121).
        // r[impl cram.record.mate_detached]
        mate_fields.pos = match ds.next_mate_pos.decode(ctx)? {
            0 => None,
            np => Some(
                Pos1::try_from(i64::from(np))
                    .map_err(|_| DecodeError::InvalidMatePosition { value: np })?
                    .to_zero_based(),
            ),
        };
        mate_fields.template_len = NonZeroI32::new(ds.template_size.decode(ctx)?);
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
    // r[impl cram.stream.md_nm]
    let mut stored_md = false;
    let mut stored_nm = false;
    if let Some(line) = tag_lines.get(tag_line_idx) {
        for (head, enc) in line {
            stored_md |= head.starts_with(b"MD");
            stored_nm |= head.starts_with(b"NM");
            // Serialize to BAM binary aux format: tag[0] tag[1] type value_bytes.
            // Decoding the value directly into `aux_buf` skips the
            // per-tag `Vec<u8>` allocation that the old `decode` API
            // forced.
            let start = aux_buf.len();
            aux_buf.extend_from_slice(head);
            enc.decode_into(ctx, aux_buf)?;
            // r[impl cram.record.cf_tag]
            if *head == *b"cFC" && aux_buf.len() == start.wrapping_add(4) {
                let flags = aux_buf.last().copied().unwrap_or(0);
                stored_md |= flags & 1 != 0;
                stored_nm |= flags & 2 != 0;
                aux_buf.truncate(start);
            }
        }
    }

    // r[impl cram.record.rg_tag]
    // Add RG tag if read group >= 0. RG IDs are pre-parsed from the header
    // once at open time (cached on `CramShared::read_group_ids`); used to
    // be a per-record header text scan + String allocation, which dominated
    // `decode_record` profiles.
    // Appended last, after any generated MD/NM, as `cram_to_bam` does.
    let rg_id = usize::try_from(read_group).ok().and_then(|i| read_group_ids.get(i));
    let append_rg = |aux_buf: &mut Vec<u8>| {
        if let Some(rg_id) = rg_id {
            aux_buf.extend_from_slice(b"RGZ");
            aux_buf.extend_from_slice(rg_id.as_bytes());
            aux_buf.push(0); // null terminator for Z-type
        }
    };

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

        // r[impl cram.index.multi_ref_slices]
        // r[impl cram.slice.multi_ref_skip]
        // Multi-ref slices interleave records from multiple references; for a
        // single-tid query the foreign-tid records are guaranteed to be
        // discarded, so don't reconstruct their bases/CIGAR. We still need to
        // drain every per-record data series in lockstep (they're independent
        // streams) and to compute `ref_consumed` for the SliceMateInfo.end_pos
        // used by `resolve_mate_tlen`.
        //
        // Reconstruction sets matching_bases / indel_bases for the kept-record
        // push later; for foreign-tid we don't push, so 0 is fine.
        let generate = GenerateMdNm {
            md: B::MD_NM && decode_md && !seq_unknown && !stored_md,
            nm: B::MD_NM && decode_md && !seq_unknown && !stored_nm,
        };
        let (ref_consumed, matching_bases, indel_bases) = if foreign_tid {
            let span =
                scan_features_for_refspan(ch, ctx, feature_count, read_length, feature_byte_buf)?;
            (span, 0, 0)
        } else {
            let (reference_seq, ref_start_0based) = reference.for_tid(record_ref_id);
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
                generate.any().then(|| {
                    md_buf.start();
                    &mut *md_buf
                }),
            )?;
            if generate.any() {
                md_buf.finish_into(generate, aux_buf);
            }
            (result.ref_consumed, result.matching_bases, result.indel_bases)
        };

        // MQ + quality come after features regardless. They must still be
        // read for foreign-tid records to keep the per-stream cursors
        // synchronized for subsequent records in the slice.
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "mapping quality is 0..=255 encoded as ITF8 i32; only lower 8 bits are valid"
        )]
        let mapq = MappingQuality::from_byte(ds.mapping_quality.decode(ctx)? as u8);

        // r[impl cram.record.quality]
        // Quality scores. Bulk decode for External (single external-block
        // lookup + one memcpy of `read_length` bytes) instead of the
        // per-byte `decode(ctx)` loop that dominated profiles.
        if quality_as_array {
            if foreign_tid {
                // Drain into a discardable scratch — qual_buf stays empty.
                feature_byte_buf.clear();
                ds.quality_score.decode_n_into(ctx, read_length, feature_byte_buf)?;
            } else {
                ds.quality_score.decode_n_into(ctx, read_length, qual_buf)?;
            }
        } else if !foreign_tid {
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
        // r[impl cram.edge.seq_unknown+3]
        // The bases and `RL` qualities were decoded only to keep the series
        // in step; SEQ and QUAL are `*`.
        if seq_unknown {
            bases_buf.clear();
            qual_buf.clear();
        }

        // r[impl cram.record.end_pos]
        let end_pos_raw = pos_0based.as_i64().wrapping_add(i64::from(ref_consumed));
        let end_pos_inclusive_raw =
            if ref_consumed == 0 { pos_0based.as_i64() } else { end_pos_raw.wrapping_sub(1) };
        let end_pos = Pos0::try_from(end_pos_inclusive_raw).map_err(|_| {
            super::reader::CramError::InvalidPosition { value: end_pos_inclusive_raw }
        })?;

        append_rg(aux_buf);
        if foreign_tid {
            return Ok(SliceMateInfo {
                pos: pos_0based.as_i64(),
                end_pos: end_pos_raw,
                mate_line,
                slot: None,
                bam_flags: raw_flags,
                tid: record_ref_id,
            });
        }

        let slot = out.push(
            pos_0based,
            end_pos,
            false,
            bam_flags,
            mapq,
            matching_bases,
            indel_bases,
            name_buf,
            cigar_buf,
            bases_buf,
            qual_buf,
            aux_buf,
            record_ref_id,
            mate_fields,
        )?;

        return Ok(SliceMateInfo {
            pos: pos_0based.as_i64(),
            end_pos: end_pos_raw,
            mate_line,
            slot: Some(slot),
            bam_flags: raw_flags,
            tid: record_ref_id,
        });
    }

    append_rg(aux_buf);
    // Unmapped read
    for _ in 0..read_length {
        let base_byte = ds.base.decode(ctx)?;
        bases_buf.push(B::from_ascii(base_byte));
    }

    // Quality. See the mapped-read branch above for the bulk-decode
    // rationale.
    if quality_as_array {
        ds.quality_score.decode_n_into(ctx, read_length, qual_buf)?;
    } else {
        qual_buf.resize(read_length, 0xFF);
    }
    // r[impl cram.edge.seq_unknown+3]
    if seq_unknown {
        bases_buf.clear();
        qual_buf.clear();
    }

    // r[impl cram.edge.unmapped_reads+2]
    // Unmapped read — kept so filter_raw can decide.
    // htslib compat: end_pos = pos (ignore CIGAR for unmapped reads).
    //
    // Apply the same foreign-tid check as the mapped path, so unmapped reads
    // placed on other references of a multi-ref slice stay out — and so do
    // unplaced ones (reference -1), which are on no reference at all.
    let slot = if foreign_tid {
        None
    } else {
        Some(out.push(
            pos_0based,
            pos_0based,
            unplaced,
            bam_flags,
            MappingQuality::from_byte(0), // unmapped reads have MAPQ=0 per spec
            0,                            // matching_bases
            0,                            // indel_bases
            name_buf,
            &[], // empty cigar — unmapped reads don't have features decoded
            bases_buf,
            qual_buf,
            aux_buf,
            record_ref_id,
            mate_fields,
        )?)
    };
    Ok(SliceMateInfo {
        pos: pos_0based.as_i64(),
        end_pos: pos_0based.as_i64(),
        mate_line,
        slot,
        bam_flags: raw_flags,
        tid: record_ref_id,
    })
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
fn resolve_mate_tlen<B: SeqBase>(infos: &[SliceMateInfo], out: &mut DecodedSlice<B>) {
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
        // Each record's mate is the next entry in the chain (last wraps to
        // first).
        let rest_tlen = first_tlen.neg();
        for (ci, &j) in chain.iter().enumerate() {
            let Some(slot) = infos[j].slot else { continue };
            let mate = &infos[chain[if ci + 1 < chain.len() { ci + 1 } else { 0 }]];
            let Some(r) = out.records.get_mut(slot as usize) else { continue };
            r.mate_fields.template_len =
                NonZeroI32::new(if ci == 0 { first_tlen } else { rest_tlen });
            r.mate_fields.tid = mate.tid;
            // An unmapped mate (`pos < 0`) has no position.
            r.mate_fields.pos = Pos0::try_from(mate.pos).ok();
        }
    }
}

// r[impl cram.stream.md_nm]
/// htslib's `MD`/`NM` as `cram_decode_seq` builds them while it applies a
/// read's features; [`Self::start`] begins a read.
#[derive(Debug, Default)]
pub(crate) struct MdNm {
    md: Vec<u8>,
    nm: u32,
    /// The match run not yet written; -1 once the read ran off the
    /// reference, after which htslib writes no more numbers.
    dist: i64,
}

impl MdNm {
    fn start(&mut self) {
        self.md.clear();
        self.nm = 0;
        self.dist = 0;
    }

    /// `BLOCK_APPEND_UINT(md_dist)` as htslib guards it.
    fn flush(&mut self) {
        if self.dist >= 0 {
            self.md.extend_from_slice(itoa::Buffer::new().format(self.dist).as_bytes());
        }
    }

    /// `add_md_char`: the run, then the reference base, unguarded.
    fn mismatch(&mut self, ref_base: u8) {
        self.md.extend_from_slice(itoa::Buffer::new().format(self.dist).as_bytes());
        self.md.push(ref_base);
        self.dist = 0;
        self.nm = self.nm.saturating_add(1);
    }

    /// `n` implicit matches over `reference` (`None` past its end); a
    /// reference `N` is a mismatch.
    fn matches<B: SeqBase>(&mut self, reference: Option<&[B::Ref]>, n: usize) {
        let n = i64::try_from(n).unwrap_or(i64::MAX);
        match reference {
            Some(bases) => {
                let mut run = 0i64;
                for &b in bases {
                    if B::ref_is_n(b) {
                        self.dist = self.dist.saturating_add(run);
                        run = 0;
                        self.mismatch(b'N');
                    } else {
                        run = run.saturating_add(1);
                    }
                }
                self.dist = self.dist.saturating_add(run);
            }
            None if self.dist >= 0 => self.dist = self.dist.saturating_add(n),
            None => {}
        }
    }

    /// `X`: a substitution of the reference base, `None` past its end.
    fn substitution<B: SeqBase>(&mut self, ref_base: Option<B::Ref>) {
        match ref_base {
            Some(b) => self.mismatch(B::ref_ascii(b)),
            None => {
                self.flush();
                self.dist = -1;
            }
        }
    }

    /// `B`: a read base stored as it is, a mismatch wherever it lies.
    fn read_base<B: SeqBase>(&mut self, ref_base: Option<B::Ref>) {
        self.flush();
        match ref_base {
            Some(b) => {
                self.md.push(B::ref_ascii(b));
                self.nm = self.nm.saturating_add(1);
                self.dist = 0;
            }
            None => self.dist = -1,
        }
    }

    /// `b`: read bases stored as they are, each a mismatch.
    fn read_bases<B: SeqBase>(&mut self, reference: &[B::Ref], start: usize, n: usize) {
        self.flush();
        for i in 0..n {
            if i > 0 {
                self.md.push(b'0');
            }
            let Some(&b) = reference.get(start.saturating_add(i)) else {
                self.dist = -1;
                break;
            };
            self.md.push(B::ref_ascii(b));
            self.nm = self.nm.saturating_add(1);
        }
        // htslib resets the run even when the bases ran off the reference.
        self.dist = 0;
    }

    fn deletion<B: SeqBase>(&mut self, reference: &[B::Ref], start: usize, n: usize) {
        self.flush();
        let Some(have) = reference.get(start..) else {
            self.dist = -1;
            return;
        };
        let deleted = have.get(..n).unwrap_or(have);
        self.md.push(b'^');
        self.md.extend(deleted.iter().map(|&b| B::ref_ascii(b)));
        self.nm = self.nm.saturating_add(u32::try_from(deleted.len()).unwrap_or(u32::MAX));
        if deleted.len() == n {
            self.dist = 0;
        } else {
            self.md.push(b'0');
            self.dist = -1;
        }
    }

    fn insertion(&mut self, n: usize) {
        self.nm = self.nm.saturating_add(u32::try_from(n).unwrap_or(u32::MAX));
    }

    /// Append the tags that are wanted, `MD:Z` then `NM` as the smallest
    /// unsigned type that holds it.
    fn finish_into(&mut self, want: GenerateMdNm, aux: &mut Vec<u8>) {
        self.flush();
        if want.md {
            aux.extend_from_slice(b"MDZ");
            aux.extend_from_slice(&self.md);
            aux.push(0);
        }
        if want.nm {
            aux.extend_from_slice(b"NM");
            if let Ok(nm) = u8::try_from(self.nm) {
                aux.extend_from_slice(&[b'C', nm]);
            } else if let Ok(nm) = u16::try_from(self.nm) {
                aux.push(b'S');
                aux.extend_from_slice(&nm.to_le_bytes());
            } else {
                aux.push(b'I');
                aux.extend_from_slice(&self.nm.to_le_bytes());
            }
        }
    }
}

/// Which of `MD` and `NM` a record gets generated.
#[derive(Debug, Clone, Copy, Default)]
struct GenerateMdNm {
    md: bool,
    nm: bool,
}

impl GenerateMdNm {
    fn any(self) -> bool {
        self.md || self.nm
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
fn ref_base_at<B: SeqBase>(reference: &[B::Ref], index: usize, warned: &mut bool) -> B::Ref {
    match reference.get(index) {
        Some(&b) => b,
        None => {
            warn_ref_short(reference, index, warned);
            B::REF_N
        }
    }
}

fn warn_ref_short<B>(reference: &[B], index: usize, warned: &mut bool) {
    if !*warned {
        warn!(
            index,
            ref_len = reference.len(),
            "reference shorter than expected during CRAM sequence reconstruction; \
             substituting N (may indicate reference/CRAM mismatch)"
        );
        *warned = true;
    }
}

/// Append the `n` reference bases from `start` as read bases — one copy of
/// the already-converted reference instead of a lookup and push per base —
/// with `N` for any past the reference's end, as [`ref_base_at`] gives.
fn extend_from_reference<B: SeqBase>(
    reference: &[B::Ref],
    start: usize,
    n: usize,
    warned: &mut bool,
    bases_buf: &mut Vec<B>,
) {
    let have = reference.get(start..).unwrap_or_default();
    let (have, missing) = if have.len() >= n {
        (have.get(..n).unwrap_or_default(), 0)
    } else {
        (have, n.wrapping_sub(have.len()))
    };
    B::extend_from_ref(have, bases_buf);
    if missing > 0 {
        warn_ref_short(reference, start.wrapping_add(have.len()), warned);
        bases_buf.resize(bases_buf.len().wrapping_add(missing), B::N);
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
fn decode_features_and_reconstruct<B: SeqBase>(
    ch: &CompressionHeader,
    ctx: &mut DecodeContext<'_>,
    feature_count: usize,
    read_length: usize,
    pos_0based: i64,
    slice_start_0based: i64,
    reference_seq: &[B::Ref],
    cigar_buf: &mut Vec<CigarOp>,
    bases_buf: &mut Vec<B>,
    feature_byte_buf: &mut Vec<u8>,
    cigar_ops_buf: &mut Vec<(u32, u8)>,
    md: Option<&mut MdNm>,
) -> Result<ReconstructResult, CramError> {
    let ds = &ch.data_series;
    let mut md = md.filter(|_| B::MD_NM);

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

    /// Emit `n` matching bases copied from the reference, advancing both
    /// `read_pos` and `ref_pos`: the fill up to the next feature and the tail.
    #[expect(clippy::too_many_arguments, reason = "shared helper for ref-match emit")]
    fn emit_ref_run<B: SeqBase>(
        reference: &[B::Ref],
        ref_offset: usize,
        n: usize,
        ref_pos: &mut usize,
        read_pos: &mut usize,
        ref_warned: &mut bool,
        bases_buf: &mut Vec<B>,
        cigar_ops: &mut Vec<(u32, u8)>,
        matching_bases: &mut u32,
        md: Option<&mut MdNm>,
    ) {
        if n == 0 {
            return;
        }
        let start = ref_offset.wrapping_add(*ref_pos);
        if let Some(md) = md {
            md.matches::<B>(reference.get(start..start.saturating_add(n)), n);
        }
        extend_from_reference::<B>(reference, start, n, ref_warned, bases_buf);
        let n32 = u32::try_from(n).unwrap_or(u32::MAX);
        push_cigar_op(cigar_ops, n32, 0); // M
        *matching_bases = matching_bases.saturating_add(n32);
        *read_pos = read_pos.wrapping_add(n);
        *ref_pos = ref_pos.wrapping_add(n);
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
        {
            emit_ref_run(
                reference_seq,
                ref_offset,
                feat_target.min(read_length).saturating_sub(read_pos),
                &mut ref_pos,
                &mut read_pos,
                &mut ref_warned,
                bases_buf,
                cigar_ops_buf,
                &mut matching_bases,
                md.as_deref_mut(),
            );
        }
        let ref_index = ref_offset.wrapping_add(ref_pos);

        // Apply the feature inline. Payload-bearing variants decode straight
        // into `feature_byte_buf` (cleared, then filled) — no per-feature alloc.
        // r[related cram.slice.ref_bounds_warning+2]
        match fc {
            // X — substitution: 1 read base, 1 ref base, M cigar.
            b'X' => {
                let bs = ds.base_sub.decode(ctx)?;
                let ref_base = ref_base_at::<B>(
                    reference_seq,
                    ref_offset.wrapping_add(ref_pos),
                    &mut ref_warned,
                );
                let read_base =
                    ch.preservation.substitution_matrix.substitute(B::ref_ascii(ref_base), bs);
                bases_buf.push(B::from_ascii(read_base));
                if let Some(md) = md.as_deref_mut() {
                    md.substitution::<B>(reference_seq.get(ref_index).copied());
                }
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
                    bases_buf.push(B::from_ascii(b));
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "insertion length is bounded by read_length (validated); fits in u32"
                )]
                let len_u32 = len as u32;
                if let Some(md) = md.as_deref_mut() {
                    md.insertion(len);
                }
                push_cigar_op(cigar_ops_buf, len_u32, 1); // I
                indel_bases = indel_bases.saturating_add(len_u32);
                read_pos = read_pos.wrapping_add(len);
            }
            // i — single-base insertion.
            b'i' => {
                let base = ds.base.decode(ctx)?;
                bases_buf.push(B::from_ascii(base));
                if let Some(md) = md.as_deref_mut() {
                    md.insertion(1);
                }
                push_cigar_op(cigar_ops_buf, 1, 1); // I
                indel_bases = indel_bases.saturating_add(1);
                read_pos = read_pos.wrapping_add(1);
            }
            // D — deletion: 0 read, N ref, D cigar.
            b'D' => {
                let len = ds.deletion_length.decode(ctx)? as u32;
                if let Some(md) = md.as_deref_mut() {
                    md.deletion::<B>(reference_seq, ref_index, len as usize);
                }
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
                    bases_buf.push(B::from_ascii(b));
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
                bases_buf.push(B::from_ascii(base));
                if let Some(md) = md.as_deref_mut() {
                    md.read_base::<B>(reference_seq.get(ref_index).copied());
                }
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
                    bases_buf.push(B::from_ascii(b));
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "bases length is bounded by read_length (validated); fits in u32"
                )]
                let len_u32 = len as u32;
                if let Some(md) = md.as_deref_mut() {
                    md.read_bases::<B>(reference_seq, ref_index, len);
                }
                push_cigar_op(cigar_ops_buf, len_u32, 0); // M
                matching_bases = matching_bases.saturating_add(len_u32);
                read_pos = read_pos.wrapping_add(len);
                ref_pos = ref_pos.wrapping_add(len);
            }
            // r[impl cram.feature.quality_only]
            // q — quality block: consumed so the codec stream stays aligned,
            // value discarded. It carries no base and moves nothing.
            // r[unimpl cram.feature.quality_block]
            b'q' => {
                feature_byte_buf.clear();
                ds.quality_block.decode_into(ctx, feature_byte_buf)?;
            }
            // r[impl cram.feature.quality_only]
            // Q — single quality score: same shape as q but no payload.
            // r[unimpl cram.feature.quality_score]
            b'Q' => {
                let _qual = ds.quality_score.decode(ctx)?;
            }
            _ => return Err(CramError::UnknownFeatureCode { feature_code: fc }),
        }
    }

    // Tail: fill reference matches until the read length is reached.
    {
        emit_ref_run(
            reference_seq,
            ref_offset,
            read_length.saturating_sub(read_pos),
            &mut ref_pos,
            &mut read_pos,
            &mut ref_warned,
            bases_buf,
            cigar_ops_buf,
            &mut matching_bases,
            md,
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

/// Scan the feature streams for a record without reconstructing its
/// sequence or CIGAR. Used for records in multi-ref slices whose
/// `record_ref_id` does not match the query `tid` — we still need to
/// advance every per-record data series in lockstep with the other
/// records in the slice, and we still need `ref_consumed` (for
/// `SliceMateInfo.end_pos` used by `resolve_mate_tlen`), but we don't
/// need bases/CIGAR/ref-lookup since the record will be discarded.
///
/// Stream consumption MUST match `decode_features_and_reconstruct`
/// exactly — same `ds.X.decode*` calls in the same order for the same
/// feature codes. Drift here desynchronises every subsequent record in
/// the slice.
/// `r[impl cram.slice.multi_ref_skip]`
fn scan_features_for_refspan(
    ch: &CompressionHeader,
    ctx: &mut DecodeContext<'_>,
    feature_count: usize,
    read_length: usize,
    feature_byte_buf: &mut Vec<u8>,
) -> Result<u32, CramError> {
    let ds = &ch.data_series;
    let mut read_pos = 0usize;
    let mut ref_pos = 0usize;
    let mut feature_read_pos = 0u32;

    for _ in 0..feature_count {
        let fc = ds.feature_code.decode(ctx)?;
        let fp = ds.feature_pos.decode(ctx)? as u32;
        feature_read_pos = feature_read_pos.wrapping_add(fp);
        let feat_target = (feature_read_pos as usize).saturating_sub(1);

        // Catch up read_pos and ref_pos to the feature anchor; matching
        // bases advance both by 1 each.
        let fill = feat_target.min(read_length).saturating_sub(read_pos);
        read_pos = read_pos.wrapping_add(fill);
        ref_pos = ref_pos.wrapping_add(fill);

        match fc {
            b'X' => {
                let _bs = ds.base_sub.decode(ctx)?;
                read_pos = read_pos.wrapping_add(1);
                ref_pos = ref_pos.wrapping_add(1);
            }
            b'I' => {
                feature_byte_buf.clear();
                ds.insertion.decode_into(ctx, feature_byte_buf)?;
                read_pos = read_pos.wrapping_add(feature_byte_buf.len());
            }
            b'i' => {
                let _base = ds.base.decode(ctx)?;
                read_pos = read_pos.wrapping_add(1);
            }
            b'D' => {
                let len = ds.deletion_length.decode(ctx)? as u32;
                ref_pos = ref_pos.wrapping_add(len as usize);
            }
            b'S' => {
                feature_byte_buf.clear();
                ds.soft_clip.decode_into(ctx, feature_byte_buf)?;
                read_pos = read_pos.wrapping_add(feature_byte_buf.len());
            }
            b'H' => {
                let _len = ds.hard_clip.decode(ctx)?;
            }
            b'N' => {
                let len = ds.ref_skip.decode(ctx)? as u32;
                ref_pos = ref_pos.wrapping_add(len as usize);
            }
            b'P' => {
                let _len = ds.padding.decode(ctx)?;
            }
            b'B' => {
                let _base = ds.base.decode(ctx)?;
                let _qual = ds.quality_score.decode(ctx)?;
                read_pos = read_pos.wrapping_add(1);
                ref_pos = ref_pos.wrapping_add(1);
            }
            b'b' => {
                feature_byte_buf.clear();
                ds.bases_block.decode_into(ctx, feature_byte_buf)?;
                let len = feature_byte_buf.len();
                read_pos = read_pos.wrapping_add(len);
                ref_pos = ref_pos.wrapping_add(len);
            }
            // r[impl cram.feature.quality_only]
            b'q' => {
                feature_byte_buf.clear();
                ds.quality_block.decode_into(ctx, feature_byte_buf)?;
            }
            // r[impl cram.feature.quality_only]
            b'Q' => {
                let _qual = ds.quality_score.decode(ctx)?;
            }
            _ => return Err(CramError::UnknownFeatureCode { feature_code: fc }),
        }
    }

    // Tail: matching bases advance read_pos and ref_pos in lockstep up
    // to read_length.
    let tail = read_length.saturating_sub(read_pos);
    ref_pos = ref_pos.wrapping_add(tail);

    #[expect(
        clippy::cast_possible_truncation,
        reason = "ref_pos bounded by reference length which fits in u32; BAM positions are capped at 2^29"
    )]
    Ok(ref_pos as u32)
}

fn read_itf8(cursor: &mut &[u8]) -> Result<u32, CramError> {
    varint::read_itf8_from(cursor)
        .ok_or_else(|| CramError::Truncated { context: "slice header itf8" })
}

/// MD5 of `bytes` after ASCII-uppercasing, computed without allocating
/// an intermediate `Vec<u8>`. Slice references for a CRAM region can be
/// 1 MB+; allocating per slice for verification was a real cost on
/// pileup-heavy workloads.
fn md5_uppercase_streaming(bytes: &[u8]) -> [u8; 16] {
    let mut ctx = md5::Context::new();
    let mut buf = [0u8; 4096];
    for chunk in bytes.chunks(buf.len()) {
        let dst = buf.get_mut(..chunk.len()).expect("chunks(N) yields chunks of length ≤ N");
        for (d, &s) in dst.iter_mut().zip(chunk) {
            *d = s.to_ascii_uppercase();
        }
        ctx.consume(&*dst);
    }
    ctx.finalize().into()
}

fn read_ltf8(cursor: &mut &[u8]) -> Result<u64, CramError> {
    varint::read_ltf8_from(cursor)
        .ok_or_else(|| CramError::Truncated { context: "slice header ltf8" })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::test_util::tid;
    use hegel::prelude::*;

    fn load_test_cram() -> Vec<u8> {
        std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram")).unwrap()
    }

    /// Reference implementation that allocates a Vec — what `md5_uppercase_streaming`
    /// must match for any input.
    fn md5_uppercase_alloc(bytes: &[u8]) -> [u8; 16] {
        let upper: Vec<u8> = bytes.iter().map(|b| b.to_ascii_uppercase()).collect();
        md5::compute(&upper).into()
    }

    // r[verify cram.edge.reference_mismatch]
    // Streaming MD5 must produce identical output to the allocating
    // form across input lengths spanning sub-block, single-block,
    // multi-block, and unaligned-tail cases.
    #[hegel::test]
    fn md5_uppercase_streaming_matches_alloc(tc: TestCase) {
        let bytes = tc.draw(gs::binary().max_size(16_383));
        assert_eq!(md5_uppercase_streaming(&bytes), md5_uppercase_alloc(&bytes));
    }

    #[test]
    fn md5_uppercase_streaming_known_vectors() {
        // Empty input.
        let empty_md5: [u8; 16] = md5::compute(b"").into();
        assert_eq!(md5_uppercase_streaming(b""), empty_md5);

        // Lowercase string uppercases to the canonical MD5 of "ACGT".
        let acgt_md5: [u8; 16] = md5::compute(b"ACGT").into();
        assert_eq!(md5_uppercase_streaming(b"acgt"), acgt_md5);

        // Mixed case + non-letters (Ns are common in references).
        let mixed_md5: [u8; 16] = md5::compute(b"ACGTNACGT").into();
        assert_eq!(md5_uppercase_streaming(b"AcGtNaCgT"), mixed_md5);
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
            SliceReference::Window { bases: Some(&fake_ref), start: ref_start },
            &bam_header,
            &[],
            Keep::Tid(tid(0)),
            false,
            &mut SliceScratch::<Base>::default(),
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

    // r[verify cram.record.ref_ids]
    /// A slice on a reference the header does not have fails, rather than
    /// handing out a `Tid` that names nothing.
    #[test]
    fn a_slice_reference_the_header_lacks_is_an_error() {
        let header = crate::bam::BamHeader::from_sam_text("@SQ\tSN:c0\tLN:10\n").unwrap();
        assert_eq!(RefSeqId::Ref(0).resolve(&header).unwrap(), RefSeqId::Ref(tid(0)));
        assert_eq!(RefSeqId::Unmapped.resolve(&header).unwrap(), RefSeqId::Unmapped);
        assert_eq!(RefSeqId::Multi.resolve(&header).unwrap(), RefSeqId::Multi);
        assert!(matches!(
            RefSeqId::Ref(1).resolve(&header),
            Err(DecodeError::InvalidTid { source: crate::reader::TidError::TidOutOfRange { .. } })
        ));
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
        let reference = &[Base::A, Base::C, Base::G, Base::T];
        let mut warned = false;

        // In-bounds access should return the actual base
        let base = ref_base_at::<Base>(reference, 0, &mut warned);
        assert_eq!(base, Base::A);
        assert!(!warned, "should not warn for in-bounds access");

        // Out-of-bounds access should return b'N' and set warned flag
        let base = ref_base_at::<Base>(reference, 10, &mut warned);
        assert_eq!(base, Base::Unknown);
        assert!(warned, "should warn for out-of-bounds access");

        // Second out-of-bounds access should still return b'N' but not re-warn
        // (warned is already true, so the warn! is skipped)
        let base = ref_base_at::<Base>(reference, 20, &mut warned);
        assert_eq!(base, Base::Unknown);
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

    fn md_of(md: &mut MdNm) -> (String, Vec<u8>) {
        let mut aux = Vec::new();
        md.finish_into(GenerateMdNm { md: true, nm: true }, &mut aux);
        (String::from_utf8(md.md.clone()).unwrap(), aux)
    }

    // r[verify cram.stream.md_nm]
    /// A reference skip leaves the match run going (htslib measured
    /// `5M100N5M` as one run); reference `N`s in a match are mismatches.
    #[test]
    fn md_runs_through_skips_and_counts_reference_n() {
        let reference = Base::from_ascii_vec(b"ACGTANNCGT".to_vec());
        let mut md = MdNm::default();
        md.start();
        md.matches::<Base>(reference.get(0..5), 5);
        // a skip: no call at all
        md.matches::<Base>(reference.get(4..10), 6);
        let (text, aux) = md_of(&mut md);
        assert_eq!(text, "6N0N3");
        assert_eq!(aux.get(aux.len() - 4..), Some(&b"NMC\x02"[..]));
    }

    // r[verify cram.stream.md_nm]
    /// `b` bases are each a mismatch, written with a `0` between them, even
    /// where they equal the reference; `D` writes `^` and the deleted bases.
    #[test]
    fn md_of_bases_blocks_and_deletions() {
        let reference = Base::from_ascii_vec(b"ACGTACGT".to_vec());
        let mut md = MdNm::default();
        md.start();
        md.matches::<Base>(reference.get(0..2), 2);
        md.read_bases::<Base>(&reference, 2, 2);
        md.deletion::<Base>(&reference, 4, 2);
        md.insertion(3);
        md.matches::<Base>(reference.get(6..8), 2);
        let (text, _) = md_of(&mut md);
        assert_eq!(text, "2G0T0^AC2");
        assert_eq!(md.nm, 7);
    }

    // r[verify cram.stream.md_nm]
    /// NM takes the smallest unsigned type that holds it.
    #[test]
    fn nm_is_the_smallest_unsigned_type() {
        for (nm, want) in [
            (255u32, &b"NMC\xff"[..]),
            (256, &b"NMS\x00\x01"[..]),
            (70_000, &b"NMI\x70\x11\x01\x00"[..]),
        ] {
            let mut md = MdNm::default();
            md.start();
            md.nm = nm;
            let mut aux = Vec::new();
            md.finish_into(GenerateMdNm { md: false, nm: true }, &mut aux);
            assert_eq!(aux, want, "{nm}");
        }
    }

    // r[verify cram.record.seq_qual_lengths]
    /// A record is pushed only when it has as many qualities as bases.
    #[hegel::test]
    fn push_requires_a_quality_per_base(tc: TestCase) {
        let n_bases = tc.draw(gs::integers::<usize>().max_value(20));
        let n_qual = tc.draw(gs::integers::<usize>().max_value(20));
        let mut slice = DecodedSlice::<IupacBase>::default();
        let result = slice.push(
            Pos0::ZERO,
            Pos0::ZERO,
            false,
            BamFlags::from(4u16),
            MappingQuality::from_byte(0),
            0,
            0,
            b"r",
            &[],
            &vec![IupacBase::A; n_bases],
            &vec![30; n_qual],
            &[],
            None,
            MateFields::default(),
        );
        if n_bases == n_qual {
            assert_eq!(result.unwrap(), 0);
            assert_eq!(slice.records.len(), 1);
        } else {
            assert!(matches!(
                result,
                Err(CramError::QualityLengthMismatch { bases, qualities })
                    if bases == n_bases && qualities == n_qual
            ));
            assert!(slice.records.is_empty());
        }
    }
}
