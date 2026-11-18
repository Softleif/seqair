//! A slice's data series as the records fill them, and the container that
//! frames them: compression header, slice header, blocks.

use super::CramWriteError;
use super::mates::{self, MateRecord};
use super::{CramVersion, rans, rans_nx16, tok3};
use crate::cram::varint::{write_itf8, write_ltf8};
use crate::io::CompressionLevel;
use crate::reader::Tid;
use rustc_hash::FxHashMap;

/// The data series, in the order a record's decode first reads them. Each is
/// stored `EXTERNAL` in a block of its own, whose content id is its index
/// plus one (the core block is 0).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Series {
    Bf,
    Cf,
    Rl,
    Ap,
    Rg,
    Rn,
    Mf,
    Ns,
    Np,
    Ts,
    Nf,
    Tl,
    Fn,
    Fc,
    Fp,
    Dl,
    Mq,
    Ba,
    Qs,
    Bs,
    In,
    Sc,
    Rs,
    Pd,
    Hc,
    Bb,
    /// `BB`'s lengths, in a block of their own.
    BbLen,
}

const SERIES: [Series; 27] = [
    Series::Bf,
    Series::Cf,
    Series::Rl,
    Series::Ap,
    Series::Rg,
    Series::Rn,
    Series::Mf,
    Series::Ns,
    Series::Np,
    Series::Ts,
    Series::Nf,
    Series::Tl,
    Series::Fn,
    Series::Fc,
    Series::Fp,
    Series::Dl,
    Series::Mq,
    Series::Ba,
    Series::Qs,
    Series::Bs,
    Series::In,
    Series::Sc,
    Series::Rs,
    Series::Pd,
    Series::Hc,
    Series::Bb,
    Series::BbLen,
];

impl Series {
    fn content_id(self) -> i32 {
        (self as i32).saturating_add(1)
    }

    // r[impl cram_writer.data_series]
    /// The series' key and encoding in the data series map; `None` for
    /// `BbLen`, which `BB`'s encoding names.
    fn encoding(self, out: &mut Vec<u8>) -> Option<()> {
        let key: &[u8; 2] = match self {
            Self::Bf => b"BF",
            Self::Cf => b"CF",
            Self::Rl => b"RL",
            Self::Ap => b"AP",
            Self::Rg => b"RG",
            Self::Rn => b"RN",
            Self::Mf => b"MF",
            Self::Ns => b"NS",
            Self::Np => b"NP",
            Self::Ts => b"TS",
            Self::Nf => b"NF",
            Self::Tl => b"TL",
            Self::Fn => b"FN",
            Self::Fc => b"FC",
            Self::Fp => b"FP",
            Self::Dl => b"DL",
            Self::Mq => b"MQ",
            Self::Ba => b"BA",
            Self::Qs => b"QS",
            Self::Bs => b"BS",
            Self::In => b"IN",
            Self::Sc => b"SC",
            Self::Rs => b"RS",
            Self::Pd => b"PD",
            Self::Hc => b"HC",
            Self::Bb => b"BB",
            Self::BbLen => return None,
        };
        out.extend_from_slice(key);
        match self {
            Self::Rn | Self::In | Self::Sc => byte_array_stop(out, 0, self.content_id()),
            Self::Bb => byte_array_len(
                out,
                |p| external(p, Self::BbLen.content_id()),
                |p| {
                    external(p, self.content_id());
                },
            ),
            _ => external(out, self.content_id()),
        }
        Some(())
    }
}

/// An encoding descriptor: its id, then its parameters' length and bytes.
fn encoding(out: &mut Vec<u8>, id: i32, params: impl FnOnce(&mut Vec<u8>)) {
    let mut p = Vec::new();
    params(&mut p);
    write_itf8(out, id);
    write_small_len(out, p.len());
    out.extend_from_slice(&p);
}

fn external(out: &mut Vec<u8>, content_id: i32) {
    encoding(out, 1, |p| write_itf8(p, content_id));
}

fn byte_array_stop(out: &mut Vec<u8>, stop: u8, content_id: i32) {
    encoding(out, 5, |p| {
        p.push(stop);
        write_itf8(p, content_id);
    });
}

fn byte_array_len(
    out: &mut Vec<u8>,
    len: impl FnOnce(&mut Vec<u8>),
    val: impl FnOnce(&mut Vec<u8>),
) {
    encoding(out, 4, |p| {
        len(p);
        val(p);
    });
}

/// A `HUFFMAN` code of one symbol, which takes no bits: a constant.
fn constant(out: &mut Vec<u8>, value: i32) {
    encoding(out, 3, |p| {
        write_itf8(p, 1);
        write_itf8(p, value);
        write_itf8(p, 1);
        write_itf8(p, 0);
    });
}

/// The length of a structure the writer itself bounds (an encoding's
/// parameters, a map of a few dozen entries).
fn write_small_len(out: &mut Vec<u8>, len: usize) {
    debug_assert!(i32::try_from(len).is_ok(), "header structure of {len} bytes");
    write_itf8(out, i32::try_from(len).unwrap_or(i32::MAX));
}

fn len_i32(what: &'static str, len: usize) -> Result<i32, CramWriteError> {
    i32::try_from(len).map_err(|_| CramWriteError::ValueTooLarge { what, value: len as u64 })
}

/// What a slice's records are on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SliceKind {
    /// One reference; `reference` when the FASTA has it and the records are
    /// stored against it.
    Mapped {
        tid: Tid,
        reference: bool,
    },
    Unplaced,
}

/// A tag's values, by its key (`name[0] << 16 | name[1] << 8 | type`).
struct TagBlock {
    key: i32,
    data: Vec<u8>,
}

/// One slice's records, as data series, until it is assembled.
#[derive(Default)]
pub(super) struct SliceBuilder {
    kind: Option<SliceKind>,
    series: [Vec<u8>; SERIES.len()],
    tags: Vec<TagBlock>,
    /// Each record's mate fields; the mate series and `CF` are written
    /// when the slice is complete, once its pairs are known.
    mates: Vec<MateRecord>,
    tag_lines: Vec<Vec<u8>>,
    line_ids: FxHashMap<Vec<u8>, i32>,
    n_records: u32,
    n_bases: u64,
    /// 1-based position of the first record; 0 for unplaced.
    start: i64,
    /// 0-based exclusive end of the furthest record.
    end: i64,
    prev_ap: i64,
}

/// htslib's `seqs_per_slice` and `bases_per_slice` defaults.
pub(super) const MAX_RECORDS: u32 = 10_000;
const MAX_BASES: u64 = 5_000_000;

impl SliceBuilder {
    pub(super) fn kind(&self) -> Option<SliceKind> {
        self.kind
    }

    pub(super) fn n_records(&self) -> u32 {
        self.n_records
    }

    /// Whether the slice must be written before a record of `kind` with
    /// `read_len` bases, holding at most `max_records`.
    pub(super) fn ends_before(&self, kind: SliceKind, read_len: i32, max_records: u32) -> bool {
        self.kind.is_some_and(|k| {
            k != kind
                || self.n_records >= max_records
                || self.n_bases.saturating_add(u64::from(read_len.unsigned_abs())) > MAX_BASES
        })
    }

    /// Start a record at 1-based `pos1` (0 unplaced) in a slice of `kind`.
    pub(super) fn begin_record(&mut self, kind: SliceKind, pos1: i64, read_len: i32) {
        if self.kind.is_none() {
            self.kind = Some(kind);
            // An unplaced slice starts at 0 whatever its records' `POS`.
            self.start = if kind == SliceKind::Unplaced { 0 } else { pos1 };
            self.prev_ap = self.start;
        }
        self.n_records = self.n_records.saturating_add(1);
        self.n_bases = self.n_bases.saturating_add(u64::from(read_len.unsigned_abs()));
        let ap = pos1.saturating_sub(self.prev_ap);
        self.prev_ap = pos1;
        // Both are 0 or a 1-based position below `i32::MAX` (`check` refuses
        // `Pos0::MAX`), so their difference fits.
        debug_assert!(i32::try_from(ap).is_ok(), "AP delta {ap}");
        self.int(Series::Ap, i32::try_from(ap).unwrap_or(i32::MAX));
    }

    // r[impl cram_writer.slice_extent]
    /// Extend the slice to cover a record ending (0-based, exclusive) at
    /// `end`.
    pub(super) fn reach(&mut self, end: i64) {
        self.end = self.end.max(end);
    }

    /// The byte offset the next value of `series` goes to.
    pub(super) fn len(&self, series: Series) -> usize {
        self.series.get(series as usize).map_or(0, Vec::len)
    }

    pub(super) fn push_mate(&mut self, mate: MateRecord) {
        self.mates.push(mate);
    }

    pub(super) fn int(&mut self, series: Series, value: i32) {
        write_itf8(self.buf(series), value);
    }

    pub(super) fn byte(&mut self, series: Series, value: u8) {
        self.buf(series).push(value);
    }

    pub(super) fn bytes(&mut self, series: Series, value: &[u8]) {
        self.buf(series).extend_from_slice(value);
    }

    fn buf(&mut self, series: Series) -> &mut Vec<u8> {
        let index = series as usize;
        debug_assert!(index < SERIES.len());
        #[allow(clippy::indexing_slicing, reason = "every Series is an index of `series`")]
        &mut self.series[index]
    }

    // r[impl cram_writer.tags]
    /// Store a tag's value: `head` is its name and type, `value` its BAM
    /// bytes (with the NUL of a `Z` or `H`).
    pub(super) fn tag(&mut self, head: [u8; 3], value: &[u8]) {
        let key = tag_key(head);
        let block = match self.tags.iter().position(|t| t.key == key) {
            Some(i) => self.tags.get_mut(i),
            None => {
                self.tags.push(TagBlock { key, data: Vec::new() });
                self.tags.last_mut()
            }
        };
        let Some(block) = block else { return };
        match head[2] {
            b'Z' | b'H' => {
                block.data.extend_from_slice(value);
                block.data.push(b'\t');
            }
            b'B' => {
                write_itf8(&mut block.data, i32::try_from(value.len()).unwrap_or(i32::MAX));
                block.data.extend_from_slice(value);
            }
            _ => block.data.extend_from_slice(value),
        }
    }

    /// The tag dictionary index of the tag line `line` (its tags' 3-byte
    /// heads, concatenated).
    pub(super) fn tag_line(&mut self, line: &[u8]) -> i32 {
        if let Some(&id) = self.line_ids.get(line) {
            return id;
        }
        let id = i32::try_from(self.tag_lines.len()).unwrap_or(i32::MAX);
        self.tag_lines.push(line.to_vec());
        self.line_ids.insert(line.to_vec(), id);
        id
    }

    /// Empty the slice, keeping its buffers.
    pub(super) fn clear(&mut self) {
        let Self {
            kind,
            series,
            tags,
            mates,
            tag_lines,
            line_ids,
            n_records,
            n_bases,
            start,
            end,
            prev_ap,
        } = self;
        *kind = None;
        series.iter_mut().for_each(Vec::clear);
        tags.clear();
        mates.clear();
        tag_lines.clear();
        line_ids.clear();
        *n_records = 0;
        *n_bases = 0;
        *start = 0;
        *end = 0;
        *prev_ap = 0;
    }

    /// The slice's 1-based start and its span, cut at `length`, its
    /// reference's length.
    fn extent(&self, length: u64) -> (i64, i64) {
        match self.kind {
            Some(SliceKind::Mapped { .. }) => {
                let end = self.end.min(i64::try_from(length).unwrap_or(i64::MAX));
                (self.start, end.saturating_sub(self.start.saturating_sub(1)).max(0))
            }
            _ => (0, 0),
        }
    }

    /// What assembling the slice needs from the writer: the records before
    /// it, its reference's length, and the MD5 of `reference` (uppercased
    /// bases) over its extent when it is stored against them.
    pub(super) fn frame(
        &self,
        record_counter: u64,
        length: u64,
        reference: Option<&[u8]>,
    ) -> Frame {
        let md5 = match (self.kind, reference) {
            (Some(SliceKind::Mapped { reference: true, .. }), Some(bases)) => {
                let (start, span) = self.extent(length);
                let from = usize::try_from(start.saturating_sub(1)).unwrap_or(usize::MAX);
                let to = from.saturating_add(usize::try_from(span).unwrap_or(0));
                *md5::compute(bases.get(from..to).unwrap_or_default())
            }
            _ => [0; 16],
        };
        Frame { record_counter, length, md5 }
    }

    // r[impl cram_writer.container_header]
    // r[impl cram_writer.slice_header]
    // r[impl cram_writer.blocks]
    /// The container of this slice.
    pub(super) fn assemble(
        &self,
        frame: &Frame,
        codecs: &mut Codecs,
    ) -> Result<Assembled, CramWriteError> {
        let &Frame { record_counter, length, md5 } = frame;
        let kind = self.kind.unwrap_or(SliceKind::Unplaced);
        let (start, span) = self.extent(length);
        let ref_id = match kind {
            SliceKind::Mapped { tid, .. } => i32::try_from(tid.as_u32()).unwrap_or(i32::MAX),
            SliceKind::Unplaced => -1,
        };
        let start_i32 = i32::try_from(start).map_err(|_| CramWriteError::ValueTooLarge {
            what: "slice start",
            value: start.unsigned_abs(),
        })?;
        let span_i32 = i32::try_from(span).map_err(|_| CramWriteError::ValueTooLarge {
            what: "slice span",
            value: span.unsigned_abs(),
        })?;
        let n_records = i32::try_from(self.n_records).unwrap_or(i32::MAX);

        // r[impl cram_writer.mates]
        let mates = mates::series(
            &self.mates,
            self.series.get(Series::Rn as usize).map_or(&[], Vec::as_slice),
        );
        // The external blocks with data, series first, then tags.
        let blocks: Vec<(i32, &[u8])> = SERIES
            .into_iter()
            .filter_map(|series| {
                let data = match series {
                    Series::Cf => &mates.cf,
                    Series::Mf => &mates.mf,
                    Series::Ns => &mates.ns,
                    Series::Np => &mates.np,
                    Series::Ts => &mates.ts,
                    Series::Nf => &mates.nf,
                    _ => self.series.get(series as usize)?,
                };
                Some((series.content_id(), data.as_slice()))
            })
            .chain(self.tags.iter().map(|t| (t.key, t.data.as_slice())))
            .filter(|(_, data)| !data.is_empty())
            .collect();
        let mut external = Vec::new();
        for &(id, data) in &blocks {
            let content =
                if id == Series::Rn.content_id() { Content::Names } else { Content::Data };
            codecs.block(&mut external, 4, id, data, content)?;
        }
        let n_external = blocks.len();

        let mut compression_header = Vec::new();
        Codecs::raw_block(
            &mut compression_header,
            1,
            0,
            &self.compression_header(matches!(kind, SliceKind::Mapped { reference: true, .. })),
        )?;

        let mut slice_header = Vec::new();
        write_itf8(&mut slice_header, ref_id);
        write_itf8(&mut slice_header, start_i32);
        write_itf8(&mut slice_header, span_i32);
        write_itf8(&mut slice_header, n_records);
        write_ltf8(&mut slice_header, record_counter);
        write_itf8(&mut slice_header, len_i32("block count", n_external.saturating_add(1))?);
        write_itf8(&mut slice_header, len_i32("block count", n_external)?);
        for &(id, _) in &blocks {
            write_itf8(&mut slice_header, id);
        }
        write_itf8(&mut slice_header, -1);
        slice_header.extend_from_slice(&md5);

        let mut slice = Vec::new();
        Codecs::raw_block(&mut slice, 2, 0, &slice_header)?;
        Codecs::raw_block(&mut slice, 5, 0, &[])?;
        slice.extend_from_slice(&external);

        let n_blocks = n_external.saturating_add(3);
        let length = compression_header.len().saturating_add(slice.len());
        let mut bytes = Vec::with_capacity(length.saturating_add(64));
        bytes.extend_from_slice(&len_i32("container length", length)?.to_le_bytes());
        write_itf8(&mut bytes, ref_id);
        write_itf8(&mut bytes, start_i32);
        write_itf8(&mut bytes, span_i32);
        write_itf8(&mut bytes, n_records);
        write_ltf8(&mut bytes, record_counter);
        write_ltf8(&mut bytes, self.n_bases);
        write_itf8(&mut bytes, len_i32("block count", n_blocks)?);
        write_itf8(&mut bytes, 1);
        write_itf8(&mut bytes, len_i32("landmark", compression_header.len())?);
        push_crc(&mut bytes, 0);
        bytes.extend_from_slice(&compression_header);
        bytes.extend_from_slice(&slice);

        Ok(Assembled {
            bytes,
            alignment_start: start,
            alignment_span: span,
            slice_offset: compression_header.len() as u64,
            slice_size: slice.len() as u64,
        })
    }

    fn compression_header(&self, reference: bool) -> Vec<u8> {
        // r[impl cram_writer.preservation_map]
        let mut map = Vec::new();
        write_itf8(&mut map, 5);
        map.extend_from_slice(b"RN\x01AP\x01RR");
        map.push(u8::from(reference));
        map.extend_from_slice(b"SM");
        // r[impl cram_writer.substitution]
        map.extend_from_slice(&[SUBSTITUTION_MATRIX_BYTE; 5]);
        map.extend_from_slice(b"TD");
        let mut td = Vec::new();
        for line in &self.tag_lines {
            td.extend_from_slice(line);
            td.push(0);
        }
        if td.is_empty() {
            td.push(0);
        }
        write_small_len(&mut map, td.len());
        map.extend_from_slice(&td);

        let mut ds = Vec::new();
        // Every series, used or not: htslib's decoder wants an encoding for
        // one it reads nothing of (`QS` of a record without bases).
        write_small_len(&mut ds, SERIES.len().saturating_sub(1));
        for series in SERIES {
            series.encoding(&mut ds);
        }

        let mut tags = Vec::new();
        write_small_len(&mut tags, self.tags.len());
        for tag in &self.tags {
            write_itf8(&mut tags, tag.key);
            #[allow(clippy::cast_sign_loss, reason = "the low byte of a key is the BAM type")]
            match (tag.key & 0xFF) as u8 {
                b'Z' | b'H' => byte_array_stop(&mut tags, b'\t', tag.key),
                b'B' => byte_array_len(
                    &mut tags,
                    |p| external(p, tag.key),
                    |p| {
                        external(p, tag.key);
                    },
                ),
                ty => byte_array_len(
                    &mut tags,
                    |p| constant(p, fixed_size(ty)),
                    |p| {
                        external(p, tag.key);
                    },
                ),
            }
        }

        let mut out = Vec::new();
        for part in [map, ds, tags] {
            write_small_len(&mut out, part.len());
            out.extend_from_slice(&part);
        }
        out
    }
}

/// The size of a fixed-size tag type's value.
fn fixed_size(ty: u8) -> i32 {
    match ty {
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        _ => 1,
    }
}

/// A tag's key from its name and type.
pub(super) fn tag_key([a, b, ty]: [u8; 3]) -> i32 {
    i32::from(a) << 16 | i32::from(b) << 8 | i32::from(ty)
}

/// Each reference base's other bases in `ACGTN` order get codes 0–3:
/// `0b00_01_10_11`.
pub(super) const SUBSTITUTION_MATRIX_BYTE: u8 = 0x1B;

/// What assembling a slice needs from the writer.
pub(super) struct Frame {
    record_counter: u64,
    length: u64,
    md5: [u8; 16],
}

/// A written container and where its slice is in it.
pub(super) struct Assembled {
    pub bytes: Vec<u8>,
    pub alignment_start: i64,
    pub alignment_span: i64,
    pub slice_offset: u64,
    pub slice_size: u64,
}

fn push_crc(out: &mut Vec<u8>, from: usize) {
    let mut crc = libdeflater::Crc::new();
    crc.update(out.get(from..).unwrap_or_default());
    out.extend_from_slice(&crc.sum().to_le_bytes());
}

/// What a block holds, for the codecs it is tried with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Content {
    /// The SAM header: gzip only, which every reader takes.
    Header,
    /// The `RN` block: the name tokeniser too, in CRAM 3.1.
    Names,
    Data,
}

/// The block codecs, with their reusable state.
pub(super) struct Codecs {
    level: CompressionLevel,
    version: CramVersion,
    gzip: libdeflater::Compressor,
    scratch: Vec<u8>,
}

/// Blocks shorter than this are stored raw: no codec's header pays for
/// itself.
const MIN_COMPRESSED: usize = 32;

impl Codecs {
    pub(super) fn new(level: CompressionLevel, version: CramVersion) -> Self {
        let gzip_level = libdeflater::CompressionLvl::new(i32::from(level.get()))
            .unwrap_or_else(|_| libdeflater::CompressionLvl::default());
        Self { level, version, gzip: libdeflater::Compressor::new(gzip_level), scratch: Vec::new() }
    }

    // r[impl cram_writer.compression]
    /// Append a block of `data`, stored by whichever codec gives the fewest
    /// bytes: gzip, and rANS 4×8 (CRAM 3.0) or rANS Nx16 with htslib's
    /// normal profile and, for names, the name tokeniser (CRAM 3.1).
    fn block(
        &mut self,
        out: &mut Vec<u8>,
        content_type: u8,
        content_id: i32,
        data: &[u8],
        content: Content,
    ) -> Result<(), CramWriteError> {
        if data.len() < MIN_COMPRESSED || self.level == CompressionLevel::Level0 {
            return Self::raw_block(out, content_type, content_id, data);
        }
        let mut best: Option<(u8, Vec<u8>)> = None;
        let mut consider = |method: u8, bytes: Vec<u8>| {
            if bytes.len() < best.as_ref().map_or(data.len(), |(_, b)| b.len()) {
                best = Some((method, bytes));
            }
        };
        self.scratch.resize(self.gzip.gzip_compress_bound(data.len()), 0);
        if let Ok(n) = self.gzip.gzip_compress(data, &mut self.scratch) {
            consider(1, self.scratch.get(..n).unwrap_or_default().to_vec());
        }
        match (content, self.version) {
            (Content::Header, _) => {}
            (_, CramVersion::V3_0) => {
                for order in [rans::Order::Zero, rans::Order::One] {
                    if let Some(bytes) = rans::encode(data, order) {
                        consider(4, bytes);
                    }
                }
            }
            (_, CramVersion::V3_1) => {
                for flags in rans_nx16::normal_profile(data.len()) {
                    if let Some(bytes) = rans_nx16::encode(data, flags) {
                        consider(5, bytes);
                    }
                }
                if content == Content::Names
                    && let Some(bytes) = tok3::encode(data)
                {
                    consider(8, bytes);
                }
            }
        }
        match best {
            Some((method, bytes)) => {
                write_block(out, method, content_type, content_id, data.len(), &bytes)
            }
            None => Self::raw_block(out, content_type, content_id, data),
        }
    }

    pub(super) fn raw_block(
        out: &mut Vec<u8>,
        content_type: u8,
        content_id: i32,
        data: &[u8],
    ) -> Result<(), CramWriteError> {
        write_block(out, 0, content_type, content_id, data.len(), data)
    }
}

fn write_block(
    out: &mut Vec<u8>,
    method: u8,
    content_type: u8,
    content_id: i32,
    raw_len: usize,
    data: &[u8],
) -> Result<(), CramWriteError> {
    let start = out.len();
    out.push(method);
    out.push(content_type);
    write_itf8(out, content_id);
    write_itf8(out, len_i32("block size", data.len())?);
    write_itf8(out, len_i32("block size", raw_len)?);
    out.extend_from_slice(data);
    push_crc(out, start);
    Ok(())
}

// r[impl cram_writer.header_container]
/// The container holding the SAM header text.
pub(super) fn header_container(text: &str, codecs: &mut Codecs) -> Result<Vec<u8>, CramWriteError> {
    let mut data = Vec::with_capacity(text.len().saturating_add(4));
    data.extend_from_slice(&len_i32("header text", text.len())?.to_le_bytes());
    data.extend_from_slice(text.as_bytes());
    let mut block = Vec::new();
    codecs.block(&mut block, 0, 0, &data, Content::Header)?;

    let mut bytes = Vec::with_capacity(block.len().saturating_add(32));
    bytes.extend_from_slice(&len_i32("container length", block.len())?.to_le_bytes());
    for _ in 0..4 {
        write_itf8(&mut bytes, 0); // reference, start, span, records
    }
    write_ltf8(&mut bytes, 0);
    write_ltf8(&mut bytes, 0);
    write_itf8(&mut bytes, 1);
    write_itf8(&mut bytes, 1);
    write_itf8(&mut bytes, 0);
    push_crc(&mut bytes, 0);
    bytes.extend_from_slice(&block);
    Ok(bytes)
}

/// htslib's CRAM 3 EOF container (`cram_check_EOF`'s template).
pub(super) const EOF_CONTAINER: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];
