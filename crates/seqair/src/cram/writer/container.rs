//! A slice's data series as the records fill them, and the container that
//! frames them: compression header, slice header, blocks.

use super::CramWriteError;
use super::embed::Held;
use super::fqzcomp::QualRecord;
use super::mates::{self, MateRecord};
use super::metrics::{Method, Metrics, Plan, Trial};
use super::{CramProfile, CramVersion, arith, fqzcomp, rans, rans_nx16, tok3};
use crate::cram::varint::{decode_itf8, write_itf8, write_ltf8};
use crate::io::CompressionLevel;
use crate::reader::Tid;
use rustc_hash::FxHashMap;
use std::sync::Arc;

/// Every entry of [`arith::PROFILE`].
const ARITH: [Method; 8] = [
    Method::Arith(0),
    Method::Arith(1),
    Method::Arith(2),
    Method::Arith(3),
    Method::Arith(4),
    Method::Arith(5),
    Method::Arith(6),
    Method::Arith(7),
];

/// Every entry of [`fqzcomp::Strategy::ALL`].
const FQZ: [Method; 4] = [Method::Fqz(0), Method::Fqz(1), Method::Fqz(2), Method::Fqz(3)];

/// Every entry of [`rans_nx16::normal_profile`].
const NX16: [Method; 6] = [
    Method::Nx16(0),
    Method::Nx16(1),
    Method::Nx16(2),
    Method::Nx16(3),
    Method::Nx16(4),
    Method::Nx16(5),
];

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
    /// Read after `CF`, but last here: only a multi-reference slice has it
    /// (`r[cram_writer.multi_ref.records]`).
    Ri,
}

const SERIES: [Series; 28] = [
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
    Series::Ri,
];

impl Series {
    fn content_id(self) -> i32 {
        (self as i32).saturating_add(1)
    }

    /// The one value a series of integers or bytes holds, if it holds no
    /// other; 0 for one that holds none. `None` for a byte-array series.
    pub(super) fn constant(self, data: &[u8]) -> Option<i32> {
        match self {
            Self::Rn | Self::In | Self::Sc | Self::Bb | Self::BbLen => None,
            Self::Fc | Self::Ba | Self::Qs | Self::Bs => {
                let first = data.first().copied().unwrap_or(0);
                data.iter().all(|&b| b == first).then_some(i32::from(first))
            }
            _ => {
                let mut rest = data;
                let mut first = None;
                while !rest.is_empty() {
                    let (value, n) = decode_itf8(rest)?;
                    rest = rest.get(n..)?;
                    if *first.get_or_insert(value) != value {
                        return None;
                    }
                }
                // ITF8 holds an `i32`'s bits.
                #[allow(clippy::cast_possible_wrap, reason = "the bits of an i32")]
                Some(first.map_or(0, |v| v as i32))
            }
        }
    }

    // r[impl cram_writer.data_series+2]
    /// The series' key and encoding in the data series map, `constant` if
    /// it holds one value; `None` for `BbLen`, which `BB`'s encoding names.
    fn encoding(self, out: &mut Vec<u8>, constant: Option<i32>) -> Option<()> {
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
            Self::Ri => b"RI",
            Self::BbLen => return None,
        };
        out.extend_from_slice(key);
        if let Some(value) = constant {
            constant_code(out, value);
            return Some(());
        }
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
fn constant_code(out: &mut Vec<u8>, value: i32) {
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
    /// Records of several references the FASTA has
    /// (`r[cram_writer.multi_ref]`).
    Multi,
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
    /// Where each record's `QS` values start, for fqzcomp.
    qs_starts: Vec<usize>,
    tag_lines: Vec<Vec<u8>>,
    line_ids: FxHashMap<Vec<u8>, i32>,
    /// The tag block after the one last written, and the last tag line:
    /// records mostly repeat their tags in order, so both are tried first.
    next_tag: usize,
    last_line: Option<i32>,
    /// Records on a reference the FASTA lacks, encoded when the slice is
    /// complete (`r[cram_writer.embedded_reference]`).
    held: Vec<Held>,
    /// The reference the slice embeds; empty for none.
    embedded: Vec<u8>,
    /// Each run of records on one reference id: the id, the run's first
    /// 1-based position and furthest end, for the index of a
    /// multi-reference slice (`r[cram_writer.crai+2]`).
    runs: Vec<(Option<Tid>, i64, i64)>,
    n_records: u32,
    n_bases: u64,
    /// 1-based position of the first record; 0 for unplaced.
    start: i64,
    /// 0-based exclusive end of the furthest record.
    end: i64,
    prev_ap: i64,
}

/// Whether a slice of `n_records` and `n_bases` holding at most
/// `max_records` is full before a record of `read_len` bases.
pub(super) fn is_full(n_records: u32, n_bases: u64, read_len: i32, max_records: u32) -> bool {
    n_records >= max_records
        || n_bases.saturating_add(u64::from(read_len.unsigned_abs()))
            > MAX_BASES.max(u64::from(max_records).saturating_mul(BASES_PER_RECORD))
}

/// htslib's `seqs_per_slice` and `bases_per_slice` defaults.
pub(super) const MAX_RECORDS: u32 = 10_000;
/// htslib's archive profile's `seqs_per_slice`.
pub(super) const ARCHIVE_RECORDS: u32 = 100_000;
const MAX_BASES: u64 = 5_000_000;
/// htslib's bases per slice when a profile sets the records per slice
/// (`seqs_per_slice * 500`); never below `MAX_BASES`, so a test's small
/// `records_per_slice` keeps the default bases cap.
const BASES_PER_RECORD: u64 = 500;

impl SliceBuilder {
    pub(super) fn kind(&self) -> Option<SliceKind> {
        self.kind
    }

    pub(super) fn n_records(&self) -> u32 {
        self.n_records
    }

    /// The runs of records on one reference id a multi-reference slice holds.
    pub(super) fn reference_runs(&self) -> usize {
        self.runs.len()
    }

    /// Whether the slice must be written before a record of `kind` with
    /// `read_len` bases, holding at most `max_records`.
    pub(super) fn ends_before(&self, kind: SliceKind, read_len: i32, max_records: u32) -> bool {
        self.kind.is_some_and(|k| {
            k != kind || is_full(self.n_records, self.n_bases, read_len, max_records)
        })
    }

    /// Hold a record of a slice of `kind` to encode when the slice is
    /// complete; it counts towards the slice's size.
    pub(super) fn hold(&mut self, kind: SliceKind, record: Held) {
        self.kind.get_or_insert(kind);
        self.n_records = self.n_records.saturating_add(1);
        let read_len = record.checked.read_len;
        self.n_bases = self.n_bases.saturating_add(u64::from(read_len.unsigned_abs()));
        self.held.push(record);
    }

    /// The held records, the slice left empty to encode them into.
    pub(super) fn take_held(&mut self) -> Vec<Held> {
        self.kind = None;
        self.n_records = 0;
        self.n_bases = 0;
        std::mem::take(&mut self.held)
    }

    /// Embed `bases`, the reference the slice's records are stored against.
    pub(super) fn embed(&mut self, bases: Vec<u8>) {
        self.embedded = bases;
    }

    /// Start a record at 1-based `pos1` (0 unplaced) in a slice of `kind`.
    pub(super) fn begin_record(&mut self, kind: SliceKind, pos1: i64, read_len: i32) {
        if self.kind.is_none() {
            self.kind = Some(kind);
            // An unplaced or multi-reference slice starts at 0 whatever its
            // records' `POS`.
            self.start = if matches!(kind, SliceKind::Mapped { .. }) { pos1 } else { 0 };
            self.prev_ap = self.start;
        }
        let qs = self.len(Series::Qs);
        self.qs_starts.push(qs);
        self.n_records = self.n_records.saturating_add(1);
        self.n_bases = self.n_bases.saturating_add(u64::from(read_len.unsigned_abs()));
        let ap = pos1.saturating_sub(self.prev_ap);
        self.prev_ap = pos1;
        // Both are 0 or a 1-based position below `i32::MAX` (`check` refuses
        // `Pos0::MAX`), so their difference fits.
        debug_assert!(i32::try_from(ap).is_ok(), "AP delta {ap}");
        self.int(Series::Ap, i32::try_from(ap).unwrap_or(i32::MAX));
    }

    // r[impl cram_writer.slice_extent+2]
    /// Extend the slice to cover a record ending (0-based, exclusive) at
    /// `end`.
    pub(super) fn reach(&mut self, end: i64) {
        self.end = self.end.max(end);
    }

    // r[impl cram_writer.multi_ref.records]
    /// A record of reference `tid` at 1-based `pos1` (0 unplaced) ending at
    /// 1-based `aend` (inclusive, as htslib's decode sets it): its `RI` in a
    /// multi-reference slice, and its run for the index.
    pub(super) fn reference_id(&mut self, tid: Option<Tid>, pos1: i64, aend: i64) {
        if self.kind != Some(SliceKind::Multi) {
            return;
        }
        let ri = tid.map_or(-1, |t| i32::try_from(t.as_u32()).unwrap_or(i32::MAX));
        self.int(Series::Ri, ri);
        match self.runs.last_mut() {
            Some((run, _, end)) if *run == tid => *end = (*end).max(aend),
            _ => self.runs.push((tid, pos1, aend)),
        }
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
        debug_assert!(
            !matches!(
                series,
                Series::Cf | Series::Mf | Series::Ns | Series::Np | Series::Ts | Series::Nf
            ),
            "{series:?} is written by mates::series when the slice is assembled"
        );
        let index = series as usize;
        debug_assert!(index < SERIES.len());
        #[allow(clippy::indexing_slicing, reason = "every Series is an index of `series`")]
        &mut self.series[index]
    }

    /// Store an `X` feature's substitution: reference row `row` (`ACGTN`),
    /// and `alt`, the read base's index among the row's other bases. `BS`
    /// holds `row * 4 + alt` until the slice's matrix codes it.
    pub(super) fn substitution(&mut self, row: u8, alt: u8) {
        debug_assert!(row < 5 && alt < 4);
        self.buf(Series::Bs).push(row.saturating_mul(4).saturating_add(alt));
    }

    // r[impl cram_writer.tags]
    /// Store a tag's value: `head` is its name and type, `value` its BAM
    /// bytes (with the NUL of a `Z` or `H`).
    pub(super) fn tag(&mut self, head: [u8; 3], value: &[u8]) {
        let key = tag_key(head);
        let hinted = self.tags.get(self.next_tag).is_some_and(|t| t.key == key);
        let found =
            if hinted { Some(self.next_tag) } else { self.tags.iter().position(|t| t.key == key) };
        let index = found.unwrap_or_else(|| {
            self.tags.push(TagBlock { key, data: Vec::new() });
            self.tags.len().saturating_sub(1)
        });
        self.next_tag = index.saturating_add(1);
        let block = self.tags.get_mut(index);
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
        if let Some(id) = self.last_line
            && usize::try_from(id)
                .ok()
                .and_then(|i| self.tag_lines.get(i))
                .is_some_and(|l| l == line)
        {
            return id;
        }
        if let Some(&id) = self.line_ids.get(line) {
            self.last_line = Some(id);
            return id;
        }
        let id = i32::try_from(self.tag_lines.len()).unwrap_or(i32::MAX);
        self.tag_lines.push(line.to_vec());
        self.line_ids.insert(line.to_vec(), id);
        self.last_line = Some(id);
        id
    }

    /// Empty the slice, keeping its buffers.
    pub(super) fn clear(&mut self) {
        let Self {
            kind,
            series,
            tags,
            mates,
            qs_starts,
            tag_lines,
            line_ids,
            next_tag,
            last_line,
            held,
            embedded,
            runs,
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
        qs_starts.clear();
        tag_lines.clear();
        line_ids.clear();
        *next_tag = 0;
        *last_line = None;
        held.clear();
        embedded.clear();
        runs.clear();
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
        slice_index: u64,
        record_counter: u64,
        length: u64,
        reference: Option<&[u8]>,
    ) -> Frame {
        let md5 = match (self.kind, reference) {
            _ if !self.embedded.is_empty() => {
                debug_assert_eq!(
                    i64::try_from(self.embedded.len()).ok(),
                    Some(self.extent(length).1)
                );
                *md5::compute(&self.embedded)
            }
            (Some(SliceKind::Mapped { reference: true, .. }), Some(bases)) => {
                let (start, span) = self.extent(length);
                let from = usize::try_from(start.saturating_sub(1)).unwrap_or(usize::MAX);
                let to = from.saturating_add(usize::try_from(span).unwrap_or(0));
                *md5::compute(bases.get(from..to).unwrap_or_default())
            }
            _ => [0; 16],
        };
        Frame { slice_index, record_counter, length, md5 }
    }

    // r[impl cram_writer.container_header+2]
    // r[impl cram_writer.slice_header+2]
    // r[impl cram_writer.blocks+2]
    /// The container of this slice.
    pub(super) fn assemble(
        &self,
        frame: &Frame,
        codecs: &mut Codecs,
    ) -> Result<Assembled, CramWriteError> {
        let &Frame { slice_index, record_counter, length, md5 } = frame;
        // First, so a trial slice that fails anywhere still releases the
        // slices waiting on it.
        let mut trial =
            TrialGuard { metrics: Arc::clone(&codecs.metrics), slice_index, trials: Vec::new() };
        let kind = self.kind.unwrap_or(SliceKind::Unplaced);
        let (start, span) = self.extent(length);
        let ref_id = match kind {
            SliceKind::Mapped { tid, .. } => i32::try_from(tid.as_u32()).unwrap_or(i32::MAX),
            SliceKind::Multi => -2,
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

        let (sm, bs) =
            substitution_matrix(self.series.get(Series::Bs as usize).map_or(&[], Vec::as_slice));
        // r[impl cram_writer.mates+2]
        let mates = mates::series(
            &self.mates,
            self.series.get(Series::Rn as usize).map_or(&[], Vec::as_slice),
        );
        let empty = Vec::new();
        let series_data = SERIES.map(|series| {
            let data = match series {
                Series::Cf => &mates.cf,
                Series::Mf => &mates.mf,
                Series::Ns => &mates.ns,
                Series::Np => &mates.np,
                Series::Ts => &mates.ts,
                Series::Nf => &mates.nf,
                Series::Bs => &bs,
                _ => self.series.get(series as usize).unwrap_or(&empty),
            };
            (series, data.as_slice())
        });
        let multi = kind == SliceKind::Multi;
        // r[impl cram_writer.data_series+2]
        let constants = series_data.map(|(series, data)| series.constant(data));
        // The external blocks, series first, then tags: every series that is
        // no constant and holds data.
        let blocks: Vec<(i32, &[u8])> = series_data
            .iter()
            .zip(&constants)
            .filter(|((series, data), constant)| {
                constant.is_none() && !data.is_empty() && (*series != Series::Ri || multi)
            })
            .map(|(&(series, data), _)| (series.content_id(), data))
            .chain(self.tags.iter().map(|t| (t.key, t.data.as_slice())))
            .chain(
                (!self.embedded.is_empty())
                    .then_some((EMBEDDED_REFERENCE, self.embedded.as_slice())),
            )
            .collect();
        let quals = if codecs.version == CramVersion::V3_1(CramProfile::Archive) {
            self.qual_records()
        } else {
            Vec::new()
        };
        let mut external = Vec::new();
        let plan = codecs.metrics.plan(slice_index);
        for &(id, data) in &blocks {
            let content = if id == Series::Rn.content_id() {
                Content::Names
            } else if id == Series::Qs.content_id() {
                Content::Qualities(&quals)
            } else {
                Content::Data
            };
            let tried = match &plan {
                Plan::Trial => Tried::Record(&mut trial.trials),
                Plan::Learned(choices) => {
                    match choices.get(&id).and_then(|c| c.methods(data.len())) {
                        Some(methods) => Tried::Only(methods),
                        None => Tried::All,
                    }
                }
            };
            codecs.block(&mut external, 4, id, data, content, tried)?;
        }
        drop(trial);
        let n_external = blocks.len();

        let mut compression_header = Vec::new();
        Codecs::raw_block(
            &mut compression_header,
            1,
            0,
            &self.compression_header(
                matches!(kind, SliceKind::Mapped { reference: true, .. } | SliceKind::Multi),
                sm,
                &constants,
            ),
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
        write_itf8(
            &mut slice_header,
            if self.embedded.is_empty() { -1 } else { EMBEDDED_REFERENCE },
        );
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

        // r[impl cram_writer.crai+2]
        let spans = match kind {
            SliceKind::Multi => self
                .runs
                .iter()
                .map(|&(tid, from, end)| (tid, from, end.saturating_sub(from).saturating_add(1)))
                .collect(),
            SliceKind::Mapped { tid, .. } => vec![(Some(tid), start, span)],
            SliceKind::Unplaced => vec![(None, start, span)],
        };
        Ok(Assembled {
            bytes,
            spans,
            slice_offset: compression_header.len() as u64,
            slice_size: slice.len() as u64,
        })
    }

    /// Each record's share of the `QS` block and the flags fqzcomp models.
    fn qual_records(&self) -> Vec<QualRecord> {
        debug_assert_eq!(self.qs_starts.len(), self.mates.len(), "one of each per record");
        let total = self.len(Series::Qs);
        let ends = self.qs_starts.iter().skip(1).copied().chain(std::iter::once(total));
        self.qs_starts
            .iter()
            .zip(ends)
            .zip(&self.mates)
            .map(|((&start, end), mate)| QualRecord {
                len: u32::try_from(end.saturating_sub(start)).unwrap_or(u32::MAX),
                reverse: mate.flags.is_reverse(),
                read2: mate.flags.is_second_in_template(),
            })
            .collect()
    }

    fn compression_header(
        &self,
        reference: bool,
        sm: [u8; 5],
        constants: &[Option<i32>; SERIES.len()],
    ) -> Vec<u8> {
        // r[impl cram_writer.preservation_map+2]
        let mut map = Vec::new();
        write_itf8(&mut map, 5);
        map.extend_from_slice(b"RN\x01AP\x01RR");
        map.push(u8::from(reference));
        map.extend_from_slice(b"SM");
        map.extend_from_slice(&sm);
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
        let multi = self.kind == Some(SliceKind::Multi);
        // `BbLen` is `BB`'s, and `RI` only in a multi-reference slice.
        write_small_len(&mut ds, SERIES.len().saturating_sub(if multi { 1 } else { 2 }));
        for (series, &constant) in
            SERIES.into_iter().zip(constants).filter(|&(s, _)| s != Series::Ri || multi)
        {
            series.encoding(&mut ds, constant);
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
                    |p| constant_code(p, fixed_size(ty)),
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

/// The content id of an embedded reference's block: no series' (1 to
/// [`SERIES`]' length) and no tag's (from `A` shifted 16 bits up).
const EMBEDDED_REFERENCE: i32 = 0x100;

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

// r[impl cram_writer.substitution+2]
/// The slice's substitution matrix from its `BS` pairs (`row * 4 + alt`):
/// the `SM` bytes, and `BS` coded by them. Each row codes its other bases by
/// how often the slice substitutes them, most often 0, ties in `ACGTN` order.
pub(super) fn substitution_matrix(pairs: &[u8]) -> ([u8; 5], Vec<u8>) {
    let mut counts = [[0u32; 4]; 5];
    for &pair in pairs {
        if let Some(c) =
            counts.get_mut(usize::from(pair / 4)).and_then(|row| row.get_mut(usize::from(pair % 4)))
        {
            *c = c.saturating_add(1);
        }
    }
    let mut codes = [[0u8; 4]; 5];
    let mut sm = [0u8; 5];
    for ((row, code), byte) in counts.iter().zip(&mut codes).zip(&mut sm) {
        let mut order = [0u8, 1, 2, 3];
        // Stable: ties keep `ACGTN` order.
        order.sort_by_key(|&alt| std::cmp::Reverse(row.get(usize::from(alt)).copied()));
        for (rank, &alt) in (0u8..).zip(&order) {
            if let Some(c) = code.get_mut(usize::from(alt)) {
                *c = rank;
            }
        }
        *byte = code.iter().fold(0, |b, &c| b << 2 | c);
    }
    let coded = pairs
        .iter()
        .map(|&p| {
            codes.get(usize::from(p / 4)).and_then(|row| row.get(usize::from(p % 4))).copied()
        })
        .map(Option::unwrap_or_default)
        .collect();
    (sm, coded)
}

/// What assembling a slice needs from the writer.
pub(super) struct Frame {
    /// The data container's index in the file, from 0.
    slice_index: u64,
    record_counter: u64,
    length: u64,
    md5: [u8; 16],
}

/// A written container and where its slice is in it.
pub(super) struct Assembled {
    pub bytes: Vec<u8>,
    /// The index entries' references, alignment starts and spans.
    pub spans: Vec<(Option<Tid>, i64, i64)>,
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
pub(super) enum Content<'a> {
    /// The SAM header: gzip only, which every reader takes.
    Header,
    /// The `RN` block: the name tokeniser too, in CRAM 3.1.
    Names,
    /// The `QS` block, with each record's share of it: fqzcomp too, in the
    /// archive profile.
    Qualities(&'a [QualRecord]),
    Data,
}

/// The block codecs, with their reusable state.
pub(super) struct Codecs {
    level: CompressionLevel,
    version: CramVersion,
    gzip: libdeflater::Compressor,
    scratch: Vec<u8>,
    metrics: Arc<Metrics>,
}

/// Which methods [`Codecs::block`] tries.
pub(super) enum Tried<'a> {
    /// Every one, recording their sizes.
    Record(&'a mut Vec<Trial>),
    /// Every one.
    All,
    /// The learned ones (raw always).
    Only([Option<Method>; 2]),
}

/// Records a trial slice's sizes when it is dropped, so a slice that fails
/// or panics part way still releases the slices that wait for it.
struct TrialGuard {
    metrics: Arc<Metrics>,
    slice_index: u64,
    trials: Vec<Trial>,
}

impl Drop for TrialGuard {
    fn drop(&mut self) {
        if self.metrics.is_trial(self.slice_index) {
            self.metrics.record(self.slice_index, std::mem::take(&mut self.trials));
        }
    }
}

/// Blocks shorter than this are stored raw: no codec's header pays for
/// itself.
const MIN_COMPRESSED: usize = 32;

impl Codecs {
    /// Codecs of the same level and version that learn with these.
    pub(super) fn fork(&self) -> Self {
        Self::new(self.level, self.version, Arc::clone(&self.metrics))
    }

    pub(super) fn new(
        level: CompressionLevel,
        version: CramVersion,
        metrics: Arc<Metrics>,
    ) -> Self {
        let gzip_level = libdeflater::CompressionLvl::new(i32::from(level.get()))
            .unwrap_or_else(|_| libdeflater::CompressionLvl::default());
        Self {
            level,
            version,
            gzip: libdeflater::Compressor::new(gzip_level),
            scratch: Vec::new(),
            metrics,
        }
    }

    // r[impl cram_writer.profile]
    /// Every method the version and profile have for `content`, in
    /// [`rank`]'s order, so a tie goes to the faster decode.
    fn methods(&self, content: Content<'_>) -> Vec<Method> {
        let mut methods = Vec::new();
        let profile = match (content, self.version) {
            (Content::Header, _) => return vec![Method::Gzip],
            (_, CramVersion::V3_0) => {
                methods.extend([Method::Rans4x8Order0, Method::Rans4x8Order1, Method::Gzip]);
                return methods;
            }
            (_, CramVersion::V3_1(profile)) => profile,
        };
        methods.extend(NX16);
        methods.push(Method::Gzip);
        if content == Content::Names {
            methods.push(Method::Tok3);
        }
        if profile == CramProfile::Archive {
            methods.extend(ARITH);
            if matches!(content, Content::Qualities(_)) {
                methods.extend(FQZ);
            }
        }
        methods
    }

    /// `data` stored by `method`, or `None` where it does not apply; `nx16`
    /// holds what the rANS Nx16 methods share on `data`.
    fn encode(
        &mut self,
        method: Method,
        data: &[u8],
        nx16: &rans_nx16::Prepared<'_>,
        content: Content<'_>,
    ) -> Option<(u8, Vec<u8>)> {
        match method {
            Method::Raw => None,
            Method::Gzip => {
                self.scratch.resize(self.gzip.gzip_compress_bound(data.len()), 0);
                let n = self.gzip.gzip_compress(data, &mut self.scratch).ok()?;
                Some((1, self.scratch.get(..n)?.to_vec()))
            }
            Method::Rans4x8Order0 => Some((4, rans::encode(data, rans::Order::Zero)?)),
            Method::Rans4x8Order1 => Some((4, rans::encode(data, rans::Order::One)?)),
            Method::Nx16(i) => {
                let flags = *rans_nx16::normal_profile(data.len()).get(usize::from(i))?;
                Some((5, rans_nx16::encode_prepared(nx16, flags)?))
            }
            Method::Tok3 => Some((8, tok3::encode(data)?)),
            Method::Arith(i) => {
                Some((6, arith::encode(data, *arith::PROFILE.get(usize::from(i))?)?))
            }
            Method::Fqz(i) => {
                let Content::Qualities(records) = content else { return None };
                let strategy = *fqzcomp::Strategy::ALL.get(usize::from(i))?;
                Some((7, fqzcomp::encode(data, records, strategy)?))
            }
        }
    }

    // r[impl cram_writer.compression+3]
    /// Append a block of `data`, stored by whichever of the methods `tried`
    /// names gives the fewest bytes, raw included.
    fn block(
        &mut self,
        out: &mut Vec<u8>,
        content_type: u8,
        content_id: i32,
        data: &[u8],
        content: Content<'_>,
        tried: Tried<'_>,
    ) -> Result<(), CramWriteError> {
        if data.len() < MIN_COMPRESSED || self.level == CompressionLevel::Level0 {
            return Self::raw_block(out, content_type, content_id, data);
        }
        let mut best: Option<(u8, Vec<u8>)> = None;
        let mut sizes = Vec::new();
        let methods: Vec<Method> = match tried {
            Tried::Only(methods) => methods.into_iter().flatten().collect(),
            Tried::All | Tried::Record(_) => self.methods(content),
        };
        let mut declined = true;
        let nx16 = rans_nx16::Prepared::new(data);
        for method in methods {
            let Some((id, bytes)) = self.encode(method, data, &nx16, content) else { continue };
            declined = false;
            sizes.push((method, bytes.len()));
            if bytes.len() < best.as_ref().map_or(data.len(), |(_, b)| b.len()) {
                best = Some((id, bytes));
            }
        }
        if declined && matches!(tried, Tried::Only([Some(m), _]) if m != Method::Raw) {
            // The learned methods did not apply to this block (as htslib,
            // which redoes the trial): try them all.
            return self.block(out, content_type, content_id, data, content, Tried::All);
        }
        if let Tried::Record(trials) = tried {
            sizes.push((Method::Raw, data.len()));
            trials.push(Trial { content_id, len: data.len(), sizes });
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
    codecs.block(&mut block, 0, 0, &data, Content::Header, Tried::All)?;

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

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    // r[verify cram_writer.slicing+2]
    /// A slice ends at its record cap, or before the read that would take it
    /// past its bases cap: 5,000,000 at 10,000 records, 500 a record above.
    #[test]
    fn a_slice_ends_at_its_records_or_bases() {
        let kind = SliceKind::Unplaced;
        let mut slice = SliceBuilder::default();
        assert!(!slice.ends_before(kind, 150, 10_000));
        for _ in 0..3 {
            slice.begin_record(kind, 0, 1_000_000);
        }
        assert!(slice.ends_before(kind, 2, 3));
        assert!(!slice.ends_before(kind, 2_000_000, 4));
        assert!(slice.ends_before(kind, 2_000_001, 4));
        // 100,000 records: 50,000,000 bases.
        assert!(!slice.ends_before(kind, 47_000_000, 100_000));
        assert!(slice.ends_before(kind, 47_000_001, 100_000));
    }
}
