//! fqzcomp quality codec (CRAM 3.1 block compression method 7).
//!
//! A context-modelling compressor for quality strings: each quality is coded
//! with one of 2^16 adaptive models, picked by a context mixed from the
//! previous qualities, the position in the record and how often the quality
//! has changed so far. The layout of that context is set by a parameter block
//! at the start of the stream.
//!
//! This follows htscodecs' `fqzcomp_qual.c` (`uncompress_block_fqz2f`), which
//! writes the files, where it differs from the `CRAMcodecs` pseudocode; the
//! `cram.codec.fqzcomp.*` rules in the spec list the differences.

// See rans.rs: lazy `ok_or_else(|| CramError::...)` keeps error construction and its
// `drop_in_place<CramError>` off the per-record path.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use super::codec_io::{Uint7Error, read_u8, read_u16_le, read_uint7, split_off};
use super::range_coder::{AdaptiveModel, RangeDecoder, SymFreq, decode_symbol};
use super::reader::{CramError, check_alloc_size};

/// The only fqzcomp stream version htscodecs reads or writes (`FQZ_VERS`).
const VERSION: u8 = 5;

/// Global flags.
const GFLAG_MULTI_PARAM: u8 = 1;
const GFLAG_HAVE_STAB: u8 = 2;
const GFLAG_DO_REV: u8 = 4;

/// Parameter block flags. Bit 1 is reserved and ignored.
const PFLAG_DO_DEDUP: u8 = 2;
/// The record length is stored once and reused (htscodecs' `fixed_len`).
const PFLAG_FIXED_LEN: u8 = 4;
const PFLAG_DO_SEL: u8 = 8;
const PFLAG_HAVE_QMAP: u8 = 16;
const PFLAG_HAVE_PTAB: u8 = 32;
const PFLAG_HAVE_DTAB: u8 = 64;
const PFLAG_HAVE_QTAB: u8 = 128;

/// The quality context is 16 bits wide.
const CONTEXTS: usize = 1 << 16;

/// Positions at or past this share the last `ptab` entry.
const PTAB_SIZE: usize = 1024;
/// Deltas at or past this share the last `dtab` entry.
const DTAB_SIZE: usize = 256;
/// `qtab`, `qmap` and `stab` cover every byte value.
const BYTE_TABLE: usize = 256;

/// The most first-level run lengths `read_array` holds (htscodecs' `R[1024]`).
const MAX_RUNS: usize = 1024;

/// A stream that claims a large output only gets this much output capacity
/// per input byte up front; the rest grows as records decode. Quality data
/// compresses 2–10×, so real blocks rarely grow; a tiny corrupt stream
/// claiming hundreds of MiB allocates little before it fails.
const OUTPUT_CAPACITY_PER_INPUT_BYTE: usize = 32;
const MIN_OUTPUT_CAPACITY: usize = 1 << 16;

/// Decode an fqzcomp block into the concatenated qualities of its records.
///
/// The output size comes from the stream; [`super::block`] checks it against
/// the block header's uncompressed size.
// r[impl cram.codec.fqzcomp]
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    let mut cur = src;
    let out_len = read_uint7(&mut cur).map_err(|e| match e {
        Uint7Error::Truncated => CramError::Truncated { context: "fqzcomp output size" },
        Uint7Error::Overflow => CramError::Uint7Overflow,
    })?;
    let out_len = out_len as usize;
    // r[impl cram.codec.fqzcomp.alloc]
    check_alloc_size(out_len, "fqzcomp output")?;

    let params = Params::parse(&mut cur)?;
    if out_len == 0 {
        // htscodecs decodes nothing and never looks at the range coder.
        return Ok(Vec::new());
    }

    let mut rc = RangeDecoder::new(cur)?;
    let capacity = out_len
        .min(src.len().saturating_mul(OUTPUT_CAPACITY_PER_INPUT_BYTE).max(MIN_OUTPUT_CAPACITY));
    let mut out = Vec::with_capacity(capacity);
    Decoder::new(&params).run(&mut rc, out_len, &mut out)?;
    Ok(out)
}

// ── Parameters ───────────────────────────────────────────────────────

/// The global parameters and every parameter block.
#[derive(Debug)]
struct Params {
    do_rev: bool,
    /// Map the selector through `stab` (else the selector is the block index).
    have_stab: bool,
    max_sel: u8,
    stab: [u16; BYTE_TABLE],
    /// At least one; the first also drives the per-quality context update.
    blocks: Vec<ParamBlock>,
    /// The largest `max_sym` over all blocks: the quality models' size.
    max_sym: u8,
}

/// One parameter block, with `ptab` and `dtab` pre-shifted into place.
#[derive(Debug)]
struct ParamBlock {
    context: u16,
    fixed_len: bool,
    do_dedup: bool,
    do_sel: bool,
    max_sym: u8,
    qshift: u32,
    qloc: u32,
    sloc: u32,
    qmask: u32,
    qmap: [u8; BYTE_TABLE],
    qtab: [u16; BYTE_TABLE],
    /// `ptab[p] << ploc`, low 16 bits — only those reach the context.
    ptab: [u16; PTAB_SIZE],
    /// `dtab[d] << dloc`, low 16 bits.
    dtab: [u16; DTAB_SIZE],
}

impl Params {
    // r[impl cram.codec.fqzcomp.params]
    fn parse(src: &mut &[u8]) -> Result<Self, CramError> {
        let truncated = || CramError::Truncated { context: "fqzcomp parameters" };
        // htscodecs' `fqz_read_parameters` wants 10 bytes before it starts.
        if src.len() < 10 {
            return Err(truncated());
        }
        let version = read_u8(src).ok_or_else(truncated)?;
        if version != VERSION {
            return Err(CramError::FqzcompVersion { version });
        }
        let gflags = read_u8(src).ok_or_else(truncated)?;
        let nparam =
            if gflags & GFLAG_MULTI_PARAM != 0 { read_u8(src).ok_or_else(truncated)? } else { 1 };
        if nparam == 0 {
            return Err(CramError::FqzcompNoParams);
        }
        let mut max_sel = if nparam > 1 { nparam } else { 0 };

        let have_stab = gflags & GFLAG_HAVE_STAB != 0;
        let mut stab = [0u16; BYTE_TABLE];
        if have_stab {
            max_sel = read_u8(src).ok_or_else(truncated)?;
            read_array(src, &mut stab, "stab")?;
        }

        let mut blocks = Vec::new();
        let mut max_sym = 0;
        for block in 0..nparam {
            let p = ParamBlock::parse(src)?;
            if p.do_sel && max_sel == 0 {
                return Err(CramError::FqzcompSelectorWithoutRange { block });
            }
            max_sym = max_sym.max(p.max_sym);
            blocks.push(p);
        }

        Ok(Self { do_rev: gflags & GFLAG_DO_REV != 0, have_stab, max_sel, stab, blocks, max_sym })
    }
}

impl ParamBlock {
    // r[impl cram.codec.fqzcomp.param_block]
    fn parse(src: &mut &[u8]) -> Result<Self, CramError> {
        let truncated = || CramError::Truncated { context: "fqzcomp parameter block" };
        if src.len() < 7 {
            return Err(truncated());
        }
        let context = read_u16_le(src).ok_or_else(truncated)?;
        let pflags = read_u8(src).ok_or_else(truncated)?;
        let max_sym = read_u8(src).ok_or_else(truncated)?;
        let [qbits, qshift] = nibbles(read_u8(src).ok_or_else(truncated)?);
        let [qloc, sloc] = nibbles(read_u8(src).ok_or_else(truncated)?);
        let [ploc, dloc] = nibbles(read_u8(src).ok_or_else(truncated)?);

        let qmap = if pflags & PFLAG_HAVE_QMAP != 0 {
            // Symbols past the stored map read htscodecs' unset `INT_MAX`
            // entries, which it truncates to 0xFF on output.
            let mut qmap = [0xFF; BYTE_TABLE];
            let stored = split_off(src, usize::from(max_sym)).ok_or_else(truncated)?;
            for (m, &q) in qmap.iter_mut().zip(stored) {
                *m = q;
            }
            qmap
        } else {
            identity_u8()
        };

        // Without `qbits` the quality sub-context is masked away, and neither
        // htscodecs' writer nor its reader store `qtab` then, whatever the flag.
        let mut qtab = identity_u16();
        if qbits > 0 && pflags & PFLAG_HAVE_QTAB != 0 {
            read_array(src, &mut qtab, "qtab")?;
        }

        let mut ptab = [0u16; PTAB_SIZE];
        if pflags & PFLAG_HAVE_PTAB != 0 {
            read_array(src, &mut ptab, "ptab")?;
            shift_into_place(&mut ptab, ploc);
        }
        let mut dtab = [0u16; DTAB_SIZE];
        if pflags & PFLAG_HAVE_DTAB != 0 {
            read_array(src, &mut dtab, "dtab")?;
            shift_into_place(&mut dtab, dloc);
        }

        Ok(Self {
            context,
            fixed_len: pflags & PFLAG_FIXED_LEN != 0,
            do_dedup: pflags & PFLAG_DO_DEDUP != 0,
            do_sel: pflags & PFLAG_DO_SEL != 0,
            max_sym,
            qshift: u32::from(qshift),
            qloc: u32::from(qloc),
            sloc: u32::from(sloc),
            // qbits <= 15
            qmask: (1u32 << qbits).wrapping_sub(1),
            qmap,
            qtab,
            ptab,
            dtab,
        })
    }
}

/// The high and low nibble of `b`.
fn nibbles(b: u8) -> [u8; 2] {
    [b >> 4, b & 0x0F]
}

fn identity_u8() -> [u8; BYTE_TABLE] {
    std::array::from_fn(|i| u8::try_from(i).unwrap_or(u8::MAX))
}

fn identity_u16() -> [u16; BYTE_TABLE] {
    std::array::from_fn(|i| u16::try_from(i).unwrap_or(u16::MAX))
}

/// The low 16 bits of `x`: all the 16-bit context sees of a table entry.
#[inline]
fn low16(x: u32) -> u16 {
    let [a, b, _, _] = x.to_le_bytes();
    u16::from_le_bytes([a, b])
}

/// Shift every entry left by `loc` (a nibble), as htscodecs does once after
/// parsing so the per-quality update is a plain add.
fn shift_into_place(table: &mut [u16], loc: u8) {
    for v in table {
        *v = low16(u32::from(*v) << loc);
    }
}

/// Read a table stored as htscodecs' `read_array`: run lengths per output
/// value (runs of 255 continue), themselves run-length coded by following a
/// repeated byte with a count of further copies.
// r[impl cram.codec.fqzcomp.array]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "sums of at most 1024 run bytes times 255 fit in usize; indices are bounds-checked"
)]
fn read_array(src: &mut &[u8], out: &mut [u16], table: &'static str) -> Result<(), CramError> {
    let malformed = || CramError::FqzcompMalformedTable { table };
    let size = out.len().min(MAX_RUNS);
    let bytes = *src;

    // Level one: expand the repeat counts into run lengths.
    let mut runs = [0u8; MAX_RUNS];
    let mut n_runs = 0usize;
    let mut total = 0usize;
    let mut last = None;
    let mut i = 0usize;
    while total < size {
        let Some(&run) = bytes.get(i) else { break };
        // n_runs < MAX_RUNS here: checked at the end of every iteration.
        *runs.get_mut(n_runs).ok_or_else(malformed)? = run;
        n_runs += 1;
        total += usize::from(run);
        if last == Some(run) {
            i += 1;
            let &copies = bytes.get(i).ok_or_else(malformed)?;
            total += usize::from(run) * usize::from(copies);
            // htscodecs stops copying once the sum passes the table size,
            // though it already added all copies to the sum.
            let mut copies = copies;
            while copies > 0 && total <= size && n_runs < MAX_RUNS {
                *runs.get_mut(n_runs).ok_or_else(malformed)? = run;
                n_runs += 1;
                copies -= 1;
            }
        }
        if n_runs >= MAX_RUNS {
            return Err(malformed());
        }
        last = Some(run);
        i += 1;
    }
    let consumed = i;
    let runs = runs.get(..n_runs).ok_or_else(malformed)?;

    // Level two: value `v` fills the next run of entries.
    let mut next_run = runs.iter();
    let mut filled = 0usize;
    let mut value = 0u16;
    while filled < size {
        let mut run_len = 0usize;
        loop {
            let &part = next_run.next().ok_or_else(malformed)?;
            run_len += usize::from(part);
            if part != 255 {
                break;
            }
        }
        let end = (filled + run_len).min(size);
        out.get_mut(filled..end).ok_or_else(malformed)?.fill(value);
        filled = end;
        // At most one value per run byte, so value < MAX_RUNS.
        value += 1;
    }

    *src = bytes.get(consumed..).ok_or_else(malformed)?;
    Ok(())
}

// ── Models ───────────────────────────────────────────────────────────

/// The 2^16 quality models, each over `nsym` symbols, created on first use.
///
/// Models live back to back in one arena, each a header entry (its `freq`
/// holds the model's total, which stays <= `MAX_FREQ` < 2^16) followed by
/// `nsym` symbol entries. A block reaches a small fraction of the contexts,
/// so the arena stays small and dense, where htscodecs allocates and
/// initialises all 2^16 models (67 MB at 256 symbols) up front.
// r[impl cram.codec.fqzcomp.models]
#[derive(Debug)]
struct QualModels {
    /// Per context, the arena offset just past its model; 0 if not created.
    ends: Vec<u32>,
    arena: Vec<SymFreq>,
    nsym: u16,
}

impl QualModels {
    fn new(nsym: u16) -> Self {
        Self { ends: vec![0; CONTEXTS], arena: Vec::new(), nsym }
    }

    /// Decode one quality with the model for `ctx`.
    #[inline]
    fn decode(&mut self, ctx: u16, rc: &mut RangeDecoder<'_>) -> Option<u16> {
        let stride = usize::from(self.nsym).wrapping_add(1);
        let end = match *self.ends.get(usize::from(ctx))? {
            0 => self.create(ctx)?,
            end => end as usize,
        };
        let model = self.arena.get_mut(end.checked_sub(stride)?..end)?;
        let (head, syms) = model.split_first_mut()?;
        let mut total = u32::from(head.freq);
        let sym = decode_symbol(rc, &mut total, syms)?;
        head.freq = u16::try_from(total).ok()?;
        Some(sym)
    }

    /// Append a fresh model for `ctx` (every symbol at frequency 1) and
    /// return its end offset.
    #[cold]
    #[inline(never)]
    fn create(&mut self, ctx: u16) -> Option<usize> {
        self.arena.push(SymFreq { freq: self.nsym, sym: 0 });
        self.arena.extend((0..self.nsym).map(|sym| SymFreq { freq: 1, sym }));
        let end = self.arena.len();
        *self.ends.get_mut(usize::from(ctx))? = u32::try_from(end).ok()?;
        Some(end)
    }
}

/// The per-record state the quality context is computed from.
#[derive(Debug, Default)]
struct QualState {
    qctx: u32,
    /// Qualities left in the record, counted before the current one is taken off.
    p: u32,
    delta: u32,
    prevq: u8,
}

impl ParamBlock {
    /// Fold quality `q` into the state and return the next context.
    // r[impl cram.codec.fqzcomp.context]
    #[inline]
    fn update(&self, st: &mut QualState, q: u8, sel: u32) -> u16 {
        let qi = usize::from(q);
        st.qctx = st
            .qctx
            .wrapping_shl(self.qshift)
            .wrapping_add(u32::from(self.qtab.get(qi).copied().unwrap_or(0)));
        let pos = (st.p as usize).min(PTAB_SIZE.wrapping_sub(1));
        let delta = (st.delta as usize).min(DTAB_SIZE.wrapping_sub(1));
        let ctx = (st.qctx & self.qmask)
            .wrapping_shl(self.qloc)
            .wrapping_add(u32::from(self.ptab.get(pos).copied().unwrap_or(0)))
            .wrapping_add(u32::from(self.dtab.get(delta).copied().unwrap_or(0)))
            .wrapping_add(sel.wrapping_shl(self.sloc));
        st.delta = st.delta.wrapping_add(u32::from(st.prevq != q));
        st.prevq = q;
        st.p = st.p.wrapping_sub(1);
        low16(ctx)
    }
}

// ── Decoding ─────────────────────────────────────────────────────────

/// All models of one stream.
struct Decoder<'p> {
    params: &'p Params,
    qual: QualModels,
    len: [AdaptiveModel<256>; 4],
    rev: AdaptiveModel<2>,
    dup: AdaptiveModel<2>,
    sel: AdaptiveModel<256>,
}

impl<'p> Decoder<'p> {
    fn new(params: &'p Params) -> Self {
        Self {
            params,
            qual: QualModels::new(u16::from(params.max_sym).wrapping_add(1)),
            len: std::array::from_fn(|_| AdaptiveModel::new(256)),
            rev: AdaptiveModel::new(2),
            dup: AdaptiveModel::new(2),
            // Only decoded from when a block has `do_sel`, which needs max_sel > 0.
            sel: AdaptiveModel::new(usize::from(params.max_sel).wrapping_add(1)),
        }
    }

    fn run(
        &mut self,
        rc: &mut RangeDecoder<'_>,
        out_len: usize,
        out: &mut Vec<u8>,
    ) -> Result<(), CramError> {
        let params = self.params;
        // htscodecs passes the first block to the context update and qmap
        // for every record, whichever block the selector chose.
        let first = params.blocks.first().ok_or_else(|| CramError::FqzcompNoParams)?;
        let corrupt = || CramError::FqzcompRangeCoder;

        let mut last_len: Option<u32> = None;
        // Reversal waits until the end: a duplicate copies the bytes as decoded.
        let mut reversed: Vec<(usize, usize)> = Vec::new();

        while out.len() < out_len {
            // r[impl cram.codec.fqzcomp.selector]
            let sel = if first.do_sel { self.sel.decode(rc).ok_or_else(corrupt)? } else { 0 };
            let index = if params.have_stab {
                params.stab.get(usize::from(sel)).copied().unwrap_or(0)
            } else {
                sel
            };
            let block = params.blocks.get(usize::from(index)).ok_or_else(|| {
                CramError::FqzcompParamOutOfRange { index, nparam: params.blocks.len() }
            })?;

            // r[impl cram.codec.fqzcomp.record]
            let len = match last_len {
                Some(len) if block.fixed_len => len,
                _ => {
                    let len = self.decode_len(rc).ok_or_else(corrupt)?;
                    last_len = Some(len);
                    len
                }
            };
            let remaining = out_len.wrapping_sub(out.len());
            let len_usize = len as usize;
            if len == 0 || len_usize > remaining {
                return Err(CramError::FqzcompRecordLength { len, remaining });
            }
            let start = out.len();

            // r[impl cram.codec.fqzcomp.reverse]
            if params.do_rev && self.rev.decode(rc).ok_or_else(corrupt)? != 0 {
                reversed.push((start, len_usize));
            }

            if block.do_dedup && self.dup.decode(rc).ok_or_else(corrupt)? != 0 {
                let from = start
                    .checked_sub(len_usize)
                    .ok_or_else(|| CramError::FqzcompDuplicateTooLong { len, available: start })?;
                out.extend_from_within(from..start);
                continue;
            }

            self.decode_qualities(rc, first, block.context, len, u32::from(sel), out)
                .ok_or_else(corrupt)?;
        }

        for (start, len) in reversed {
            if let Some(rec) = out.get_mut(start..start.wrapping_add(len)) {
                rec.reverse();
            }
        }
        Ok(())
    }

    /// A record length: four bytes, least significant first, one model each.
    fn decode_len(&mut self, rc: &mut RangeDecoder<'_>) -> Option<u32> {
        let mut len = 0u32;
        for (shift, model) in (0u32..).step_by(8).zip(&mut self.len) {
            len |= u32::from(model.decode(rc)?).wrapping_shl(shift);
        }
        Some(len)
    }

    /// Decode `len` qualities of one record, starting from `context`.
    #[inline]
    fn decode_qualities(
        &mut self,
        rc: &mut RangeDecoder<'_>,
        first: &ParamBlock,
        context: u16,
        len: u32,
        sel: u32,
        out: &mut Vec<u8>,
    ) -> Option<()> {
        let mut st = QualState { p: len, ..QualState::default() };
        let mut ctx = context;
        for _ in 0..len {
            // Quality models have max_sym + 1 <= 256 symbols.
            let q = u8::try_from(self.qual.decode(ctx, rc)?).ok()?;
            ctx = first.update(&mut st, q, sel);
            out.push(first.qmap.get(usize::from(q)).copied().unwrap_or(0xFF));
        }
        Some(())
    }
}

#[cfg(test)]
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code with fixed-size inputs"
)]
mod tests {
    use super::*;

    fn array(bytes: &[u8], size: usize) -> Result<(Vec<u16>, usize), CramError> {
        let mut src = bytes;
        let mut out = vec![0u16; size];
        read_array(&mut src, &mut out, "test")?;
        Ok((out, bytes.len() - src.len()))
    }

    // r[verify cram.codec.fqzcomp.array]
    #[test]
    fn read_array_spec_example() {
        // CRAMcodecs §6 ReadArray: A = {0,1,3,4,5,6,7,7,7,7} has run lengths
        // R = {1,1,0,1,1,1,1,4}, stored as R2 = {1,1,+0,0,1,1,+2,4}.
        let (a, used) = array(&[1, 1, 0, 0, 1, 1, 2, 4, 99], 10).unwrap();
        assert_eq!(a, [0, 1, 3, 4, 5, 6, 7, 7, 7, 7]);
        assert_eq!(used, 8, "the trailing byte belongs to what follows");
    }

    #[test]
    fn read_array_long_runs_continue_past_255() {
        // CRAMcodecs: a run of 600 is 255 255 90, stored as 255 255 +0 90.
        // Then 424 of value 1: 255 169.
        let (a, used) = array(&[255, 255, 0, 90, 255, 169], 1024).unwrap();
        assert!(a[..600].iter().all(|&v| v == 0));
        assert!(a[600..].iter().all(|&v| v == 1));
        assert_eq!(used, 6);
    }

    #[test]
    fn read_array_htscodecs_all_zero_stab() {
        // The `stab` of the hts-specs q4.0 vector: one run of 256 zeros.
        let (a, used) = array(&[0xff, 0x01], 256).unwrap();
        assert!(a.iter().all(|&v| v == 0));
        assert_eq!(used, 2);
    }

    #[test]
    fn read_array_rejects_malformed() {
        // A repeated byte with no count after it.
        assert!(matches!(array(&[1, 1], 256), Err(CramError::FqzcompMalformedTable { .. })));
        // Run lengths end on a 255 continuation.
        assert!(matches!(array(&[255], 256), Err(CramError::FqzcompMalformedTable { .. })));
        // Input ends before the table is full.
        assert!(matches!(array(&[3, 4], 256), Err(CramError::FqzcompMalformedTable { .. })));
        // 1024 runs of zero never fill the table: too many first-level runs.
        assert!(matches!(
            array(&[0, 0, 255, 0, 255, 0, 255, 0, 255], 256),
            Err(CramError::FqzcompMalformedTable { .. })
        ));
    }

    // r[verify cram.codec.fqzcomp.params]
    // r[verify cram.codec.fqzcomp.param_block]
    #[test]
    fn parses_the_q4_vector_parameters() {
        let src = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/codecs/fqzcomp/q4.0"
        ))
        .unwrap();
        let mut cur = src.as_slice();
        // 0x89 0x9b 0x58 = (9 << 14) | (0x1b << 7) | 0x58
        assert_eq!(read_uint7(&mut cur), Ok(151_000));
        let p = Params::parse(&mut cur).unwrap();
        // Read off the bytes: 05 | 02 | 03 ff 01 | 00 00 7c 04 22 0d af | 02 0c 12 24 | ...
        assert!(p.have_stab && !p.do_rev);
        assert_eq!(p.max_sel, 3);
        assert_eq!(p.stab, [0; 256]);
        assert_eq!(p.blocks.len(), 1);
        let b = &p.blocks[0];
        assert_eq!(b.context, 0);
        // pflags 0x7c: dtab, ptab, qmap, sel, fixed length; no dedup, no qtab.
        assert!(b.fixed_len && b.do_sel && !b.do_dedup);
        assert_eq!(b.max_sym, 4);
        assert_eq!((b.qmask, b.qshift, b.qloc, b.sloc), (0b11, 2, 0, 13));
        assert_eq!(&b.qmap[..5], &[2, 12, 18, 36, 0xFF]);
        assert_eq!(b.qtab, identity_u16());
        // ploc = 10, dloc = 15: every entry is a multiple of 2^10 / 2^15.
        assert!(b.ptab.iter().all(|&v| v % (1 << 10) == 0));
        assert!(b.ptab.iter().any(|&v| v != 0));
        assert!(b.dtab.iter().all(|&v| v % (1 << 15) == 0));
    }

    #[test]
    fn rejects_other_versions_and_empty_blocks() {
        let mut v4: &[u8] = &[4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        assert!(matches!(Params::parse(&mut v4), Err(CramError::FqzcompVersion { version: 4 })));
        let mut none: &[u8] = &[5, GFLAG_MULTI_PARAM, 0, 0, 0, 0, 0, 0, 0, 0];
        assert!(matches!(Params::parse(&mut none), Err(CramError::FqzcompNoParams)));
        // do_sel with no selector range.
        let mut sel: &[u8] = &[5, 0, 0, 0, PFLAG_DO_SEL, 1, 0, 0, 0, 0, 0];
        assert!(matches!(
            Params::parse(&mut sel),
            Err(CramError::FqzcompSelectorWithoutRange { block: 0 })
        ));
    }

    #[test]
    fn empty_output_needs_no_range_coder() {
        // Size 0, then ten bytes of a minimal parameter block and nothing else.
        let src = [0, 5, 0, 0, 0, 0, 1, 0, 0, 0, 0];
        assert_eq!(decode(&src).unwrap(), Vec::<u8>::new());
    }

    // r[verify cram.codec.fqzcomp.alloc]
    #[test]
    fn output_size_is_bounded() {
        // 300 000 000 as uint7: past the allocation limit.
        let src = [0x81, 0x8f, 0x86, 0xc6, 0x00, 5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        assert!(matches!(decode(&src), Err(CramError::AllocationTooLarge { .. })));

        // 200 MiB claimed by a 20-byte stream: only a small buffer up front.
        let src = [0xe4, 0x80, 0x80, 0x00, 5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let mut cur = &src[..];
        assert_eq!(read_uint7(&mut cur), Ok(200 << 20));
        let params = Params::parse(&mut cur).unwrap();
        let mut rc = RangeDecoder::new(cur).unwrap();
        let mut out = Vec::with_capacity(MIN_OUTPUT_CAPACITY);
        assert!(Decoder::new(&params).run(&mut rc, 200 << 20, &mut out).is_err());
        assert!(out.capacity() < 1 << 20, "grew to {}", out.capacity());
    }

    // r[verify cram.codec.fqzcomp.models]
    #[test]
    fn quality_models_are_created_on_first_use() {
        let src = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/codecs/fqzcomp/q40+dir.0"
        ))
        .unwrap();
        let mut cur = src.as_slice();
        let out_len = read_uint7(&mut cur).unwrap() as usize;
        let params = Params::parse(&mut cur).unwrap();
        let mut rc = RangeDecoder::new(cur).unwrap();
        let mut decoder = Decoder::new(&params);
        let mut out = Vec::new();
        decoder.run(&mut rc, out_len, &mut out).unwrap();

        let stride = usize::from(params.max_sym) + 2;
        let created = decoder.qual.ends.iter().filter(|&&e| e != 0).count();
        assert_eq!(decoder.qual.arena.len(), created * stride);
        assert!(created < CONTEXTS / 2, "{created} of {CONTEXTS} contexts");
    }
}
