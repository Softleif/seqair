//! A plain fqzcomp decoder, transliterated from the `CRAMcodecs` §6
//! pseudocode, for testing [`super::fqzcomp`] against.
//!
//! It shares no code with the production decoder — not even the range coder
//! ([`RangeDecoder`] and [`Model`] here follow `CRAMcodecs` §4 "Range coding"
//! and "Adaptive Modelling" line by line) — and trades every bit of speed for
//! being obviously right: models in a `HashMap`, tables kept as read and
//! shifted at use, state in plain variables named as in the pseudocode.
//!
//! Where the pseudocode disagrees with htscodecs' `fqzcomp_qual.c`, which
//! writes the files, this follows htscodecs; every such place says so in a
//! `htscodecs:` comment. Beyond the format it applies one crate policy, the
//! output size limit, so that it fails exactly where the production decoder
//! does.
//!
//! Every failure is `None`: the reference only has to say *whether* a stream
//! decodes, and to what.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "a transliteration of the pseudocode; each index and sum is bounded by the checks \
              beside it, and the differential tests and fuzz target run it on arbitrary input"
)]

use std::collections::HashMap;

use super::reader::MAX_ALLOC_SIZE;

/// Reads bytes off the front of the input.
struct Input<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Input<'_> {
    /// `ReadUint8`.
    fn read_u8(&mut self) -> Option<u8> {
        let b = *self.bytes.get(self.pos)?;
        self.pos = self.pos.checked_add(1)?;
        Some(b)
    }

    /// `ReadUint16`, little-endian.
    fn read_u16(&mut self) -> Option<u16> {
        let lo = self.read_u8()?;
        let hi = self.read_u8()?;
        Some(u16::from(lo) | (u16::from(hi) << 8))
    }

    /// `ReadUint7`: 7 bits per byte, most significant group first, the top
    /// bit set on every byte but the last; at most 5 bytes for 32 bits.
    fn read_uint7(&mut self) -> Option<u32> {
        let mut value: u64 = 0;
        for _ in 0..5 {
            let b = self.read_u8()?;
            value = (value << 7) | u64::from(b & 0x7F);
            if b & 0x80 == 0 {
                return u32::try_from(value & 0xFFFF_FFFF).ok();
            }
        }
        None
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.pos)
    }
}

// ── CRAMcodecs §4 "Range coding" ─────────────────────────────────────

/// The range decoder: `range` and `code`, 32-bit unsigned.
pub struct RangeDecoder<'a> {
    input: Input<'a>,
    range: u32,
    code: u32,
}

impl<'a> RangeDecoder<'a> {
    /// `RangeDecodeCreate`.
    pub fn new(bytes: &'a [u8]) -> Option<Self> {
        let mut input = Input { bytes, pos: 0 };
        let mut code: u64 = 0;
        for _ in 0..5 {
            code = (code << 8) + u64::from(input.read_u8()?);
        }
        let code = u32::try_from(code & 0xFFFF_FFFF).ok()?;
        Some(Self { input, range: u32::MAX, code })
    }

    /// `RangeGetFreq`.
    fn get_freq(&mut self, tot_freq: u32) -> Option<u32> {
        self.range = self.range.checked_div(tot_freq)?;
        self.code.checked_div(self.range)
    }

    /// `RangeDecode`. 32-bit unsigned arithmetic wraps, as in C.
    fn decode(&mut self, sym_low: u32, sym_freq: u32) -> Option<()> {
        self.code = self.code.wrapping_sub(sym_low.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(sym_freq);
        while self.range < (1 << 24) {
            self.range <<= 8;
            self.code = (self.code << 8).wrapping_add(u32::from(self.input.read_u8()?));
        }
        Some(())
    }
}

/// An adaptive model: symbols `S` and frequencies `F`, kept roughly most
/// frequent first.
#[derive(Clone)]
pub struct Model {
    total_freq: u32,
    s: Vec<u16>,
    f: Vec<u32>,
}

impl Model {
    /// `ModelCreate`.
    pub fn new(num_sym: usize) -> Self {
        let mut s = Vec::new();
        let mut f = Vec::new();
        for i in 0..num_sym {
            s.push(u16::try_from(i).unwrap_or(u16::MAX));
            f.push(1);
        }
        Self { total_freq: u32::try_from(num_sym).unwrap_or(u32::MAX), s, f }
    }

    /// `ModelDecode`. A coded value past every symbol is corrupt input.
    pub fn decode(&mut self, rc: &mut RangeDecoder<'_>) -> Option<u16> {
        let freq = rc.get_freq(self.total_freq)?;
        let mut x = 0;
        let mut acc: u32 = 0;
        while acc.checked_add(*self.f.get(x)?)? <= freq {
            acc = acc.checked_add(*self.f.get(x)?)?;
            x += 1;
        }
        rc.decode(acc, *self.f.get(x)?)?;
        *self.f.get_mut(x)? += 16;
        self.total_freq += 16;
        if self.total_freq > (1 << 16) - 17 {
            self.renormalise();
        }
        let sym = *self.s.get(x)?;
        if x > 0 && self.f.get(x)? > self.f.get(x - 1)? {
            self.f.swap(x, x - 1);
            self.s.swap(x, x - 1);
        }
        Some(sym)
    }

    /// `ModelRenormalise`.
    fn renormalise(&mut self) {
        self.total_freq = 0;
        for f in &mut self.f {
            *f -= *f / 2;
            self.total_freq += *f;
        }
    }
}

// ── CRAMcodecs §6 "FQZComp quality codec" ────────────────────────────

/// One parameter block, as `FQZDecodeSingleParam` reads it.
struct Param {
    context: u16,
    flags: u8,
    max_sym: u8,
    qbits: u32,
    qshift: u32,
    qloc: u32,
    sloc: u32,
    ploc: u32,
    dloc: u32,
    /// `None` without `have_qmap`.
    qmap: Option<Vec<u8>>,
    qtab: Vec<u32>,
    ptab: Vec<u32>,
    dtab: Vec<u32>,
}

impl Param {
    fn do_dedup(&self) -> bool {
        self.flags & 2 != 0
    }
    /// htscodecs: flag 4 means the length is stored once (`fixed_len`); the
    /// pseudocode calls it `do_len` and reads it the other way round.
    fn fixed_len(&self) -> bool {
        self.flags & 4 != 0
    }
    fn do_sel(&self) -> bool {
        self.flags & 8 != 0
    }
    fn have_qmap(&self) -> bool {
        self.flags & 16 != 0
    }
    fn have_ptab(&self) -> bool {
        self.flags & 32 != 0
    }
    fn have_dtab(&self) -> bool {
        self.flags & 64 != 0
    }
    fn have_qtab(&self) -> bool {
        self.flags & 128 != 0
    }
}

/// `ReadArray(n)`.
fn read_array(input: &mut Input<'_>, n: usize) -> Option<Vec<u32>> {
    // Convert R2 to R.
    let mut r: Vec<u8> = Vec::new();
    let mut z: usize = 0;
    let mut last: Option<u8> = None;
    while z < n {
        let run = input.read_u8()?;
        r.push(run);
        z += usize::from(run);
        if Some(run) == last {
            let copy = input.read_u8()?;
            // htscodecs adds all copies to z first, then pushes copies only
            // while z <= n (so none if they overshoot) and R has room.
            z += usize::from(run) * usize::from(copy);
            for _ in 0..copy {
                if z > n || r.len() >= 1024 {
                    break;
                }
                r.push(run);
            }
        }
        // htscodecs: R holds at most 1023 entries.
        if r.len() >= 1024 {
            return None;
        }
        last = Some(run);
    }

    // Convert R to A.
    let mut a = vec![0u32; n];
    let mut i: u32 = 0;
    let mut j: usize = 0;
    let mut z: usize = 0;
    while z < n {
        let mut run_len: usize = 0;
        loop {
            let part = *r.get(j)?;
            j += 1;
            run_len += usize::from(part);
            if part != 255 {
                break;
            }
        }
        for _ in 0..run_len {
            // htscodecs stops at the table's end.
            if z >= n {
                break;
            }
            a[z] = i;
            z += 1;
        }
        i += 1;
    }
    Some(a)
}

/// `FQZDecodeSingleParam`.
fn decode_single_param(input: &mut Input<'_>) -> Option<Param> {
    // htscodecs: fewer than 7 bytes left is an error up front.
    if input.remaining() < 7 {
        return None;
    }
    let context = input.read_u16()?;
    let flags = input.read_u8()?;
    let max_sym = input.read_u8()?;
    let x = input.read_u8()?;
    let (qbits, qshift) = (u32::from(x / 16), u32::from(x % 16));
    let x = input.read_u8()?;
    let (qloc, sloc) = (u32::from(x / 16), u32::from(x % 16));
    let x = input.read_u8()?;
    let (ploc, dloc) = (u32::from(x / 16), u32::from(x % 16));

    let mut p = Param {
        context,
        flags,
        max_sym,
        qbits,
        qshift,
        qloc,
        sloc,
        ploc,
        dloc,
        qmap: None,
        qtab: (0..256).collect(),
        ptab: vec![0; 1024],
        dtab: vec![0; 256],
    };
    if p.have_qmap() {
        let mut qmap = Vec::new();
        for _ in 0..max_sym {
            qmap.push(input.read_u8()?);
        }
        p.qmap = Some(qmap);
    }
    // htscodecs: qtab is only stored when qbits > 0.
    if p.have_qtab() && p.qbits > 0 {
        p.qtab = read_array(input, 256)?;
    }
    // htscodecs: without the flag the table is all zeros, which adds
    // nothing — the same as the pseudocode leaving the term out.
    if p.have_ptab() {
        p.ptab = read_array(input, 1024)?;
    }
    if p.have_dtab() {
        p.dtab = read_array(input, 256)?;
    }
    Some(p)
}

/// The global parameters, as `FQZDecodeParams` reads them.
struct Params {
    gflags: u8,
    max_sel: u8,
    stab: Vec<u32>,
    params: Vec<Param>,
    max_sym: u8,
}

/// `FQZDecodeParams`.
fn decode_params(input: &mut Input<'_>) -> Option<Params> {
    // htscodecs: fewer than 10 bytes left is an error up front.
    if input.remaining() < 10 {
        return None;
    }
    let vers = input.read_u8()?;
    if vers != 5 {
        return None;
    }
    let gflags = input.read_u8()?;
    let (nparam, mut max_sel) = if gflags & 1 != 0 {
        let nparam = input.read_u8()?;
        // htscodecs: max_sel is nparam only when there is more than one.
        (nparam, if nparam > 1 { nparam } else { 0 })
    } else {
        (1, 0)
    };
    // htscodecs: no parameter blocks is an error.
    if nparam == 0 {
        return None;
    }
    let mut stab = vec![0; 256];
    if gflags & 2 != 0 {
        max_sel = input.read_u8()?;
        stab = read_array(input, 256)?;
    }
    let mut max_sym = 0;
    let mut params = Vec::new();
    for _ in 0..nparam {
        let p = decode_single_param(input)?;
        // htscodecs: a selector with no selector values is an error.
        if p.do_sel() && max_sel == 0 {
            return None;
        }
        if max_sym < p.max_sym {
            max_sym = p.max_sym;
        }
        params.push(p);
    }
    Some(Params { gflags, max_sel, stab, params, max_sym })
}

/// The models of `FQZCreateModels`; quality models are created on first
/// use, which is the same as creating all 2^16 up front.
struct Models {
    qual: HashMap<u16, Model>,
    qual_syms: usize,
    len: [Model; 4],
    dup: Model,
    rev: Model,
    sel: Model,
}

/// `DecodeLength`.
fn decode_length(models: &mut Models, rc: &mut RangeDecoder<'_>) -> Option<u32> {
    let [l0, l1, l2, l3] = &mut models.len;
    let mut rec_len = u32::from(l0.decode(rc)?);
    rec_len += u32::from(l1.decode(rc)?) << 8;
    rec_len += u32::from(l2.decode(rc)?) << 16;
    rec_len += u32::from(l3.decode(rc)?) << 24;
    Some(rec_len)
}

/// The variables `FQZUpdateContext` updates.
struct QualState {
    qctx: u32,
    pos: u32,
    delta: u32,
    prevq: u16,
    sel: u32,
}

/// `FQZUpdateContext`, with 32-bit unsigned arithmetic wrapping as in C.
fn update_context(params: &Param, st: &mut QualState, q: u16) -> u16 {
    // htscodecs starts from 0 (`last = 0; // pm->context`): the starting
    // context is only the record's first context. The pseudocode starts
    // every update from `params.context`.
    let mut ctx: u32 = 0;
    let qtab_q = params.qtab.get(usize::from(q)).copied().unwrap_or(0);
    st.qctx = (st.qctx << params.qshift).wrapping_add(qtab_q);
    let qmask = (1u32 << params.qbits) - 1;
    ctx = ctx.wrapping_add((st.qctx & qmask) << params.qloc);
    if params.have_ptab() {
        let p = st.pos.min(1023) as usize;
        ctx = ctx.wrapping_add(params.ptab.get(p).copied().unwrap_or(0) << params.ploc);
    }
    if params.have_dtab() {
        let d = st.delta.min(255) as usize;
        ctx = ctx.wrapping_add(params.dtab.get(d).copied().unwrap_or(0) << params.dloc);
        if st.prevq != q {
            st.delta += 1;
        }
        st.prevq = q;
    }
    if params.do_sel() {
        ctx = ctx.wrapping_add(st.sel << params.sloc);
    }
    u16::try_from(ctx & 0xFFFF).unwrap_or(0)
}

/// `FQZDecode`: decode an fqzcomp stream, or `None` if it is invalid.
pub fn decode(src: &[u8]) -> Option<Vec<u8>> {
    let mut input = Input { bytes: src, pos: 0 };
    let buf_len = input.read_uint7()? as usize;
    // Crate policy, not the format: the production decoder refuses this too.
    if buf_len > MAX_ALLOC_SIZE {
        return None;
    }
    let params = decode_params(&mut input)?;
    // htscodecs: an empty block decodes without starting the range coder.
    if buf_len == 0 {
        return Some(Vec::new());
    }

    // FQZCreateModels.
    let mut rc = RangeDecoder::new(input.bytes.get(input.pos..)?)?;
    let mut models = Models {
        qual: HashMap::new(),
        qual_syms: usize::from(params.max_sym) + 1,
        len: std::array::from_fn(|_| Model::new(256)),
        dup: Model::new(2),
        rev: Model::new(2),
        sel: Model::new(usize::from(params.max_sel) + 1),
    };

    // htscodecs updates the context with, and maps qualities through, the
    // first parameter block for every record; the selected block only gives
    // the starting context and the fixed-length and dedup flags.
    let first = params.params.first()?;
    let do_rev = params.gflags & 4 != 0;
    let have_stab = params.gflags & 2 != 0;

    let mut output: Vec<u8> = Vec::new();
    // htscodecs: `first_len` and `last_len` are global, not per block.
    let mut first_len = true;
    let mut last_len: u32 = 0;
    let mut rev: Vec<bool> = Vec::new();
    let mut len: Vec<usize> = Vec::new();
    let mut st = QualState { qctx: 0, pos: 0, delta: 0, prevq: 0, sel: 0 };

    let mut i: usize = 0;
    while i < buf_len {
        // FQZNewRecord.
        // htscodecs: the selector is read if the first block has do_sel
        // (the pseudocode: if max_sel > 0).
        st.sel = 0;
        if first.do_sel() {
            st.sel = u32::from(models.sel.decode(&mut rc)?);
        }
        // htscodecs: without stab the selector is the block index itself
        // (the pseudocode: block 0).
        let x = if have_stab { *params.stab.get(st.sel as usize)? } else { st.sel };
        let param = params.params.get(x as usize)?;

        let rec_len = if !param.fixed_len() || first_len {
            let rec_len = decode_length(&mut models, &mut rc)?;
            last_len = rec_len;
            first_len = false;
            rec_len
        } else {
            last_len
        };
        // htscodecs: a record must be non-empty and fit in the output.
        if rec_len == 0 || rec_len as usize > buf_len - i {
            return None;
        }
        let rec_len = rec_len as usize;
        st.pos = u32::try_from(rec_len).ok()?;

        // htscodecs: do_rev is a global flag.
        if do_rev {
            rev.push(models.rev.decode(&mut rc)? != 0);
            len.push(rec_len);
        }

        let mut is_dup = false;
        if param.do_dedup() && models.dup.decode(&mut rc)? > 0 {
            is_dup = true;
        }
        if is_dup {
            // htscodecs: there must be rec_len bytes to copy.
            if rec_len > i {
                return None;
            }
            for j in 0..rec_len {
                output.push(output[i + j - rec_len]);
            }
            i += rec_len;
            continue;
        }
        st.qctx = 0;
        st.delta = 0;
        st.prevq = 0;
        let mut ctx = param.context;

        for _ in 0..rec_len {
            let qual_syms = models.qual_syms;
            let model = models.qual.entry(ctx).or_insert_with(|| Model::new(qual_syms));
            let q = model.decode(&mut rc)?;
            // htscodecs: symbols past the stored qmap come out as 0xFF.
            let out = match &first.qmap {
                Some(qmap) => qmap.get(usize::from(q)).copied().unwrap_or(0xFF),
                None => u8::try_from(q).ok()?,
            };
            output.push(out);
            ctx = update_context(first, &mut st, q);
            i += 1;
            st.pos -= 1;
        }
    }

    // ReverseQualities. htscodecs walks every record; the pseudocode only
    // moves on past reversed ones.
    if do_rev {
        let mut i = 0;
        for (rec_len, rev) in len.into_iter().zip(rev) {
            if rev {
                output.get_mut(i..i + rec_len)?.reverse();
            }
            i += rec_len;
        }
    }
    Some(output)
}

#[cfg(test)]
#[path = "../../tests/support/htscodecs_fqz.rs"]
mod htscodecs_fqz;

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "test code; indices bounded by the generators"
)]
mod tests {
    use std::{fs, io::Read, path::PathBuf};

    use super::htscodecs_fqz::{builtin_stream, custom_stream};
    use crate::cram::fqzcomp;
    use hegel::TestCase;
    use hegel::generators as gs;

    fn codecs_dir() -> PathBuf {
        PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/codecs"
        ))
    }

    /// First column of every line of the gzipped original, as binary
    /// qualities (the vectors were made from ASCII phred+33).
    fn original(name: &str) -> Vec<u8> {
        let gz = fs::read(codecs_dir().join("gzip").join(format!("{name}.gz"))).unwrap();
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(gz.as_slice()).read_to_end(&mut text).unwrap();
        text.split(|&b| b == b'\n')
            .flat_map(|line| line.split(|b| b.is_ascii_whitespace()).next().unwrap_or_default())
            .map(|&ascii| ascii.checked_sub(33).expect("phred+33"))
            .collect()
    }

    // r[verify cram.codec.fqzcomp]
    #[test]
    fn reference_decodes_every_hts_specs_vector() {
        let mut n = 0;
        for entry in fs::read_dir(codecs_dir().join("fqzcomp")).unwrap() {
            let path = entry.unwrap().path();
            let file_name = path.file_name().unwrap().to_str().unwrap().to_owned();
            let (name, _) = file_name.rsplit_once('.').unwrap();
            let got = super::decode(&fs::read(&path).unwrap());
            assert!(got.as_deref() == Some(original(name).as_slice()), "{file_name}");
            n += 1;
        }
        assert_eq!(n, 16);
    }

    /// Both decoders on `bytes`: the same output, or both an error.
    fn assert_agree(bytes: &[u8]) {
        let production = fqzcomp::decode(bytes);
        let reference = super::decode(bytes);
        match (&production, &reference) {
            (Ok(p), Some(r)) => assert!(p == r, "outputs differ"),
            (Err(_), None) => {}
            _ => panic!(
                "production {:?} but reference {}",
                production.as_ref().map(Vec::len),
                if reference.is_some() { "decodes" } else { "fails" }
            ),
        }
    }

    // r[verify cram.codec.fqzcomp.params]
    // r[verify cram.codec.fqzcomp.param_block]
    // r[verify cram.codec.fqzcomp.array]
    // r[verify cram.codec.fqzcomp.models]
    // r[verify cram.codec.fqzcomp.selector]
    // r[verify cram.codec.fqzcomp.record]
    // r[verify cram.codec.fqzcomp.context]
    // r[verify cram.codec.fqzcomp.reverse]
    #[hegel::test(test_cases = 200)]
    fn both_decoders_give_back_htscodecs_input(tc: TestCase) {
        let stream = if tc.draw(gs::booleans()) { builtin_stream(&tc) } else { custom_stream(&tc) };
        tc.assume(stream.is_some());
        let (records, stream) = stream.unwrap();
        assert!(super::decode(&stream).as_ref() == Some(&records.quals));
        assert!(fqzcomp::decode(&stream).unwrap() == records.quals);
    }

    /// Damaged streams: the decoders agree on whether each decodes, and to
    /// what. They share no code, so this pins the production decoder's error
    /// paths (and its output for streams htscodecs would never write) to the
    /// pseudocode's.
    #[hegel::test(test_cases = 500)]
    fn decoders_agree_on_damaged_streams(tc: TestCase) {
        let stream = if tc.draw(gs::booleans()) { builtin_stream(&tc) } else { custom_stream(&tc) };
        tc.assume(stream.is_some());
        let (_, mut bytes) = stream.unwrap();
        match tc.draw(gs::integers::<u8>().max_value(2)) {
            // Flip bytes, mostly in the parameter block.
            0 => {
                let flips: Vec<(usize, u8)> = tc.draw(
                    gs::vecs(gs::tuples!(
                        gs::integers::<usize>().max_value(bytes.len() - 1),
                        gs::integers::<u8>().min_value(1)
                    ))
                    .min_size(1)
                    .max_size(4),
                );
                let near_start = tc.draw(gs::booleans());
                for (at, x) in flips {
                    let at = if near_start { at % 64 % bytes.len() } else { at };
                    bytes[at] ^= x;
                }
            }
            // Cut the stream short.
            1 => {
                let cut = tc.draw(gs::integers::<usize>().max_value(bytes.len()));
                bytes.truncate(cut);
            }
            // Keep a prefix and replace the rest with noise.
            _ => {
                let cut = tc.draw(gs::integers::<usize>().max_value(bytes.len()));
                bytes.truncate(cut);
                bytes.extend(tc.draw(gs::binary().max_size(256)));
            }
        }
        assert_agree(&bytes);
    }

    #[hegel::test(test_cases = 500)]
    fn decoders_agree_on_arbitrary_bytes(tc: TestCase) {
        // Half the time past the version check: a small output size, then 5.
        let mut bytes = Vec::new();
        if tc.draw(gs::booleans()) {
            bytes.extend([tc.draw(gs::integers::<u8>().max_value(0x7F)), 5]);
        }
        bytes.extend(tc.draw(gs::binary().max_size(512)));
        assert_agree(&bytes);
    }
}
