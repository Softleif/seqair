//! A reference decoder for the arithmetic coder: the pseudocode of
//! `CRAMcodecs` §4 "Range coding" (`RangeDecodeCreate`, `RangeGetFreq`,
//! `RangeDecode`, `ModelCreate`, `ModelDecode`, `ModelRenormalise`,
//! `DecodeOrder0/1`, `DecodeRLE0/1`, `DecodeStripe`, `DecodeExt`,
//! `ArithDecode`) and the rANS Nx16 `DecodePackMeta` / `DecodePack` it reuses,
//! transliterated line by line.
//!
//! It exists to be obviously right, not fast, and shares no code with
//! [`super::decode`]: tests and the fuzzer compare the two. It departs from
//! the pseudocode only where the pseudocode is wrong or cannot run, each
//! place marked `Correction`: the NoSize test (inverted in the spec), the
//! stripe length (always stored, as htscodecs writes and reads it), the
//! ORDER field (htscodecs reads only 1 as order-1), the RLE loop's `i ← run+1`
//! (`i ← i+run+1`) and its writes past the end, and a nesting bound on
//! stripes so a crafted stream cannot exhaust the stack.
//!
//! Where the input is invalid it returns `None`, which is also where the
//! pseudocode would divide by zero or read out of bounds. On such input the
//! production decoder's verdict can differ; where both succeed they must
//! produce the same bytes. The known places where one succeeds and the other
//! does not (all on streams no encoder writes):
//! - an empty body: production decodes it to nothing (htscodecs), the
//!   pseudocode still reads `max_sym` and the range coder's first bytes;
//! - a run that reaches the output length on a continuation part:
//!   production stops reading parts there (htscodecs), the pseudocode reads
//!   on;
//! - a range smaller than the model total: production decodes as htscodecs
//!   does (`RC_GetFreq` returns 0), the pseudocode divides by zero;
//! - PACK: production rejects a packed length above the unpacked one, which
//!   the pseudocode accepts, and decodes an index past the symbol map to 0,
//!   where the pseudocode has no symbol;
//! - lengths: production requires the result, and every stripe part, to have
//!   exactly the declared length, and a uint7 of at most five bytes; the
//!   pseudocode only fills what it indexes.

#![allow(
    clippy::arithmetic_side_effects,
    reason = "the pseudocode's arithmetic: counters bounded by the output length, divisors checked \
              non-zero first, and `wrapping_*` where the spec's 32-bit registers wrap"
)]

use std::io::Read;

const ORDER: u8 = 0x01;
const EXT: u8 = 0x04;
const STRIPE: u8 = 0x08;
const NO_SIZE: u8 = 0x10;
const CAT: u8 = 0x20;
const RLE: u8 = 0x40;
const PACK: u8 = 0x80;

/// Not in the spec: resource bounds, the same as the production decoder's.
const MAX_STRIPE_DEPTH: u32 = 4;
const MAX_LEN: usize = 256 * 1024 * 1024;

/// A length the decoder may allocate for.
fn bounded(len: usize) -> Option<usize> {
    (len <= MAX_LEN).then_some(len)
}

/// Decode an arithmetic-coder stream; `len` is used when the stream does
/// not store its length (`NoSize`). `None` for an invalid stream.
pub fn decode(src: &[u8], len: usize) -> Option<Vec<u8>> {
    let mut input = Input { bytes: src, pos: 0 };
    arith_decode(&mut input, bounded(len)?, 0)
}

/// The input stream the pseudocode's `ReadUint8` reads from.
struct Input<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Input<'a> {
    fn read_uint8(&mut self) -> Option<u8> {
        let b = *self.bytes.get(self.pos)?;
        self.pos = self.pos.checked_add(1)?;
        Some(b)
    }

    /// `ReadUint7`: 7 bits a byte, most significant group first, the top
    /// bit set on all but the last byte.
    fn read_uint7(&mut self) -> Option<usize> {
        let mut value: u64 = 0;
        loop {
            let c = self.read_uint8()?;
            value = (value << 7) | u64::from(c & 127);
            if value > u64::from(u32::MAX) {
                return None;
            }
            if c < 128 {
                return usize::try_from(value).ok();
            }
        }
    }

    /// `ReadData(n)`: the next `n` bytes.
    fn read_data(&mut self, n: usize) -> Option<&'a [u8]> {
        let end = self.pos.checked_add(n)?;
        let data = self.bytes.get(self.pos..end)?;
        self.pos = end;
        Some(data)
    }

    fn rest(&self) -> &'a [u8] {
        self.bytes.get(self.pos..).unwrap_or_default()
    }
}

/// `RangeDecodeCreate` and the decoder's state.
struct RangeCoder {
    range: u32,
    code: u32,
}

impl RangeCoder {
    fn create(input: &mut Input<'_>) -> Option<Self> {
        let range = u32::MAX;
        let mut code: u64 = 0;
        for _ in 0..=4 {
            code = (code << 8) + u64::from(input.read_uint8()?);
        }
        code &= 0xFFFF_FFFF;
        Some(Self { range, code: u32::try_from(code).ok()? })
    }

    /// `RangeGetFreq`. A zero `range` is where the pseudocode divides by
    /// zero.
    fn get_freq(&mut self, tot_freq: u32) -> Option<u32> {
        self.range = self.range.checked_div(tot_freq)?;
        self.code.checked_div(self.range)
    }

    /// `RangeDecode`; `code` and `range` are 32-bit registers.
    fn decode(&mut self, input: &mut Input<'_>, sym_low: u32, sym_freq: u32) -> Option<()> {
        self.code = self.code.wrapping_sub(sym_low.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(sym_freq);
        while self.range < (1 << 24) {
            self.range <<= 8;
            self.code = (self.code << 8).wrapping_add(u32::from(input.read_uint8()?));
        }
        Some(())
    }
}

/// `ModelCreate` / `ModelDecode` / `ModelRenormalise`: symbols `S` and
/// frequencies `F`, as `(S_i, F_i)` pairs.
struct Model {
    total_freq: u32,
    entries: Vec<(u16, u32)>,
}

impl Model {
    fn create(num_sym: u16) -> Self {
        Self { total_freq: u32::from(num_sym), entries: (0..num_sym).map(|i| (i, 1)).collect() }
    }

    fn decode(&mut self, rc: &mut RangeCoder, input: &mut Input<'_>) -> Option<u16> {
        let freq = rc.get_freq(self.total_freq)?;
        let mut x = 0usize;
        let mut acc = 0u32;
        while acc.checked_add(self.entries.get(x)?.1)? <= freq {
            acc = acc.checked_add(self.entries.get(x)?.1)?;
            x = x.checked_add(1)?;
        }
        let f_x = self.entries.get(x)?.1;
        rc.decode(input, acc, f_x)?;
        self.entries.get_mut(x)?.1 = f_x.checked_add(16)?;
        self.total_freq = self.total_freq.checked_add(16)?;
        if self.total_freq > (1 << 16) - 17 {
            self.renormalise();
        }
        let sym = self.entries.get(x)?.0;
        if let Some(prev) = x.checked_sub(1)
            && self.entries.get(x)?.1 > self.entries.get(prev)?.1
        {
            self.entries.swap(x, prev);
        }
        Some(sym)
    }

    fn renormalise(&mut self) {
        self.total_freq = 0;
        for (_, f) in &mut self.entries {
            *f -= *f / 2;
            self.total_freq = self.total_freq.saturating_add(*f);
        }
    }
}

/// The `max_sym` byte that starts every range-coded body (0 meaning 256).
fn read_max_sym(input: &mut Input<'_>) -> Option<u16> {
    let max_sym = input.read_uint8()?;
    Some(if max_sym == 0 { 256 } else { u16::from(max_sym) })
}

fn byte(sym: u16) -> Option<u8> {
    u8::try_from(sym).ok()
}

fn decode_order0(input: &mut Input<'_>, len: usize) -> Option<Vec<u8>> {
    let max_sym = read_max_sym(input)?;
    let mut model_lit = Model::create(max_sym);
    let mut rc = RangeCoder::create(input)?;
    let mut out = Vec::new();
    for _ in 0..len {
        out.push(byte(model_lit.decode(&mut rc, input)?)?);
    }
    Some(out)
}

fn decode_order1(input: &mut Input<'_>, len: usize) -> Option<Vec<u8>> {
    let max_sym = read_max_sym(input)?;
    let mut model_lit: Vec<Model> = (0..max_sym).map(|_| Model::create(max_sym)).collect();
    let mut rc = RangeCoder::create(input)?;
    let mut out = Vec::new();
    let mut last = 0usize;
    for _ in 0..len {
        let sym = model_lit.get_mut(last)?.decode(&mut rc, input)?;
        out.push(byte(sym)?);
        last = usize::from(sym);
    }
    Some(out)
}

/// `DecodeRLE0` (`ORDER1 == false`) and `DecodeRLE1`, which differ only in
/// the literal's context.
fn decode_rle<const ORDER1: bool>(input: &mut Input<'_>, len: usize) -> Option<Vec<u8>> {
    let max_sym = read_max_sym(input)?;
    let literal_models = if ORDER1 { max_sym } else { 1 };
    let mut model_lit: Vec<Model> = (0..literal_models).map(|_| Model::create(max_sym)).collect();
    let mut model_run: Vec<Model> = (0..=257).map(|_| Model::create(4)).collect();
    let mut rc = RangeCoder::create(input)?;
    let mut out = vec![0u8; len];
    let mut last = 0usize;
    let mut i = 0usize;
    while i < len {
        let ctx = if ORDER1 { last } else { 0 };
        let lit = model_lit.get_mut(ctx)?.decode(&mut rc, input)?;
        *out.get_mut(i)? = byte(lit)?;
        last = usize::from(lit);
        let mut part = model_run.get_mut(last)?.decode(&mut rc, input)?;
        let mut run = usize::from(part);
        let mut rctx = 256;
        while part == 3 {
            part = model_run.get_mut(rctx)?.decode(&mut rc, input)?;
            rctx = 257;
            run = run.checked_add(usize::from(part))?;
        }
        for j in 1..=run {
            // Correction: the pseudocode writes the run past the output's end.
            if let Some(o) = out.get_mut(i.checked_add(j)?) {
                *o = byte(lit)?;
            }
        }
        // Correction: the pseudocode's `i ← run+1` means `i ← i+run+1`.
        i = i.checked_add(run)?.checked_add(1)?;
    }
    Some(out)
}

fn decode_stripe(input: &mut Input<'_>, len: usize, depth: u32) -> Option<Vec<u8>> {
    if depth >= MAX_STRIPE_DEPTH {
        return None;
    }
    let n = usize::from(input.read_uint8()?);
    if n == 0 {
        // The pseudocode divides by zero.
        return None;
    }
    let ulen = |j: usize| (len / n) + usize::from((len % n) > j);
    let mut clen = Vec::new();
    for _ in 0..n {
        clen.push(input.read_uint7()?);
    }
    let mut t = Vec::new();
    for (j, &c) in clen.iter().enumerate() {
        let mut sub = Input { bytes: input.read_data(c)?, pos: 0 };
        t.push(arith_decode(&mut sub, ulen(j), depth + 1)?);
    }
    let mut out = vec![0u8; len];
    for (j, t_j) in t.iter().enumerate() {
        for i in 0..ulen(j) {
            *out.get_mut(i * n + j)? = *t_j.get(i)?;
        }
    }
    Some(out)
}

fn decode_ext(input: &mut Input<'_>, len: usize) -> Option<Vec<u8>> {
    let data = input.rest();
    if !data.starts_with(b"BZh") {
        return None;
    }
    let mut out = Vec::new();
    let limit = u64::try_from(len).ok()?.checked_add(1)?;
    bzip2::read::BzDecoder::new(data).take(limit).read_to_end(&mut out).ok()?;
    (out.len() == len).then_some(out)
}

/// `DecodePackMeta`: the symbol map and the packed length.
fn decode_pack_meta(input: &mut Input<'_>) -> Option<(Vec<u8>, usize, usize)> {
    let nsym = usize::from(input.read_uint8()?);
    let p = input.read_data(nsym)?.to_vec();
    let len = input.read_uint7()?;
    Some((p, nsym, len))
}

fn decode_pack(data: &[u8], p: &[u8], nsym: usize, len: usize) -> Option<Vec<u8>> {
    let (per_byte, mask) = match nsym {
        0..=1 => return Some(vec![*p.first()?; len]),
        2 => (8, 1),
        3..=4 => (4, 3),
        5..=16 => (2, 15),
        _ => return None,
    };
    let bits = 8 / per_byte;
    let mut out = Vec::new();
    let mut j = 0usize;
    let mut v = 0u8;
    for i in 0..len {
        if i % per_byte == 0 {
            v = *data.get(j)?;
            j += 1;
        }
        out.push(*p.get(usize::from(v & mask))?);
        v >>= bits;
    }
    Some(out)
}

fn arith_decode(input: &mut Input<'_>, mut len: usize, depth: u32) -> Option<Vec<u8>> {
    let flags = input.read_uint8()?;
    // Correction: the pseudocode reads the length when NoSize *is* set. And
    // a stripe stores it whatever NoSize says (htscodecs).
    if flags & NO_SIZE == 0 || flags & STRIPE != 0 {
        len = bounded(input.read_uint7()?)?;
    }
    if flags & STRIPE != 0 {
        return decode_stripe(input, len, depth);
    }
    let mut pack = None;
    if flags & PACK != 0 {
        let pack_len = len;
        let (p, nsym, packed_len) = decode_pack_meta(input)?;
        len = bounded(packed_len)?;
        pack = Some((p, nsym, pack_len));
    }
    // Correction: order-1 only when the ORDER field is exactly 1 (htscodecs
    // reads the reserved 2 and 3 as order-0).
    let order1 = flags & 0x03 == ORDER;
    let data = if flags & CAT != 0 {
        input.read_data(len)?.to_vec()
    } else if flags & EXT != 0 {
        decode_ext(input, len)?
    } else if flags & RLE != 0 {
        if order1 { decode_rle::<true>(input, len)? } else { decode_rle::<false>(input, len)? }
    } else if order1 {
        decode_order1(input, len)?
    } else {
        decode_order0(input, len)?
    };
    match pack {
        Some((p, nsym, pack_len)) => decode_pack(&data, &p, nsym, pack_len),
        None => Some(data),
    }
}

#[cfg(test)]
#[allow(clippy::indexing_slicing, reason = "test code")]
mod tests {
    use std::path::Path;

    use super::*;

    const CODECS: &str =
        concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs");

    /// The vectors' README: the u32 file as-is, the quality files' first
    /// whitespace-separated column of each line with newlines dropped.
    fn expected(stem: &str) -> Vec<u8> {
        let path = Path::new(CODECS).join("gzip").join(format!("{stem}.gz"));
        let mut original = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_end(&mut original)
            .unwrap();
        if stem == "u32" {
            return original;
        }
        original
            .split(|&b| b == b'\n')
            .filter_map(|line| line.split(u8::is_ascii_whitespace).find(|f| !f.is_empty()))
            .flatten()
            .copied()
            .collect()
    }

    /// Pins the reference itself to real data, independently of the
    /// production decoder it is compared with.
    // r[verify cram.codec.arith+2]
    #[test]
    fn decodes_the_hts_specs_range_vectors() {
        let mut count = 0;
        for entry in std::fs::read_dir(Path::new(CODECS).join("range")).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_owned();
            let (stem, _) = name.rsplit_once('.').unwrap();
            let src = std::fs::read(&path).unwrap();
            let decoded = decode(&src, 0).unwrap_or_else(|| panic!("{name}: rejected"));
            assert!(decoded == expected(stem), "{name}: content differs");
            count += 1;
        }
        assert!(count >= 32, "expected the 32 hts-specs range vectors, found {count}");
    }
}
