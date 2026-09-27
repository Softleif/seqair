//! The WGSL source of the kernel, generated per band class.
//!
//! WGSL has no private arrays of run-time size, so the row buffer is sized by
//! a `const` and there is one shader per class of half-width. The source is
//! assembled here rather than kept in a `.wgsl` file because two of its
//! variants are not a constant away from each other: the unrolled style
//! replaces the three row arrays with `3 * (2W + 2)` named scalars, which is
//! the one way to be sure a driver keeps them in registers instead of indexing
//! them in scratch memory.
//!
//! Both styles perform the batch kernel's per-lane arithmetic, in its order:
//! the recurrence of `batch::batch_kernel`, the flush of subnormal cells, the
//! power-of-two renormalisation of the row every eight rows from the maximum
//! of the row that crosses into the strip, and `total + (m + i)` on the last
//! row. The row buffer is indexed by *band position* rather than by column --
//! position `p` of row `i` is column `i + offset - half_width + p` -- so a
//! cell reads its diagonal at `p` and the cell above at `p + 1` of the
//! previous row and its left neighbour at `p - 1` of its own. Updated in place
//! in ascending `p`, `p` and `p + 1` are still the previous row's when they are
//! read and `p - 1` is already this row's.

use core::fmt::Write as _;

use super::plan::{GAP_CONTINUATION, INDEL_TO_MATCH, MATCH_TO_DELETION, MATCH_TO_INSERTION};

/// How the row buffer is held.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Style {
    /// Three private arrays indexed by the band position, and a loop over it.
    Loop,
    /// The band loop unrolled at generation time, the arrays replaced by named
    /// scalars: nothing is indexed at run time, so nothing has to live in
    /// memory that a register could hold.
    Unrolled,
}

/// Whether the compiler may fuse a multiply into the add after it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Contraction {
    /// As written, which lets Metal's fast-math and SPIR-V's default both
    /// contract `a * b + c` into a fused multiply-add.
    Allowed,
    /// Every product and the last row's pair sum pass through an `or` with a
    /// uniform zero the compiler cannot see through, so each is rounded to
    /// `f32` before anything consumes it. One integer op per product.
    Blocked,
}

/// Everything a shader is specialised on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Variant {
    /// The widest half-width this shader's row buffer holds; a pair with a
    /// narrower band masks the positions past its own.
    pub half_width: u32,
    pub style: Style,
    pub contraction: Contraction,
    pub workgroup_size: u32,
}

/// The shader's source for `variant`. Writing into a `String` only fails if a
/// `Display` impl does, and none here can.
pub fn source(variant: Variant) -> Result<String, core::fmt::Error> {
    let mut out = String::with_capacity(8192);
    let Variant { half_width, style, contraction, workgroup_size } = variant;
    out.push_str(HEADER);
    writeln!(out, "const HALF: u32 = {half_width}u;")?;
    writeln!(out, "const WORKGROUP: u32 = {workgroup_size}u;")?;
    writeln!(out, "const MATCH_TO_INSERTION: u32 = {MATCH_TO_INSERTION}u;")?;
    writeln!(out, "const MATCH_TO_DELETION: u32 = {MATCH_TO_DELETION}u;")?;
    writeln!(out, "const INDEL_TO_MATCH: u32 = {INDEL_TO_MATCH}u;")?;
    writeln!(out, "const GAP_CONTINUATION: u32 = {GAP_CONTINUATION}u;")?;
    out.push_str(
        match contraction {
            Contraction::Allowed => {
                "fn product(a: f32, b: f32) -> f32 { return a * b; }\n\
                 fn rounded(x: f32) -> f32 { return x; }\n"
            }
            Contraction::Blocked => {
                "fn product(a: f32, b: f32) -> f32 { return rounded(a * b); }\n\
                 fn rounded(x: f32) -> f32 { return bitcast<f32>(bitcast<u32>(x) | params.zero); }\n"
            }
        },
    );
    out.push_str(HELPERS);
    match style {
        Style::Loop => out.push_str(LOOP_MAIN),
        Style::Unrolled => unrolled_main(&mut out, half_width)?,
    }
    Ok(out)
}

const HEADER: &str = r"
// One (read, haplotype) pair per invocation. See `compair::gpu::shader`.

struct Pair {
    rows: u32,
    read_len: u32,
    weights: u32,
    hap_len: u32,
    offset: i32,
    half_width: u32,
    init: f32,
    pad: u32,
}

struct RowRecord {
    packed: u32,
    spread: f32,
    mismatched: f32,
}

struct Params {
    first: u32,
    count: u32,
    zero: u32,
    pad: u32,
}

@group(0) @binding(0) var<storage, read> pairs: array<Pair>;
@group(0) @binding(1) var<storage, read> rows: array<RowRecord>;
@group(0) @binding(2) var<storage, read> weights: array<f32>;
@group(0) @binding(3) var<storage, read_write> scores: array<vec2<u32>>;
@group(0) @binding(4) var<uniform> params: Params;
// `plan::TRANSITIONS`: match-to-match by (insertion, deletion) quality, then
// the four transitions of one quality each.
@group(0) @binding(5) var<storage, read> transitions: array<f32>;

const STRIP_ROWS: u32 = 8u;
";

/// The cell, the renormalisation and the per-row inputs, shared by both
/// styles so that the arithmetic is written once.
const HELPERS: &str = r"
const SPAN: u32 = 2u * HALF + 1u;
// `f32::MIN_POSITIVE`, as a hex float so no decimal parse can round it.
const TINY: f32 = 0x1p-126f;

fn flush(x: f32) -> f32 {
    return select(x, 0.0, x < TINY);
}

// `scaling::normalising_shift_f32`: the shift that puts a normal `max` into
// [1, 2), zero for anything non-positive, subnormal or not finite.
fn normalising_shift(max: f32) -> i32 {
    if !(max > 0.0) {
        return 0;
    }
    let exponent = i32((bitcast<u32>(max) >> 23u) & 0xffu) - 127;
    if exponent == -127 || exponent == 128 {
        return 0;
    }
    return -exponent;
}

// `scaling::exp2_f32` over the range `normalising_shift` produces.
fn exp2_exact(shift: i32) -> f32 {
    if shift < -126 {
        return 0.0;
    }
    return bitcast<f32>(u32(shift + 127) << 23u);
}

struct Row {
    // The first weight of this row's base's track, minus one column, as a
    // wrapping index: position `p` reads `track + p` (see `row_inputs`).
    track: u32,
    unknown: bool,
    spread: f32,
    mismatched: f32,
    match_to_match: f32,
    match_to_insertion: f32,
    match_to_deletion: f32,
    indel_to_match: f32,
    gap_continuation: f32,
    // Live positions are `lo..=hi`: inside the pair's own band and inside
    // columns `1..=hap_len`.
    lo: i32,
    hi: i32,
}

fn row_inputs(pair: Pair, row: u32, first: i32) -> Row {
    let record = rows[pair.rows + row - 1u];
    let base = record.packed & 0xffu;
    let insertion = (record.packed >> 8u) & 0xffu;
    let deletion = (record.packed >> 16u) & 0xffu;
    let gap = record.packed >> 24u;
    var out: Row;
    out.unknown = base >= 4u;
    // Column `first + p` reads site `first + p - 1`. Where that is outside the
    // haplotype the cell is masked, so the index may wrap: an out-of-bounds
    // storage read in WGSL returns some in-bounds value or zero, never
    // undefined behaviour, and the value is discarded.
    out.track = u32(i32(pair.weights + min(base, 3u) * pair.hap_len) + first - 1);
    out.spread = record.spread;
    out.mismatched = record.mismatched;
    out.match_to_match = transitions[insertion * 256u + deletion];
    out.match_to_insertion = transitions[MATCH_TO_INSERTION + insertion];
    out.match_to_deletion = transitions[MATCH_TO_DELETION + deletion];
    out.indel_to_match = transitions[INDEL_TO_MATCH + gap];
    out.gap_continuation = transitions[GAP_CONTINUATION + gap];
    out.lo = max(0, 1 - first);
    out.hi = min(i32(2u * pair.half_width), i32(pair.hap_len) - first);
    return out;
}

struct Cell {
    m: f32,
    i: f32,
    d: f32,
}

// One cell of `batch::batch_kernel`, in its order of operations.
fn cell(
    r: Row,
    p: i32,
    diag_m: f32,
    diag_indel: f32,
    up_m: f32,
    up_i: f32,
    left_m: f32,
    left_d: f32,
) -> Cell {
    let weight = select(weights[r.track + u32(p)], 1.0, r.unknown);
    let prior = r.mismatched + product(weight, r.spread);
    let m = prior * (product(diag_m, r.match_to_match) + product(diag_indel, r.indel_to_match));
    let i = product(up_m, r.match_to_insertion) + product(up_i, r.gap_continuation);
    let d = product(left_m, r.match_to_deletion) + product(left_d, r.gap_continuation);
    let keep = p >= r.lo && p <= r.hi;
    var out: Cell;
    out.m = select(0.0, flush(m), keep);
    out.i = select(0.0, flush(i), keep);
    out.d = select(0.0, flush(d), keep);
    return out;
}

fn row0_cell(pair: Pair, p: i32) -> f32 {
    let column = pair.offset - i32(pair.half_width) + p;
    let live = p <= i32(2u * pair.half_width) && column >= 0 && column <= i32(pair.hap_len);
    return select(0.0, pair.init, live);
}
";

const LOOP_MAIN: &str = r"
@compute @workgroup_size(WORKGROUP)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    if gid.x >= params.count {
        return;
    }
    let index = params.first + gid.x;
    let pair = pairs[index];

    // One slot past the band, always zero: the `up` of the last position.
    var m: array<f32, SPAN + 1u>;
    var ins: array<f32, SPAN + 1u>;
    var del: array<f32, SPAN + 1u>;

    var crossing = 0.0;
    for (var p = 0u; p < SPAN; p++) {
        let cell = row0_cell(pair, i32(p));
        del[p] = cell;
        crossing = max(crossing, cell);
    }

    var exponent = 0i;
    var total = 0.0;
    let origin = pair.offset - i32(pair.half_width);
    for (var row = 1u; row <= pair.read_len; row++) {
        if (row - 1u) % STRIP_ROWS == 0u {
            let shift = normalising_shift(crossing);
            exponent += shift;
            if shift != 0 {
                let lift = exp2_exact(shift);
                for (var p = 0u; p < SPAN; p++) {
                    m[p] *= lift;
                    ins[p] *= lift;
                    del[p] *= lift;
                }
            }
        }
        let r = row_inputs(pair, row, origin + i32(row));
        let summing = row == pair.read_len;
        var left_m = 0.0;
        var left_d = 0.0;
        var running = 0.0;
        for (var p = 0u; p < SPAN; p++) {
            let c = cell(r, i32(p), m[p], ins[p] + del[p], m[p + 1u], ins[p + 1u], left_m, left_d);
            m[p] = c.m;
            ins[p] = c.i;
            del[p] = c.d;
            running = max(max(max(running, c.m), c.i), c.d);
            if summing {
                total = total + rounded(c.m + c.i);
            }
            left_m = c.m;
            left_d = c.d;
        }
        if row % STRIP_ROWS == 0u {
            crossing = running;
        }
    }
    scores[index] = vec2<u32>(bitcast<u32>(total), bitcast<u32>(exponent));
}
";

/// `LOOP_MAIN` with every band loop unrolled and `m[p]` spelled `m_p`.
fn unrolled_main(out: &mut String, half_width: u32) -> core::fmt::Result {
    let span = 2 * half_width + 1;
    out.push_str(
        r"
@compute @workgroup_size(WORKGROUP)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    if gid.x >= params.count {
        return;
    }
    let index = params.first + gid.x;
    let pair = pairs[index];
",
    );
    for p in 0..=span {
        writeln!(out, "    var m_{p} = 0.0; var i_{p} = 0.0; var d_{p} = 0.0;")?;
    }
    out.push_str("    var crossing = 0.0;\n");
    for p in 0..span {
        writeln!(out, "    d_{p} = row0_cell(pair, {p}); crossing = max(crossing, d_{p});")?;
    }
    out.push_str(
        r"
    var exponent = 0i;
    var total = 0.0;
    let origin = pair.offset - i32(pair.half_width);
    for (var row = 1u; row <= pair.read_len; row++) {
        if (row - 1u) % STRIP_ROWS == 0u {
            let shift = normalising_shift(crossing);
            exponent += shift;
            if shift != 0 {
                let lift = exp2_exact(shift);
",
    );
    for p in 0..span {
        writeln!(out, "                m_{p} *= lift; i_{p} *= lift; d_{p} *= lift;")?;
    }
    out.push_str(
        r"            }
        }
        let r = row_inputs(pair, row, origin + i32(row));
        let summing = row == pair.read_len;
        var left_m = 0.0;
        var left_d = 0.0;
        var running = 0.0;
        var c: Cell;
",
    );
    for p in 0..span {
        let up = p + 1;
        writeln!(
            out,
            "        c = cell(r, {p}, m_{p}, i_{p} + d_{p}, m_{up}, i_{up}, left_m, left_d);\n        \
             m_{p} = c.m; i_{p} = c.i; d_{p} = c.d;\n        \
             running = max(max(max(running, c.m), c.i), c.d);\n        \
             if summing {{ total = total + rounded(c.m + c.i); }}\n        \
             left_m = c.m; left_d = c.d;"
        )?;
    }
    out.push_str(
        r"        if row % STRIP_ROWS == 0u {
            crossing = running;
        }
    }
    scores[index] = vec2<u32>(bitcast<u32>(total), bitcast<u32>(exponent));
}
",
    );
    Ok(())
}
