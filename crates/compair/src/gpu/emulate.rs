//! The shader, line for line, in Rust.
//!
//! Two uses. It is the test that the *algorithm* the shader runs -- the band-
//! relative row buffer updated in place, the masks, the renormalisation --
//! is the batch kernel's, bit for bit, which a machine without a GPU can run;
//! and with `flush` set it is the shader as a GPU that flushes subnormal
//! intermediates executes it, which is how the notes attribute the scores a
//! GPU gets wrong in the last bits.

use super::{
    device::finish,
    plan::{
        GAP_CONTINUATION, GpuPairs, INDEL_TO_MATCH, MATCH_TO_DELETION, MATCH_TO_INSERTION,
        PairRecord, RowRecord, TRANSITIONS,
    },
};
use crate::types::Log10Likelihood;

/// What the GPU does to a subnormal result of an arithmetic operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub enum Subnormals {
    /// Kept, as IEEE and the CPU kernels do.
    Kept,
    /// Flushed to zero after every multiply and add, as Metal's fast-math
    /// and RADV's default float mode do.
    Flushed,
}

impl GpuPairs {
    /// Every pair scored by the Rust transcription of the shader.
    #[doc(hidden)]
    #[must_use]
    pub fn emulate(&self, subnormals: Subnormals) -> Vec<Log10Likelihood> {
        self.pairs
            .iter()
            .map(|pending| match pending.class {
                Some(_) => {
                    let (sum, exponent) =
                        run(&pending.record, &self.rows, &self.weights, subnormals);
                    finish(sum, exponent)
                }
                None => Log10Likelihood::IMPOSSIBLE,
            })
            .collect()
    }
}

#[allow(
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    reason = "the shader's own i32/u32 arithmetic, reproduced"
)]
fn run(
    pair: &PairRecord,
    rows: &[RowRecord],
    weights: &[f32],
    subnormals: Subnormals,
) -> (f32, i32) {
    let ftz = |x: f32| match subnormals {
        Subnormals::Kept => x,
        Subnormals::Flushed if x.abs() < f32::MIN_POSITIVE => 0.0,
        Subnormals::Flushed => x,
    };
    let mul = |a: f32, b: f32| ftz(a * b);
    let add = |a: f32, b: f32| ftz(a + b);
    let flush = |x: f32| if x < f32::MIN_POSITIVE { 0.0 } else { x };

    let half = pair.half_width as i32;
    let span = (2 * pair.half_width + 1) as usize;
    let mut m = vec![0.0f32; span + 1];
    let mut ins = vec![0.0f32; span + 1];
    let mut del = vec![0.0f32; span + 1];

    let mut crossing = 0.0f32;
    for p in 0..span {
        let column = pair.offset - half + p as i32;
        let live = p as i32 <= 2 * half && column >= 0 && column <= pair.hap_len as i32;
        let cell = if live { pair.init } else { 0.0 };
        if let Some(slot) = del.get_mut(p) {
            *slot = cell;
        }
        crossing = crossing.max(cell);
    }

    let mut exponent = 0i32;
    let mut total = 0.0f32;
    let origin = pair.offset - half;
    for row in 1..=pair.read_len {
        if (row - 1) % 8 == 0 {
            let shift = crate::scaling::normalising_shift_f32(crossing);
            exponent += shift;
            if shift != 0 {
                let lift = crate::scaling::exp2_f32(shift);
                for track in [&mut m, &mut ins, &mut del] {
                    for value in track.iter_mut() {
                        *value = mul(*value, lift);
                    }
                }
            }
        }
        let Some(record) = rows.get((pair.rows + row - 1) as usize) else {
            return (0.0, 0);
        };
        let (spread, mismatched) = (record.spread, record.mismatched);
        let base = record.packed & 0xff;
        let [insertion, deletion, gap] =
            [8, 16, 24].map(|shift| ((record.packed >> shift) & 0xff) as usize);
        let lookup = |index: usize| TRANSITIONS.get(index).copied().unwrap_or(0.0);
        let mm = lookup(insertion * 256 + deletion);
        let mti = lookup(MATCH_TO_INSERTION + insertion);
        let mtd = lookup(MATCH_TO_DELETION + deletion);
        let itm = lookup(INDEL_TO_MATCH + gap);
        let gc = lookup(GAP_CONTINUATION + gap);
        let unknown = base >= 4;
        let first = origin + row as i32;
        let track = (pair.weights + base.min(3) * pair.hap_len) as i32 + first - 1;
        let lo = (1 - first).max(0);
        let hi = (2 * half).min(pair.hap_len as i32 - first);
        let summing = row == pair.read_len;
        let (mut left_m, mut left_d, mut running) = (0.0f32, 0.0f32, 0.0f32);
        for p in 0..span {
            let at = |track: &[f32], index: usize| track.get(index).copied().unwrap_or(0.0);
            let (diag_m, diag_indel) = (at(&m, p), add(at(&ins, p), at(&del, p)));
            let (up_m, up_i) = (at(&m, p + 1), at(&ins, p + 1));
            let keep = p as i32 >= lo && p as i32 <= hi;
            let weight = if unknown {
                1.0
            } else if keep {
                weights.get((track + p as i32) as usize).copied().unwrap_or(0.0)
            } else {
                0.0
            };
            let prior = add(mismatched, mul(weight, spread));
            let new_m = mul(prior, add(mul(diag_m, mm), mul(diag_indel, itm)));
            let new_i = add(mul(up_m, mti), mul(up_i, gc));
            let new_d = add(mul(left_m, mtd), mul(left_d, gc));
            let (new_m, new_i, new_d) =
                if keep { (flush(new_m), flush(new_i), flush(new_d)) } else { (0.0, 0.0, 0.0) };
            for (track, value) in [(&mut m, new_m), (&mut ins, new_i), (&mut del, new_d)] {
                if let Some(slot) = track.get_mut(p) {
                    *slot = value;
                }
            }
            running = running.max(new_m).max(new_i).max(new_d);
            if summing {
                total = add(total, add(new_m, new_i));
            }
            left_m = new_m;
            left_d = new_d;
        }
        if row % 8 == 0 {
            crossing = running;
        }
    }
    (total, exponent)
}
