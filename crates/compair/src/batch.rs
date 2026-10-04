//! One alignment per lane: a batch of haplotypes scored against one read in
//! lockstep, row by row.
//!
//! The two banded kernels put eight *cells of one alignment* in a vector, so
//! every lane of every step belongs to the same pair -- and pays for it: the
//! neighbours are lane shifts, the band's edges cross a vector diagonally and
//! have to be masked lane by lane, and ~23 % of the cells a strip computes are
//! outside the band. This kernel puts eight *alignments* in a vector instead.
//! Lane `k` holds candidate haplotype `k` at the same `(row, column)` as every
//! other lane, so:
//!
//! - the three neighbours are the same lane of the previous column and of the
//!   row buffer -- no shifts at all;
//! - the band is a property of the row, not of the lane, so the column loop
//!   *is* the band and there is no per-lane edge mask. What stays is one
//!   lanewise mask for haplotypes shorter than the longest in the batch;
//! - the read's row is the same for every lane, so its base, its two emission
//!   terms and its five transitions are splats held for the whole row;
//! - the haplotype columns are blocks of lanes, so lane `k`'s column `j` is
//!   lane `k` of block `j` and a column is five contiguous loads, as in the
//!   strip kernel.
//!
//! Below its plan it is the pairs kernel's sweep (`lanes`), with no column
//! shift and the read's row as splats: a batch is a group of pairs that share
//! a read and a band.
//!
//! The arithmetic, the order of the operations, the flush to zero, the
//! power-of-two renormalisation every eight rows and the free start
//! are the strip kernel's, cell for cell. So a batch is bit-identical to
//! running [`align_strips`] once per haplotype, which is the gate this kernel
//! is worth nothing without.
//!
//! What it does *not* do is make a lane free. A batch of `n < BATCH`
//! haplotypes computes the full eight lanes and throws `8 - n` away, so the
//! per-alignment cost falls only with the fill.
//!
//! [`align_strips`]: crate::align_strips

use std::borrow::Borrow;

use fearless_simd::Level;
use seqair_types::Strand;

use crate::{
    banded::{Band, LANE_MAX, Lane, RowTracks, Window, Workspace, code},
    emission::Emission,
    haplotype::Haplotype,
    lanes::{Cells, ColumnBlock, LanesView, SharedRows, lanes_kernel},
    read::Read,
    reference::trusted,
    types::Log10Likelihood,
};

/// Alignments per batch: the lane count, which is what the whole kernel is.
pub const BATCH: usize = LANE_MAX;

/// How many haplotypes a batch needs before it beats scoring them one at a
/// time through the strip kernel.
///
/// A batch computes all [`BATCH`] lanes whatever the caller asked for, so its
/// cost is nearly flat in the group size and the per-alignment cost falls only
/// with the fill. Measured against the strip kernel *on the same lane*, over a
/// 200 bp haplotype through the default band, at both read lengths (the two
/// agreed to within 0.02x at every point):
///
/// | haplotypes | 1 | 2 | 3 | 4 | 5 | **6** | 7 | 8 |
/// |---|---|---|---|---|---|---|---|---|
/// | M4 Pro (NEON) | 0.18x | 0.36x | 0.53x | 0.70x | 0.86x | **1.03x** | 1.18x | 1.33x |
/// | 3950X (AVX2) | 0.28x | 0.54x | 0.79x | **1.02x** | 1.25x | 1.48x | 1.69x | 1.90x |
///
/// **The two machines disagree -- the M4 crosses at six, the 3950X at four --
/// and six is the number that is safe on both.** Four would put the M4 into the
/// batch kernel at half fill, where its strip kernel is 1.4x faster; six leaves
/// the 3950X's half-full groups about 1.15x on the table. The number is the
/// ratio between the two kernels' per-cell costs and nothing else, so a
/// per-target constant is defensible -- but not on two machines.
///
/// `examples/batchfill.rs` regenerates the table, racing the two kernels on
/// the same lane, so the crossover it reports is the cost of an empty lane and
/// not the gap between two lanes. (The table was measured on the hand-written
/// intrinsics lane that `fearless_simd` replaced, which measured at parity
/// with it on both kernels.)
///
/// Whatever a machine says, the number is bounded below: the batch kernel
/// cannot be worth running below `BATCH / 2` lanes, because at that fill it is
/// doing more than twice the work for the same answers.
pub const BATCH_BREAK_EVEN: usize = 6;

// A break-even above `BATCH` would mean the batch kernel never runs, and one at
// or below `BATCH / 2` would mean it runs where it is doing more than twice the
// work for the same answers.
const _: () = assert!(
    BATCH_BREAK_EVEN > BATCH / 2 && BATCH_BREAK_EVEN <= BATCH,
    "the break-even has to be reachable and above half fill"
);

/// One batch of haplotypes scored, plus how many lanes that kernel fills --
/// which is how many haplotypes a group may hold, so it travels with it.
///
/// A function pointer rather than a type parameter because the eight-lane
/// kernel is not one type: it is an instance per SIMD level, chosen at run
/// time by `dispatch!` inside the function this points to. One indirect call
/// per batch is nothing, and every caller passes a constant.
#[derive(Clone, Copy)]
pub(crate) struct BatchKernel {
    lanes: usize,
    run: fn(BatchView<'_>, &mut BatchBuffer, usize, Band) -> [Log10Likelihood; BATCH],
}

impl BatchKernel {
    /// [`batch_kernel`] over one lane type, for the paths that can name it.
    const fn over<L: Lane<Token = ()>>() -> Self {
        Self {
            lanes: if L::LANES < BATCH { L::LANES } else { BATCH },
            run: |plan, buffer, read_len, band| batch_kernel::<L>((), plan, buffer, read_len, band),
        }
    }
}

/// The batch kernel at the best SIMD level this CPU has.
const SIMD_BATCH_KERNEL: BatchKernel = BatchKernel { lanes: BATCH, run: batch_kernel_simd };

/// [`SIMD_BATCH_KERNEL`]'s function, in [`BatchKernel`]'s shape.
fn batch_kernel_simd(
    plan: BatchView<'_>,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    crate::simd::batch_kernel_at(Level::new(), plan, buffer, read_len, band)
}

/// The three matrices of one read row, a block of lanes per column: column
/// `j` is `cells[j]`, one-based, and `j = width + 1` is the pad the last row's
/// `up` reads.
#[derive(Debug, Default)]
pub(crate) struct BatchBuffer {
    cells: Vec<Cells>,
}

/// The half of what one batch hoists out of its loops that depends on the
/// haplotypes, the emission and the read's strand, and on nothing else about
/// the read. The other half is the strip kernel's own `RowTracks`: every
/// lane reads the same row entry, so the batch kernel splats the same
/// numbers the strip kernel loads, and one fill per read serves both.
#[derive(Debug, Default)]
pub(crate) struct BatchPlan {
    /// Per haplotype column, a block of lanes: lane `k`'s column `j`,
    /// one-based, is lane `k` of `columns[j]`, so `j = 0` is the free-start
    /// column and `j = width + 1` a pad.
    columns: Vec<ColumnBlock>,
    /// `1 / h` per lane, the read's free start; zero for an absent or empty
    /// haplotype, which then scores [`Log10Likelihood::IMPOSSIBLE`].
    init: Window,
    /// `h + 1` per lane, so `column < lengths[k]` is `column <= h`.
    past_end: Window,
    /// The longest haplotype in the batch.
    width: usize,
}

/// A batch plan and a read's rows, borrowed the way the kernel reads them.
#[derive(Clone, Copy)]
pub(crate) struct BatchView<'a> {
    rows: &'a RowTracks,
    columns: &'a [ColumnBlock],
    init: &'a Window,
    past_end: &'a Window,
    width: usize,
}

impl BatchPlan {
    /// Rebuilds the column half for this group of haplotypes on `strand`,
    /// deriving lane `k`'s columns `columns(h_k)` and leaving the rest at
    /// their sentinels. `None` only if a haplotype cannot be addressed, which
    /// its constructor rules out.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        reason = "the f32 narrowing is the point of this kernel, and a haplotype is a few hundred bases"
    )]
    pub(crate) fn fill<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        emission: &E,
        strand: Strand,
        columns: impl Fn(usize) -> Option<(usize, usize)>,
    ) -> Option<()> {
        self.width = haplotypes.iter().map(|h| h.borrow().len()).max().unwrap_or(0);
        let Self { columns: blocks, init, past_end, width } = self;
        blocks.clear();
        blocks.resize(*width + 2, ColumnBlock::UNVISITED);
        *init = [0.0; BATCH];
        *past_end = [0.0; BATCH];

        for (lane, haplotype) in haplotypes.iter().enumerate().take(BATCH) {
            let haplotype = haplotype.borrow();
            let h = haplotype.len();
            if h == 0 {
                continue;
            }
            *init.get_mut(lane)? = 1.0 / h as f32;
            *past_end.get_mut(lane)? = (h + 1) as f32;
            // Only the columns asked for; for one read, the columns the band
            // can reach. Everything else keeps the sentinels just written,
            // which is what an unvisited column holds anyway -- the kernel's
            // column loop is the band, so it never reads them. Deriving them
            // cost ~29% of a call at 200 bp haplotypes and grew without bound
            // with haplotype length, because the band's span is `read +
            // width` however long the haplotype is.
            let Some((lo, hi)) = columns(h) else { continue };
            for (index, block) in (lo..=hi).zip(blocks.get_mut(lo + 1..=hi + 1)?) {
                let weights = emission.site_weights(haplotype.site(index)?, strand);
                let site_base = code(weights.base);
                *block.base.get_mut(lane)? = site_base;
                match weights.converted {
                    Some(converted_base) => {
                        *block.converted.get_mut(lane)? = code(converted_base);
                        *block.rate.get_mut(lane)? = weights.rate as f32;
                        *block.unconverted.get_mut(lane)? = (1.0 - weights.rate) as f32;
                    }
                    None => *block.plain.get_mut(lane)? = site_base,
                }
            }
        }
        Some(())
    }

    /// This plan with a read's rows, as the kernel reads them.
    pub(crate) fn view<'a>(&'a self, rows: &'a RowTracks) -> BatchView<'a> {
        BatchView {
            rows,
            columns: &self.columns,
            init: &self.init,
            past_end: &self.past_end,
            width: self.width,
        }
    }
}

impl Workspace {
    /// One read against its candidate haplotypes, through whichever kernel is
    /// faster for the number of them. **This is the entry point a caller
    /// wants**; [`Workspace::align_batch`] and [`Workspace::align_strips_simd`]
    /// are the two kernels it chooses between.
    ///
    /// `haplotypes` may hold the haplotypes or references to them. The scores
    /// replace `out`'s contents, in the haplotypes' order.
    ///
    /// The choice is per batch, not per call, and it has to be: a caller with
    /// nine haplotypes gets one full batch and a remainder of one, and running
    /// the remainder through eight lanes to score a single haplotype is what
    /// made nine haplotypes *slower* than eight before this existed. Full
    /// batches go through the batch kernel, a remainder shorter than
    /// [`BATCH_BREAK_EVEN`] goes through the strip kernel one haplotype at a
    /// time, and because the two kernels are bit-identical the split cannot
    /// change a score.
    pub fn align_candidates<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        for group in haplotypes.chunks(BATCH) {
            if group.len() >= BATCH_BREAK_EVEN {
                self.batch_group::<H, E>(group, read, emission, band, out, SIMD_BATCH_KERNEL);
            } else {
                for haplotype in group {
                    out.push(self.align_strips_simd(haplotype.borrow(), read, emission, band));
                }
            }
        }
    }

    /// Up to [`BATCH`] haplotypes against one read, one haplotype per lane.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order. More than
    /// [`BATCH`] haplotypes are scored in successive batches; a batch short of
    /// [`BATCH`] still computes every lane, so the cost per alignment falls
    /// only with the fill. Prefer [`Workspace::align_candidates`], which is
    /// this kernel where it wins and the strip kernel where it does not.
    pub fn align_batch<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        self.batch::<H, E>(haplotypes, read, emission, band, out, SIMD_BATCH_KERNEL);
    }

    /// [`Workspace::align_batch`] at a given `fearless_simd` level, so the
    /// parity tests can hold every level the CPU has to the scalar kernel.
    #[doc(hidden)]
    pub fn align_batch_at<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        level: Level,
        haplotypes: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        let rows = !read.is_empty() && self.plan.rows.fill(read, emission).is_some();
        for group in haplotypes.chunks(BATCH) {
            let scores = if rows && self.fill_batch(group, read, emission, band).is_some() {
                crate::simd::batch_kernel_at(
                    level,
                    self.batch_plan.view(&self.plan.rows),
                    &mut self.batch_rows,
                    read.len(),
                    band,
                )
            } else {
                [Log10Likelihood::IMPOSSIBLE; BATCH]
            };
            out.extend(scores.iter().zip(group).map(|(score, haplotype)| {
                trusted(*score, haplotype.borrow(), read, emission, band)
            }));
        }
    }

    /// [`Workspace::align_batch`] with one lane, the bit-parity oracle.
    #[doc(hidden)]
    pub fn align_batch_scalar<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        self.batch::<H, E>(haplotypes, read, emission, band, out, BatchKernel::over::<f32>());
    }

    fn batch<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        kernel: BatchKernel,
    ) {
        out.clear();
        for group in haplotypes.chunks(kernel.lanes) {
            self.batch_group::<H, E>(group, read, emission, band, out, kernel);
        }
    }

    /// One batch of at most [`BATCH`] haplotypes, appended to `out`.
    fn batch_group<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        group: &[H],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        kernel: BatchKernel,
    ) {
        let scores = if !read.is_empty()
            && self.plan.rows.fill(read, emission).is_some()
            && self.fill_batch(group, read, emission, band).is_some()
        {
            (kernel.run)(
                self.batch_plan.view(&self.plan.rows),
                &mut self.batch_rows,
                read.len(),
                band,
            )
        } else {
            [Log10Likelihood::IMPOSSIBLE; BATCH]
        };
        out.extend(
            scores.iter().zip(group).map(|(score, haplotype)| {
                trusted(*score, haplotype.borrow(), read, emission, band)
            }),
        );
    }

    /// The batch plan's column half for one read's band.
    fn fill_batch<H: Borrow<Haplotype>, E: Emission>(
        &mut self,
        group: &[H],
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Option<()> {
        let r = read.len();
        self.batch_plan.fill(group, emission, read.strand(), |h| band.columns(h, r))
    }
}

/// Every lane's own haplotype against the one read whose rows `plan` holds:
/// the shared sweep with no shift (`delta_k = 0`), every lane's read the same
/// length, and the read's row as splats.
#[inline(always)]
pub(crate) fn batch_kernel<L: Lane>(
    token: L::Token,
    plan: BatchView<'_>,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    let view = LanesView {
        columns: plan.columns,
        origin: 0,
        init: plan.init,
        front: &[0.0; BATCH],
        past_end: plan.past_end,
        lengths: &[read_len; BATCH],
        rows_len: read_len,
        width: plan.width,
        offset: band.offset,
        half_width: band.half_width,
    };
    lanes_kernel::<L, _>(token, view, &SharedRows(plan.rows), &mut buffer.cells)
}

#[cfg(test)]
mod tests {
    use fearless_simd::Level;
    use hegel::TestCase;
    use hegel::generators as gs;

    use super::{BATCH, BatchBuffer, Cells};
    use crate::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};

    /// One of a cell's three matrices.
    type Matrix = fn(&Cells) -> &[f32; 8];

    /// Every cell of lane `lane` past column `h` in one matrix, as `(column,
    /// value)` for the ones that are not zero.
    fn survivors(cells: &[Cells], matrix: Matrix, lane: usize, h: usize) -> Vec<(usize, f32)> {
        cells
            .iter()
            .enumerate()
            .skip(h + 1)
            .filter_map(|(column, cell)| {
                matrix(cell).get(lane).filter(|&&value| value != 0.0).map(|&value| (column, value))
            })
            .collect()
    }

    fn lanes_past_their_ends(buffer: &BatchBuffer, lengths: &[usize]) -> Vec<String> {
        let matrices: [(&str, Matrix); 3] =
            [("m", |cell| &cell.m), ("i", |cell| &cell.i), ("d", |cell| &cell.d)];
        let mut found = Vec::new();
        for (lane, &h) in lengths.iter().enumerate() {
            for (name, matrix) in matrices {
                let cells = survivors(&buffer.cells, matrix, lane, h);
                if !cells.is_empty() {
                    found.push(format!("lane {lane} (h = {h}) {name}: {cells:?}"));
                }
            }
        }
        found
    }

    /// Past a lane's haplotype every cell the batch kernel leaves in its
    /// row buffer is zero, in all three matrices, at every level.
    ///
    /// No score shows this: a deletion leaking rightward past a lane's
    /// end never reaches that lane's total, only its running maximum,
    /// and so only which cells flush. The column loop
    /// masks only the columns past the batch's shortest haplotype, so
    /// this is what pins where that split falls -- one column late and
    /// the shortest lane's `d` survives at `h + 1`.
    #[hegel::test]
    fn nothing_survives_past_a_lanes_haplotype(tc: TestCase) {
        let read_len = tc.draw(gs::integers::<usize>().min_value(6).max_value(39));
        let cuts: [usize; 8] = tc.draw(gs::arrays(gs::integers::<usize>().max_value(9)));
        let lanes = tc.draw(gs::integers::<usize>().min_value(1).max_value(BATCH));
        let half_width = tc.draw(gs::integers::<u32>().min_value(3).max_value(11));
        let seed = tc.draw(gs::integers::<u64>());
        let base = |index: usize| {
            let bits = seed.rotate_left(u32::try_from(index % 64).unwrap_or(0)) >> 62;
            [Base::A, Base::C, Base::G, Base::T][usize::try_from(bits).unwrap_or(0)]
        };
        let full = read_len + 8;
        let haplotypes: Vec<Haplotype> = cuts
            .iter()
            .take(lanes)
            .map(|&cut| Haplotype::new((0..full - cut).map(base).collect::<Vec<_>>()))
            .collect();
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let lengths: Vec<usize> = haplotypes.iter().map(Haplotype::len).collect();
        let quals = vec![BaseQuality::from_byte(30); read_len];
        let q = BaseQuality::from_byte(40);
        let read = Read::uniform(
            (0..read_len).map(|index| base(index + 3)).collect::<Vec<_>>(),
            &quals,
            q,
            q,
            BaseQuality::from_byte(10),
            Strand::OT,
        )
        .expect("valid read");
        let band = Band::new(2 * half_width, 2).expect("valid band");

        for level in [Level::new(), Level::fallback()] {
            let mut workspace = Workspace::new();
            let mut out = Vec::new();
            workspace.align_batch_at(
                level,
                &refs,
                &read,
                &StandardEmission::default(),
                band,
                &mut out,
            );
            let found = lanes_past_their_ends(&workspace.batch_rows, &lengths);
            assert!(found.is_empty(), "{level:?}: {}", found.join("; "));
        }
    }
}
