//! Candidate haplotypes prepared once and scored against many reads.
//!
//! A plan has two halves with different lifetimes. The row tracks are a
//! function of the read and the emission; the column tracks are a function of
//! the haplotype, the emission and the read's *strand*, and of nothing else
//! about the read. A caller scoring R reads against H haplotypes -- a variant
//! caller's shape, and rastair's -- rebuilt both halves R * H times through
//! [`Workspace::align_strips_simd`], where R row fills and 2 * H column fills
//! are all there is. [`Candidates`] does that many.
//!
//! The column tracks are derived over the **whole** haplotype, not the band's
//! columns, so one fill serves every read whatever its band. The kernels read
//! the band and nothing else, so the extra columns cannot change a score, and
//! the tests hold [`Candidates`] to the per-pair path bit for bit.

use std::borrow::Borrow;

use fearless_simd::Level;
use seqair_types::Strand;

use crate::{
    banded::{Band, ColumnTracks, PlanView, Shape, Workspace},
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    strips::{RowBuffer, strip_kernel},
    types::Log10Likelihood,
};

/// One haplotype's column tracks, one set per strand, each derived the first
/// time a read on that strand asks for it.
#[derive(Debug, Default)]
pub(crate) struct PreparedColumns {
    /// `[OT, OB]`.
    strands: [ColumnTracks; 2],
    ready: [bool; 2],
}

impl PreparedColumns {
    /// This haplotype's tracks for a read on `strand`, derived if they are
    /// not yet. `None` for an unknown strand, which `Read` rejects, or a
    /// haplotype the plan cannot address.
    fn for_strand<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        emission: &E,
        strand: Strand,
    ) -> Option<&ColumnTracks> {
        let slot = match strand {
            Strand::OT => 0,
            Strand::OB => 1,
            Strand::Unknown => return None,
        };
        let ready = self.ready.get_mut(slot)?;
        let tracks = self.strands.get_mut(slot)?;
        if !*ready {
            let last = haplotype.len().checked_sub(1)?;
            tracks.fill(haplotype, emission, strand, (0, last))?;
            *ready = true;
        }
        Some(tracks)
    }
}

impl Workspace {
    /// These haplotypes under this emission, ready to score read after read
    /// against.
    ///
    /// Every score is bit-identical to [`Workspace::align_strips_simd`] on
    /// the same pair; what changes is the setup. A read's row tracks are
    /// derived once for all the haplotypes, and a haplotype's column tracks
    /// once per strand for all the reads, rather than both once per pair.
    ///
    /// The haplotypes and the emission are borrowed for as long as the
    /// [`Candidates`] lives, so nothing it has derived can go stale under it,
    /// and every call here starts from nothing: the workspace keeps the
    /// allocations between calls, never the contents.
    ///
    /// `emission` is taken by value; pass `&emission` to borrow one, as
    /// `Emission` is implemented for references.
    pub fn candidates<'w, 'h, H: Borrow<Haplotype>, E: Emission>(
        &'w mut self,
        haplotypes: &'h [H],
        emission: E,
    ) -> Candidates<'w, 'h, H, E> {
        if self.prepared.len() < haplotypes.len() {
            self.prepared.resize_with(haplotypes.len(), PreparedColumns::default);
        }
        for prepared in &mut self.prepared {
            prepared.ready = [false; 2];
        }
        Candidates { workspace: self, haplotypes, emission }
    }
}

/// Candidate haplotypes under one emission, scored against one read at a
/// time; see [`Workspace::candidates`].
pub struct Candidates<'w, 'h, H, E> {
    workspace: &'w mut Workspace,
    haplotypes: &'h [H],
    emission: E,
}

impl<H, E> core::fmt::Debug for Candidates<'_, '_, H, E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Candidates").field("haplotypes", &self.haplotypes.len()).finish()
    }
}

impl<H: Borrow<Haplotype>, E: Emission> Candidates<'_, '_, H, E> {
    #[must_use]
    pub fn len(&self) -> usize {
        self.haplotypes.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.haplotypes.is_empty()
    }

    /// `read` against every haplotype through the eight-lane strip kernel, at
    /// the best SIMD level this CPU has.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order.
    pub fn align_strips_simd(&mut self, read: &Read, band: Band, out: &mut Vec<Log10Likelihood>) {
        self.align_strips_simd_at(Level::new(), read, band, out);
    }

    /// [`Candidates::align_strips_simd`] at a given `fearless_simd` level, so
    /// the parity tests can hold every level the CPU has to the scalar kernel.
    #[doc(hidden)]
    pub fn align_strips_simd_at(
        &mut self,
        level: Level,
        read: &Read,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        self.each(read, band, out, |plan, rows, shape, band| {
            crate::simd::strip_kernel_at(level, plan, rows, shape, band)
        });
    }

    /// `read` against every haplotype through the scalar strip kernel.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order.
    pub fn align_strips(&mut self, read: &Read, band: Band, out: &mut Vec<Log10Likelihood>) {
        self.each(read, band, out, |plan, rows, shape, band| {
            strip_kernel::<f32>((), plan, rows, shape, band)
        });
    }

    /// Fills the read's rows once, then runs `kernel` per haplotype on its
    /// prepared columns. Where [`Workspace::fill_plan`] would have returned
    /// `None` -- an empty read or haplotype, a band that misses the
    /// haplotype -- the score is [`Log10Likelihood::IMPOSSIBLE`], as it is
    /// there.
    fn each(
        &mut self,
        read: &Read,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        mut kernel: impl FnMut(PlanView<'_>, &mut RowBuffer, Shape, Band) -> Log10Likelihood,
    ) {
        out.clear();
        let Workspace { plan, rows, prepared, .. } = &mut *self.workspace;
        let r = read.len();
        let rows_ready = r > 0 && plan.rows.fill(read, &self.emission).is_some();
        for (haplotype, columns) in self.haplotypes.iter().zip(prepared.iter_mut()) {
            let haplotype = haplotype.borrow();
            let h = haplotype.len();
            let score = if rows_ready && band.columns(h, r).is_some() {
                columns.for_strand(haplotype, &self.emission, read.strand()).map(|columns| {
                    let view = PlanView { rows: &plan.rows, columns };
                    kernel(view, rows, Shape { haplotype: h, read: r }, band)
                })
            } else {
                None
            };
            out.push(score.unwrap_or(Log10Likelihood::IMPOSSIBLE));
        }
    }
}
