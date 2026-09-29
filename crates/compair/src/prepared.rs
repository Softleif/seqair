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
//! The same holds for the batch kernel's interleaved columns, which depend on
//! a group of eight haplotypes, the emission and the strand:
//! [`Candidates::align`] derives them once per group and strand, and the
//! batch kernel reads the same row tracks the strip kernel does.
//!
//! The column tracks are derived over the **whole** haplotype, not the band's
//! columns, so one fill serves every read whatever its band. The kernels read
//! the band and nothing else, so the extra columns cannot change a score, and
//! the tests hold [`Candidates`] to the per-pair path bit for bit.

use std::borrow::Borrow;

use fearless_simd::Level;
use seqair_types::Strand;

use crate::{
    banded::{Band, ColumnTracks, PlanView, RowTracks, Shape, Workspace},
    batch::{BATCH, BATCH_BREAK_EVEN},
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    reference::trusted,
    strips::{RowBuffer, strip_kernel},
    types::Log10Likelihood,
};

/// Something derived per strand, `[OT, OB]`, each the first time a read on
/// that strand asks for it.
#[derive(Debug, Default)]
pub(crate) struct Prepared<T> {
    strands: [T; 2],
    ready: [bool; 2],
}

impl<T> Prepared<T> {
    /// The one for `strand`, derived by `fill` if it is not yet. `None` for
    /// an unknown strand, which `Read` rejects, or where `fill` fails.
    fn for_strand(
        &mut self,
        strand: Strand,
        fill: impl FnOnce(&mut T) -> Option<()>,
    ) -> Option<&T> {
        let slot = match strand {
            Strand::OT => 0,
            Strand::OB => 1,
            Strand::Unknown => return None,
        };
        let ready = self.ready.get_mut(slot)?;
        let prepared = self.strands.get_mut(slot)?;
        if !*ready {
            fill(prepared)?;
            *ready = true;
        }
        Some(prepared)
    }

    /// Everything derived so far is stale; the allocations are kept.
    fn forget(&mut self) {
        self.ready = [false; 2];
    }
}

/// Keeps `prepared` at least `len` long and forgets all of it.
fn forget_all<T: Default>(prepared: &mut Vec<Prepared<T>>, len: usize) {
    if prepared.len() < len {
        prepared.resize_with(len, Prepared::default);
    }
    prepared.iter_mut().for_each(Prepared::forget);
}

/// Every column of a haplotype of `h` bases, as `Band::columns` spells a
/// range: `None` for an empty one.
fn whole(h: usize) -> Option<(usize, usize)> {
    Some((0, h.checked_sub(1)?))
}

impl Workspace {
    /// These haplotypes under this emission, ready to score read after read
    /// against.
    ///
    /// Every score is bit-identical to the per-pair entry points on the same
    /// pair; what changes is the setup. A read's row tracks are derived once
    /// for all the haplotypes, and the haplotypes' column tracks once per
    /// strand for all the reads, rather than both once per pair.
    ///
    /// The haplotypes and the emission are borrowed for as long as the
    /// [`Candidates`] lives, so nothing it has derived can go stale under it,
    /// and every call here starts from nothing: the workspace keeps the
    /// allocations between calls, never the contents.
    ///
    /// `emission` is taken by value; pass `&emission` to borrow one, as
    /// `Emission` is implemented for references.
    ///
    /// ```
    /// use compair::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};
    ///
    /// // One haplotype per length of a `CA` repeat.
    /// let haplotypes: Vec<Haplotype> = (3..=9)
    ///     .map(|units| format!("CCGTAATGCC{}AGAGTTTTTC", "CA".repeat(units)))
    ///     .map(|sequence| Haplotype::from_ascii(sequence.as_bytes()))
    ///     .collect();
    /// let mut workspace = Workspace::new();
    /// let mut candidates = workspace.candidates(&haplotypes, StandardEmission::default());
    ///
    /// let (mut scores, mut calls) = (Vec::new(), Vec::new());
    /// // Reads of 5 and 6 units, both starting at haplotype position 4.
    /// for seq in [&b"AATGCCCACACACACAAGAGTT"[..], b"AATGCCCACACACACACAAGAGTT"] {
    ///     let read = Read::uniform(
    ///         Base::from_ascii_vec(seq.to_vec()),
    ///         &vec![BaseQuality::from_byte(30); seq.len()],
    ///         BaseQuality::from_byte(45),
    ///         BaseQuality::from_byte(45),
    ///         BaseQuality::from_byte(10),
    ///         Strand::OT,
    ///     )?;
    ///     // One score per haplotype, in their order.
    ///     candidates.align(&read, Band::anchored(4), &mut scores);
    ///     let best = scores.iter().enumerate().max_by(|a, b| a.1.get().total_cmp(&b.1.get()));
    ///     calls.extend(best.map(|(index, _)| index + 3));
    /// }
    /// assert_eq!(calls, [5, 6]);
    /// # Ok::<(), compair::Error>(())
    /// ```
    pub fn candidates<'w, 'h, H: Borrow<Haplotype>, E: Emission>(
        &'w mut self,
        haplotypes: &'h [H],
        emission: E,
    ) -> Candidates<'w, 'h, H, E> {
        forget_all(&mut self.prepared, haplotypes.len());
        forget_all(&mut self.prepared_batches, haplotypes.len().div_ceil(BATCH));
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

/// The pieces of a strip alignment on prepared columns.
struct StripParts<'a, E> {
    rows: &'a RowTracks,
    buffer: &'a mut RowBuffer,
    emission: &'a E,
    read: &'a Read,
    band: Band,
}

impl<E: Emission> StripParts<'_, E> {
    /// One haplotype through `kernel` on its prepared columns. Where the
    /// per-pair path's plan would have failed -- an empty haplotype, a band
    /// that misses it -- the score is [`Log10Likelihood::IMPOSSIBLE`], as it
    /// is there.
    fn score(
        &mut self,
        haplotype: &Haplotype,
        prepared: &mut Prepared<ColumnTracks>,
        kernel: &mut impl FnMut(PlanView<'_>, &mut RowBuffer, Shape, Band) -> Log10Likelihood,
    ) -> Log10Likelihood {
        let (h, r) = (haplotype.len(), self.read.len());
        let emission = self.emission;
        let strand = self.read.strand();
        let score = self.band.columns(h, r).and_then(|_| {
            let columns = prepared
                .for_strand(strand, |tracks| tracks.fill(haplotype, emission, strand, whole(h)?))?;
            let view = PlanView { rows: self.rows, columns };
            let score = kernel(view, self.buffer, Shape { haplotype: h, read: r }, self.band);
            Some(trusted(score, haplotype, self.read, emission, self.band))
        });
        score.unwrap_or(Log10Likelihood::IMPOSSIBLE)
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

    /// `read` against every haplotype, through whichever kernel is faster for
    /// the number of them: [`Workspace::align_candidates`]'s dispatch, on
    /// prepared columns. **This is the entry point a caller wants.**
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order.
    pub fn align(&mut self, read: &Read, band: Band, out: &mut Vec<Log10Likelihood>) {
        self.align_at(Level::new(), read, band, out);
    }

    /// [`Candidates::align`] at a given `fearless_simd` level.
    #[doc(hidden)]
    pub fn align_at(
        &mut self,
        level: Level,
        read: &Read,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        let Self { workspace, haplotypes, emission } = self;
        let Workspace { plan, rows, batch_rows, prepared, prepared_batches, .. } = &mut **workspace;
        let r = read.len();
        if r == 0 || plan.rows.fill(read, &*emission).is_none() {
            out.resize(haplotypes.len(), Log10Likelihood::IMPOSSIBLE);
            return;
        }
        let mut kernel = |view: PlanView<'_>, buffer: &mut RowBuffer, shape: Shape, band: Band| {
            crate::simd::strip_kernel_at(level, view, buffer, shape, band)
        };
        let mut strips =
            StripParts { rows: &plan.rows, buffer: rows, emission: &*emission, read, band };
        let strand = read.strand();
        for ((group, batch), columns) in haplotypes
            .chunks(BATCH)
            .zip(prepared_batches.iter_mut())
            .zip(prepared.chunks_mut(BATCH))
        {
            if group.len() >= BATCH_BREAK_EVEN {
                let scores = batch
                    .for_strand(strand, |plan| plan.fill(group, &*emission, strand, whole))
                    .map(|batch| {
                        let view = batch.view(strips.rows);
                        crate::simd::batch_kernel_at(level, view, batch_rows, r, band)
                    })
                    .unwrap_or([Log10Likelihood::IMPOSSIBLE; BATCH]);
                out.extend(scores.iter().zip(group.iter()).map(|(score, haplotype)| {
                    trusted(*score, haplotype.borrow(), read, &*emission, band)
                }));
            } else {
                for (haplotype, prepared) in group.iter().zip(columns.iter_mut()) {
                    out.push(strips.score(haplotype.borrow(), prepared, &mut kernel));
                }
            }
        }
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
        self.strips(read, band, out, |view, buffer, shape, band| {
            crate::simd::strip_kernel_at(level, view, buffer, shape, band)
        });
    }

    /// `read` against every haplotype through the scalar strip kernel.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order.
    pub fn align_strips(&mut self, read: &Read, band: Band, out: &mut Vec<Log10Likelihood>) {
        self.strips(read, band, out, |view, buffer, shape, band| {
            strip_kernel::<f32>((), view, buffer, shape, band)
        });
    }

    /// Fills the read's rows once, then runs `kernel` per haplotype on its
    /// prepared columns.
    fn strips(
        &mut self,
        read: &Read,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        mut kernel: impl FnMut(PlanView<'_>, &mut RowBuffer, Shape, Band) -> Log10Likelihood,
    ) {
        out.clear();
        let Self { workspace, haplotypes, emission } = self;
        let Workspace { plan, rows, prepared, .. } = &mut **workspace;
        if read.is_empty() || plan.rows.fill(read, &*emission).is_none() {
            out.resize(haplotypes.len(), Log10Likelihood::IMPOSSIBLE);
            return;
        }
        let mut strips =
            StripParts { rows: &plan.rows, buffer: rows, emission: &*emission, read, band };
        for (haplotype, prepared) in haplotypes.iter().zip(prepared.iter_mut()) {
            out.push(strips.score(haplotype.borrow(), prepared, &mut kernel));
        }
    }
}
