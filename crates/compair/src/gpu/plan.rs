//! The pairs of one launch, as the flat buffers the shader reads.
//!
//! A read's row tracks and a haplotype's column weights are uploaded once
//! however many pairs use them, and a pair is three indices and a band. That
//! is the shape a variant caller has -- a few haplotypes and many reads per
//! locus -- and it is why the plan takes reads and haplotypes separately
//! rather than pairs of them.

use std::borrow::Borrow;

use bytemuck::{Pod, Zeroable};
use seqair_types::{Base, BaseQuality, Strand};

use super::GpuError;
use crate::{
    banded::Band, emission::Emission, emission::SiteWeights, haplotype::Haplotype, read::Read,
    transitions::Transition, types::ErrorTable,
};

/// A read in a [`GpuPairs`], returned by [`GpuPairs::push_read`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadSlot(u32);

/// A haplotype in a [`GpuPairs`], returned by [`GpuPairs::push_haplotype`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HaplotypeSlot(u32);

/// One pair as the shader reads it. The layout is the WGSL `Pair`'s.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Pod, Zeroable)]
pub(crate) struct PairRecord {
    /// The read's first entry in the row buffer, in rows.
    pub(crate) rows: u32,
    pub(crate) read_len: u32,
    /// The haplotype's first entry in the weight buffer.
    pub(crate) weights: u32,
    pub(crate) hap_len: u32,
    pub(crate) offset: i32,
    pub(crate) half_width: u32,
    /// `1 / h`, the read's free start, divided in `f32` as the CPU kernels do.
    pub(crate) init: f32,
    pub(crate) pad: u32,
}

/// One read row as the shader reads it, twelve bytes.
///
/// The five transitions are not here: each is a function of the row's
/// insertion, deletion and gap qualities, so the row carries those three
/// bytes and the shader looks the transitions up in [`TRANSITIONS`]. Only the
/// emission's two terms, which an `Emission` may derive from anything, are
/// carried as values.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Pod, Zeroable)]
pub(crate) struct RowRecord {
    /// The base selector (`0..4` for a known base, `4` for `N`) in the low
    /// byte, then the insertion, deletion and gap-continuation quality bytes.
    pub(crate) packed: u32,
    /// `(1 - eps) - eps / 3`, as in `banded::RowTracks`.
    pub(crate) spread: f32,
    pub(crate) mismatched: f32,
}

/// Offsets into [`TRANSITIONS`]: `match_to_match` by insertion and deletion
/// quality, then the four one-quality transitions.
pub(crate) const MATCH_TO_INSERTION: usize = 256 * 256;
pub(crate) const MATCH_TO_DELETION: usize = MATCH_TO_INSERTION + 256;
pub(crate) const INDEL_TO_MATCH: usize = MATCH_TO_DELETION + 256;
pub(crate) const GAP_CONTINUATION: usize = INDEL_TO_MATCH + 256;
const TRANSITION_TABLE: usize = GAP_CONTINUATION + 256;

/// Every transition a read row can carry, narrowed to `f32` exactly as the
/// CPU plans narrow `Read`'s: the same `Transition::from_qualities`, so a
/// lookup is the value the CPU kernels use, bit for bit.
pub(crate) static TRANSITIONS: std::sync::LazyLock<Box<[f32]>> = std::sync::LazyLock::new(|| {
    let mut table = vec![0.0f32; TRANSITION_TABLE];
    let errors = ErrorTable::get();
    for insertion in 0..=255u8 {
        for deletion in 0..=255u8 {
            let t = Transition::from_qualities(
                errors,
                BaseQuality::from_byte(insertion),
                BaseQuality::from_byte(deletion),
                BaseQuality::from_byte(0),
            );
            if let Some(slot) = table.get_mut(usize::from(insertion) * 256 + usize::from(deletion))
            {
                *slot = narrow(t.match_to_match);
            }
        }
    }
    for quality in 0..=255u8 {
        let t = Transition::from_qualities(
            errors,
            BaseQuality::from_byte(quality),
            BaseQuality::from_byte(quality),
            BaseQuality::from_byte(quality),
        );
        let at = usize::from(quality);
        for (offset, value) in [
            (MATCH_TO_INSERTION, t.match_to_insertion),
            (MATCH_TO_DELETION, t.match_to_deletion),
            (INDEL_TO_MATCH, t.indel_to_match),
            (GAP_CONTINUATION, t.gap_continuation),
        ] {
            if let Some(slot) = table.get_mut(offset + at) {
                *slot = narrow(value);
            }
        }
    }
    table.into_boxed_slice()
});

#[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point of this kernel")]
fn narrow(value: f64) -> f32 {
    value as f32
}

/// Read bases have this many weight tracks per haplotype, one per known base.
const TRACKS: usize = 4;

/// The selector for a read `N`, which matches every column at full weight.
const SELECT_N: u32 = 4;

#[derive(Debug, Clone, Copy)]
struct ReadEntry {
    first_row: u32,
    len: u32,
    strand: Strand,
}

#[derive(Debug, Clone, Copy)]
struct HaplotypeEntry {
    first_weight: u32,
    len: u32,
    strand: Strand,
}

/// A pair and the band class it runs in; `None` for a pair no alignment can
/// relate, which never reaches the GPU.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Pending {
    pub(crate) record: PairRecord,
    pub(crate) class: Option<u32>,
}

/// Band half-widths are rounded up to a multiple of this to pick a shader, so
/// that a caller whose widths vary by a few columns compiles a few shaders
/// rather than one per width. A pair narrower than its class masks the
/// positions past its own band, which costs their arithmetic and changes no
/// value.
pub(crate) const CLASS_STEP: u32 = 4;

/// Reads, haplotypes and the pairs between them, for one
/// [`GpuAligner::submit`](super::GpuAligner::submit).
///
/// The emission is folded in when a read or haplotype is pushed, the way the
/// CPU kernels fold it into their plan: a read's rows carry the emission's
/// `epsilon`, a haplotype's weights carry its `site_weights` for one strand.
/// A pair is only meaningful between a read and a haplotype pushed with the
/// **same** emission, which the plan cannot check; the strand it does check.
#[derive(Debug, Default)]
pub struct GpuPairs {
    pub(crate) rows: Vec<RowRecord>,
    pub(crate) weights: Vec<f32>,
    reads: Vec<ReadEntry>,
    haplotypes: Vec<HaplotypeEntry>,
    pub(crate) pairs: Vec<Pending>,
}

impl GpuPairs {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Empties the plan, keeping its allocations.
    pub fn clear(&mut self) {
        self.rows.clear();
        self.weights.clear();
        self.reads.clear();
        self.haplotypes.clear();
        self.pairs.clear();
    }

    /// The number of pairs pushed.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// The bytes a launch of this plan uploads.
    #[must_use]
    pub fn upload_bytes(&self) -> usize {
        size_of_val(self.rows.as_slice())
            + size_of_val(self.weights.as_slice())
            + self.pairs.len() * size_of::<PairRecord>()
    }

    /// Adds every read, haplotype and pair of `other`, and returns where its
    /// pairs landed: the scores of `other`'s pairs are that range of the
    /// launch's scores, in `other`'s order.
    ///
    /// This is how several callers' plans become one launch: each fills its
    /// own `GpuPairs`, and one of them -- or a thread that owns the aligner --
    /// appends the rest. `other`'s slots name its own reads and haplotypes,
    /// not the copies here, so pairs between the two plans are made before
    /// appending, not after.
    pub fn append(&mut self, other: &GpuPairs) -> Result<core::ops::Range<usize>, GpuError> {
        let rows = index_u32(self.rows.len(), "read rows")?;
        let weights = index_u32(self.weights.len(), "haplotype weights")?;
        index_u32(self.rows.len() + other.rows.len(), "read rows")?;
        index_u32(self.reads.len() + other.reads.len(), "reads")?;
        index_u32(self.haplotypes.len() + other.haplotypes.len(), "haplotypes")?;
        // The shader forms its weight index in `i32`, as `push_haplotype` checks.
        if i32::try_from(self.weights.len() + other.weights.len()).is_err() {
            return Err(GpuError::TooLarge { what: "haplotype weights" });
        }
        let first = self.pairs.len();
        self.rows.extend_from_slice(&other.rows);
        self.weights.extend_from_slice(&other.weights);
        self.reads.extend(
            other
                .reads
                .iter()
                .map(|entry| ReadEntry { first_row: entry.first_row + rows, ..*entry }),
        );
        self.haplotypes.extend(
            other.haplotypes.iter().map(|entry| HaplotypeEntry {
                first_weight: entry.first_weight + weights,
                ..*entry
            }),
        );
        self.pairs.extend(other.pairs.iter().map(|pending| Pending {
            record: PairRecord {
                rows: pending.record.rows + rows,
                weights: pending.record.weights + weights,
                ..pending.record
            },
            class: pending.class,
        }));
        Ok(first..self.pairs.len())
    }

    /// Every read against every haplotype, the shape of
    /// [`Workspace::align_reads`](crate::Workspace::align_reads), and where
    /// their scores will land: read `r` against haplotype `h` is
    /// `range.start + r * haplotypes.len() + h` of the launch's scores, as
    /// `align_reads` lays them out.
    ///
    /// Each read is uploaded once, and each haplotype once per strand its
    /// reads are on, so this is [`push_read`], [`push_haplotype`] and
    /// [`push_pair`] with the strand bookkeeping done. An error (a plan past
    /// the kernel's 32-bit indices) leaves part of the call pushed.
    ///
    /// [`push_read`]: Self::push_read
    /// [`push_haplotype`]: Self::push_haplotype
    /// [`push_pair`]: Self::push_pair
    pub fn push_reads<H: Borrow<Haplotype>, R: Borrow<Read>, E: Emission>(
        &mut self,
        haplotypes: &[H],
        reads: &[(R, Band)],
        emission: &E,
    ) -> Result<core::ops::Range<usize>, GpuError> {
        let start = self.pairs.len();
        // A strand's haplotypes are pushed together, so their slots are
        // consecutive and the first one names them all.
        let (mut ot, mut ob) = (None, None);
        for (read, band) in reads {
            let read = read.borrow();
            let strand = read.strand();
            let first = if strand == Strand::OB { &mut ob } else { &mut ot };
            let first = match *first {
                Some(slot) => slot,
                None => {
                    let slot = index_u32(self.haplotypes.len(), "haplotypes")?;
                    for haplotype in haplotypes {
                        self.push_haplotype(haplotype.borrow(), strand, emission)?;
                    }
                    *first.insert(slot)
                }
            };
            let read_slot = self.push_read(read, emission)?;
            for index in 0..haplotypes.len() {
                let slot = first
                    .checked_add(index_u32(index, "haplotypes")?)
                    .ok_or(GpuError::TooLarge { what: "haplotypes" })?;
                self.push_pair(read_slot, HaplotypeSlot(slot), *band)?;
            }
        }
        Ok(start..self.pairs.len())
    }

    /// Adds a read's row tracks: `Plan::fill`'s per-row arithmetic, row for
    /// row.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the f32 narrowing is the point of this kernel"
    )]
    pub fn push_read<E: Emission>(
        &mut self,
        read: &Read,
        emission: &E,
    ) -> Result<ReadSlot, GpuError> {
        let slot = ReadSlot(index_u32(self.reads.len(), "reads")?);
        let first_row = index_u32(self.rows.len(), "read rows")?;
        let len = index_u32(read.len(), "read length")?;
        index_u32(self.rows.len() + read.len(), "read rows")?;
        self.rows.reserve(read.len());
        let qualities =
            read.insertion_quals().iter().zip(read.deletion_quals()).zip(read.gap_quals());
        for (index, ((insertion, deletion), gap)) in qualities.enumerate() {
            let Some(observation) = read.observation(index) else {
                return Err(GpuError::TooLarge { what: "read length" });
            };
            let eps = crate::emission::epsilon(emission, observation);
            let select = observation
                .base
                .known_index()
                .and_then(|k| u32::try_from(k).ok())
                .unwrap_or(SELECT_N);
            let packed = select
                | u32::from(insertion.as_byte()) << 8
                | u32::from(deletion.as_byte()) << 16
                | u32::from(gap.as_byte()) << 24;
            self.rows.push(RowRecord {
                packed,
                spread: ((1.0 - eps) - eps / 3.0) as f32,
                mismatched: (eps / 3.0) as f32,
            });
        }
        self.reads.push(ReadEntry { first_row, len, strand: read.strand() });
        Ok(slot)
    }

    /// Adds a haplotype's column weights for reads on `strand`: one track per
    /// known read base, each holding what `banded::prior` selects for that
    /// base at every column.
    pub fn push_haplotype<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        strand: Strand,
        emission: &E,
    ) -> Result<HaplotypeSlot, GpuError> {
        let slot = HaplotypeSlot(index_u32(self.haplotypes.len(), "haplotypes")?);
        let first_weight = index_u32(self.weights.len(), "haplotype weights")?;
        let h = haplotype.len();
        let len = index_u32(h, "haplotype length")?;
        let end = self
            .weights
            .len()
            .checked_add(
                h.checked_mul(TRACKS).ok_or(GpuError::TooLarge { what: "haplotype weights" })?,
            )
            .ok_or(GpuError::TooLarge { what: "haplotype weights" })?;
        // The shader forms its weight index in `i32`.
        if i32::try_from(end).is_err() {
            return Err(GpuError::TooLarge { what: "haplotype weights" });
        }
        let start = self.weights.len();
        self.weights.resize(end, 0.0);
        let Some(tracks) = self.weights.get_mut(start..end) else {
            return Err(GpuError::TooLarge { what: "haplotype weights" });
        };
        for index in 0..h {
            let Some(site) = haplotype.site(index) else {
                return Err(GpuError::TooLarge { what: "haplotype length" });
            };
            let site = emission.site_weights(site, strand);
            for observed in Base::KNOWN {
                let Some(track) = observed.known_index() else { continue };
                if let Some(slot) = tracks.get_mut(track * h + index) {
                    *slot = weight(site, observed);
                }
            }
        }
        self.haplotypes.push(HaplotypeEntry { first_weight, len, strand });
        Ok(slot)
    }

    /// Adds the pair of `read` and `haplotype` under `band`. Scores come back
    /// in the order pairs were pushed.
    pub fn push_pair(
        &mut self,
        read: ReadSlot,
        haplotype: HaplotypeSlot,
        band: Band,
    ) -> Result<(), GpuError> {
        let read_entry = self.reads.get(read.0 as usize).copied().ok_or(GpuError::UnknownSlot)?;
        let hap =
            self.haplotypes.get(haplotype.0 as usize).copied().ok_or(GpuError::UnknownSlot)?;
        if read_entry.strand != hap.strand {
            return Err(GpuError::StrandMismatch {
                read: read_entry.strand,
                haplotype: hap.strand,
            });
        }
        let half_width =
            u32::try_from(band.half_width).map_err(|_| GpuError::TooLarge { what: "band" })?;
        let (h, r) = (hap.len as usize, read_entry.len as usize);
        // No read, no haplotype, or a band that misses the haplotype: the CPU
        // kernels return `IMPOSSIBLE` without running, and so does this one.
        // Keeping these off the GPU is also what bounds every offset the
        // shader adds in `i32`.
        let class = band.columns(h, r).map(|_| half_width.next_multiple_of(CLASS_STEP));
        #[allow(
            clippy::cast_precision_loss,
            reason = "a haplotype is at most a few thousand bases"
        )]
        let init = if h == 0 { 0.0 } else { 1.0 / h as f32 };
        let record = PairRecord {
            rows: read_entry.first_row,
            read_len: read_entry.len,
            weights: hap.first_weight,
            hap_len: hap.len,
            offset: band.offset(),
            half_width,
            init,
            pad: 0,
        };
        self.pairs.push(Pending { record, class });
        Ok(())
    }
}

/// What `banded::prior` selects as the weight of a read base `observed` at a
/// column with these site weights, in `f32`: the same four outcomes in the
/// same precedence, narrowed the way `Plan::fill` narrows them.
#[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point of this kernel")]
fn weight(site: SiteWeights, observed: Base) -> f32 {
    if site.base == Base::Unknown {
        return 1.0;
    }
    match site.converted {
        None if observed == site.base => 1.0,
        None => 0.0,
        Some(converted) if observed == converted => site.rate as f32,
        Some(_) if observed == site.base => (1.0 - site.rate) as f32,
        Some(_) => 0.0,
    }
}

fn index_u32(value: usize, what: &'static str) -> Result<u32, GpuError> {
    u32::try_from(value).map_err(|_| GpuError::TooLarge { what })
}

#[cfg(test)]
mod tests {
    use super::GpuPairs;
    use crate::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand};

    /// `push_reads` is the loop every caller wrote by hand: each haplotype
    /// pushed once per strand, then each read and its pairs, read-major.
    #[test]
    fn push_reads_is_the_hand_written_loop() -> Result<(), Box<dyn std::error::Error>> {
        let emission = StandardEmission::default();
        let haplotypes =
            [Haplotype::from_ascii(b"ACGTTGCAACGT"), Haplotype::from_ascii(b"ACGTTGCCAACGT")];
        let read = |seq: &[u8], strand| {
            Read::uniform(
                Base::from_ascii_vec(seq.to_vec()),
                &vec![BaseQuality::from_byte(30); seq.len()],
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(10),
                strand,
            )
        };
        let reads = [
            (read(b"GTTGCA", Strand::OB)?, Band::anchored(2)),
            (read(b"ACGTTG", Strand::OT)?, Band::anchored(0)),
            (read(b"TGCCAA", Strand::OB)?, Band::new(8, 4)?),
        ];

        let mut by_hand = GpuPairs::new();
        let mut slots = Vec::new();
        for strand in [Strand::OB, Strand::OT] {
            for haplotype in &haplotypes {
                slots.push((strand, by_hand.push_haplotype(haplotype, strand, &emission)?));
            }
        }
        for (read, band) in &reads {
            let read_slot = by_hand.push_read(read, &emission)?;
            for &(strand, haplotype) in &slots {
                if strand == read.strand() {
                    by_hand.push_pair(read_slot, haplotype, *band)?;
                }
            }
        }

        let mut plan = GpuPairs::new();
        plan.push_read(&reads[0].0, &emission)?;
        let range = plan.push_reads(&haplotypes, &reads, &emission)?;
        assert_eq!(range, 0..reads.len() * haplotypes.len());

        let records = |plan: &GpuPairs| {
            plan.pairs
                .iter()
                .map(|p| {
                    (
                        p.record.read_len,
                        p.record.hap_len,
                        p.record.offset,
                        p.record.half_width,
                        p.class,
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(records(&plan), records(&by_hand));
        // The leading read shifts every row index by its length, nothing else.
        assert_eq!(plan.rows.get(reads[0].0.len()..), Some(&by_hand.rows[..]));
        assert_eq!(plan.weights, by_hand.weights);
        Ok(())
    }
}
