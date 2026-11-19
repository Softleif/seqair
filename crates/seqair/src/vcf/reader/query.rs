//! Region queries and partitions: which chunks a query reads, which of
//! their records it yields (`docs/spec/2-vcf-1-reader.md`, "Region queries",
//! "Partitions").

use super::error::VcfReadError;
use super::file::{FileShared, Input};
use super::region::{Match, Region, RegionContig, RegionEnd};
use crate::bam::bgzf::{BgzfError, VirtualOffset};
use crate::bam::index::{Chunk, merge_overlapping_chunks};
use crate::bam::region_buf::RegionBuf;
use seqair_types::Pos0;
use std::num::NonZeroUsize;

/// What a query does with a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Verdict {
    Take,
    Skip,
    /// The record starts past the last range of its group: no later record
    /// of the group's chunks can be in it.
    Stop,
}

/// One contig's merged ranges, 0-based half-open, sorted and disjoint.
#[derive(Debug)]
struct Group {
    tid: u32,
    ranges: Vec<(i64, i64)>,
}

/// The ranges a query yields records of, and how it matches them.
#[derive(Debug)]
pub(crate) struct Query {
    /// In file order: chunks are tagged with their group's index.
    groups: Vec<Group>,
    matching: Match,
}

impl Query {
    pub(crate) fn matching(&self) -> Match {
        self.matching
    }

    // r[impl vcf_reader.query.overlap+2]
    // r[impl vcf_reader.query.many+2]
    /// Whether the record on contig `tid` at 0-based `pos0` (−1 for the
    /// telomere), read from a chunk of `group`, is in the query. `rlen` is
    /// asked for only for a record that starts before a range.
    pub(crate) fn verdict(
        &self,
        group: usize,
        tid: u32,
        pos0: i64,
        rlen: impl FnOnce() -> i64,
    ) -> Verdict {
        let Some(g) = self.groups.get(group) else { return Verdict::Skip };
        if g.tid != tid {
            return Verdict::Skip;
        }
        // As tabix: the telomere counts as the first base.
        let start = pos0.max(0);
        let i = g.ranges.partition_point(|&(_, end)| end <= start);
        let Some(&(range_start, _)) = g.ranges.get(i) else { return Verdict::Stop };
        if range_start <= start {
            return Verdict::Take;
        }
        if self.matching == Match::StartsIn {
            return Verdict::Skip;
        }
        let end = pos0.saturating_add(rlen()).max(start.saturating_add(1));
        if end > range_start { Verdict::Take } else { Verdict::Skip }
    }
}

/// The contig-table tid `contig` names in `file`.
fn resolve(file: &FileShared, contig: &RegionContig) -> Result<u32, VcfReadError> {
    match contig {
        RegionContig::Name(name) => file
            .contigs
            .find(&file.header, name)
            .map(|c| c.tid())
            .ok_or_else(|| VcfReadError::UnknownContig { name: name.clone() }),
        RegionContig::Id(id) if id.header != file.header.id() => {
            Err(VcfReadError::ForeignKey { name: id.name.clone() })
        }
        RegionContig::Id(id) => Ok(id.tid()),
    }
}

// r[impl vcf_reader.query+2]
// r[impl vcf_reader.query.contigs+2]
// r[impl vcf_reader.query.many+2]
/// Plan a query of `regions`: its ranges, merged per contig, and the
/// chunks to read, in file order, each tagged with its group.
///
/// # Errors
/// [`VcfReadError::NoIndex`], [`VcfReadError::UnknownContig`],
/// [`VcfReadError::ForeignKey`].
pub(crate) fn plan(
    file: &FileShared,
    regions: impl IntoIterator<Item = Region>,
    matching: Match,
) -> Result<(Query, Vec<(Chunk, usize)>), VcfReadError> {
    let index = file.index.as_ref().ok_or(VcfReadError::NoIndex)?;
    let mut by_tid: Vec<(u32, (i64, i64))> = Vec::new();
    for region in regions {
        let tid = resolve(file, region.target())?;
        if !region.is_empty() {
            by_tid.push((tid, region.span()));
        }
    }
    by_tid.sort_unstable();
    let mut groups: Vec<(Group, Vec<Chunk>)> = Vec::new();
    for (tid, (start, end)) in by_tid {
        match groups.last_mut() {
            Some((g, _)) if g.tid == tid => match g.ranges.last_mut() {
                // Overlapping or adjacent: one range.
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => g.ranges.push((start, end)),
            },
            _ => groups.push((Group { tid, ranges: vec![(start, end)] }, Vec::new())),
        }
    }
    // r[impl vcf_reader.query.contigs+2]
    // A contig the index does not list has no records: no chunks, no group.
    groups.retain_mut(|(g, chunks)| {
        if let Some(reference) = index.reference(g.tid) {
            for &(start, end) in &g.ranges {
                chunks.extend(index.chunks(reference, start, end));
            }
            chunks.sort_by_key(|c| c.begin);
            merge_overlapping_chunks(chunks);
        }
        !chunks.is_empty()
    });
    groups.sort_by_key(|(_, chunks)| chunks.first().map(|c| c.begin));
    let mut tagged = Vec::new();
    let mut query = Query { groups: Vec::with_capacity(groups.len()), matching };
    for (i, (group, chunks)) in groups.into_iter().enumerate() {
        tagged.extend(chunks.into_iter().map(|c| (c, i)));
        query.groups.push(group);
    }
    Ok((query, tagged))
}

/// A query's input: its chunks, read through one [`RegionBuf`] that the
/// reader keeps across [`requery`](super::Reader::requery)s.
pub(crate) struct RegionInput {
    buf: RegionBuf<'static, Input>,
    chunks: Vec<(Chunk, usize)>,
    /// The next chunk to seek to.
    next: usize,
    /// Where the chunk being read ends; `None` between chunks.
    end: Option<VirtualOffset>,
}

impl RegionInput {
    pub(crate) fn new(input: Input, chunks: Vec<(Chunk, usize)>) -> Result<Self, BgzfError> {
        let plain: Vec<Chunk> = chunks.iter().map(|&(c, _)| c).collect();
        let buf = RegionBuf::owning(input, &plain)?;
        Ok(Self { buf, chunks, next: 0, end: None })
    }

    /// Read `chunks` instead, keeping the file handle and block cache.
    pub(crate) fn replan(&mut self, chunks: Vec<(Chunk, usize)>) {
        let plain: Vec<Chunk> = chunks.iter().map(|&(c, _)| c).collect();
        self.buf.replan(&plain);
        self.chunks = chunks;
        self.next = 0;
        self.end = None;
    }

    /// Whether a record starts here: inside the current chunk, else at the
    /// start of the next one.
    pub(crate) fn has_record(&mut self) -> Result<bool, BgzfError> {
        loop {
            if let Some(end) = self.end {
                // `fill_buf` may load the next block, so the end is checked
                // against where the returned bytes start.
                if self.buf.virtual_offset() < end {
                    let (at, bytes) = self.buf.fill_buf()?;
                    if !bytes.is_empty() && at < end {
                        return Ok(true);
                    }
                }
                self.end = None;
            }
            let Some(&(chunk, _)) = self.chunks.get(self.next) else { return Ok(false) };
            self.next = self.next.saturating_add(1);
            self.buf.seek_virtual(chunk.begin)?;
            self.end = Some(chunk.end);
        }
    }

    /// The group of the chunk being read.
    pub(crate) fn group(&self) -> Option<usize> {
        self.chunks.get(self.next.checked_sub(1)?).map(|&(_, g)| g)
    }

    /// Leave `group`'s chunks.
    pub(crate) fn skip_group(&mut self, group: usize) {
        if self.group() == Some(group) {
            self.end = None;
            while self.chunks.get(self.next).is_some_and(|&(_, g)| g == group) {
                self.next = self.next.saturating_add(1);
            }
        }
    }

    pub(crate) fn fill(&mut self) -> Result<&[u8], BgzfError> {
        self.buf.fill_buf().map(|(_, bytes)| bytes)
    }

    pub(crate) fn advance(&mut self, n: usize) {
        self.buf.consume(n);
    }

    pub(crate) fn virtual_offset(&self) -> VirtualOffset {
        self.buf.virtual_offset()
    }
}

/// A partition of a file: one region of one indexed contig, read
/// [`Match::StartsIn`] so that every record is in exactly one partition
/// (`r[vcf_reader.partitions+2]`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Partition {
    region: Region,
}

impl Partition {
    /// The partition's region; its contig is a [`ContigId`](crate::vcf::ContigId)
    /// of the file it was made from.
    pub fn region(&self) -> &Region {
        &self.region
    }
}

// r[impl vcf_reader.partitions+2]
// r[impl vcf_reader.partitions.complete+2]
/// Tile every indexed contig of `file` into about `n` partitions in all, of
/// roughly equal compressed bytes, in file order.
pub(crate) fn partitions(
    file: &FileShared,
    n: NonZeroUsize,
) -> Result<Vec<Partition>, VcfReadError> {
    let index = file.index.as_ref().ok_or(VcfReadError::NoIndex)?;
    // (tid, reference, first byte, bytes) of every contig with records.
    let mut contigs: Vec<(u32, u32, u64, u64)> = index
        .contigs()
        .filter_map(|tid| {
            let reference = index.reference(tid)?;
            let chunks = index.chunks(reference, 0, i64::MAX);
            let first = chunks.iter().map(|c| c.begin.block_offset()).min()?;
            let last = chunks.iter().map(|c| c.end.block_offset()).max()?;
            Some((tid, reference, first, last.saturating_sub(first).max(1)))
        })
        .collect();
    contigs.sort_by_key(|&(_, _, first, _)| first);
    let total: u128 = contigs.iter().map(|&(.., bytes)| u128::from(bytes)).sum::<u128>().max(1);
    let n = u128::try_from(n.get()).unwrap_or(u128::MAX);
    let mut partitions = Vec::new();
    for (tid, reference, first, bytes) in contigs {
        let Some(contig) = file.contigs.get(tid) else { continue };
        // This contig's share of `n`, at least one.
        let k = n
            .saturating_mul(u128::from(bytes))
            .saturating_add(total / 2)
            .checked_div(total)
            .unwrap_or(1)
            .max(1);
        let mut boundaries: Vec<Pos0> = Vec::new();
        let mut leaves = index.leaf_offsets(reference).into_iter().peekable();
        let mut reached = 0u64;
        for j in 1..k {
            let share = u128::from(bytes).saturating_mul(j).checked_div(k).unwrap_or(0);
            let target = u128::from(first).saturating_add(share);
            // The first leaf window past the last boundary whose records start
            // at or after the target byte.
            let found = leaves.by_ref().find(|&(pos, offset)| {
                reached = reached.max(offset.block_offset());
                u128::from(reached) >= target && pos > boundaries.last().map_or(0, |b| b.as_u64())
            });
            let Some(pos) = found.and_then(|(pos, _)| Pos0::try_from(pos).ok()) else { break };
            boundaries.push(pos);
        }
        let mut start = Pos0::ZERO;
        let contig = RegionContig::Id(contig);
        for end in boundaries.into_iter().map(RegionEnd::Before).chain([RegionEnd::ContigEnd]) {
            let region = match end {
                RegionEnd::Before(end) => Region::new(contig.clone(), start..end),
                RegionEnd::ContigEnd => Region::new(contig.clone(), start..),
            };
            if let RegionEnd::Before(end) = end {
                start = end;
            }
            // Boundaries strictly increase from zero: never reversed.
            if let Ok(region) = region {
                partitions.push(Partition { region });
            }
        }
    }
    Ok(partitions)
}
