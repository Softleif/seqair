//! The SAM header text a CRAM file carries, and the reference behind each of
//! its targets.

use super::CramWriteError;
use crate::bam::header::BamHeader;
use crate::cram::reader::parse_read_group_ids;
use crate::fasta::IndexedFastaReader;
use crate::reader::Tid;
use rustc_hash::FxHashMap;
use seqair_types::{Pos0, SmolStr};

/// A reference's id and uppercased bases.
pub(super) type Hashed = (Tid, Vec<u8>);

/// The header the writer writes, and what it needs of it per record.
pub(super) struct CramHeader {
    text: String,
    targets: Vec<Target>,
    /// `@RG` IDs by the index the `RG` data series stores.
    read_groups: FxHashMap<Box<[u8]>, i32>,
}

struct Target {
    name: SmolStr,
    length: u64,
    /// `None` when the FASTA lacks the reference.
    reference: Option<TargetReference>,
}

/// A target the FASTA holds.
pub(super) struct TargetReference {
    m5: M5,
}

/// Where the header's `M5` stands against the FASTA's bases.
enum M5 {
    /// Computed from them, or checked.
    Verified,
    /// The header's, not checked yet: that costs reading the reference,
    /// which the writer does when the first record on it arrives.
    Unchecked(SmolStr),
    /// Not the bases': every record on the reference is refused.
    Mismatch { header: SmolStr, fasta: SmolStr },
}

impl CramHeader {
    // r[impl cram_writer.header_text]
    // r[impl cram_writer.header_m5]
    // r[impl cram_writer.header_m5.reuse]
    /// The header, and the uppercased bases of the first reference whose
    /// `M5` it computed, for the first slice on it.
    pub(super) fn new(
        header: &BamHeader,
        fasta: &mut IndexedFastaReader,
    ) -> Result<(Self, Option<Hashed>), CramWriteError> {
        let text = header.header_text();
        let mut targets: Vec<Target> = header
            .targets()
            .map(|t| Target { name: t.name, length: t.length, reference: None })
            .collect();
        let has_sq = text.lines().any(|l| l.starts_with("@SQ"));
        let mut out =
            String::with_capacity(text.len().saturating_add(targets.len().saturating_mul(80)));
        let mut bases = Vec::new();
        let mut kept = None;
        let mut sq_lines = 0usize;
        // Without `@SQ` lines, they go after a leading `@HD`, else first.
        let mut synthesized = has_sq;
        if !synthesized && !text.starts_with("@HD") {
            synthesize(&mut out, &mut targets, fasta, &mut bases, &mut kept)?;
            synthesized = true;
        }
        for line in text.lines() {
            out.push_str(line);
            if line.starts_with("@SQ") {
                let sq_name = field(line, "SN");
                let target = targets.get_mut(sq_lines).filter(|t| Some(t.name.as_str()) == sq_name);
                let Some(target) = target else {
                    return Err(CramWriteError::HeaderTargetsMismatch {
                        index: sq_lines,
                        sq_name: sq_name.map(SmolStr::new),
                        target: header_target_name(header, sq_lines),
                    });
                };
                let (reference, m5) =
                    describe(target, sq_lines, field(line, "M5"), fasta, &mut bases, &mut kept)?;
                target.reference = reference;
                push_m5(&mut out, m5);
                sq_lines = sq_lines.saturating_add(1);
            }
            out.push('\n');
            if !synthesized {
                synthesize(&mut out, &mut targets, fasta, &mut bases, &mut kept)?;
                synthesized = true;
            }
        }
        if has_sq && sq_lines != targets.len() {
            return Err(CramWriteError::HeaderTargetsMismatch {
                index: sq_lines,
                sq_name: None,
                target: header_target_name(header, sq_lines),
            });
        }

        // r[impl cram_writer.record.read_group]
        let ids = parse_read_group_ids(&out);
        let java_unknown = ids.last().is_some_and(|id| id == "UNKNOWN");
        let keep = ids.len().saturating_sub(usize::from(java_unknown));
        let mut read_groups = FxHashMap::default();
        for (i, id) in ids.iter().take(keep).enumerate() {
            let i = i32::try_from(i).map_err(|_| CramWriteError::ValueTooLarge {
                what: "@RG count",
                value: i as u64,
            })?;
            read_groups.entry(Box::from(id.as_bytes())).or_insert(i);
        }
        let kept = kept.and_then(|(index, bases)| {
            Some((Tid::of_header(u32::try_from(index).ok()?, header)?, bases))
        });
        Ok((Self { text: out, targets, read_groups }, kept))
    }

    pub(super) fn text(&self) -> &str {
        &self.text
    }

    pub(super) fn n_targets(&self) -> usize {
        self.targets.len()
    }

    pub(super) fn target_len(&self, tid: Tid) -> u64 {
        self.targets.get(tid.as_usize()).map_or(0, |t| t.length)
    }

    /// The reference behind `tid`, `None` if the FASTA lacks it.
    pub(super) fn reference(&self, tid: Tid) -> Option<&TargetReference> {
        self.targets.get(tid.as_usize())?.reference.as_ref()
    }

    /// The `RG` series value of an `RG:Z` tag's value, `None` if it must
    /// stay a tag.
    pub(super) fn read_group(&self, id: &[u8]) -> Option<i32> {
        self.read_groups.get(id).copied()
    }

    /// The uppercased bases of `tid` into `bases`, checked against the
    /// header's `M5` the first time.
    pub(super) fn fetch(
        &mut self,
        tid: Tid,
        fasta: &mut IndexedFastaReader,
        bases: &mut Vec<u8>,
    ) -> Result<(), CramWriteError> {
        let Some(target) = self.targets.get_mut(tid.as_usize()) else { return Ok(()) };
        let mismatch = |name: &SmolStr, header: &SmolStr, fasta: &SmolStr| {
            CramWriteError::ReferenceMd5Mismatch {
                contig: name.clone(),
                header: header.clone(),
                fasta: fasta.clone(),
            }
        };
        if let Some(TargetReference { m5: M5::Mismatch { header, fasta } }) = &target.reference {
            return Err(mismatch(&target.name, header, fasta));
        }
        fetch_upper(&target.name, 0..target.length, fasta, bases)?;
        if let Some(reference) = &mut target.reference
            && let M5::Unchecked(header) = &reference.m5
        {
            let found = md5_hex(bases);
            if !header.eq_ignore_ascii_case(&found) {
                let err = mismatch(&target.name, header, &found);
                reference.m5 = M5::Mismatch { header: header.clone(), fasta: found };
                return Err(err);
            }
            reference.m5 = M5::Verified;
        }
        Ok(())
    }
}

/// `@SQ` lines for every target, for a header text without any.
fn synthesize(
    out: &mut String,
    targets: &mut [Target],
    fasta: &mut IndexedFastaReader,
    bases: &mut Vec<u8>,
    kept: &mut Option<(usize, Vec<u8>)>,
) -> Result<(), CramWriteError> {
    for (index, target) in targets.iter_mut().enumerate() {
        out.push_str("@SQ\tSN:");
        out.push_str(&target.name);
        out.push_str("\tLN:");
        out.push_str(itoa::Buffer::new().format(target.length));
        let (reference, m5) = describe(target, index, None, fasta, bases, kept)?;
        target.reference = reference;
        push_m5(out, m5);
        out.push('\n');
    }
    Ok(())
}

/// The reference behind `target` (target `index`), and the `M5` to add to
/// its line: none if the FASTA lacks it or the line has one, which is
/// checked when the bases are first read. The first reference it hashes goes
/// to `kept` whole; the others are hashed a chunk at a time in `scratch`.
fn describe(
    target: &Target,
    index: usize,
    m5: Option<&str>,
    fasta: &mut IndexedFastaReader,
    scratch: &mut Vec<u8>,
    kept: &mut Option<(usize, Vec<u8>)>,
) -> Result<(Option<TargetReference>, Option<SmolStr>), CramWriteError> {
    let Some(entry) = fasta.index().get(&target.name) else { return Ok((None, None)) };
    if entry.length != target.length {
        return Err(CramWriteError::ReferenceLengthMismatch {
            contig: target.name.clone(),
            header: target.length,
            fasta: entry.length,
        });
    }
    match m5 {
        Some(m5) => Ok((Some(TargetReference { m5: M5::Unchecked(SmolStr::new(m5)) }), None)),
        None => {
            // r[impl cram_writer.header_m5.reuse]
            let m5 = if kept.is_none() {
                let mut whole = Vec::new();
                fetch_upper(&target.name, 0..target.length, fasta, &mut whole)?;
                let m5 = md5_hex(&whole);
                *kept = Some((index, whole));
                m5
            } else {
                let mut md5 = md5::Context::new();
                let mut start = 0;
                while start < target.length {
                    let end = start.saturating_add(HASH_CHUNK).min(target.length);
                    fetch_upper(&target.name, start..end, fasta, scratch)?;
                    md5.consume(&scratch);
                    start = end;
                }
                SmolStr::new(format!("{:x}", md5.finalize()))
            };
            Ok((Some(TargetReference { m5: M5::Verified }), Some(m5)))
        }
    }
}

/// The bases of a reference not kept are hashed this many at a time.
const HASH_CHUNK: u64 = if cfg!(test) { 7 } else { 1 << 24 };

/// The uppercased bases of `name` over `range` (0-based) into `bases`.
fn fetch_upper(
    name: &str,
    range: std::ops::Range<u64>,
    fasta: &mut IndexedFastaReader,
    bases: &mut Vec<u8>,
) -> Result<(), CramWriteError> {
    bases.clear();
    let Some(last) = range.end.checked_sub(1).filter(|&l| l >= range.start) else {
        return Ok(());
    };
    let pos = |p: u64| {
        u32::try_from(p)
            .ok()
            .and_then(Pos0::new)
            .ok_or(CramWriteError::ValueTooLarge { what: "reference length", value: range.end })
    };
    fasta.fetch_seq_into(name, (pos(range.start)?..=pos(last)?).into(), bases)?;
    bases.make_ascii_uppercase();
    Ok(())
}

pub(super) fn md5_hex(bases: &[u8]) -> SmolStr {
    SmolStr::new(format!("{:x}", md5::compute(bases)))
}

fn push_m5(out: &mut String, m5: Option<SmolStr>) {
    if let Some(m5) = m5 {
        out.push_str("\tM5:");
        out.push_str(&m5);
    }
}

/// The value of `tag` on a header line.
fn field<'a>(line: &'a str, tag: &str) -> Option<&'a str> {
    line.split('\t').skip(1).find_map(|f| f.strip_prefix(tag)?.strip_prefix(':'))
}

fn header_target_name(header: &BamHeader, index: usize) -> Option<SmolStr> {
    header.target_names().nth(index).map(SmolStr::new)
}
