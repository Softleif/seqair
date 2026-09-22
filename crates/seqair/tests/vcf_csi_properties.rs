//! Generated CSI output from the VCF writer, refereed by bcftools.
//!
//! `vcf_index_comparison.rs` generates for TBI. CSI is the index the writer
//! actually hands back from `Writer::finish` for both `.vcf.gz` and `.bcf`, and
//! it was covered by four hand-written records on two contigs.
//!
//! Byte equality against `bcftools index -c` is not a property that can hold:
//! htslib walks its bins out of a hash table, so their order is arbitrary,
//! it ends the last chunk at the BGZF EOF rather than at the last record, and
//! it writes an `n_no_coor` trailer seqair omits. What *can* be compared is
//! behaviour — the two indexes must send bcftools to the same records — and the
//! parts of the binary layout that have a single right answer: the header, the
//! reference count, and the pseudo-bin's mapped/unmapped tallies.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::vcf::{Alleles, ContigDef, ContigId, OutputFormat, VcfHeader, Writer};
use seqair_types::{Base, Pos1};
use std::path::Path;
use std::sync::Arc;

/// The largest contig a `depth=5`, `min_shift=14` index can address: bins run
/// out at `2^(14 + 3*5)`, so a longer contig forces a deeper index — see
/// `a_contig_past_the_depth_5_bin_limit_is_still_reachable`.
const BIN_LIMIT: u32 = 1 << 29;

const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

/// A record, kept as the coordinates the index has to get right.
#[derive(Debug, Clone)]
struct Rec {
    contig: usize,
    pos: u32,
    /// REF length: 1 for an SNV, more for a deletion, which is how a record
    /// comes to span more than the bin its start position falls in.
    ref_len: u32,
    deleted: Vec<Base>,
}

impl Rec {
    /// The closed interval of reference positions the record covers.
    fn end(&self) -> u32 {
        self.pos.saturating_add(self.ref_len).saturating_sub(1)
    }

    fn overlaps(&self, contig: usize, start: u32, end: u32) -> bool {
        self.contig == contig && self.pos <= end && self.end() >= start
    }
}

#[derive(Debug, Clone)]
struct Layout {
    /// `(name, length)` for each contig in header order. Some may hold no
    /// records: a reference with nothing in it is its own case in CSI output.
    contigs: Vec<(String, u32)>,
    /// Sorted by `(contig, pos)` — the writer rejects anything else.
    records: Vec<Rec>,
}

#[hegel::composite]
fn arb_layout(tc: &TestCase) -> Layout {
    let n_contigs = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(3));
    let contigs: Vec<(String, u32)> = (0..n_contigs)
        .map(|i| {
            // Lengths spread across the bin levels a record can land on, and
            // over the depth-5 limit as well: past `BIN_LIMIT` the index has
            // to deepen, and that is the band this file exists for. The top
            // band straddles the boundary rather than sitting clear of it,
            // because the off-by-one there is real — a contig of exactly
            // `BIN_LIMIT - 256` still needs depth 6, since a record may end
            // past the declared length.
            let len = tc.draw_silent(hegel::one_of!(
                gs::integers::<u32>().min_value(100_000).max_value(5_000_000),
                gs::integers::<u32>().min_value(5_000_000).max_value(BIN_LIMIT - 1),
                gs::integers::<u32>().min_value(BIN_LIMIT - 1_000).max_value(700_000_000),
            ));
            (format!("chr{}", i.saturating_add(1)), len)
        })
        .collect();

    let n_records = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(12));
    let mut raw: Vec<(usize, u32)> = (0..n_records)
        .map(|_| {
            let c = tc.draw_silent(gs::integers::<usize>().max_value(n_contigs.saturating_sub(1)));
            let len = contigs[c].1;
            // Leave room for the longest REF a deletion can produce.
            let pos =
                tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(len.saturating_sub(8)));
            (c, pos)
        })
        .collect();
    raw.sort_unstable();

    let records = raw
        .into_iter()
        .map(|(contig, pos)| {
            let deleted = tc.draw_silent(gs::vecs(gs::sampled_from(&ACGT)).max_size(4));
            let ref_len = u32::try_from(deleted.len()).unwrap_or(0).saturating_add(1);
            Rec { contig, pos, ref_len, deleted }
        })
        .collect();

    Layout { contigs, records }
}

/// A region to query: a contig, and an interval biased towards the records —
/// anchoring on a real position rather than drawing uniformly is what keeps
/// most generated queries from being trivially empty.
#[hegel::composite]
fn arb_region(tc: &TestCase, layout: Layout) -> (usize, u32, u32) {
    let contig =
        tc.draw_silent(gs::integers::<usize>().max_value(layout.contigs.len().saturating_sub(1)));
    let len = layout.contigs[contig].1;
    let anchor = if layout.records.is_empty() {
        tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(len))
    } else {
        let pick = tc.draw_silent(gs::sampled_from(&layout.records));
        pick.pos
    };
    let back = tc.draw_silent(gs::integers::<u32>().max_value(2000));
    let span = tc.draw_silent(hegel::one_of!(
        gs::integers::<u32>().max_value(20),
        gs::integers::<u32>().max_value(200_000),
    ));
    let start = anchor.saturating_sub(back).max(1);
    let end = start.saturating_add(span).min(len);
    (contig, start, end)
}

// ── Writing ────────────────────────────────────────────────────────────

fn write_indexed(dir: &Path, layout: &Layout, format: OutputFormat) -> std::path::PathBuf {
    let mut builder = VcfHeader::builder();
    let ids: Vec<ContigId> = layout
        .contigs
        .iter()
        .map(|(name, len)| {
            builder.register_contig(name, ContigDef { length: Some(u64::from(*len)) }).unwrap()
        })
        .collect();
    let header = Arc::new(builder.samples().build().unwrap());

    let path = dir.join(match format {
        OutputFormat::Bcf => "generated.bcf",
        _ => "generated.vcf.gz",
    });
    let file = std::fs::File::create(&path).unwrap();
    let writer = Writer::new(file, format);
    let mut writer = writer.write_header(&header).unwrap();
    for rec in &layout.records {
        let alleles = if rec.deleted.is_empty() {
            Alleles::snv(Base::A, Base::T).unwrap()
        } else {
            Alleles::deletion(Base::A, &rec.deleted).unwrap()
        };
        writer
            .begin_record(&ids[rec.contig], Pos1::new(rec.pos).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass()
            .emit()
            .unwrap();
    }
    let (_inner, index) = writer.finish().unwrap();
    let index = index.expect("a compressed writer hands back an index");
    index.write_to_path(csi_path(&path)).unwrap();
    path
}

/// Where htslib looks for the index of `path`.
fn csi_path(path: &Path) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".csi");
    s.into()
}

/// `bcftools view -H -r <region>`, as `(CHROM, POS)` pairs in file order.
fn query(path: &Path, region: &str) -> Vec<(String, u32)> {
    let out = std::process::Command::new("bcftools")
        .args(["view", "-H", "-r", region])
        .arg(path)
        .output()
        .expect("bcftools not found");
    assert!(
        out.status.success(),
        "bcftools view -r {region} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut f = l.split('\t');
            let chrom = f.next().unwrap().to_owned();
            let pos = f.next().unwrap().parse().unwrap();
            (chrom, pos)
        })
        .collect()
}

fn expected(layout: &Layout, contig: usize, start: u32, end: u32) -> Vec<(String, u32)> {
    layout
        .records
        .iter()
        .filter(|r| r.overlaps(contig, start, end))
        .map(|r| (layout.contigs[contig].0.clone(), r.pos))
        .collect()
}

fn decompress(path: &Path) -> Vec<u8> {
    let file = std::fs::File::open(path).unwrap();
    let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(file);
    let mut data = Vec::new();
    reader.read_to_end(&mut data).unwrap();
    data
}

// ── Properties ─────────────────────────────────────────────────────────

// r[verify csi.write_format]
// r[verify csi.write_tabix_aux]
// r[verify index_builder.csi_format]
/// The CSI `Writer::finish` returns for a `.vcf.gz` must make `bcftools view
/// -r` return exactly the records overlapping the region — no more, and in
/// particular no fewer.
#[hegel::test(test_cases = 25)]
fn seqair_csi_makes_bcftools_return_exactly_the_overlapping_records(tc: TestCase) {
    let layout = tc.draw(arb_layout().print_as_debug());
    let (contig, start, end) = tc.draw(arb_region(layout.clone()).print_as_debug());

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::VcfGz);
    let region = format!("{}:{start}-{end}", layout.contigs[contig].0);

    assert_eq!(query(&path, &region), expected(&layout, contig, start, end), "region {region}");

    let hits = expected(&layout, contig, start, end).len();
    tc.event(if hits == 0 { "empty region" } else { "region with hits" });
    if layout.records.iter().any(|r| r.ref_len > 1) {
        tc.event("a record spans more than one base");
    }
}

// r[verify csi.write_format]
// r[verify index_builder.csi_format]
/// The same for BCF, whose CSI carries no tabix aux block and lists every
/// reference rather than only the ones with records.
#[hegel::test(test_cases = 25)]
fn seqair_bcf_csi_makes_bcftools_return_exactly_the_overlapping_records(tc: TestCase) {
    let layout = tc.draw(arb_layout().print_as_debug());
    let (contig, start, end) = tc.draw(arb_region(layout.clone()).print_as_debug());

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::Bcf);
    let region = format!("{}:{start}-{end}", layout.contigs[contig].0);

    assert_eq!(query(&path, &region), expected(&layout, contig, start, end), "region {region}");
}

// r[verify csi.write]
// r[verify csi.write_loffset]
/// seqair's CSI and the one `bcftools index -c` builds for the same file must
/// send bcftools to the same records.
///
/// This is the comparison that stands in for byte equality, which htslib's
/// arbitrary bin order rules out. It is also the stricter half of the two:
/// `loffset` and the chunk boundaries are where an index can be subtly
/// conservative or subtly lossy, and only a differential query shows it.
#[hegel::test(test_cases = 25)]
fn seqair_csi_and_bcftools_csi_answer_the_same_queries(tc: TestCase) {
    let layout = tc.draw(arb_layout().print_as_debug());
    let (contig, start, end) = tc.draw(arb_region(layout.clone()).print_as_debug());

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::VcfGz);
    let region = format!("{}:{start}-{end}", layout.contigs[contig].0);

    let with_seqair = query(&path, &region);

    let theirs = dir.path().join("bcftools.csi");
    let out = std::process::Command::new("bcftools")
        .args(["index", "-c", "-f"])
        .arg(&path)
        .arg("-o")
        .arg(&theirs)
        .output()
        .expect("bcftools not found");
    assert!(out.status.success(), "bcftools index: {}", String::from_utf8_lossy(&out.stderr));
    std::fs::copy(&theirs, csi_path(&path)).unwrap();

    let with_bcftools = query(&path, &region);
    assert_eq!(with_seqair, with_bcftools, "region {region}");
    assert_eq!(with_seqair, expected(&layout, contig, start, end), "region {region}");
}

// r[verify index_builder.pseudo_bin]
// r[verify csi.pseudo_bin]
// r[verify csi.write_format]
/// The CSI binary layout, checked against the generated record set rather than
/// against another index: the header fields, the number of references listed,
/// and each reference's pseudo-bin tally of mapped and unmapped records.
///
/// `n_unmapped` is always zero here — a VCF record is always placed — so the
/// pseudo-bin's second chunk must read `(records on that contig, 0)`.
#[hegel::test(test_cases = 25)]
fn csi_header_and_pseudo_bin_match_the_record_set(tc: TestCase) {
    let layout = tc.draw(arb_layout().print_as_debug());

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::VcfGz);
    let data = decompress(&csi_path(&path));

    let i32_at = |off: usize| i32::from_le_bytes(data[off..off + 4].try_into().unwrap());
    let u64_at = |off: usize| u64::from_le_bytes(data[off..off + 8].try_into().unwrap());

    assert_eq!(&data[..4], b"CSI\x01", "magic");
    let min_shift = i32_at(4);
    let depth = i32_at(8);
    assert_eq!(min_shift, 14, "min_shift");
    // Not a fixed number any more, and not re-derived with a copy of the
    // builder's formula either: what the rule actually demands is that the
    // depth *reaches* the longest contig. Bins at depth `d` run out at
    // `2^(14 + 3d)`, so that bound must clear the longest reference, and the
    // depth must never drop below BAI's 5.
    let depth = u32::try_from(depth).expect("depth is small and positive");
    assert!(depth >= 5, "depth {depth} is shallower than BAI's");
    let longest = u64::from(layout.contigs.iter().map(|(_, len)| *len).max().unwrap());
    let reach = 1u64 << (14 + 3 * depth);
    assert!(
        reach > longest,
        "depth {depth} reaches {reach}, short of the longest contig {longest}"
    );
    // The deep band is the one this file exists for, so make it visible that
    // the generator gets there rather than assuming it does.
    tc.event_value("csi depth", f64::from(depth));
    let l_aux = usize::try_from(i32_at(12)).unwrap();

    // A tabix aux block: format=2 (VCF), col_seq=1, col_beg=2, col_end=0,
    // meta='#', skip=0, then the length-prefixed NUL-separated names.
    assert_eq!(
        &data[16..16 + 24],
        &[2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 35, 0, 0, 0, 0, 0, 0, 0],
        "tabix aux header"
    );

    let mut with_records: Vec<usize> = layout.records.iter().map(|r| r.contig).collect();
    with_records.sort_unstable();
    with_records.dedup();

    let names: Vec<String> = with_records.iter().map(|&c| layout.contigs[c].0.clone()).collect();
    let mut want_names = Vec::new();
    for name in &names {
        want_names.extend_from_slice(name.as_bytes());
        want_names.push(0);
    }
    assert_eq!(i32_at(40), i32::try_from(want_names.len()).unwrap(), "l_nm");
    assert_eq!(&data[44..44 + want_names.len()], &want_names[..], "sequence names");

    let mut off = 16 + l_aux;
    assert_eq!(
        i32_at(off),
        i32::try_from(names.len()).unwrap(),
        "n_ref counts only contigs with records"
    );
    off += 4;

    // One past the last real bin: ((1 << 3(depth+1)) - 1) / 7. At depth 5 that
    // is 37450, the number BAI uses.
    let pseudo_bin: u32 = ((1u32 << (3 * (depth + 1))) - 1) / 7 + 1;
    for &contig in &with_records {
        let n_bin = i32_at(off);
        off += 4;
        let mut seen_pseudo = false;
        for _ in 0..n_bin {
            let bin = u32::from_le_bytes(data[off..off + 4].try_into().unwrap());
            off += 4;
            let _loffset = u64_at(off);
            off += 8;
            let n_chunk = i32_at(off);
            off += 4;
            if bin == pseudo_bin {
                assert_eq!(n_chunk, 2, "pseudo-bin holds two chunks");
                let n_mapped = u64_at(off + 16);
                let n_unmapped = u64_at(off + 24);
                let want = layout.records.iter().filter(|r| r.contig == contig).count();
                assert_eq!(
                    n_mapped,
                    u64::try_from(want).unwrap(),
                    "n_mapped on {}",
                    layout.contigs[contig].0
                );
                assert_eq!(n_unmapped, 0, "a VCF record is always placed");
                seen_pseudo = true;
            }
            off += usize::try_from(n_chunk).unwrap() * 16;
        }
        assert!(seen_pseudo, "no pseudo-bin for {}", layout.contigs[contig].0);
    }
    assert_eq!(off, data.len(), "the whole index was accounted for");
}

// r[verify index_builder.csi_depth]
/// The case CSI exists for: a contig past `2^29`.
///
/// `Writer` used to build its index with `min_shift=14, depth=5` whatever the
/// header said, so bins ran out at 512 Mbp. A record above that was written to
/// the file and pushed to the index, and then no region query could reach it —
/// while the same query answered through the CSI `bcftools index -c` builds
/// returned it. The record was never lost; the index was. Nothing errored.
///
/// This is the acceptance test for `r[index_builder.csi_depth]`, kept as a
/// fixture with hard-coded coordinates because the boundary is the point: a
/// generated case that happens to land below `2^29` proves nothing here.
#[test]
fn a_contig_past_the_depth_5_bin_limit_is_still_reachable() {
    let layout = Layout {
        contigs: vec![("big".to_owned(), 600_000_000)],
        records: vec![
            Rec { contig: 0, pos: 1, ref_len: 1, deleted: Vec::new() },
            Rec { contig: 0, pos: BIN_LIMIT + 1, ref_len: 1, deleted: Vec::new() },
        ],
    };

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::VcfGz);
    let region = format!("big:{}-600000000", BIN_LIMIT);

    // Sanity: the record really is in the file.
    let all = std::process::Command::new("bcftools")
        .args(["view", "-H"])
        .arg(&path)
        .output()
        .expect("bcftools not found");
    assert_eq!(String::from_utf8(all.stdout).unwrap().lines().count(), 2, "both records written");

    assert_eq!(
        query(&path, &region),
        vec![("big".to_owned(), BIN_LIMIT + 1)],
        "a record above 2^29 must still be reachable through seqair's CSI"
    );
}

// ── The tabix aux block, read back by tabix ─────────────────────────────

/// `tabix -l <file>` — the sequence names, listed straight out of the index's
/// aux block.
fn tabix_list(path: &Path) -> Vec<String> {
    let out =
        std::process::Command::new("tabix").arg("-l").arg(path).output().expect("tabix not found");
    assert!(out.status.success(), "tabix -l failed: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty())
        .map(ToOwned::to_owned)
        .collect()
}

/// `tabix <file> <region>`, as `(CHROM, POS)` pairs in file order.
fn tabix_query(path: &Path, region: &str) -> Vec<(String, u32)> {
    let out = std::process::Command::new("tabix")
        .arg(path)
        .arg(region)
        .output()
        .expect("tabix not found");
    assert!(
        out.status.success(),
        "tabix {region} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let mut f = l.split('\t');
            let chrom = f.next().unwrap().to_owned();
            let pos = f.next().unwrap().parse().unwrap();
            (chrom, pos)
        })
        .collect()
}

// r[verify csi.write_tabix_aux]
/// The aux block, read back by the tool it exists for.
///
/// Everything else here checks the aux block by looking at the bytes seqair
/// wrote, which cannot catch a layout that is self-consistent and wrong. This
/// hands the index to `tabix` instead, which is where the format's name comes
/// from and which reaches it through htslib's own reader rather than ours.
///
/// `tabix -l` prints the sequence-name dictionary directly out of the aux
/// block, so it is the read-back the rule asks for; the region queries then
/// need the column configuration (`format=2`, `col_seq=1`, `col_beg=2`) to be
/// right as well, since tabix uses it to decide which fields of a line are the
/// coordinates. A name dictionary in the wrong order still lists the right
/// names, and a wrong `col_beg` still lists them too — only the queries catch
/// those.
#[hegel::test(test_cases = 25)]
fn tabix_reads_the_regions_through_seqairs_aux_block(tc: TestCase) {
    let layout = tc.draw(arb_layout().print_as_debug());

    let dir = tempfile::tempdir().unwrap();
    let path = write_indexed(dir.path(), &layout, OutputFormat::VcfGz);

    let mut with_records: Vec<usize> = layout.records.iter().map(|r| r.contig).collect();
    with_records.sort_unstable();
    with_records.dedup();
    let want_names: Vec<String> =
        with_records.iter().map(|&c| layout.contigs[c].0.clone()).collect();
    assert_eq!(tabix_list(&path), want_names, "tabix -l reads the aux block's name dictionary");

    // One whole-contig query per named contig: the records tabix returns must
    // be exactly the ones generated for it, in position order.
    for (&contig, name) in with_records.iter().zip(&want_names) {
        let len = layout.contigs[contig].1;
        let region = format!("{name}:1-{len}");
        let mut want: Vec<(String, u32)> = layout
            .records
            .iter()
            .filter(|r| r.contig == contig)
            .map(|r| (name.clone(), r.pos))
            .collect();
        want.sort_by_key(|&(_, pos)| pos);
        assert_eq!(tabix_query(&path, &region), want, "tabix {region}");
    }
    tc.event_value("named contigs", with_records.len() as f64);
}
