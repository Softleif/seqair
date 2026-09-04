//! Shared lazy generation of benchmark inputs.
//!
//! Included via `#[path]` by the `fasta` and `fastq` bench targets. Everything
//! lands outside the repo, never in `tests/data/`, because the inputs run to
//! hundreds of megabytes. Generation is lazy and cached: each function is a
//! no-op once its output exists.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Where generated inputs live. Override with `SEQAIR_BENCH_DATA`.
///
/// Defaults to a stable directory under the system temp dir so the (large,
/// slow-to-build) inputs survive across runs and machines.
pub fn scratch_root() -> PathBuf {
    std::env::var_os("SEQAIR_BENCH_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("seqair-benchdata"))
}

pub fn scratch(name: &str) -> PathBuf {
    scratch_root().join(name)
}

/// xorshift64* — deterministic sequence generation without a dev-dependency.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    pub fn base(&mut self) -> u8 {
        const BASES: [u8; 4] = *b"ACGT";
        BASES.get((self.next_u64() >> 33) as usize % 4).copied().unwrap_or(b'N')
    }

    /// A base at ~41% GC, the human genome-wide figure, rather than 50%.
    pub fn gc_base(&mut self) -> u8 {
        match (self.next_u64() >> 33) % 200 {
            0..=40 => b'C',
            41..=81 => b'G',
            82..=140 => b'A',
            _ => b'T',
        }
    }

    /// Quality from the 4-bin scheme current Illumina instruments emit, with
    /// the low bins growing toward the 3' end of the read. Returns the ASCII
    /// (Phred+33) character.
    pub fn illumina_qual(&mut self, pos: usize, read_len: usize) -> u8 {
        // 0 at the 5' end, 100 at the 3' end.
        let decay = pos.saturating_mul(100).checked_div(read_len.max(1)).unwrap_or(0) as u64;
        // Chance of dropping out of the top bin rises along the read.
        let roll = (self.next_u64() >> 40) % 1000;
        let threshold = 20 + decay.saturating_mul(2);
        33 + if roll >= threshold {
            37
        } else if roll >= threshold / 4 {
            23
        } else if roll >= threshold / 16 {
            12
        } else {
            2
        }
    }
}

fn samtools_faidx(path: &Path) {
    let status = Command::new("samtools")
        .arg("faidx")
        .arg(path)
        .status()
        .expect("samtools must be on PATH to generate bench indexes");
    assert!(status.success(), "samtools faidx failed for {}", path.display());
}

/// Plain uncompressed copy of the BGZF `test.fasta.gz`, with its `.fai`.
///
/// The committed `tests/data/test.fasta.fai` already describes this exact
/// layout (offsets are identical — BGZF only changes the container), so it is
/// copied rather than rebuilt.
pub fn plain_fasta() -> PathBuf {
    let out = scratch("test.fasta");
    if out.exists() {
        return out;
    }
    std::fs::create_dir_all(scratch_root()).unwrap();
    let src = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");
    let dest = std::fs::File::create(&out).unwrap();
    let status = Command::new("bgzip")
        .arg("-dc")
        .arg(src)
        .stdout(dest)
        .status()
        .expect("bgzip must be on PATH");
    assert!(status.success(), "bgzip decompress failed");
    std::fs::copy(
        concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.fai"),
        scratch("test.fasta.fai"),
    )
    .unwrap();
    out
}

/// Many-contig FASTA: 3000 contigs of 5 kb, 60-base lines.
///
/// Exposes `noodles_fasta::fai::Index::query`, which is a linear `iter().find()`
/// over records, against seqair's `FxHashMap` name lookup. Queries target the
/// last contig, the worst case for the scan.
pub fn many_contig_fasta() -> PathBuf {
    let out = scratch("many_contig.fasta");
    if out.exists() {
        return out;
    }
    std::fs::create_dir_all(scratch_root()).unwrap();
    let mut rng = Rng::new(0x5EED);
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
    let mut line = Vec::with_capacity(60);
    for c in 0..3000u32 {
        writeln!(w, ">contig_{c}").unwrap();
        for _ in 0..(5000 / 60) {
            line.clear();
            for _ in 0..60 {
                line.push(rng.base());
            }
            w.write_all(&line).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }
    w.into_inner().unwrap();
    samtools_faidx(&out);
    out
}

/// Reference-like FASTQ: 3 long contigs (20 Mb / 2 Mb / 0.5 Mb), 60-base lines.
///
/// Indexed by `samtools faidx` into a 6-column FAI. A sibling symlink with a
/// hand-written 5-column FAI is produced by [`ref_fastq_seqair_view`].
pub fn ref_fastq() -> PathBuf {
    let out = scratch("bench_ref.fq");
    if out.exists() {
        return out;
    }
    std::fs::create_dir_all(scratch_root()).unwrap();
    let mut rng = Rng::new(0xC0FFEE);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&out).unwrap());
    for (name, len) in [("ctg_a", 20_000_000usize), ("ctg_b", 2_000_000), ("ctg_c", 500_000)] {
        let seq: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        writeln!(w, "@{name}").unwrap();
        for chunk in seq.chunks(60) {
            w.write_all(chunk).unwrap();
            w.write_all(b"\n").unwrap();
        }
        writeln!(w, "+").unwrap();
        // Quality lines must mirror the sequence line structure for faidx.
        let qual: Vec<u8> = (0..len).map(|_| b'!' + (rng.next_u64() >> 40) as u8 % 40).collect();
        for chunk in qual.chunks(60) {
            w.write_all(chunk).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }
    w.into_inner().unwrap();
    samtools_faidx(&out);
    out
}

/// A symlink to [`ref_fastq`] carrying a 5-column FAI. See [`five_column_view`].
pub fn ref_fastq_seqair_view() -> PathBuf {
    let real = ref_fastq();
    five_column_view(&real, "bench_ref_seqair.fq")
}

/// Read-like FASTQ: 200k records of 150 bp, unwrapped. ~68 MB plain.
///
/// Returned as `(plain, gzip, bgzf)`.
///
/// Deliberately not uniform random. Uniform ACGT with flat quality is the
/// worst input a compression benchmark can use: it is close to incompressible,
/// so gzip ratios land near 2x where real FASTQ reaches 4-5x, and the whole
/// compressed-input picture shifts. Instead:
///
/// * Reads are 150 bp substrings of a 5 Mb synthetic genome, half of them
///   reverse-complemented, with a ~0.3% per-base error rate. That gives the
///   redundancy between records that real short-read data has.
/// * Quality uses the 4-bin scheme current Illumina instruments emit (Q2/Q12/
///   Q23/Q37), with the low bins weighted toward the 3' end. Real quality is
///   binned and highly skewed, which is most of why FASTQ compresses at all.
/// * Names follow the Illumina convention with varying tile and x/y fields.
pub fn reads_fastq() -> (PathBuf, PathBuf, PathBuf) {
    let plain = scratch("reads_150bp.fq");
    let gz = scratch("reads_150bp.fq.gz");
    let bgz = scratch("reads_150bp.fq.bgz");
    if plain.exists() && gz.exists() && bgz.exists() {
        return (plain, gz, bgz);
    }
    std::fs::create_dir_all(scratch_root()).unwrap();
    let mut rng = Rng::new(0xBEEF);

    const GENOME: usize = 5_000_000;
    const READ_LEN: usize = 150;
    const N_READS: u32 = 200_000;

    let genome: Vec<u8> = (0..GENOME).map(|_| rng.gc_base()).collect();

    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&plain).unwrap());
    let mut seq = [0u8; READ_LEN];
    let mut qual = [0u8; READ_LEN];
    for i in 0..N_READS {
        let start = (rng.next_u64() >> 20) as usize % (GENOME - READ_LEN);
        let revcomp = rng.next_u64() >> 63 == 1;
        for (j, b) in seq.iter_mut().enumerate() {
            let src = if revcomp { start + READ_LEN - 1 - j } else { start + j };
            let g = genome.get(src).copied().unwrap_or(b'N');
            let g = if revcomp { complement(g) } else { g };
            // ~0.3% sequencing error.
            *b = if (rng.next_u64() >> 40) % 1000 < 3 { rng.gc_base() } else { g };
        }
        for (j, q) in qual.iter_mut().enumerate() {
            *q = rng.illumina_qual(j, READ_LEN);
        }
        let tile = 1101 + (rng.next_u64() >> 40) % 78;
        let x = (rng.next_u64() >> 40) % 30_000;
        let y = (rng.next_u64() >> 40) % 30_000;
        writeln!(w, "@A00519:247:HG7YTDSXX:1:{tile}:{x}:{y} 1:N:0:ATTACTCG+GGCTCTGA").unwrap();
        w.write_all(&seq).unwrap();
        writeln!(w, "\n+").unwrap();
        w.write_all(&qual).unwrap();
        w.write_all(b"\n").unwrap();
        let _ = i;
    }
    w.into_inner().unwrap();
    for (dst, prog) in [(&gz, "gzip"), (&bgz, "bgzip")] {
        let dest = std::fs::File::create(dst).unwrap();
        let status = Command::new(prog)
            .arg("-c")
            .arg(&plain)
            .stdout(dest)
            .status()
            .expect("gzip/bgzip must be on PATH");
        assert!(status.success(), "{prog} failed");
    }
    (plain, gz, bgz)
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Build a 5-column view of a 6-column FASTQ index.
///
/// Produces `<stem>` (a symlink to `real`) plus `<stem>.fai` carrying only the
/// FASTA columns. seqair's FAI parser rejects the 6-column FASTQ form, but its
/// byte-span math is already correct for FASTQ — the sequence block is bounded
/// by `length`, so the reader never walks into the `+` or quality lines.
/// Dropping `qual_offset` therefore yields exactly the index a 6-column-aware
/// parser would build, measuring the performance seqair *would* deliver once
/// the parser accepts the extra column, without touching `src/`.
pub fn five_column_view(real: &Path, stem: &str) -> PathBuf {
    let out = scratch(stem);
    let fai = scratch(&format!("{stem}.fai"));
    if !out.exists() {
        // A hard link, not a symlink: both paths sit in the same scratch
        // directory, and `std::fs::hard_link` is portable where
        // `std::os::unix::fs::symlink` would leave this bench target unable
        // to compile on Windows.
        std::fs::hard_link(real, &out).unwrap();
    }
    if !fai.exists() {
        let mut six_path = real.as_os_str().to_owned();
        six_path.push(".fai");
        let six = std::fs::read_to_string(PathBuf::from(six_path)).unwrap();
        let mut five = String::new();
        for l in six.lines() {
            // name, length, seq_offset, linebases, linewidth — drop qual_offset.
            let f: Vec<&str> = l.split('\t').take(5).collect();
            assert_eq!(f.len(), 5, "malformed FAI line: {l:?}");
            five.push_str(&f.join("\t"));
            five.push('\n');
        }
        std::fs::write(&fai, five).unwrap();
    }
    out
}

/// Real reference data as FASTQ: `tests/data/test.fasta.gz` (chr19 + two spike-ins)
/// re-emitted with synthetic quality lines, preserving the source line structure
/// (50 bases/line). ~125 MB.
///
/// The synthetic [`ref_fastq`] is uniform random ACGT; this one carries real
/// base composition including the long `N` runs of the chr19 centromere, which
/// is what a caller would actually hand to `Readers::open`.
pub fn chr19_fastq() -> PathBuf {
    let out = scratch("chr19.fq");
    if out.exists() {
        return out;
    }
    let src = plain_fasta();
    let reader = std::io::BufReader::with_capacity(1 << 20, std::fs::File::open(src).unwrap());
    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&out).unwrap());
    let mut rng = Rng::new(0xFA57_0FEE);
    // Line lengths of the current record, so the quality block can mirror them.
    let mut line_lens: Vec<usize> = Vec::new();
    let mut qual_line = Vec::new();

    let flush = |w: &mut std::io::BufWriter<std::fs::File>,
                 line_lens: &mut Vec<usize>,
                 qual_line: &mut Vec<u8>,
                 rng: &mut Rng| {
        if line_lens.is_empty() {
            return;
        }
        w.write_all(b"+\n").unwrap();
        for len in line_lens.iter() {
            qual_line.clear();
            for _ in 0..*len {
                qual_line.push(b'!' + (rng.next_u64() >> 40) as u8 % 40);
            }
            w.write_all(qual_line).unwrap();
            w.write_all(b"\n").unwrap();
        }
        line_lens.clear();
    };

    for line in std::io::BufRead::lines(reader) {
        let line = line.unwrap();
        if let Some(name) = line.strip_prefix('>') {
            flush(&mut w, &mut line_lens, &mut qual_line, &mut rng);
            writeln!(w, "@{name}").unwrap();
        } else {
            line_lens.push(line.len());
            w.write_all(line.as_bytes()).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }
    flush(&mut w, &mut line_lens, &mut qual_line, &mut rng);
    w.into_inner().unwrap();
    samtools_faidx(&out);
    out
}
