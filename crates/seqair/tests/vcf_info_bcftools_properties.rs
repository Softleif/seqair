//! Generated INFO fields, written by seqair, read back by bcftools.
//!
//! The fixtures in `vcf_bcftools_roundtrip.rs` show bcftools thirteen
//! hand-written records. That is enough to catch a broken magic number and not
//! much else: every INFO type, every `Number` cardinality, and the way `A`/`R`/
//! `G` widths follow the ALT count are all pinned by a single example apiece.
//!
//! Here the records are generated and bcftools referees. The value of that is
//! the combinations — a `Number=R` array on a four-allele site, a `Character`
//! field next to a `Flag`, an `int8` array whose width is decided by one value
//! in the middle of it. Those are where this writer's bugs have been.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

use hegel::prelude::*;
use seqair::vcf::record_encoder::{
    InfoFieldDef, InfoFlag, InfoFloat, InfoFloats, InfoInt, InfoInts, InfoString,
};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;
use std::sync::Arc;

/// Every contig position a generated record may take. Small enough that a whole
/// record set lands in one BGZF block, large enough that positions are distinct.
const CONTIG_LEN: u64 = 250_000_000;
const MAX_POS: u32 = 1_000_000;

struct Setup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    /// `Number=1`, one per scalar type the VCF spec defines.
    int1: InfoInt,
    float1: InfoFloat,
    flag: InfoFlag,
    char1: InfoString,
    string1: InfoString,
    /// `Number=3` — a fixed count above one.
    fixed3: InfoInts,
    /// The cardinalities that are derived from the ALT count.
    per_alt: InfoInts,
    per_allele: InfoInts,
    per_genotype: InfoInts,
    unknown: InfoInts,
    per_alt_float: InfoFloats,
}

fn setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(CONTIG_LEN) }).unwrap();
    let mut b = builder.infos();
    let int1 = b
        .register_info(&InfoFieldDef::new("II", Number::Count(1), ValueType::Integer, "int"))
        .unwrap();
    let float1 = b
        .register_info(&InfoFieldDef::new("FF", Number::Count(1), ValueType::Float, "float"))
        .unwrap();
    let flag = b
        .register_info(&InfoFieldDef::new("FG", Number::Count(0), ValueType::Flag, "flag"))
        .unwrap();
    let char1 = b
        .register_info(&InfoFieldDef::new("CH", Number::Count(1), ValueType::Character, "char"))
        .unwrap();
    let string1 = b
        .register_info(&InfoFieldDef::new("ST", Number::Count(1), ValueType::String, "string"))
        .unwrap();
    let fixed3 = b
        .register_info(&InfoFieldDef::new("P3", Number::Count(3), ValueType::Integer, "three"))
        .unwrap();
    let per_alt = b
        .register_info(&InfoFieldDef::new(
            "AI",
            Number::AlternateBases,
            ValueType::Integer,
            "one per ALT",
        ))
        .unwrap();
    let per_allele = b
        .register_info(&InfoFieldDef::new(
            "RI",
            Number::ReferenceAlternateBases,
            ValueType::Integer,
            "one per allele",
        ))
        .unwrap();
    let per_genotype = b
        .register_info(&InfoFieldDef::new(
            "GI",
            Number::Genotypes,
            ValueType::Integer,
            "one per genotype",
        ))
        .unwrap();
    let unknown = b
        .register_info(&InfoFieldDef::new("UI", Number::Unknown, ValueType::Integer, "any count"))
        .unwrap();
    let per_alt_float = b
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "one per ALT",
        ))
        .unwrap();
    let header = Arc::new(b.samples().build().unwrap());
    Setup {
        header,
        contig,
        int1,
        float1,
        flag,
        char1,
        string1,
        fixed3,
        per_alt,
        per_allele,
        per_genotype,
        unknown,
        per_alt_float,
    }
}

// ── Generated records ──────────────────────────────────────────────────

const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

/// What a record's alleles are, kept as the *inputs* rather than as an
/// [`Alleles`]: the expected REF/ALT text is then built from the draw and not
/// from the same accessor the writer uses.
#[derive(Debug, Clone)]
enum AlleleSpec {
    /// REF only, no ALT (gVCF reference block).
    Reference(Base),
    /// One REF base, 1..=3 distinct ALT bases.
    Snv(Base, Vec<Base>),
    Insertion(Base, Vec<Base>),
    Deletion(Base, Vec<Base>),
    Complex(String, Vec<String>),
}

impl AlleleSpec {
    fn to_alleles(&self) -> Alleles {
        match self {
            Self::Reference(r) => Alleles::reference(*r),
            Self::Snv(r, alts) => Alleles::snv_multi(*r, alts).unwrap(),
            Self::Insertion(a, ins) => Alleles::insertion(*a, ins).unwrap(),
            Self::Deletion(a, del) => Alleles::deletion(*a, del).unwrap(),
            Self::Complex(r, alts) => Alleles::complex(
                r.as_str().into(),
                alts.iter().map(|a| a.as_str().into()).collect(),
            ),
        }
    }

    /// The REF column, spelled out from the draw rather than read back out of
    /// the [`Alleles`] the writer was handed.
    fn ref_text(&self) -> String {
        let s = |b: &Base| b.as_str().to_owned();
        match self {
            Self::Reference(r) | Self::Snv(r, _) | Self::Insertion(r, _) => s(r),
            Self::Deletion(anchor, deleted) => {
                let mut t = s(anchor);
                for b in deleted {
                    t.push_str(b.as_str());
                }
                t
            }
            Self::Complex(r, _) => r.clone(),
        }
    }

    /// The ALT column, spelled out from the draw. Empty means `.`.
    fn alt_texts(&self) -> Vec<String> {
        let s = |b: &Base| b.as_str().to_owned();
        match self {
            Self::Reference(_) => Vec::new(),
            Self::Snv(_, alts) => alts.iter().map(s).collect(),
            Self::Insertion(anchor, inserted) => {
                let mut t = s(anchor);
                for b in inserted {
                    t.push_str(b.as_str());
                }
                vec![t]
            }
            Self::Deletion(anchor, _) => vec![s(anchor)],
            Self::Complex(_, alts) => alts.clone(),
        }
    }

    fn n_alt(&self) -> usize {
        self.alt_texts().len()
    }
}

#[hegel::composite]
fn arb_alleles(tc: &TestCase) -> AlleleSpec {
    // Drawn by construction: picking a REF and then a subsequence of the other
    // three bases can never violate `snv_multi`'s distinctness rule, so no
    // generated case is ever rejected.
    let kind = tc.draw_silent(gs::integers::<u8>().max_value(4));
    let ref_base = tc.draw_silent(gs::sampled_from(&ACGT));
    let others: Vec<Base> = ACGT.iter().copied().filter(|b| *b != ref_base).collect();
    match kind {
        0 => AlleleSpec::Reference(ref_base),
        1 => {
            let alts = tc.draw_silent(gs::subsequences(&others).min_size(1));
            AlleleSpec::Snv(ref_base, alts)
        }
        2 => {
            let ins = tc.draw_silent(gs::vecs(gs::sampled_from(&ACGT)).min_size(1).max_size(4));
            AlleleSpec::Insertion(ref_base, ins)
        }
        3 => {
            let del = tc.draw_silent(gs::vecs(gs::sampled_from(&ACGT)).min_size(1).max_size(4));
            AlleleSpec::Deletion(ref_base, del)
        }
        _ => {
            let seq = || gs::from_regex("[ACGT]{1,4}");
            let r: String = tc.draw_silent(seq());
            let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(3));
            // Distinct from REF and from each other, by construction: a unique
            // suffix is appended rather than redrawing until the draw differs.
            let alts = (0..n)
                .map(|i| {
                    let a: String = tc.draw_silent(seq());
                    format!("{a}{}", "A".repeat(i.saturating_add(1)))
                })
                .collect();
            AlleleSpec::Complex(r, alts)
        }
    }
}

/// One record: alleles plus one value for every INFO field in the header.
#[derive(Debug, Clone)]
struct Record {
    pos: u32,
    alleles: AlleleSpec,
    int1: i32,
    float1: f32,
    flag: bool,
    char1: char,
    string1: String,
    fixed3: [i32; 3],
    /// `n_alt` values.
    per_alt: Vec<i32>,
    /// `n_alt + 1` values.
    per_allele: Vec<i32>,
    /// `(n_allele * (n_allele + 1)) / 2` values — diploid genotype count.
    per_genotype: Vec<i32>,
    /// Any count at all, including none.
    unknown: Vec<i32>,
    per_alt_float: Vec<f32>,
}

/// Integers that exercise all three BCF widths, biased towards the boundaries
/// where `r[bcf_writer.smallest_int_type]` switches. The eight most negative
/// values of each width are sentinels, so the ranges stop short of them.
fn arb_int() -> impl PrintableGenerator<i32> {
    hegel::one_of!(
        gs::integers::<i32>().min_value(-120).max_value(127),
        gs::integers::<i32>().min_value(-32_760).max_value(32_767),
        gs::integers::<i32>().min_value(i32::MIN + 8).max_value(i32::MAX),
        gs::sampled_from(vec![
            -120i32,
            127,
            128,
            -121,
            32_767,
            -32_760,
            32_768,
            -32_761,
            0,
            i32::MAX,
            i32::MIN + 8,
        ]),
    )
}

/// Finite floats only: bcftools renders a non-finite float as `nan`/`inf`,
/// which has its own coverage in the writer's unit tests, and parsing it back
/// here would only re-test the printer.
fn arb_float() -> impl PrintableGenerator<f32> {
    gs::floats::<f32>()
        .min_value(-100_000.0)
        .max_value(100_000.0)
        .allow_nan(false)
        .allow_infinity(false)
}

#[hegel::composite]
fn arb_records(tc: &TestCase) -> Vec<Record> {
    let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(4));
    let mut positions: Vec<u32> = (0..n)
        .map(|_| tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(MAX_POS)))
        .collect();
    positions.sort_unstable();

    positions
        .into_iter()
        .map(|pos| {
            let alleles = tc.draw_silent(arb_alleles());
            let n_alt = alleles.n_alt();
            let n_allele = n_alt.saturating_add(1);
            let n_genotype = n_allele.saturating_mul(n_allele.saturating_add(1)) / 2;
            let ints = |k: usize| (0..k).map(|_| tc.draw_silent(arb_int())).collect::<Vec<_>>();
            Record {
                pos,
                int1: tc.draw_silent(arb_int()),
                float1: tc.draw_silent(arb_float()),
                flag: tc.draw_silent(gs::booleans()),
                char1: tc.draw_silent(gs::sampled_from(vec!['A', 'c', '9', '_', '+', '-'])),
                string1: tc.draw_silent(gs::from_regex("[A-Za-z0-9_.+-]{1,8}")),
                fixed3: [
                    tc.draw_silent(arb_int()),
                    tc.draw_silent(arb_int()),
                    tc.draw_silent(arb_int()),
                ],
                per_alt: ints(n_alt),
                per_allele: ints(n_allele),
                per_genotype: ints(n_genotype),
                unknown: tc.draw_silent(gs::vecs(arb_int()).max_size(5)),
                per_alt_float: (0..n_alt).map(|_| tc.draw_silent(arb_float())).collect(),
                alleles,
            }
        })
        .collect()
}

// ── Writing ────────────────────────────────────────────────────────────

fn write(s: &Setup, records: &[Record], format: OutputFormat) -> Vec<u8> {
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, format);
    let mut writer = writer.write_header(&s.header).unwrap();
    for rec in records {
        let alleles = rec.alleles.to_alleles();
        let mut enc = writer
            .begin_record(&s.contig, Pos1::new(rec.pos).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        s.int1.encode(&mut enc, rec.int1);
        s.float1.encode(&mut enc, rec.float1);
        if rec.flag {
            s.flag.encode(&mut enc);
        }
        s.char1.encode(&mut enc, rec.char1.to_string().as_str());
        s.string1.encode(&mut enc, &rec.string1);
        s.fixed3.encode(&mut enc, &rec.fixed3);
        // An `A`/`R`/`G` field with no values has no legal VCF spelling, so a
        // reference-only site simply omits the per-ALT ones.
        if !rec.per_alt.is_empty() {
            s.per_alt.encode(&mut enc, &rec.per_alt);
            s.per_alt_float.encode(&mut enc, &rec.per_alt_float);
        }
        s.per_allele.encode(&mut enc, &rec.per_allele);
        s.per_genotype.encode(&mut enc, &rec.per_genotype);
        if !rec.unknown.is_empty() {
            s.unknown.encode(&mut enc, &rec.unknown);
        }
        enc.emit().unwrap();
    }
    writer.finish().unwrap();
    out
}

/// `bcftools query -f <fmt>` over BCF bytes, one line per record.
fn bcftools_query(bcf: &[u8], format: &str) -> Vec<String> {
    let mut file = tempfile::NamedTempFile::with_suffix(".bcf").unwrap();
    std::io::Write::write_all(&mut file, bcf).unwrap();
    std::io::Write::flush(&mut file).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", format])
        .arg(file.path())
        .output()
        .expect("bcftools not found");
    assert!(
        out.status.success(),
        "bcftools query failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap().lines().map(str::to_owned).collect()
}

fn join_ints(values: &[i32]) -> String {
    values.iter().map(i32::to_string).collect::<Vec<_>>().join(",")
}

/// Compare bcftools' rendering of a float list against the drawn values. BCF
/// stores `f32` exactly; bcftools prints six significant digits.
fn assert_floats_match(printed: &str, want: &[f32], what: &str) {
    let got: Vec<f32> =
        printed.split(',').map(|v| v.parse().expect("bcftools prints parseable floats")).collect();
    assert_eq!(got.len(), want.len(), "{what}: arity");
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() <= w.abs() * 1e-5, "{what}: {g} != {w}");
    }
}

// ── Properties ─────────────────────────────────────────────────────────

// r[verify record_encoder.key_encode]
// r[verify bcf_writer.typed_values]
// r[verify bcf_writer.string_encoding]
// r[verify bcf_writer.flag_encoding]
/// Every scalar INFO type the VCF spec defines — Integer, Float, Flag,
/// Character, String — comes back from bcftools as the value that went in.
#[hegel::test(test_cases = 50)]
fn info_scalar_types_round_trip_through_bcftools(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());
    let s = setup();
    let bcf = write(&s, &records, OutputFormat::Bcf);
    let lines = bcftools_query(&bcf, "%INFO/II\\t%INFO/FF\\t%INFO/FG\\t%INFO/CH\\t%INFO/ST\\n");

    assert_eq!(lines.len(), records.len(), "record count");
    for (line, rec) in lines.iter().zip(&records) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 5, "field count in {line:?}");
        assert_eq!(f[0], rec.int1.to_string(), "II at pos {}", rec.pos);
        assert_floats_match(f[1], &[rec.float1], &format!("FF at pos {}", rec.pos));
        // A Flag prints as `1` when set; bcftools writes `.` when it is absent.
        assert_eq!(f[2], if rec.flag { "1" } else { "." }, "FG at pos {}", rec.pos);
        assert_eq!(f[3], rec.char1.to_string(), "CH at pos {}", rec.pos);
        assert_eq!(f[4], rec.string1, "ST at pos {}", rec.pos);
    }

    tc.event(if records.iter().any(|r| r.flag) { "a flag is set" } else { "no flag set" });
}

// r[verify record_encoder.key_encode]
// r[verify bcf_writer.shared_variable]
// r[verify vcf_header.info_def]
/// The `Number` cardinalities that are derived from the ALT count — `A`, `R`,
/// `G` — keep exactly the width the record's allele count implies, and a fixed
/// `Number=n` and an unbounded `Number=.` keep theirs, all on the same record.
#[hegel::test(test_cases = 50)]
fn info_number_cardinalities_survive_bcftools(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());
    let s = setup();
    let bcf = write(&s, &records, OutputFormat::Bcf);
    let lines =
        bcftools_query(&bcf, "%INFO/P3\\t%INFO/AI\\t%INFO/RI\\t%INFO/GI\\t%INFO/UI\\t%INFO/AF\\n");

    assert_eq!(lines.len(), records.len(), "record count");
    for (line, rec) in lines.iter().zip(&records) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 6, "field count in {line:?}");
        let at = |what: &str| format!("{what} at pos {} ({} ALTs)", rec.pos, rec.alleles.n_alt());
        assert_eq!(f[0], join_ints(&rec.fixed3), "{}", at("P3"));
        if rec.per_alt.is_empty() {
            assert_eq!(f[1], ".", "{}", at("AI"));
            assert_eq!(f[5], ".", "{}", at("AF"));
        } else {
            assert_eq!(f[1], join_ints(&rec.per_alt), "{}", at("AI"));
            assert_floats_match(f[5], &rec.per_alt_float, &at("AF"));
        }
        assert_eq!(f[2], join_ints(&rec.per_allele), "{}", at("RI"));
        assert_eq!(f[3], join_ints(&rec.per_genotype), "{}", at("GI"));
        assert_eq!(
            f[4],
            if rec.unknown.is_empty() { ".".into() } else { join_ints(&rec.unknown) },
            "{}",
            at("UI")
        );
    }

    let max_alt = records.iter().map(|r| r.alleles.n_alt()).max().unwrap_or(0);
    tc.event(match max_alt {
        0 => "reference-only site",
        1 => "bi-allelic only",
        _ => "multi-allelic",
    });
}

// r[verify vcf_record.alleles_serialization]
// r[verify bcf_writer.shared_variable]
/// Every allele kind the type-safe `Alleles` API can build — a reference-only
/// site, a 1..3-ALT SNV, an insertion, a deletion and a complex site — spells
/// the REF and ALT columns that the draw implies, as bcftools reads them back
/// out of the BCF.
///
/// The expectation is built from the drawn bases, not from `Alleles::ref_text`,
/// so this is not the writer checked against its own accessor.
#[hegel::test(test_cases = 50)]
fn allele_kinds_spell_ref_and_alt_for_bcftools(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());
    let s = setup();
    let bcf = write(&s, &records, OutputFormat::Bcf);
    let lines = bcftools_query(&bcf, "%CHROM\\t%POS\\t%REF\\t%ALT\\n");

    assert_eq!(lines.len(), records.len(), "record count");
    for (line, rec) in lines.iter().zip(&records) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 4, "field count in {line:?}");
        assert_eq!(f[0], "chr1", "CHROM");
        assert_eq!(f[1], rec.pos.to_string(), "POS");
        assert_eq!(f[2], rec.alleles.ref_text(), "REF at pos {}", rec.pos);
        let alts = rec.alleles.alt_texts();
        let want = if alts.is_empty() { ".".to_owned() } else { alts.join(",") };
        assert_eq!(f[3], want, "ALT at pos {}", rec.pos);
    }

    let kinds = records
        .iter()
        .map(|r| match r.alleles {
            AlleleSpec::Reference(_) => "reference-only",
            AlleleSpec::Snv(_, ref a) if a.len() > 1 => "multi-allelic SNV",
            AlleleSpec::Snv(..) => "bi-allelic SNV",
            AlleleSpec::Insertion(..) => "insertion",
            AlleleSpec::Deletion(..) => "deletion",
            AlleleSpec::Complex(..) => "complex",
        })
        .collect::<std::collections::BTreeSet<_>>();
    for kind in kinds {
        tc.event(kind);
    }
}

// r[verify record_encoder.vcf_bcf_equivalence]
// r[verify vcf_writer.info_serialization]
/// seqair writes each record twice over, through two encoders that share no
/// serialization code: straight to VCF text, and to BCF that bcftools then
/// renders back as VCF text. Every column of the two renderings must agree.
///
/// This is the three-way agreement the writer's contract asks for — if the two
/// seqair paths can disagree about the same logical record, one of them is
/// wrong, and bcftools says which.
#[hegel::test(test_cases = 50)]
fn vcf_text_and_bcf_render_the_same_record(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());
    let s = setup();

    let direct: Vec<String> = String::from_utf8(write(&s, &records, OutputFormat::Vcf))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_owned)
        .collect();

    let bcf = write(&s, &records, OutputFormat::Bcf);
    let via_bcf = bcftools_query(&bcf, "%LINE\\n");

    assert_eq!(direct.len(), records.len(), "seqair wrote a record per input");
    assert_eq!(direct, via_bcf, "VCF text and BCF disagree");
}
