//! The writer checks an array's length against its declared `Number`, and a
//! rejected field leaves the record emittable: bcftools reads it back with the
//! field absent and everything else as written.
//!
//! Records are generated with arrays of every length against `A`, `R`, `G`,
//! a fixed count and `.`, on sites of one to three ALTs, with and without a GT
//! of ploidy one to three. The oracle counts are written out per ploidy rather
//! than through the writer's binomial.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
use hegel::generators as gs;
use hegel::prelude::*;
use seqair::vcf::record_encoder::FormatEncoder as _;
use seqair::vcf::{
    Alleles, Arr, ContigDef, ContigId, FormatFieldDef, FormatFloatOpts, FormatFloats, FormatGt,
    FormatInts, Genotype, Gt, InfoFieldDef, InfoFloats, InfoInt, InfoIntOpts, InfoInts, Number,
    OptArr, OutputFormat, Scalar, VcfError, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::collections::BTreeMap;
use std::process::Command;
use std::sync::Arc;

#[derive(Debug, Clone, Copy)]
enum Card {
    A,
    R,
    G,
    Count(u32),
    Any,
}

impl Card {
    fn number(self) -> Number {
        match self {
            Self::A => Number::AlternateBases,
            Self::R => Number::ReferenceAlternateBases,
            Self::G => Number::Genotypes,
            Self::Count(n) => Number::Count(n),
            Self::Any => Number::Unknown,
        }
    }

    /// The count the writer must insist on, `None` where it checks nothing.
    fn expected(self, n_alt: usize, ploidy: Option<usize>) -> Option<usize> {
        let n = n_alt + 1;
        match self {
            Self::A => Some(n_alt),
            Self::R => Some(n),
            Self::G => ploidy.map(|p| match p {
                1 => n,
                2 => n * (n + 1) / 2,
                3 => n * (n + 1) * (n + 2) / 6,
                _ => unreachable!("ploidy is drawn from 1..=3"),
            }),
            Self::Count(c) => Some(usize::try_from(c).unwrap()),
            Self::Any => None,
        }
    }
}

const CARDS: [Card; 5] = [Card::A, Card::R, Card::G, Card::Count(3), Card::Any];
const N_SAMPLES: usize = 2;

struct Setup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    dp: InfoInt,
    info_ints: Vec<(String, Card, InfoInts)>,
    info_opt: (Card, InfoIntOpts),
    info_floats: (Card, InfoFloats),
    gt: FormatGt,
    format_ints: Vec<(String, Card, FormatInts)>,
    format_floats: (Card, FormatFloats),
    format_opt: (Card, FormatFloatOpts),
}

fn setup() -> Setup {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
    let mut b = b.infos();
    let dp = b.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "depth")).unwrap();
    let info_ints = CARDS
        .iter()
        .enumerate()
        .map(|(i, &card)| {
            let name = format!("I{i}");
            let def = InfoFieldDef::<Arr<i32>>::new(leak(&name), card.number(), "x");
            (name, card, b.register_info(&def).unwrap())
        })
        .collect();
    let info_opt = (
        Card::R,
        b.register_info(&InfoFieldDef::<OptArr<i32>>::new("IO", Card::R.number(), "x")).unwrap(),
    );
    let info_floats = (
        Card::A,
        b.register_info(&InfoFieldDef::<Arr<f32>>::new("IF", Card::A.number(), "x")).unwrap(),
    );
    let mut b = b.formats();
    let gt = b.register_format(&FormatFieldDef::<Gt>::new("GT", "genotype")).unwrap();
    let format_ints = CARDS
        .iter()
        .enumerate()
        .map(|(i, &card)| {
            let name = format!("F{i}");
            let def = FormatFieldDef::<Arr<i32>>::new(leak(&name), card.number(), "x");
            (name, card, b.register_format(&def).unwrap())
        })
        .collect();
    let format_floats = (
        Card::G,
        b.register_format(&FormatFieldDef::<Arr<f32>>::new("FF", Card::G.number(), "x")).unwrap(),
    );
    let format_opt = (
        Card::R,
        b.register_format(&FormatFieldDef::<OptArr<f32>>::new("FO", Card::R.number(), "x"))
            .unwrap(),
    );
    let mut b = b.samples();
    b.add_sample("S1").unwrap();
    b.add_sample("S2").unwrap();
    let header = Arc::new(b.build().unwrap());
    Setup {
        header,
        contig,
        dp,
        info_ints,
        info_opt,
        info_floats,
        gt,
        format_ints,
        format_floats,
        format_opt,
    }
}

/// Definitions take `&'static str` names; a test leaks a handful.
fn leak(name: &str) -> &'static str {
    Box::leak(name.to_owned().into_boxed_str())
}

/// A length: the expected one when there is one and `exact` is drawn, any
/// other in `min..=7` otherwise.
fn draw_len(tc: &TestCase, expected: Option<usize>, min: usize) -> usize {
    match expected {
        Some(e) if tc.draw(gs::booleans()) => e,
        _ => tc.draw(gs::integers::<usize>().min_value(min).max_value(7)),
    }
}

fn draw_ints(tc: &TestCase, len: usize) -> Vec<i32> {
    (0..len).map(|_| tc.draw(gs::integers::<i32>().min_value(0).max_value(99))).collect()
}

fn draw_opts(tc: &TestCase, len: usize) -> Vec<Option<i32>> {
    (0..len)
        .map(|_| {
            if tc.draw(gs::booleans()) {
                Some(tc.draw(gs::integers::<i32>().min_value(0).max_value(99)))
            } else {
                None
            }
        })
        .collect()
}

/// Per-sample lengths for a FORMAT field: each sample missing (empty), exact
/// or arbitrary.
fn draw_sample_lens(tc: &TestCase, expected: Option<usize>) -> Vec<usize> {
    (0..N_SAMPLES)
        .map(|_| if tc.draw(gs::booleans()) { 0 } else { draw_len(tc, expected, 1) })
        .collect()
}

fn join<T: ToString>(values: &[T]) -> String {
    values.iter().map(ToString::to_string).collect::<Vec<_>>().join(",")
}

fn join_opt<T: ToString>(values: &[Option<T>]) -> String {
    values
        .iter()
        .map(|v| v.as_ref().map_or_else(|| ".".to_owned(), ToString::to_string))
        .collect::<Vec<_>>()
        .join(",")
}

/// Assert that an encode result is what the oracle says: `Ok` where every
/// length fits, `CountMismatch` naming the first length that does not.
fn check_result(
    result: Result<(), VcfError>,
    name: &str,
    card: Card,
    expected: Option<usize>,
    lens: &[usize],
    missing_ok: bool,
) -> bool {
    let bad = lens
        .iter()
        .copied()
        .find(|&len| !(missing_ok && len == 0) && expected.is_some_and(|e| e != len));
    match (bad, result) {
        (None, Ok(())) => true,
        (Some(found), Err(VcfError::CountMismatch { field, number, expected: e, found: f })) => {
            assert_eq!(field, name);
            assert_eq!(number, card.number(), "{name}");
            assert_eq!(Some(e), expected, "{name}");
            assert_eq!(f, found, "{name}");
            false
        }
        (bad, result) => {
            panic!("{name} ({card:?}) lens {lens:?}: expected bad={bad:?}, got {result:?}")
        }
    }
}

/// One data line of `bcftools view -H`.
fn bcftools_line(bytes: &[u8], suffix: &str) -> String {
    let mut file = tempfile::NamedTempFile::with_suffix(suffix).unwrap();
    std::io::Write::write_all(&mut file, bytes).unwrap();
    std::io::Write::flush(&mut file).unwrap();
    let out =
        Command::new("bcftools").args(["view", "-H"]).arg(file.path()).output().expect("bcftools");
    assert!(out.status.success(), "bcftools view: {}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let mut lines = text.lines();
    let line = lines.next().expect("one record").to_owned();
    assert!(lines.next().is_none());
    line
}

// r[verify record_encoder.count_check+2]
#[hegel::test(test_cases = pinned::cases(60))]
fn the_writer_accepts_exactly_the_declared_counts(tc: TestCase) {
    let s = setup();
    let n_alt = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    let alts = &[Base::C, Base::G, Base::T][..n_alt];
    let alleles = Alleles::snv_multi(Base::A, alts).unwrap();
    let ploidy =
        tc.draw(gs::booleans()).then(|| tc.draw(gs::integers::<usize>().min_value(1).max_value(3)));
    let format = if tc.draw(gs::booleans()) { OutputFormat::Bcf } else { OutputFormat::Vcf };

    // INFO: what bcftools must print, field by field.
    let mut want_info: BTreeMap<String, String> = BTreeMap::new();
    want_info.insert("DP".into(), "42".into());
    // FORMAT: key order, and per key the sample columns.
    let mut want_format: Vec<(String, Vec<String>)> = Vec::new();

    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&s.header).unwrap();
    let mut enc = writer
        .begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass();
    s.dp.encode(&mut enc, 42).unwrap();
    for (name, card, key) in &s.info_ints {
        let expected = card.expected(n_alt, None);
        let values = draw_ints(&tc, draw_len(&tc, expected, 1));
        let result = key.encode(&mut enc, &values);
        if check_result(result, name, *card, expected, &[values.len()], false) {
            want_info.insert(name.clone(), join(&values));
        }
    }
    {
        let (card, key) = &s.info_opt;
        let expected = card.expected(n_alt, None);
        let values = draw_opts(&tc, draw_len(&tc, expected, 1));
        if check_result(
            key.encode(&mut enc, &values),
            "IO",
            *card,
            expected,
            &[values.len()],
            false,
        ) {
            want_info.insert("IO".into(), join_opt(&values));
        }
    }
    {
        let (card, key) = &s.info_floats;
        let expected = card.expected(n_alt, None);
        let values: Vec<f32> =
            draw_ints(&tc, draw_len(&tc, expected, 1)).iter().map(|&v| v as f32).collect();
        if check_result(
            key.encode(&mut enc, &values),
            "IF",
            *card,
            expected,
            &[values.len()],
            false,
        ) {
            want_info.insert("IF".into(), join(&values));
        }
    }

    let mut enc = enc.begin_samples();
    if let Some(p) = ploidy {
        let gts: Vec<Genotype> = (0..N_SAMPLES)
            .map(|_| Genotype {
                alleles: (0..p)
                    .map(|_| {
                        Some(
                            tc.draw(gs::integers::<u16>().max_value(u16::try_from(n_alt).unwrap())),
                        )
                    })
                    .collect(),
                phased: (1..p).map(|_| false).collect(),
            })
            .collect();
        s.gt.encode(&mut enc, &gts).unwrap();
        let columns = gts
            .iter()
            .map(|g| g.alleles.iter().map(|a| a.unwrap().to_string()).collect::<Vec<_>>().join("/"))
            .collect();
        want_format.push(("GT".into(), columns));
    }
    let missing = |s: String| if s.is_empty() { ".".to_owned() } else { s };
    for (name, card, key) in &s.format_ints {
        let expected = card.expected(n_alt, ploidy);
        let per_sample: Vec<Vec<i32>> =
            draw_sample_lens(&tc, expected).into_iter().map(|len| draw_ints(&tc, len)).collect();
        let slices: Vec<&[i32]> = per_sample.iter().map(Vec::as_slice).collect();
        let lens: Vec<usize> = per_sample.iter().map(Vec::len).collect();
        if check_result(key.encode(&mut enc, &slices), name, *card, expected, &lens, true) {
            want_format.push((name.clone(), per_sample.iter().map(|v| missing(join(v))).collect()));
        }
    }
    {
        let (card, key) = &s.format_floats;
        let expected = card.expected(n_alt, ploidy);
        let per_sample: Vec<Vec<f32>> = draw_sample_lens(&tc, expected)
            .into_iter()
            .map(|len| draw_ints(&tc, len).iter().map(|&v| v as f32).collect())
            .collect();
        let slices: Vec<&[f32]> = per_sample.iter().map(Vec::as_slice).collect();
        let lens: Vec<usize> = per_sample.iter().map(Vec::len).collect();
        if check_result(key.encode(&mut enc, &slices), "FF", *card, expected, &lens, true) {
            want_format.push(("FF".into(), per_sample.iter().map(|v| missing(join(v))).collect()));
        }
    }
    {
        let (card, key) = &s.format_opt;
        let expected = card.expected(n_alt, ploidy);
        let per_sample: Vec<Vec<Option<f32>>> = draw_sample_lens(&tc, expected)
            .into_iter()
            .map(|len| draw_opts(&tc, len).iter().map(|v| v.map(|v| v as f32)).collect())
            .collect();
        let slices: Vec<&[Option<f32>]> = per_sample.iter().map(Vec::as_slice).collect();
        let lens: Vec<usize> = per_sample.iter().map(Vec::len).collect();
        // The multi-sample form is the trait method; the key offers only one sample.
        let result = enc.format_float_opts(key.id(), &slices);
        if check_result(result, "FO", *card, expected, &lens, true) {
            want_format
                .push(("FO".into(), per_sample.iter().map(|v| missing(join_opt(v))).collect()));
        }
    }
    enc.emit().unwrap();
    writer.finish().unwrap();

    let suffix = if format == OutputFormat::Bcf { ".bcf" } else { ".vcf" };
    let line = bcftools_line(&out, suffix);
    let cols: Vec<&str> = line.split('\t').collect();
    let got_info: BTreeMap<String, String> = cols[7]
        .split(';')
        .map(|kv| {
            let (k, v) = kv.split_once('=').expect("no flags are written");
            (k.to_owned(), v.to_owned())
        })
        .collect();
    assert_eq!(got_info, want_info, "INFO of {line}");

    if want_format.is_empty() {
        assert_eq!(cols.len(), 8, "no FORMAT column in {line}");
        return;
    }
    let keys: Vec<&str> = cols[8].split(':').collect();
    let want_keys: Vec<&str> = want_format.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, want_keys, "FORMAT keys of {line}");
    for sample in 0..N_SAMPLES {
        let got: Vec<&str> = cols[9 + sample].split(':').collect();
        for (i, (key, columns)) in want_format.iter().enumerate() {
            // htslib drops trailing missing subfields of a sample.
            let value = got.get(i).copied().unwrap_or(".");
            assert_eq!(value, columns[sample], "sample {sample} {key} of {line}");
        }
    }
}
