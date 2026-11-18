# Reading VCF and BCF

<!--
These examples are the acceptance target of the reader (docs/spec/2-vcf-*.md).
Each starts as `rust,ignore`; the phase that implements it switches it to a plain
`rust` doctest and makes it pass. `tests/vcf_reading_guide.rs` fails while any
`ignore` remains. Fixtures live in tests/data/vcf/ (see its README).
-->

Open a file once, resolve the fields you need from its header, and read
records into a type that names what you want decoded. A field the type does not
hold is a compile error, a field the header declares differently is an error at
open, and nothing on disk that means "missing" reaches you as a magic value.

## Stream records and read INFO

rastair's truth loading needs FILTER, the alleles and two INFO flags — so it
reads `(Locus, Site, Info)` and never touches the samples.

```rust
use seqair::vcf::{FilterStatus, Info, InfoFlag, InfoInt, Locus, Site, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
let file = VcfFile::open(path)?;

// Resolved once, checked against the file: a key the header declares as
// something else fails here, not on some record halfway through the file.
let cpg: InfoFlag = file.info_key("CPG")?;
let dp: InfoInt = file.info_key("DP")?;

let mut reader = file.reader::<(Locus, Site, Info)>().open()?;
let (mut kept, mut failed, mut cpg_sites) = (0, 0, 0);
while let Some((locus, site, info)) = reader.next_record()? {
    match site.filters().status() {
        FilterStatus::Pass => {}
        // `.`: no filter was applied. Some files never write PASS (GIAB
        // v5.0q, the 1000 Genomes panel): keeping only `Pass` keeps nothing.
        FilterStatus::Missing => {}
        FilterStatus::Failed(_) => {
            failed += 1;
            continue; // INFO is indexed lazily: a rejected record never pays for it
        }
    }
    kept += 1;
    if info.get(&cpg)? {
        cpg_sites += 1;
        let depth: Option<i32> = info.get(&dp)?;
        assert!(depth.is_some(), "{locus} has no DP"); // a locus displays as chr1:10437
    }
}
assert!(kept > 0 && failed > 0 && cpg_sites > 0);
# Ok(()) }
```

`FilterStatus` has three states — `Pass`, `Missing` (`.`) and `Failed` — so
deciding whether `.` counts as passing is the caller's choice, not a silent
library rule. A FILTER of `PASS;q10` is `Failed` and lists `PASS` among the
failed filters; htslib's `bcf_has_filter(…, "PASS")` calls it passing.

## One sample: no sample handles

A single-sample file reads its FORMAT fields directly. `OnlySample` checks at
open that there is exactly one sample.

```rust
use seqair::vcf::{FormatGt, FormatInt, Locus, OnlySample, Site, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let gt: FormatGt = file.format_key("GT")?;
let gq: FormatInt = file.format_key("GQ")?;

let mut reader = file.reader::<(Locus, Site, OnlySample)>().open()?;
let mut hom_alt = 0;
while let Some((_, site, sample)) = reader.next_record()? {
    let Some(genotype) = sample.get(&gt)? else { continue };
    let alleles: Vec<Option<u16>> = genotype.alleles().map(|a| a.allele).collect();
    if alleles == [Some(1), Some(1)] && site.alleles().alts().count() == 1 {
        hom_alt += 1;
        assert!(sample.get(&gq)?.is_some());
    }
}
assert!(hom_alt > 0);
# Ok(()) }
```

On a multi-sample file, `OnlySample` refuses to open rather than silently
reading the first sample; pick one by name with `Sample`:

```rust
use seqair::vcf::{FormatGt, Locus, OnlySample, PlanError, Sample, VcfFile, VcfReadError};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/panel.bcf");
let file = VcfFile::open(path)?;
let err = file.reader::<(Locus, OnlySample)>().open().unwrap_err();
assert!(matches!(err, VcfReadError::Plan(PlanError::NotSingleSample { n: 4 })));

let gt: FormatGt = file.format_key("GT")?;
let mut reader = file.reader::<(Locus, Sample)>().sample("NA12878").open()?;
while let Some((_, sample)) = reader.next_record()? {
    assert_eq!(sample.name(), "NA12878");
    assert!(sample.get(&gt)?.is_some());
}
# Ok(()) }
```

## All samples, by column or by sample

BCF stores FORMAT field by field, so reading one field across all samples is
the natural shape; reading one sample's fields is just as direct.

```rust
use seqair::vcf::{FormatGt, Locus, Samples, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/panel.bcf");
let file = VcfFile::open(path)?;
let gt: FormatGt = file.format_key("GT")?;
let na12878 = file.header().sample("NA12878")?; // a SampleIdx only a header hands out

let mut reader = file.reader::<(Locus, Samples)>().open()?;
let (_, samples) = reader.next_record()?.ok_or("empty file")?;
let column = samples.column(&gt)?.ok_or("no GT")?;
assert_eq!(column.len(), 4);
let carriers = column.filter(|g| g.as_ref().is_some_and(|g| g.alleles().any(|a| a.allele == Some(1)))).count();
assert!(carriers > 0);
assert!(samples.sample(na12878).get(&gt)?.is_some());
# Ok(()) }
```

## Region queries

A query reuses the file's parsed header and index; each call opens its own file
handle, so queries from different threads never contend.

```rust,ignore
use seqair::vcf::{Locus, VcfFile, VcfPos, VcfReadError};
use seqair_types::{Pos1, RegionString};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let region: RegionString = "chr1:600000-800000".parse()?;
let mut query = file.reader::<Locus>().query(&region)?;
let mut n = 0;
while let Some(locus) = query.next_record()? {
    // every record overlaps the region, a deletion starting before it included
    assert!(locus.end() >= VcfPos::At(Pos1::new(600_000).ok_or("pos")?));
    n += 1;
}
assert!(n > 0);

// A contig the header declares but no record uses is an empty query; one the
// file has never heard of is an error you can match on.
let unknown: RegionString = "chrUn_nope".parse()?;
match file.reader::<Locus>().query(&unknown) {
    Err(VcfReadError::UnknownContig { name }) => assert_eq!(name, "chrUn_nope"),
    other => panic!("expected UnknownContig, got {other:?}"),
}
# Ok(()) }
```

Visiting regions one at a time returns a record that spans two of them twice;
`query_many` merges the regions and yields every record once:

```rust,ignore
use seqair::vcf::{Locus, VcfFile};
use seqair_types::RegionString;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let regions: Vec<RegionString> =
    ["chr1:600000-700000", "chr1:700001-800000", "chr2:1-400000"].iter().map(|r| r.parse()).collect::<Result<_, _>>()?;
let mut query = file.reader::<Locus>().query_many(&regions)?;
let mut seen = std::collections::HashSet::new();
while let Some(locus) = query.next_record()? {
    let key = (locus.contig().name().to_owned(), locus.start());
    assert!(seen.insert(key), "a record came back twice");
}
# Ok(()) }
```

## Many threads

`VcfFile` is `Send + Sync + Clone`. Split a file into partitions and read each on
its own thread; every record belongs to exactly one partition, so the union is
the file, with nothing visited twice.

```rust,ignore
use seqair::vcf::{Locus, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
let file = VcfFile::open(path)?;
let partitions = file.partitions(4)?;
let counts = std::thread::scope(|s| {
    let handles: Vec<_> = partitions
        .iter()
        .map(|p| {
            let file = &file;
            s.spawn(move || -> Result<usize, seqair::vcf::VcfReadError> {
                let mut reader = file.reader::<Locus>().partition(p)?;
                let mut n = 0;
                while reader.next_record()?.is_some() {
                    n += 1;
                }
                Ok(n)
            })
        })
        .collect();
    handles.into_iter().map(|h| h.join().expect("thread")).collect::<Result<Vec<_>, _>>()
})?;
assert_eq!(counts.iter().sum::<usize>(), 400);
# Ok(()) }
```

## Owning records: `read_into`, `records()`

`next_record` lends a record until the next call. To keep records — one per
input in a merge, a window of neighbours — own them: `read_into` reuses your
buffer, `records()` hands out a fresh record each time.

```rust
use seqair::vcf::{Locus, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf");
// Two inputs merged by position, one owned record per input (HiPhase's pattern).
let files = [VcfFile::open(format!("{dir}/rastair.bcf"))?, VcfFile::open(format!("{dir}/rastair.vcf.gz"))?];
let mut readers = files.iter().map(|f| f.reader::<Locus>().open()).collect::<Result<Vec<_>, _>>()?;
let mut heads: Vec<Locus> = Vec::new();
for reader in &mut readers {
    // `read` hands out an owned record; later `read_into(&mut rec)` reuses its buffers.
    heads.push(reader.read()?.ok_or("empty input")?);
}
// Same record from BCF and VCF text. (The loci themselves differ: each names
// a contig of its own file's header, and a contig id never crosses headers.)
assert_eq!(heads[0].contig().name(), heads[1].contig().name());
assert_eq!((heads[0].start(), heads[0].end()), (heads[1].start(), heads[1].end()));

// Or an ordinary iterator of owned records, for peekable(), channels, collect:
let all: Vec<Locus> = files[0].reader::<Locus>().open()?.records().collect::<Result<_, _>>()?;
assert_eq!(all.len(), 400);
# Ok(()) }
```

## Positions only

riker masks the reference span of every variant. `Locus` alone skips everything
past the fixed fields; the end already accounts for INFO END, SVLEN and gVCF
blocks, as htslib computes it.

```rust
use seqair::vcf::{Locus, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let mut reader = file.reader::<Locus>().open()?;
let mut longest = 0;
while let Some(locus) = reader.next_record()? {
    // The reference bases, as 1-based inclusive positions; `None` for a
    // record on the telomere before base 1 only.
    if let Some(bases) = locus.bases() {
        longest = longest.max(bases.end().as_u32() - bases.start().as_u32() + 1);
    }
    let start0 = locus.start().pos0(); // 0-based, `None` for the telomere: no -1
    assert!(start0.is_some());
}
assert!(longest > 1); // the deletions span more than their first base
# Ok(()) }
```

## Fields a header leaves out, and fields a record lacks

Not every file declares what it uses, and not every record holds every
field. A VCF text file resolves keys against its header *and* the keys VCF
4.5 reserves (`SVTYPE`, `SVLEN`, `END`, `PS`, …), typed as VCF defines them;
any other key it does not declare reads as a String, with a warning. An
optional field is an `Option` of a key, decided once. And `has` tells a
field the record does not hold from one it holds as `.`.

```rust
use seqair::vcf::{FormatInt, Info, InfoAny, InfoFlag, InfoString, Locus, OnlySample, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let text = "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n\
chr1\t100\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50;TRID=t1\tGT:PS\t0|1:100\n\
chr1\t200\t.\tC\tT\t.\tPASS\tTRID=.\tGT:PS\t0/1:.\n\
chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t1/1\n";
let file = VcfFile::from_bytes(text.as_bytes().to_vec(), None)?;

// Reserved by VCF, so typed as VCF says, though the header is silent.
let svtype: InfoString = file.info_key("SVTYPE")?;
let ps: FormatInt = file.format_key("PS")?;
// A tool's own key: present or not, whatever its type.
let trid: InfoAny = file.info_key("TRID")?;
// An annotation this file may not have: `None` here, an error only if
// the file declares it as something a flag cannot read.
let cpg: Option<InfoFlag> = file.info_key_opt("CPG")?;
assert!(cpg.is_none());

let mut reader = file.reader::<(Locus, Info, OnlySample)>().open()?;
let mut seen = Vec::new();
while let Some((locus, info, sample)) = reader.next_record()? {
    let tandem_repeat = info.get(&trid)?.is_some_and(|v| !v.is_missing());
    let phase = match (sample.has(&ps)?, sample.get(&ps)?) {
        (false, _) => "no PS",  // FORMAT does not name PS
        (true, None) => "PS .", // named, but `.` for this sample
        (true, Some(_)) => "phased",
    };
    let svtype = info.get(&svtype)?.map(str::to_owned);
    seen.push((locus.to_string(), svtype, tandem_repeat, phase));
}
assert_eq!(
    seen,
    [
        // The undeclared SVLEN sets the end, as tabix reads it.
        ("chr1:100-150".to_owned(), Some("DEL".to_owned()), true, "phased"),
        ("chr1:200".to_owned(), None, false, "PS ."),
        ("chr1:300".to_owned(), None, false, "no PS"),
    ]
);
# Ok(()) }
```

A key the header does not declare reads only this file: a writer refuses it,
since the header it writes has no id for the field. BCF names fields by the
header's dictionary, so a BCF file has exactly the header's fields.
`OptArr<T>` (`InfoIntOpts`, `FormatFloatOpts`, …) reads a field of any
`Number` and lets elements be `.`; `Arr<T>` (`FormatFloats`, …) fails on an
element missing inside the vector.

## Errors you can act on

```rust
use seqair::vcf::{InfoInt, KeyError, KeyKind, Scalar, VcfFile, VcfReadError, ValueType};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
let file = VcfFile::open(path)?;
// BQ is declared Float: asking for an integer is caught before any record is read.
let err = file.info_key::<Scalar<i32>>("BQ").unwrap_err();
assert!(matches!(err, KeyError::TypeMismatch { declared: ValueType::Float, requested: ValueType::Integer, .. }));
let missing: Result<InfoInt, _> = file.info_key("NOPE");
assert!(matches!(missing, Err(KeyError::UnknownField { kind: KeyKind::Info, .. })));

// A malformed record says where it is; what is wrong is its source.
let text = "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t5\t.\tA\tC\t.\t.\tDP=lots\n";
let file = VcfFile::from_bytes(text.as_bytes().to_vec(), None)?;
let dp: InfoInt = file.info_key("DP")?;
let info = file.reader::<seqair::vcf::Info>().open()?.read()?.ok_or("empty")?;
let err = info.get(&dp).unwrap_err();
assert_eq!(err.to_string(), "malformed record 1 (line 5) at chr1:5");
let VcfReadError::Record { at, source } = &err else { unreachable!() };
assert_eq!((at.number, at.line), (1, Some(5)));
assert!(std::error::Error::source(source).is_some_and(|e| e.to_string().contains("\"lots\"")));
# Ok(()) }
```

## Your own record type

Tuples cover the common shapes; for anything else, implement `Fields` and
decode straight into your own struct. Everything that depends only on the file
happens once, in `plan`, through its `PlanContext`.

```rust
use seqair::vcf::{
    Allele, Fields, FilterStatus, InfoFlag, PlanContext, PlanError, RawRecord, Sections, VcfFile,
    VcfPos, VcfReadError,
};

/// What rastair's `verify` keeps of a truth record: passing SNVs only.
#[derive(Clone, Debug)]
struct TruthSnv {
    pos: VcfPos,
    ref_base: u8,
    alt_base: u8,
    cpg: bool,
}

#[derive(Clone, Debug)]
struct Truth(Option<TruthSnv>);

struct Plan {
    cpg: InfoFlag,
}

impl Fields for Truth {
    type Plan = Plan;
    fn plan(cx: &PlanContext<'_>) -> Result<Plan, PlanError> {
        Ok(Plan { cpg: cx.info_key("CPG")? })
    }
    fn needs(_: &Plan) -> Sections {
        Sections::SITE | Sections::INFO
    }
    fn new(plan: &Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let site = raw.site()?;
        let mut alts = site.alleles().alts();
        let snv = match (site.filters().status(), site.alleles().reference(), alts.next(), alts.next()) {
            (FilterStatus::Pass, &[r], Some(Allele::Bases(&[a])), None) => Some(TruthSnv {
                pos: raw.locus()?.start(),
                ref_base: r,
                alt_base: a,
                cpg: raw.info()?.get(&plan.cpg)?,
            }),
            _ => None,
        };
        Ok(Truth(snv))
    }
}

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
let file = VcfFile::open(path)?;
let mut reader = file.reader::<Truth>().open()?;
let mut snvs = Vec::new();
while let Some(Truth(snv)) = reader.next_record()? {
    snvs.extend(snv.clone());
}
assert!(!snvs.is_empty());
# Ok(()) }
```

`Fields` has no error type of its own: `Err` is a malformed record. An
outcome of your domain — a call without a genotype, a site you skip — is
part of the value, as `Truth(Option<TruthSnv>)` above or `Call` below. A type
reading one sample's FORMAT finds it with `selected_sample()`, the selection
`Sample` itself makes, so `.sample(name)` and "the only sample" work for it
too:

```rust
use seqair::vcf::{
    Fields, FormatGt, PlanContext, PlanError, RawRecord, SampleIdx, Sections, VcfFile, VcfReadError,
};

/// One sample's genotype, or why it has none.
#[derive(Clone, Debug)]
struct Call(Result<Vec<Option<u16>>, NoCall>);

#[derive(Clone, Debug, PartialEq)]
enum NoCall {
    /// The record's FORMAT has no GT.
    NoGenotype,
    /// GT is `.` for this sample.
    Missing,
}

struct CallPlan {
    sample: SampleIdx,
    gt: FormatGt,
}

impl Fields for Call {
    type Plan = CallPlan;
    fn plan(cx: &PlanContext<'_>) -> Result<CallPlan, PlanError> {
        Ok(CallPlan { sample: cx.selected_sample()?, gt: cx.format_key("GT")? })
    }
    fn needs(plan: &CallPlan) -> Sections {
        Sections::sample(plan.sample) // that sample's bytes only
    }
    fn new(plan: &CallPlan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let sample = raw.sample(plan.sample)?;
        Ok(Call(match (sample.has(&plan.gt)?, sample.get(&plan.gt)?) {
            (false, _) => Err(NoCall::NoGenotype),
            (true, None) => Err(NoCall::Missing),
            (true, Some(gt)) => Ok(gt.alleles().map(|a| a.allele).collect()),
        }))
    }
}

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/panel.bcf");
let file = VcfFile::open(path)?;
let mut reader = file.reader::<Call>().sample("NA12878").open()?;
let mut called = 0;
while let Some(Call(call)) = reader.next_record()? {
    called += usize::from(call.is_ok());
}
assert!(called > 0);

// Runtime context — a sample chosen elsewhere, a threshold, a region — lives
// in a plan you build; `with_plan` reads with it instead of calling `plan`.
let plan = CallPlan { sample: file.header().sample("HG00096")?, gt: file.format_key("GT")? };
let calls: Vec<Call> = file.reader::<Call>().with_plan(plan).open()?.into_iter().collect::<Result<_, _>>()?;
assert_eq!(calls.len(), 150);
# Ok(()) }
```

A record that makes zero or several values of yours (a multi-allelic site
split per ALT, a symbolic allele skipped) does not fit one `Fields` value:
read a tuple of components and turn each record into your values after.

```rust
use seqair::vcf::{Allele, Locus, Site, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let mut reader = file.reader::<(Locus, Site)>().open()?;
let mut variants: Vec<(Locus, Vec<u8>)> = Vec::new();
while let Some((locus, site)) = reader.next_record()? {
    for alt in site.alleles().alts() {
        if let Allele::Bases(bases) = alt {
            variants.push((locus.clone(), bases.to_vec()));
        }
    }
}
assert!(variants.len() >= 1237);
# Ok(()) }
```
