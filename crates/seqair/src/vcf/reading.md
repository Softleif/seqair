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

// Resolved once, checked against the header: a key the header declares as
// something else fails here, not on some record halfway through the file.
let cpg: InfoFlag = file.header().info_key("CPG")?;
let dp: InfoInt = file.header().info_key("DP")?;

let mut reader = file.reader::<(Locus, Site, Info)>().open()?;
let (mut passing, mut cpg_sites) = (0, 0);
while let Some((locus, site, info)) = reader.next_record()? {
    if site.filters().status() != FilterStatus::Pass {
        continue; // INFO is indexed lazily: a rejected record never pays for it
    }
    passing += 1;
    if info.get(&cpg)? {
        cpg_sites += 1;
        let depth: Option<i32> = info.get(&dp)?;
        assert!(depth.is_some(), "{}:{:?} has no DP", locus.contig().name(), locus.start());
    }
}
assert!(passing > 0 && cpg_sites > 0);
# Ok(()) }
```

`FilterStatus` has three states — `Pass`, `Missing` (`.`) and `Failed` — so
deciding whether `.` counts as passing is the caller's choice, not a silent
library rule.

## One sample: no sample handles

A single-sample file reads its FORMAT fields directly. `OnlySample` checks at
open that there is exactly one sample.

```rust,ignore
use seqair::vcf::{FormatGt, FormatInt, Locus, OnlySample, Site, VcfFile};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let gt: FormatGt = file.header().format_key("GT")?;
let gq: FormatInt = file.header().format_key("GQ")?;

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

let gt: FormatGt = file.header().format_key("GT")?;
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
let header = file.header();
let gt: FormatGt = header.format_key("GT")?;
let na12878 = header.sample("NA12878")?; // a SampleIdx only a header hands out

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

```rust,ignore
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

```rust,ignore
use seqair::vcf::{Locus, VcfFile, VcfPos};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/truth.vcf.gz");
let file = VcfFile::open(path)?;
let mut reader = file.reader::<Locus>().open()?;
let mut longest = 0;
while let Some(locus) = reader.next_record()? {
    if let (VcfPos::At(start), VcfPos::At(end)) = (locus.start(), locus.end()) {
        longest = longest.max(end.get() - start.get() + 1);
    }
}
assert!(longest > 1); // the deletions span more than their first base
# Ok(()) }
```

## Errors you can act on

```rust
use seqair::vcf::{InfoInt, KeyError, VcfFile, ValueType};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
# let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
let file = VcfFile::open(path)?;
// BQ is declared Float: asking for an integer is caught before any record is read.
let err = file.header().info_key::<seqair::vcf::Scalar<i32>>("BQ").unwrap_err();
assert!(matches!(err, KeyError::TypeMismatch { declared: ValueType::Float, requested: ValueType::Integer, .. }));
let missing: Result<InfoInt, _> = file.header().info_key("NOPE");
assert!(matches!(missing, Err(KeyError::UnknownField { .. })));
# Ok(()) }
```

## Your own record type

Tuples cover the common shapes; for anything else, implement `Fields` and
decode straight into your own struct. Everything that depends only on the header
happens once, in `plan`.

```rust,ignore
use seqair::vcf::{
    Allele, Fields, FilterStatus, InfoFlag, PlanError, RawRecord, ReadOptions, Sections, VcfFile, VcfHeader,
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
    fn plan(header: &VcfHeader, _: &ReadOptions) -> Result<Plan, PlanError> {
        Ok(Plan { cpg: header.info_key("CPG")? })
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
