# VCF/BCF Fields: choosing what a reader decodes

A reader is parameterised by the type it reads into. The type names what is copied and decoded per record, and what the caller can call: a field the type does not hold is a compile error, not a runtime check.

```rust
type Truth = (Locus, Site, Info);
let mut r = file.reader::<Truth>().open()?;
while let Some((locus, site, info)) = r.next_record()? { … }

let mut r = file.reader::<(Locus, Site, OnlySample)>().open()?;   // Err(NotSingleSample { n })
let mut r = file.reader::<(Locus, Sample)>().sample("HG002").open()?;
let mut r = file.reader::<FullRecord>().open()?;
```

The access patterns this covers: riker reads only positions (`Locus`); aardvark reads alleles and one sample's GT (`(Locus, Site, Sample)`); rastair `verify`/`train` read FILTER, alleles and two INFO flags (`(Locus, Site, Info)`); `convert` and HiPhase's writer read everything (`FullRecord`).

## The trait

> r[vcf_fields.trait+3]
> `Fields` MUST be a public trait implemented by every component, every tuple of components, `FullRecord`, and any user type:
>
> ```rust
> pub trait Fields: Sized + Send + Sync + Clone {
>     type Plan: Send + Sync;
>     fn plan(cx: &PlanContext<'_>) -> Result<Self::Plan, PlanError>;
>     fn needs(plan: &Self::Plan) -> Sections;
>     fn new(plan: &Self::Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError>;
>     fn fill(&mut self, plan: &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
>         *self = Self::new(plan, raw)?;
>         Ok(())
>     }
> }
> ```
>
> `Fields` has no error type of its own: a domain outcome that is not a malformed record (a call without GT, a site the caller skips) belongs in the value — `struct Call(Result<VariantCall, CallError>)`, `Truth(Option<TruthSnv>)` — and a record that yields 0..n values is read as a tuple of components and post-processed. Both patterns MUST be shown in the reading guide.
>
> A value of a `Fields` type always holds a record: there is no empty or placeholder state (no contig "0 with an empty name" before the first read). `new` builds one from a record; `fill` overwrites one in place, and components override it to reuse their buffers. A user type implements `plan`, `needs` and `new`, and overrides `fill` only for speed.
>
> `RawRecord` is the reader's undecoded view of one record (the BCF sections or the text line with its column offsets); its decoding helpers are the only way a `fill` reads bytes. `Sections` says which parts of a record `fill` will read: the site columns, INFO, and which samples (none, one by index, all).

r[vcf_fields.plan]
A reader MUST call `F::plan` once, when it opens, and keep the plan for every record. Everything that depends only on the header — resolved keys, the selected sample's index, type checks — MUST happen in `plan`, so a type mismatch fails at open rather than on some record. Planning MUST be cheap enough to repeat per `file.query` (it runs once per reader, not per file).

r[vcf_fields.plan_context]
`plan` MUST receive a `PlanContext` rather than the header and the options: `header()`; `info_key`, `format_key` and their `_opt` forms, resolved as the file resolves them (`r[vcf_values.file_keys+2]`), and `lenient()` (`r[vcf_values.lenient_keys]`); `filter(name)`; and `selected_sample() -> Result<SampleIdx, PlanError>`, the sample selection `Sample` itself uses (`r[vcf_fields.sample+2]`), so a user type reading one sample's FORMAT needs no logic of its own.

r[vcf_fields.with_plan]
`builder.with_plan(plan)` MUST open the reader with a plan the caller built instead of calling `F::plan`. That is the way to give a type runtime context (a length band, a region, a trimming switch): it lives in the caller's `Plan`, whose keys come from `file.info_key` and friends.

r[vcf_fields.sections]
The reader MUST copy and index only the union of what `F::needs` asks for. Bytes of sections no component needs MUST be skipped in the decompressed stream without being copied: a BCF record's `l_indiv` bytes when no sample is needed, every sample but the selected one, and, for `Locus` alone, everything after the fixed fields (the BCF `rlen` is a fixed field; VCF text computes the end while reading, `r[vcf_values.rlen]`). Decompression itself cannot be skipped.

r[vcf_fields.tuples]
`Fields` MUST be implemented for tuples of 1 to 8 elements that each implement `Fields`, by macro. A tuple's plan is the tuple of its elements' plans, its `needs` the union, and its `fill` fills each element from the same `RawRecord`. Element order has no effect on what is read. A duplicated or redundant element (`(Info, Info)`, `(Sample, Samples)`) is legal and costs only the redundant copy.

r[vcf_fields.unused_options+2]
The builder's options are consumed through the `PlanContext` (the sample selection, by `selected_sample()`). Opening MUST fail with `PlanError::UnusedOption { option }` when an option was set that nothing consumed, so `.sample("HG002")` on a reader without a `Sample` component — or with a plan given by `with_plan` — is an error, not silently ignored.

r[vcf_fields.owned]
Every component MUST be an owned value — `Clone + Send + Sync`, no lifetimes — holding its own reusable buffers. `fill` MUST reuse them, so `read_into` allocates nothing per record once the buffers have grown. A record read through one reader stays valid after the reader moves on or is dropped.

r[vcf_fields.lazy+2]
A component that indexes its bytes (`Info`, the sample components) MUST build the index on first access, not in `fill`, so a record the caller rejects on its position, FILTER or alleles never pays for it. The component MUST stay `Sync`: concurrent first accesses MUST build once, the others waiting for that build. `fill` MUST keep the previous index's buffers for the next build. A failed build (malformed data) MUST be stored and returned as the same typed error on every access. The index lives in place behind one atomic state byte (empty, building, done), so a record costs one atomic exchange; a `std::sync::OnceLock` plus a lock for the spare buffers cost four locked atomics per record. The `unsafe` this takes MUST be covered by a multi-threaded test of racing first accesses.

## Components

r[vcf_fields.locus+3]
`Locus` MUST hold the record's `ContigId`, its start and its inclusive end, and be `Clone + Eq + Hash`, cheap to clone (a contig name is a `SmolStr`). Start and end are `VcfPos`, an enum of `Telomere` (VCF POS 0, BCF POS −1) and `At(Pos1)`; POS N+1, the other telomere, is an ordinary `At`. `VcfPos::pos1()`/`pos0()` MUST return the position as `Option<Pos1>`/`Option<Pos0>` (`None` for the telomere), and `Locus::bases()` the reference bases the record spans as `Option<RangeInclusive<Pos1>>` (from base 1 when it starts at the telomere; `None` when it spans the telomere only), so htslib ports need no sentinel. `Locus` needs no section beyond the fixed fields.

r[vcf_fields.locus.bases0]
`Locus::bases0()` MUST return the same bases as `bases()` as a 0-based half-open `Option<Range<Pos0>>`, so BED-style callers (riker's masks) need no hand conversion between inclusive 1-based and exclusive 0-based ends.

r[vcf_fields.site+2]
`Site` MUST hold the ID, REF, ALT, QUAL and FILTER columns and expose:

- `ids()`: an iterator over the `;`-separated IDs as `&str`, empty for `.`;
- `alleles()`: the REF bytes, an iterator over ALT alleles as `Allele` — `Bases(&[u8])`, `Symbolic(&str)` (`<DEL>`, `<*>`, …), `Breakend(&str)`, `Overlapping` (`*`), `Missing` (a `.` inside an ALT list, kept so allele indices stay aligned) — with no ALT yielding nothing, `get(i: u16)` by GT allele index (0 is REF), `iter()` over REF and the ALTs, and `to_alleles() -> Result<Alleles, AllelesError>`, the writer's typed alleles. `Allele::as_bytes()` MUST return the allele as written. No ALT is VCF text's `.` and BCF's single `.` allele after REF (`n_allele` 2, as htslib writes a text `.`); a lone `.` is never `Bases`;
- `qual() -> Option<f32>` (`r[vcf_values.qual]`);
- `filters()`, whose `status()` is `FilterStatus::Pass` (exactly `PASS`), `Missing` (`.`, or BCF's empty vector) or `Failed(…)` iterating the `FilterId`s, plus `contains(&FilterId)`. htslib's `bcf_has_filter` counts `.` as PASS; an enum makes each caller decide. A FILTER that lists `PASS` with other filters (`PASS;q10`) is `Failed`, listing `PASS` among them; htslib's `bcf_has_filter(…, "PASS")` answers yes for it. A file may never write `PASS` (GIAB v5.0q), so code that keeps only `Pass` keeps nothing there: the guide's first example MUST match all three states.

r[vcf_fields.info+2]
`Info` MUST hold the INFO bytes and its lazy index, and expose `get(&InfoKey<V>)` and `has(&InfoKey<V>)` (`r[vcf_values.access+2]`) and `iter()` over the fields present (`r[vcf_values.dynamic+2]`).

r[vcf_fields.only_sample]
`OnlySample` MUST fail to plan with `PlanError::NotSingleSample { n }` unless the header has exactly one sample, and expose `get(&FormatKey<V>)` returning that sample's value with no sample handle. It needs only that sample's bytes. There MUST NOT be a "first sample" component: on a multi-sample file it would silently pick sample 0.

r[vcf_fields.sample+2]
`Sample` MUST read the sample the builder's `.sample(sel)` selects — a name (`&str`, `SmolStr`) or a `SampleIdx` — failing to plan with `PlanError::Key(KeyError::UnknownSample { name })` for a name the header lacks, `ForeignKey` for a `SampleIdx` of another header, and `SampleNotSelected { n }` when nothing was selected and the header has more than one sample (with exactly one, it reads that one). It exposes the same `get`/`has` as `OnlySample` and `name()`.

r[vcf_fields.samples]
`Samples` MUST hold every sample and expose both views of the field-major FORMAT data: `column(&FormatKey<V>)`, an `ExactSizeIterator` over all samples' values with `get(SampleIdx)`, and `sample(SampleIdx)`, one sample's view with the same `get` as `OnlySample`. For VCF text the per-sample column offsets are built lazily (`r[vcf_fields.lazy]`) in one vectorised pass.

r[vcf_fields.sample_idx]
`SampleIdx` MUST be obtainable only from a header (`header.sample(name)`, `header.samples()`) and carry the header's `HeaderId`; using it with a record from another header MUST fail with `ForeignKey`. A raw `usize` cannot name a sample.

r[vcf_fields.full_record]
`FullRecord` MUST hold every section (`Locus`, `Site`, `Info`, `Samples`) plus the raw record bytes and the `HeaderId`, so a later writer can copy a record unchanged when the output header is the same, and expose each part by reference.

## Reading guide

r[vcf_fields.guide]
`crates/seqair/src/vcf/reading.md` is the user-facing guide to reading, included as documentation of the public reader module, and every example in it MUST run as a doctest. It is the acceptance target for the API: examples change only when the API is improved, never to fit an implementation. A test MUST fail while any example in it is still marked `ignore`.

## User types

r[vcf_fields.user_types+2]
A user type MUST be able to implement `Fields` by hand with only public API: resolving keys and its sample in `plan` through the `PlanContext`, declaring `needs`, and decoding in `fill` through `RawRecord`'s helpers. This is the extension point for "deserialize into my struct". A `#[derive(Fields)]` macro that maps struct fields to `#[vcf(info = "…")]`, `#[vcf(format = "…")]`, `#[vcf(filter)]` and components SHOULD follow once a caller needs it; it is out of scope for the first version.
