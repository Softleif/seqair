# VCF/BCF Reader

A pure-Rust reader for VCF (plain or BGZF-compressed) and BCF, built on seqair's existing BGZF, `RegionBuf`, TBI and CSI layers. It replaces rust-htslib in rastair's `convert`, `verify` and `train`, and is shaped by the access patterns of the other Rust tools that read VCFs (riker, HiPhase, aardvark): stream and extract, k-way merge across files, read–modify–write, and one rayon task per region.

> **Sources:** [VCF45] §1.6 "Data lines", §6 "BCF specification". [TABIX], [CSI] for region queries. [htslib] `vcf.c` (`get_rlen`, `bcf_has_filter`), `tbx.c` (VCF end column) for behaviour the specs leave open. See [References](./99-references.md).

What a reader decodes is chosen by a `Fields` type ([VCF/BCF Fields](./2-vcf-2-fields.md)); how values come out is in [VCF/BCF Values](./2-vcf-3-values.md).

## Opening

r[vcf_reader.formats]
The reader MUST read plain-text VCF, BGZF-compressed VCF and BCF 2.2 (BGZF-framed, at any compression level including 0). The format MUST be detected from the content — BGZF magic, then the decompressed `BCF\x02\x02` magic or `##fileformat=` — never from the file extension. A gzip stream that is not BGZF MUST be rejected with a typed error that suggests `bgzip`. BCF 2.1 (`BCF\x02\x01`) MUST be rejected with a typed error.

r[vcf_reader.file+2]
`VcfFile::open(path)` MUST parse the header (`r[vcf_header.parse]`) and locate an index by the names of `r[index.locate]` (VCF: `csi`, then `tbi`; BCF: `csi`). `VcfFile::from_bytes(data, index)` takes the index's bytes and MUST tell TBI from CSI by content (BGZF-compressed or not). A missing index MUST NOT fail the open: sequential reading works without one, `has_index()` is `false`, and only a query or `partitions` fails, with `VcfReadError::NoIndex`. An index that exists but does not parse MUST fail the open.

r[vcf_reader.file.in_memory_index]
`VcfFile::from_bytes_indexed(data)` MUST open bytes in memory with a CSI it builds by reading every record once: each record's contig, its span (`r[vcf_values.rlen+3]`; a record at the telomere counts from the first base and covers at least it, `r[vcf_reader.query.overlap+2]`) and the virtual offset it starts at go through the writer's `IndexBuilder` (`r[index_builder.single_pass]`), as a tabix CSI naming the indexed contigs for VCF text and a CSI numbering the contig dictionary for BCF. Plain VCF text MUST first be BGZF-compressed in memory, so such a file's `format()` is `VcfGz`: the format reports what is read. Records out of order MUST fail with `OpenError::IndexRecords`, as `bcftools index` refuses them; a malformed record fails with its record error. `from_bytes` MUST NOT index: the cost is the caller's explicit choice. riker and aardvark unit-test their query code on VCF text built in the test, which had no index. A property test MUST check that queries of a file indexed so equal queries of the same data on disk indexed by bcftools, on generated records that span region boundaries through REF, END and SVLEN.

r[vcf_reader.file.header_only]
Opening a BGZF file MUST inflate no block past the one its header ends in: where records start is the header's end in that block or, for a header ending at a block's end (bcftools and htslib write it so), the next block's offset, which the block's length gives. HiPhase opens a file per block of work, and inflating the read-ahead batch after the header cost ~340 µs per open.

r[vcf_reader.file.shared+2]
`VcfFile` MUST be `Send + Sync + Clone`, with its parsed state — header, index, the field table (`r[vcf_reader.field_table]`) and its text-key lookup (`r[vcf_values.text_resolution+2]`), the contig table (`r[vcf_reader.contigs]`) and the warn-once state (`r[vcf_reader.undeclared+4]`) — behind one `Arc`. Cloning it MUST NOT re-read the file. This is what lets a rayon task query a file without reopening and re-parsing it (aardvark does that today because noodles' readers are not `Send`).

## Readers

r[vcf_reader.reader]
`file.reader::<F: Fields>()` MUST return a builder whose `open()` (or `?` on the builder's terminal method) plans `F` against the header (`r[vcf_fields.plan]`) and yields a `Reader<F>` positioned at the first record. A `Reader<F>` owns its I/O state — its own file handle, BGZF state and block cache — and shares everything else with its `VcfFile`. `Reader<F>` MUST be `Send`; it need not be `Sync`. Opening a reader MUST cost one file open and no header or index parse.

r[vcf_reader.reader.threads+2]
The builder MUST offer `threads(n)` and `threads_in(n, pool)` with `n: NonZeroUsize`, passed to the existing `BgzfStream` (`r[bgzf.stream.threads]`, `r[bgzf.stream.shared_pool]`) for sequential reading, so a caller can bound total decompression threads across many readers; `n` above 255 is capped at 255. Queries and partitions accept the setting and decompress on the calling thread.

r[vcf_reader.read_into+2]
`Reader::read_into(&mut self, rec: &mut F) -> Result<bool, VcfReadError>` MUST be the primitive: it overwrites `rec` with the next record (`Fields::fill`) and returns `false` at end of input, leaving `rec` as it was. It MUST reuse `rec`'s allocations, so steady-state reading allocates nothing per record. The first record a caller owns comes from `Reader::read(&mut self) -> Result<Option<F>, VcfReadError>` (`Fields::new`); there is no empty record to read into, since a value of `F` always holds a record.

r[vcf_reader.next_record]
`Reader::next_record(&mut self) -> Result<Option<&F>, VcfReadError>` MUST read into a record the reader owns and lend it out until the next call. It is the zero-setup streaming form; the borrow stops a caller from keeping the record, and `clone()` makes an owned copy.

r[vcf_reader.records]
`Reader::records(self)` (and `&mut self`) MUST return an `Iterator<Item = Result<F, VcfReadError>>` yielding owned records, allocating per record. It exists for `peekable()`, channels, `collect` and `par_bridge`; the docs MUST point hot loops at `read_into` or `next_record`.

r[vcf_reader.fork_equivalence]
Any number of readers opened from one `VcfFile`, on any threads, MUST yield the same records for the same input. Tests MUST read a file from several threads at once and compare each thread's records with a single-threaded read.

## Contigs

r[vcf_reader.contigs]
Records name their contig with a `ContigId` from the file's contig table. The table MUST start as the header's contigs in dictionary order, followed by contigs only the index names (a tabix index lists its own sequence names). A VCF text record whose CHROM is in neither MUST be added to the table at the end, with a single `warn!` per file naming the first such contig. ClinVar ships without `##contig` lines, so rejecting them would refuse a common file. The table is append-only and shared by all readers of the file; a reader MUST look a CHROM up only when it differs from the previous record's, since sorted files change contig rarely. A BCF record whose contig index is not in the header's contig dictionary MUST fail with a typed error.

r[vcf_reader.contig_lookup]
`VcfFile::contig(name) -> Option<ContigId>` MUST look a name up in the file's contig table — the header's contigs, the index's, and those records named — so a caller can ask whether a query of that contig would resolve before making it. `header().contig(name)` misses the contigs only the index names, and riker queries every contig its BAM has, most of which a VCF lacks.

## Regions

r[vcf_reader.region]
A `Region` MUST be a contig (`RegionContig::Name` or `RegionContig::Id`, from `&str`, `SmolStr` or `&ContigId`), a 0-based start, and an end that is either the contig's end (`RegionEnd::ContigEnd`) or an exclusive 0-based position (`RegionEnd::Before`). Its fields are private: `Region::contig(c)` names a whole contig, and `Region::new(c, range)` takes any `RangeBounds<P>` over a sealed `Coordinate` — `Pos0` ranges are 0-based (`a..b` half-open), `Pos1` ranges 1-based (`a..=b` inclusive, as VCF and region strings count). An unbounded start is the contig's first base, an unbounded end its last. `a..a` is an empty region; a start after the end MUST fail at construction with `RegionError::Reversed`. `Region` MUST convert from `&RegionString` (`chr1:100-200`, `chr1:100`, `chr1`) and parse from a string through `FromStr`.

## Region queries

r[vcf_reader.query+2]
The builder's `query(region)` (anything `Into<Region>`) MUST plan `F` like `open()` and yield a `Reader<F>` — the same type sequential reading returns — that reads only the region's records. `Reader::requery(region)` MUST point an existing reader at a new region, reusing its file handle and block cache, for callers that issue many small queries on one thread. The query MUST stream through the existing `RegionBuf` ([RegionBuf](./3-region_buf.md)) with the TBI or CSI chunks of the region, and MUST yield records in file order. A query on a file without an index MUST fail with `VcfReadError::NoIndex`.

r[vcf_reader.query.overlap+2]
Under `Match::Overlapping` (the default) a record MUST be yielded iff its reference span `[start, end]` (`r[vcf_values.rlen]`) overlaps the queried range; under `Match::StartsIn`, set with the builder's `matching(..)`, iff its start lies in the range. Since every record that starts inside the range overlaps it, the reader MUST compute the end only for records that start before the range — for VCF text that is the only case that needs an INFO/FORMAT scan. As tabix does, a record at POS 0 (the telomere) counts as starting at the first base, and a span MUST cover at least that base. The first record that starts at or past the last range of its contig ends the reading of that contig. `bcftools view -r` selects the same records as `Overlapping`, `-t` as `StartsIn`.

r[vcf_reader.query.contigs+2]
A contig the file declares but the index does not list (it has no records) MUST give an empty query, not an error. A contig name the file does not know MUST fail with `VcfReadError::UnknownContig { name }`, a variant the caller can match to skip it; a `ContigId` of another file's header MUST fail with `VcfReadError::ForeignKey`. rastair and aardvark both carry workarounds for the first case today.

r[vcf_reader.query.many+2]
`query_many(regions)` and `Reader::requery_many(regions)` MUST sort the regions, merge overlapping and adjacent ones per contig, and yield every record at most once, in file order; under `StartsIn` a record is yielded iff it starts in the union. Visiting regions one by one returns a record that spans two adjacent regions twice; `query_many` is the way to avoid that. An empty list of regions yields no record.

r[vcf_reader.query.options]
`matching(..)` set on a builder that ends in `open()` or `partition(..)` MUST fail with `PlanError::UnusedOption`: sequential reading has no region, and a partition always reads `StartsIn`.

## Partitions

r[vcf_reader.partitions+2]
`file.partitions(n)` with `n: NonZeroUsize` MUST tile every contig the index lists into disjoint `Partition`s — whole regions from the contig's first base to its end, so a contig gets at least one — about `n` in all, with boundaries at the index's leaf-bin positions chosen so that partitions hold roughly equal compressed bytes. Partitions MUST come in file order, be `Send + Sync + Clone`, and name their contig by `ContigId`, so a partition carries its file's identity. `Partition::region()` returns its region.

r[vcf_reader.partitions.ownership+2]
The builder's `partition(&partition)` MUST yield exactly the records whose start lies inside the partition (`Match::StartsIn`) — not those that merely overlap it — so a record spanning a boundary belongs to one partition only. HiPhase guards against the double visit by hand ("we have already written though, so don't write it again"); partition reads make that guard unnecessary.

r[vcf_reader.partitions.complete+2]
The records of all partitions of a file, concatenated in partition order, MUST equal one sequential read of the file's records on indexed contigs. A property test MUST check this — no record missing, none duplicated, order preserved — against a sequential read as the oracle, for random files (VCF.gz with CSI and with TBI, BCF with CSI) and random `n`.

## Streams

r[vcf_reader.stream]
`VcfStream::new(input)` (any `Read + Send + 'static`; `VcfStream::stdin()` for standard input) MUST read VCF, BGZF VCF or BCF from input that cannot seek — a pipe, standard input — detecting the format and parsing the header as `VcfFile` does without seeking or reading the input twice. `stream.reader::<F>()` MUST return a builder that ends in `open()` only — a stream has no index, so no query method exists on it — and yields the same `Reader<F>` a file gives; the stream is consumed by it.

## Undeclared fields

r[vcf_reader.field_table]
A VCF text file's records and a file's lenient keys (`r[vcf_values.lenient_keys]`) MUST resolve INFO and FORMAT keys against a field table: the header's string dictionary, extended with

- every VCF 4.5 reserved INFO and FORMAT key (the specification's tables of reserved INFO keys and genotype keys, and its structural-variant INFO and FORMAT keys; the `M[0-9]+[ACGTUN]` patterns excepted) that the header does not declare in that section, with its reserved `Type` and `Number`;
- every other key a record uses or a lenient key resolves, as `Type=String, Number=.`.

Keys new to the dictionary get ids past the header's (the reserved ones first, then others in order of first use); a name the header declares in the other section keeps its id. The table is append-only and shared by every reader of the file, so a record read on one thread and a key resolved on another agree on the id. htslib adds an undeclared key to its header as `Type=String` too (with `Number=1`; `Str` reads either) but knows no reserved definitions; noodles applies the reserved ones. Without either, aardvark's undeclared `SVTYPE` and `TRID` silently turned breakends and tandem repeats into insertions. A BCF record names fields by dictionary index, so an index the dictionary lacks is a malformed record and MUST fail with a typed error; a BCF file's table adds only the names lenient keys resolve, at ids no record can name.

r[vcf_reader.undeclared+4]
A VCF text record that uses an INFO or FORMAT key the header does not declare — reserved or not, resolved through a lenient key or not — MUST log a single `warn!` per file, naming the first such key and the file's path. Resolving a key logs at most a `debug!`: a lenient key for a field no record carries is the common case (aardvark asks for `SVTYPE` and `TRID` on every input, and warned four times per ordinary run under a warn-on-resolve rule, BCF inputs included, which cannot hold such a key at all). A warning only on access never fires for the caller who could not resolve the key and so never reads it: aardvark's undeclared `SVTYPE` would otherwise vanish silently behind an `info_key_opt` that returns `None`.

r[vcf_reader.undeclared.filter]
A VCF text FILTER name the header does not declare MUST still fail the record (`FilterStatus::Failed`), so a record never passes for a filter the header forgot; no `FilterId` can name it, so `Failed` does not list it, and it counts towards the same single `warn!` per file. htslib instead adds the name to its header; seqair's header is shared and immutable. Dropping the name would turn a record whose only filter is undeclared into `Missing`, which `bcf_has_filter` and many callers count as passing.

## Limits and errors

r[vcf_reader.limits]
Before allocating, the reader MUST reject: a VCF text line longer than 256 MiB; a BCF record whose `l_shared + l_indiv` exceeds 256 MiB; more alleles than BCF's 16-bit `n_allele` allows; a genotype ploidy above 256; a BCF typed vector whose declared length runs past its section. Each is a distinct typed error (`r[io.fuzz.alloc_limits]`).

r[vcf_reader.errors+3]
Errors MUST be a `VcfReadError` enum (wrapping `BgzfError`, `VcfHeaderError`, the index errors) with typed variants only (`r[io.errors.typed_variants]`). Every per-record error MUST say where the record is, as far as the reader knows: its 1-based number in the stream when reading from the first record, its 1-based line for VCF text read from the start, its virtual offset when the input is BGZF, the file's path, and its `CHROM:POS` once those decoded — for VCF text whenever the first two columns parse, whatever is wrong with the rest of the line. A query or partition reads from where the index's chunks start, so a count of the records it read depends on the index's layout and points nowhere a user can find: a record error inside one MUST carry no number (`RecordPlace::number` is `None`), and its virtual offset and `CHROM:POS` say where it is. A malformed text value MUST carry the offending bytes (the first 64). An error's `Display` MUST NOT repeat its source: error reporters (anyhow, color_eyre) print the chain, and a message that embeds its cause prints it twice. A key error MUST say which section (`KeyKind`: INFO, FORMAT, FILTER) it looked in; an unknown sample MUST be one variant (`KeyError::UnknownSample`) whether a component or a user type resolved it; a BCF contig index outside the header's contig dictionary is `RecordError::ContigIndexOutOfRange`, distinct from `VcfReadError::UnknownContig { name }`.

r[vcf_reader.errors.open]
Everything that fails while opening a file — detecting its format, reading and parsing its header, reading and parsing its index, (re)opening it for a reader — MUST be one variant, `VcfReadError::Open { path, source: OpenError }`, with the file's path (`None` for bytes in memory and streams) and a typed `OpenError` saying what failed (I/O, BGZF, `VcfHeaderError`, the formats `r[vcf_reader.formats]` rejects, the index). A caller that wraps errors with the path showed it twice for I/O errors and not at all for header errors such as an unknown `Number`; with the path in every open error it wraps nothing. An index read from its own file names that file in its `OpenError` variant, as it is not the one in `path`. No error MUST display a path twice.

r[vcf_reader.debug]
`Debug` of a record component, of a record view and of `VcfFile` MUST NOT print the file's shared state (header, dictionaries, tables): it prints what the record holds — contig and positions, IDs, alleles, QUAL, FILTER names, the fields present with their values — so a record in a log line or a panic message stays short. `Locus` MUST implement `Display` as `CHROM:POS`, `CHROM:START-END` when it spans more than one position.

r[vcf_reader.fuzz]
There MUST be fuzz targets for BCF records, VCF text records, header parsing, and a full-stack target that opens a VCF or BCF from bytes and reads every field of every record through `FullRecord` (`r[io.fuzz.readers]`).

## Implementation

r[vcf_reader.reuse_io+2]
Decompression MUST go through the existing `BgzfStream` (sequential) and `RegionBuf` (queries); index parsing through `BamIndex`'s tabix parser and `CsiIndex`. A tabix index numbers only the sequences it lists, by name (the TBI header, or a CSI's tabix aux block), and those names MUST map to the file's contig table; a CSI without names (BCF) numbers the header's contig dictionary. The reader MUST NOT add a second BGZF or index implementation.

r[vcf_reader.features]
The `cram` cargo feature (default on) MUST gate the CRAM reader and its codec dependencies; without it the crate MUST build, the unified reader MUST refuse a CRAM file with a typed error, and every other format MUST behave as with it. The `log` feature MUST forward `tracing` events to the `log` crate (`tracing/log`), for applications that install a `log` logger and no `tracing` subscriber.

r[vcf_reader.text_scan]
Line, tab, `;`, `=`, `,` and `:` scanning in VCF text MUST be vectorised with `fearless_simd` (`r[io.simd_portable]`), building one bitmask of separators per vector and walking its set bits rather than calling `memchr` once per field: per-entry `memchr` measured 2–3× slower on real INFO columns, whose entries average 10–25 bytes. A scalar version MUST exist and property tests MUST check both produce identical splits (`r[io.platform_optimizations]`).

r[vcf_reader.passthrough_header]
`reader.header()` MUST return the file's `VcfHeader`, the type `Writer::write_header` takes, so a reader's header can be written unchanged.

## Testing

r[vcf_reader.oracles]
Reading MUST be checked against independent oracles:

- seqair's own writer, as a generator: random records written as VCF, VCF.gz and BCF MUST read back with the values that were written (`r[vcf_values.equivalence]`);
- bcftools: `bcftools view` of a file and seqair's reading of it MUST agree on every field of every record, for VCF text and BCF, on generated files and on real ones (GIAB, rastair output, a gnomAD slice, a multi-sample panel slice);
- noodles: the same comparison on VCF text, as a second implementation.
