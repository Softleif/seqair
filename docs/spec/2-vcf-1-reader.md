# VCF/BCF Reader

A pure-Rust reader for VCF (plain or BGZF-compressed) and BCF, built on seqair's existing BGZF, `RegionBuf`, TBI and CSI layers. It replaces rust-htslib in rastair's `convert`, `verify` and `train`, and is shaped by the access patterns of the other Rust tools that read VCFs (riker, HiPhase, aardvark): stream and extract, k-way merge across files, read–modify–write, and one rayon task per region.

> **Sources:** [VCF45] §1.6 "Data lines", §6 "BCF specification". [TABIX], [CSI] for region queries. [htslib] `vcf.c` (`get_rlen`, `bcf_has_filter`), `tbx.c` (VCF end column) for behaviour the specs leave open. See [References](./99-references.md).

What a reader decodes is chosen by a `Fields` type ([VCF/BCF Fields](./2-vcf-2-fields.md)); how values come out is in [VCF/BCF Values](./2-vcf-3-values.md).

## Opening

r[vcf_reader.formats]
The reader MUST read plain-text VCF, BGZF-compressed VCF and BCF 2.2 (BGZF-framed, at any compression level including 0). The format MUST be detected from the content — BGZF magic, then the decompressed `BCF\x02\x02` magic or `##fileformat=` — never from the file extension. A gzip stream that is not BGZF MUST be rejected with a typed error that suggests `bgzip`. BCF 2.1 (`BCF\x02\x01`) MUST be rejected with a typed error.

r[vcf_reader.file]
`VcfFile::open(path)` MUST parse the header (`r[vcf_header.parse]`) and locate an index by the names of `r[index.locate]` (VCF: `csi`, then `tbi`; BCF: `csi`). A missing index MUST NOT fail the open: sequential reading works without one, and only a query fails, with `VcfReadError::NoIndex`. An index that exists but does not parse MUST fail the open.

r[vcf_reader.file.shared+2]
`VcfFile` MUST be `Send + Sync + Clone`, with its parsed state — header, index, the field table (`r[vcf_reader.field_table]`) and its text-key lookup (`r[vcf_values.text_resolution+2]`), the contig table (`r[vcf_reader.contigs]`) and the warn-once state (`r[vcf_reader.undeclared+2]`) — behind one `Arc`. Cloning it MUST NOT re-read the file. This is what lets a rayon task query a file without reopening and re-parsing it (aardvark does that today because noodles' readers are not `Send`).

## Readers

r[vcf_reader.reader]
`file.reader::<F: Fields>()` MUST return a builder whose `open()` (or `?` on the builder's terminal method) plans `F` against the header (`r[vcf_fields.plan]`) and yields a `Reader<F>` positioned at the first record. A `Reader<F>` owns its I/O state — its own file handle, BGZF state and block cache — and shares everything else with its `VcfFile`. `Reader<F>` MUST be `Send`; it need not be `Sync`. Opening a reader MUST cost one file open and no header or index parse.

r[vcf_reader.reader.threads]
The builder MUST offer `threads(n)` and `threads_in(n, pool)` for sequential reading, passed to the existing `BgzfStream` (`r[bgzf.stream.threads]`, `r[bgzf.stream.shared_pool]`), so a caller can bound total decompression threads across many readers. Region queries decompress on the calling thread.

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

## Region queries

r[vcf_reader.query]
`file.query::<F>(region)` MUST open a reader (`r[vcf_reader.reader]`) and position it on the region; `reader.query(region)` MUST do the same on an existing reader, reusing its file handle and block cache, for callers that issue many small queries on one thread. A region is a contig (`ContigId` or name) and an optional inclusive `Pos1` range. The query MUST stream through the existing `RegionBuf` ([RegionBuf](./3-region_buf.md)) with TBI or CSI chunks, and MUST yield records in file order.

r[vcf_reader.query.overlap]
A record MUST be yielded iff its reference span `[start, end]` (`r[vcf_values.rlen]`) overlaps the queried range. Since every record that starts inside the range overlaps it, the reader MUST compute the end only for records that start before the range — for VCF text that is the only case that needs an INFO/FORMAT scan.

r[vcf_reader.query.contigs]
A contig the file declares but the index does not list (it has no records) MUST give an empty query, not an error. A contig the file does not know MUST fail with `VcfReadError::UnknownContig { name }`, a variant the caller can match to skip it. rastair and aardvark both carry workarounds for the first case today.

r[vcf_reader.query.many]
`query_many(regions)` MUST sort the regions, merge overlapping and adjacent ones per contig, and yield every record at most once, in file order. Visiting regions one by one returns a record that spans two adjacent regions twice; `query_many` is the way to avoid that.

## Partitions

r[vcf_reader.partitions]
`file.partitions(n)` MUST tile every contig the index lists into about `n` disjoint `Partition`s with boundaries taken from the index, roughly equal in compressed bytes. Partitions are `Send + Sync + Clone` and carry their file's identity.

r[vcf_reader.partitions.ownership]
`file.query_partition::<F>(&partition)` MUST yield exactly the records whose start lies inside the partition — not those that merely overlap it — so a record spanning a boundary belongs to one partition only. HiPhase guards against the double visit by hand ("we have already written though, so don't write it again"); partition reads make that guard unnecessary.

r[vcf_reader.partitions.complete]
The records of all partitions of a file, concatenated in partition order, MUST equal one sequential read of the file's records on indexed contigs. A property test MUST check this — no record missing, none duplicated, order preserved — against a sequential read as the oracle, for random files and random `n`.

## Undeclared fields

r[vcf_reader.field_table]
A VCF text file's INFO and FORMAT keys MUST resolve against a field table: the header's string dictionary, extended with

- every VCF 4.5 reserved INFO and FORMAT key (the specification's tables of reserved INFO keys and genotype keys, and its structural-variant INFO and FORMAT keys; the `M[0-9]+[ACGTUN]` patterns excepted) that the header does not declare in that section, with its reserved `Type` and `Number`;
- every other key a record uses or a caller resolves, as `Type=String, Number=.`.

Keys new to the dictionary get ids past the header's (the reserved ones first, then others in order of first use); a name the header declares in the other section keeps its id. The table is append-only and shared by every reader of the file, so a record read on one thread and a key resolved on another agree on the id. htslib adds an undeclared key to its header as `Type=String` too (with `Number=1`; `Str` reads either) but knows no reserved definitions; noodles applies the reserved ones. Without either, aardvark's undeclared `SVTYPE` and `TRID` silently turned breakends and tandem repeats into insertions. A BCF file's table is its header's dictionary alone: a BCF record names fields by dictionary index, so an index the dictionary lacks is a malformed record and MUST fail with a typed error.

r[vcf_reader.undeclared+2]
Reading an undeclared key MUST be visible in the log even when no record is ever indexed: resolving a reserved key the header does not declare logs a `debug!`, resolving any other undeclared key logs one `warn!` per name ("not declared in the header, read as String"), and a record that uses an undeclared key that VCF does not reserve and nobody resolved logs a single `warn!` per file naming the first such key. A warning only on access never fires for the caller who could not resolve the key and so never reads it.

r[vcf_reader.undeclared.filter]
A VCF text FILTER name the header does not declare MUST still fail the record (`FilterStatus::Failed`), so a record never passes for a filter the header forgot; no `FilterId` can name it, so `Failed` does not list it, and it counts towards the same single `warn!` per file. htslib instead adds the name to its header; seqair's header is shared and immutable. Dropping the name would turn a record whose only filter is undeclared into `Missing`, which `bcf_has_filter` and many callers count as passing.

## Limits and errors

r[vcf_reader.limits]
Before allocating, the reader MUST reject: a VCF text line longer than 256 MiB; a BCF record whose `l_shared + l_indiv` exceeds 256 MiB; more alleles than BCF's 16-bit `n_allele` allows; a genotype ploidy above 256; a BCF typed vector whose declared length runs past its section. Each is a distinct typed error (`r[io.fuzz.alloc_limits]`).

r[vcf_reader.errors+2]
Errors MUST be a `VcfReadError` enum (wrapping `BgzfError`, `VcfHeaderError`, the index errors) with typed variants only (`r[io.errors.typed_variants]`). Every per-record error MUST say where the record is, as far as the reader knows: its 1-based number in the stream, its 1-based line for VCF text read from the start, its virtual offset when the input is BGZF, the file's path, and its `CHROM:POS` once those decoded. A malformed text value MUST carry the offending bytes (the first 64). An error's `Display` MUST NOT repeat its source: error reporters (anyhow, color_eyre) print the chain, and a message that embeds its cause prints it twice. A key error MUST say which section (`KeyKind`: INFO, FORMAT, FILTER) it looked in; an unknown sample MUST be one variant (`KeyError::UnknownSample`) whether a component or a user type resolved it; a BCF contig index outside the header's contig dictionary is `RecordError::ContigIndexOutOfRange`, distinct from `VcfReadError::UnknownContig { name }`.

r[vcf_reader.debug]
`Debug` of a record component, of a record view and of `VcfFile` MUST NOT print the file's shared state (header, dictionaries, tables): it prints what the record holds — contig and positions, IDs, alleles, QUAL, FILTER names, the fields present with their values — so a record in a log line or a panic message stays short. `Locus` MUST implement `Display` as `CHROM:POS`, `CHROM:START-END` when it spans more than one position.

r[vcf_reader.fuzz]
There MUST be fuzz targets for BCF records, VCF text records, header parsing, and a full-stack target that opens a VCF or BCF from bytes and reads every field of every record through `FullRecord` (`r[io.fuzz.readers]`).

## Implementation

r[vcf_reader.reuse_io]
Decompression MUST go through the existing `BgzfStream` (sequential) and `RegionBuf` (queries); index parsing through `BamIndex::from_tabix_*` (generalised from a `BamHeader` to a contig-name lookup) and `CsiIndex`. The reader MUST NOT add a second BGZF or index implementation.

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
