//! Pure-Rust BAM/SAM/CRAM/FASTA reader, pileup engine, and VCF/BCF writer.
//!
//! seqair provides indexed random access to alignment files (BAM, SAM, CRAM),
//! reference sequences (FASTA), and a column-based pileup engine. It also
//! includes VCF/BCF writers with type-safe allele representation and single-pass
//! index co-production.
//!
//! # Quick start: pileup at a genomic region
//!
//! The most common workflow is opening an alignment + reference, fetching a
//! region, and iterating pileup columns:
//!
//! ```no_run
//! use seqair::reader::Readers;
//! use seqair::bam::pileup::PileupOp;
//! use seqair_types::Pos0;
//! use std::path::Path;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Open BAM + FASTA together — auto-detects format (BAM/SAM/CRAM)
//! let mut readers = Readers::open(Path::new("sample.bam"), Path::new("reference.fa"))?;
//!
//! // Fetch a region and get a pileup iterator
//! let tid = readers.header().tid("chr19").expect("contig not found");
//! let start = Pos0::new(6_100_000).unwrap();
//! let end = Pos0::new(6_200_000).unwrap();
//! let mut pileup = readers.pileup(tid, start, end)?;
//!
//! // Iterate over columns (one per reference position)
//! for column in pileup.by_ref() {
//!     let _pos = column.pos();
//!     let _ref_base = column.reference_base();
//!     let _depth = column.depth();
//!
//!     // Inspect each read's contribution at this position
//!     for aln in column.alignments() {
//!         match &aln.op {
//!             PileupOp::Match { base, qual, .. } => {
//!                 // Read shows `base` with quality `qual`
//!                 let _ = (base, qual);
//!             }
//!             PileupOp::Deletion { del_len } => {
//!                 // Read has a deletion of `del_len` bases spanning this position
//!                 let _ = del_len;
//!             }
//!             PileupOp::Insertion { base, qual, insert_len, .. } => {
//!                 // Read has an insertion after the previous match
//!                 let _ = (base, qual, insert_len);
//!             }
//!             PileupOp::RefSkip => {
//!                 // Spliced alignment skips this position (RNA-seq)
//!             }
//!         }
//!     }
//! }
//!
//! // Recover the RecordStore for reuse in the next region
//! readers.recover_store(&mut pileup);
//! # Ok(())
//! # }
//! ```
//!
//! # Multi-threaded pileup
//!
//! Both alignment and FASTA readers support [`fork()`](reader::Readers::fork)
//! for cheap per-thread copies that share the parsed index and header via `Arc`:
//!
//! ```no_run
//! # use seqair::reader::Readers;
//! # use seqair_types::Pos0;
//! # use std::path::Path;
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let readers = Readers::open(Path::new("sample.bam"), Path::new("ref.fa"))?;
//!
//! // Each thread forks its own reader — no locking, shared index
//! std::thread::scope(|s| {
//!     for _ in 0..4 {
//!         let mut forked = readers.fork().unwrap();
//!         s.spawn(move || {
//!             let pileup = forked.pileup(0,
//!                 Pos0::new(0).unwrap(),
//!                 Pos0::new(100_000).unwrap(),
//!             ).unwrap();
//!             for _column in pileup {
//!                 // process in parallel...
//!             }
//!         });
//!     }
//! });
//! # Ok(())
//! # }
//! ```
//!
//! # Writing VCF/BCF from pileup results
//!
//! After computing variant calls from pileup data, write them to VCF or BCF.
//! Use [`VcfHeaderBuilder::from_bam_header`](vcf::VcfHeaderBuilder::from_bam_header)
//! to copy contig names and lengths from the alignment header:
//!
//! ```
//! use seqair::vcf::{
//!     Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
//!     FormatGt, FormatInt, InfoInt,
//! };
//! use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
//! use seqair_types::{Base, Pos1};
//! use std::sync::Arc;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // 1. Build a header with typed field keys resolved at setup time.
//! //    In a real app, use VcfHeaderBuilder::from_bam_header() to copy contigs.
//! let mut builder = VcfHeader::builder();
//! let chr1 = builder.register_contig("chr1", ContigDef { length: Some(248_956_422) })?;
//! let dp_info: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Total read depth")
//! )?;
//! let gt_fmt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let dp_fmt: FormatInt = builder.register_format(
//!     &FormatFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Sample read depth")
//! )?;
//! let header = Arc::new(builder.add_sample("sample1")?.build()?);
//!
//! // 2. Create a writer (VCF text to buffer, or BCF to file)
//! let mut output = Vec::new();
//! let writer = Writer::new(&mut output, OutputFormat::Vcf);
//! let mut writer = writer.write_header(&header)?;
//!
//! // 3. Encode records with type-safe alleles and pre-resolved field keys
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let pos = Pos1::new(12345).unwrap();
//! let enc = writer.begin_record(&chr1, pos, &alleles, Some(30.0))?;
//! let mut enc = enc.filter_pass();
//! dp_info.encode(&mut enc, 50);
//! let mut enc = enc.begin_samples(1);
//! gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
//! dp_fmt.encode(&mut enc, &[45]);
//! enc.emit()?;
//!
//! writer.finish()?;
//!
//! let vcf_text = String::from_utf8(output)?;
//! assert!(vcf_text.contains("chr1\t12345\t.\tA\tT\t30\tPASS\tDP=50"));
//! # Ok(())
//! # }
//! ```
//!
//! # Type-safe alleles
//!
//! The [`Alleles`](vcf::Alleles) enum prevents invalid REF/ALT combinations
//! at construction time:
//!
//! ```
//! use seqair::vcf::alleles::Alleles;
//! use seqair_types::Base;
//!
//! // SNV: single ref base, single alt base
//! let snv = Alleles::snv(Base::A, Base::T).unwrap();
//! assert_eq!(snv.ref_text(), "A");
//! assert_eq!(snv.rlen(), 1);
//!
//! // Insertion: anchor + inserted bases → REF=A, ALT=ACGT
//! let ins = Alleles::insertion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
//! assert_eq!(ins.ref_text(), "A");
//!
//! // Deletion: anchor + deleted bases → REF=ACGT, ALT=A
//! let del = Alleles::deletion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
//! assert_eq!(del.ref_text(), "ACGT");
//! assert_eq!(del.rlen(), 4);
//!
//! // Can't create SNV where ref == alt
//! assert!(Alleles::snv(Base::A, Base::A).is_err());
//! // Can't create empty insertion
//! assert!(Alleles::insertion(Base::A, &[]).is_err());
//! ```
//!
//! # High-performance BCF encoding
//!
//! For hot paths (millions of records), use [`Writer`](vcf::Writer) with
//! pre-resolved typed keys from [`register_info`](vcf::VcfHeaderBuilder::register_info)
//! and [`register_format`](vcf::VcfHeaderBuilder::register_format) for zero-allocation
//! direct BCF encoding:
//!
//! ```
//! use seqair::vcf::{
//!     Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
//!     FormatGt, InfoInt,
//! };
//! use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
//! use seqair_types::{Base, Pos1};
//! use std::sync::Arc;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Resolve keys once at setup — no string lookups per record
//! let mut builder = VcfHeader::builder();
//! let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) })?;
//! let dp: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Depth")
//! )?;
//! let gt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let header = Arc::new(builder.add_sample("sample1")?.build()?);
//!
//! let mut buf = Vec::new();
//! let writer = Writer::new(&mut buf, OutputFormat::Bcf);
//! let mut writer = writer.write_header(&header)?;
//!
//! // Hot loop: zero allocations, pre-resolved dict indices
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let enc = writer.begin_record(&contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))?;
//! let mut enc = enc.filter_pass();
//! dp.encode(&mut enc, 50);
//! let mut enc = enc.begin_samples(1);
//! gt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
//! enc.emit()?;
//!
//! writer.finish()?;
//! assert!(!buf.is_empty()); // valid BCF bytes
//! # Ok(())
//! # }
//! ```
//!
//! # Reading FASTA reference sequences
//!
//! ```no_run
//! use seqair::fasta::IndexedFastaReader;
//! use seqair_types::{Base, Pos0};
//! use std::path::Path;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut reader = IndexedFastaReader::open(Path::new("reference.fa"))?;
//!
//! // Fetch a region as raw ASCII bytes
//! let seq = reader.fetch_seq("chr1", Pos0::new(1000).unwrap(), Pos0::new(2000).unwrap())?;
//!
//! // Convert to Base enum for downstream use
//! let bases = Base::from_ascii_vec(seq);
//! // bases is now Vec<Base> with A/C/G/T/Unknown variants
//! let _ = bases;
//! # Ok(())
//! # }
//! ```

pub mod bam;
pub mod cram;
pub mod fasta;
pub mod reader;
pub mod sam;
pub(crate) mod utils;
pub mod vcf;

pub use reader::{IndexedReader, Readers};
