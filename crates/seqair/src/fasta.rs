//! Reference sequence access. Open a FASTA file with [`IndexedFastaReader`] to fetch
//! subsequences by name and coordinate range. Supports plain and bgzf-compressed FASTA.
//!
//! An indexed **FASTQ** reference works through the same reader. A FASTQ index is
//! not a separate format — `samtools faidx` and `samtools fqidx` both write
//! `<file>.fai` with the five FASTA columns plus `qual_offset` — and the
//! sequence-block geometry means the same thing in both, so [`FaiEntry`] simply
//! carries the extra column as [`FaiEntry::qual_offset`]. Because a fetch is
//! bounded by the record's `length`, it never reaches the `+` separator or the
//! quality block, which makes coordinate queries immune to the `@`/`+`-in-quality
//! ambiguity that defeats line-oriented FASTQ parsers.

mod gzi;
mod index;
mod reader;
mod strip;

pub use gzi::{BlockLocation, GziError, GziIndex};
pub use index::{FaiEntry, FaiEntryError, FaiError, FastaIndex};
pub use reader::{FastaError, IndexedFastaReader};
