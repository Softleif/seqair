//! The INFO and FORMAT keys VCF 4.5 reserves, with their definitions, for
//! the keys a VCF text header leaves out (`r[vcf_reader.field_table]`):
//! the tables of reserved INFO keys and genotype keys (§1.6.1), and the
//! structural-variant INFO and FORMAT keys (§3, §4). The base-modification
//! patterns (`M[0-9]+[ACGTUN]` and their `DPM`/`ADM` forms) are left out:
//! a table cannot hold a pattern.

use crate::vcf::header::{FormatDef, InfoDef, Number, ValueType};
use seqair_types::SmolStr;

macro_rules! defs {
    ($def:ident: $(($name:literal, $number:expr, $typ:ident, $description:literal)),* $(,)?) => {
        [$((
            $name,
            $def {
                number: $number,
                typ: ValueType::$typ,
                description: SmolStr::new_static($description),
                attributes: Vec::new(),
            },
        )),*]
    };
}

use Number::{
    AlternateBases as A, BaseModification as M, Count, Genotypes as G, LocalAlternateBases as LA,
    LocalGenotypes as LG, LocalReferenceAlternateBases as LR, Ploidy as P,
    ReferenceAlternateBases as R, Unknown as Dot,
};

const ONE: Number = Count(1);

/// Reserved INFO keys: the reserved table, then structural variants and
/// tandem repeats.
pub(crate) static INFO: &[(&str, InfoDef)] = &defs!(InfoDef:
    ("AA", ONE, String, "Ancestral allele"),
    ("AC", A, Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"),
    ("AD", R, Integer, "Total read depth for each allele"),
    ("ADF", R, Integer, "Read depth for each allele on the forward strand"),
    ("ADR", R, Integer, "Read depth for each allele on the reverse strand"),
    ("AF", A, Float, "Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)"),
    ("AN", ONE, Integer, "Total number of alleles in called genotypes"),
    ("BQ", ONE, Float, "RMS base quality"),
    ("CIGAR", A, String, "Cigar string describing how to align an alternate allele to the reference allele"),
    ("DB", Count(0), Flag, "dbSNP membership"),
    ("DP", ONE, Integer, "Combined depth across samples"),
    ("END", ONE, Integer, "Deprecated. Present for backwards compatibility with earlier versions of VCF."),
    ("H2", Count(0), Flag, "HapMap2 membership"),
    ("H3", Count(0), Flag, "HapMap3 membership"),
    ("MQ", ONE, Float, "RMS mapping quality"),
    ("MQ0", ONE, Integer, "Number of MAPQ == 0 reads"),
    ("NS", ONE, Integer, "Number of samples with data"),
    ("SB", Count(4), Integer, "Strand bias"),
    ("SOMATIC", Count(0), Flag, "Somatic mutation (for cancer genomics)"),
    ("VALIDATED", Count(0), Flag, "Validated by follow-up experiment"),
    ("1000G", Count(0), Flag, "1000 Genomes membership"),
    ("IMPRECISE", Count(0), Flag, "Imprecise structural variation"),
    ("NOVEL", Count(0), Flag, "Indicates a novel structural variation"),
    ("SVTYPE", ONE, String, "Type of structural variant"),
    ("SVLEN", A, Integer, "Length of structural variant"),
    ("CIPOS", Dot, Integer, "Confidence interval around POS for symbolic structural variants"),
    ("CIEND", Dot, Integer, "Confidence interval around the inferred END for symbolic structural variants"),
    ("HOMLEN", A, Integer, "Length of base pair identical micro-homology at breakpoints"),
    ("HOMSEQ", A, String, "Sequence of base pair identical micro-homology at breakpoints"),
    ("BKPTID", A, String, "ID of the assembled alternate allele in the assembly file"),
    ("MEINFO", Dot, String, "Mobile element info of the form NAME,START,END,POLARITY"),
    ("METRANS", Dot, String, "Mobile element transduction info of the form CHR,START,END,POLARITY"),
    ("DGVID", A, String, "ID of this element in Database of Genomic Variation"),
    ("DBVARID", A, String, "ID of this element in DBVAR"),
    ("DBRIPID", A, String, "ID of this element in DBRIP"),
    ("MATEID", A, String, "ID of mate breakend"),
    ("PARID", A, String, "ID of partner breakend"),
    ("EVENT", A, String, "ID of associated event"),
    ("EVENTTYPE", A, String, "Type of associated event"),
    ("CILEN", Dot, Integer, "Confidence interval for the SVLEN field"),
    ("CN", A, Float, "Copy number of CNV/breakpoint"),
    ("CICN", Dot, Float, "Confidence interval around copy number"),
    ("SVCLAIM", A, String, "Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively"),
    ("RN", A, Integer, "Total number of repeat sequences in this allele"),
    ("RUS", Dot, String, "Repeat unit sequence of the corresponding repeat sequence"),
    ("RUL", Dot, Integer, "Repeat unit length of the corresponding repeat sequence"),
    ("RUC", Dot, Float, "Repeat unit count of corresponding repeat sequence"),
    ("RB", Dot, Integer, "Total number of bases in the corresponding repeat sequence"),
    ("CIRUC", Dot, Float, "Confidence interval around RUC"),
    ("CIRB", Dot, Integer, "Confidence interval around RB"),
    ("RUB", Dot, Integer, "Number of bases in each individual repeat unit"),
);

/// Reserved FORMAT keys: the genotype table, then structural variants.
pub(crate) static FORMAT: &[(&str, FormatDef)] = &defs!(FormatDef:
    ("AD", R, Integer, "Read depth for each allele"),
    ("ADF", R, Integer, "Read depth for each allele on the forward strand"),
    ("ADR", R, Integer, "Read depth for each allele on the reverse strand"),
    ("DP", ONE, Integer, "Read depth"),
    ("EC", A, Integer, "Expected alternate allele counts"),
    ("LEN", ONE, Integer, "Length of <*> reference block"),
    ("FT", ONE, String, "Filter indicating if this genotype was \"called\""),
    ("GL", G, Float, "Genotype likelihoods"),
    ("GP", G, Float, "Genotype posterior probabilities"),
    ("GQ", ONE, Integer, "Conditional genotype quality"),
    ("GT", ONE, String, "Genotype"),
    ("HQ", Count(2), Integer, "Haplotype quality"),
    ("LA", Dot, Integer, "Reserved"),
    ("LAA", Dot, Integer, "1-based indices into ALT, indicating which alleles are relevant (local) for the current sample"),
    ("LAD", LR, Integer, "Local-allele representation of AD"),
    ("LADF", LR, Integer, "Local-allele representation of ADF"),
    ("LADR", LR, Integer, "Local-allele representation of ADR"),
    ("LEC", LA, Integer, "Local-allele representation of EC"),
    ("LGL", LG, Float, "Local-allele representation of GL"),
    ("LGP", LG, Float, "Local-allele representation of GP"),
    ("LPL", LG, Integer, "Local-allele representation of PL"),
    ("LPP", LG, Integer, "Local-allele representation of PP"),
    ("M5mC", M, Float, "Alias for M27551C 5-Methylcytosine"),
    ("DPM5mC", M, Integer, "Alias for DPM27551C"),
    ("ADM5mC", M, Integer, "Alias for ADM27551C"),
    ("M5hmC", M, Float, "Alias for M76792C 5-Hydroxymethylcytosine"),
    ("DPM5hmC", M, Integer, "Alias for DPM76792C"),
    ("ADM5hmC", M, Integer, "Alias for ADM76792C"),
    ("M5fC", M, Float, "Alias for M76794C 5-Formylcytosine"),
    ("DPM5fC", M, Integer, "Alias for DPM76794C"),
    ("ADM5fC", M, Integer, "Alias for ADM76794C"),
    ("M5caC", M, Float, "Alias for M76793C 5-Carboxylcytosine"),
    ("DPM5caC", M, Integer, "Alias for DPM76793C"),
    ("ADM5caC", M, Integer, "Alias for ADM76793C"),
    ("M5hmU", M, Float, "Alias for M16964T 5-Hydroxymethyluracil"),
    ("DPM5hmU", M, Integer, "Alias for DPM16964T"),
    ("ADM5hmU", M, Integer, "Alias for ADM16964T"),
    ("M5fU", M, Float, "Alias for M80961T 5-Formyluracil"),
    ("DPM5fU", M, Integer, "Alias for DPM80961T"),
    ("ADM5fU", M, Integer, "Alias for ADM80961T"),
    ("M5caU", M, Float, "Alias for M17477T 5-Carboxyluracil"),
    ("DPM5caU", M, Integer, "Alias for DPM17477T"),
    ("ADM5caU", M, Integer, "Alias for ADM17477T"),
    ("M6mA", M, Float, "Alias for M28871A 6-Methyladenine"),
    ("DPM6mA", M, Integer, "Alias for DPM28871A"),
    ("ADM6mA", M, Integer, "Alias for ADM28871A"),
    ("M8oxoG", M, Float, "Alias for M44605G 8-Oxoguanine"),
    ("DPM8oxoG", M, Integer, "Alias for DPM44605G"),
    ("ADM8oxoG", M, Integer, "Alias for ADM44605G"),
    ("MXaoN", M, Float, "Alias for M18107N Xanthosine"),
    ("DPMXaoN", M, Integer, "Alias for DPM18107N"),
    ("ADMXaoN", M, Integer, "Alias for ADM18107N"),
    ("MQ", ONE, Integer, "RMS mapping quality"),
    ("PL", G, Integer, "Phred-scaled genotype likelihoods rounded to the closest integer"),
    ("PP", G, Integer, "Phred-scaled genotype posterior probabilities rounded to the closest integer"),
    ("PQ", ONE, Integer, "Phasing quality"),
    ("PS", ONE, Integer, "Phase set"),
    ("PSL", P, String, "Phase set list"),
    ("PSO", P, Integer, "Phase set list ordinal"),
    ("PSQ", P, Integer, "Phase set list quality"),
    ("CN", ONE, Float, "Copy number"),
    ("CICN", Count(2), Float, "Confidence interval around copy number"),
    ("CNQ", ONE, Float, "Copy number genotype quality"),
    ("CNL", G, Float, "Copy number genotype likelihood"),
    ("CNP", G, Float, "Copy number posterior probabilities"),
    ("NQ", ONE, Integer, "Phred style probability score that the variant is novel"),
    ("HAP", ONE, Integer, "Unique haplotype identifier"),
    ("AHAP", ONE, Integer, "Unique identifier of ancestral haplotype"),
);

/// The definition an undeclared key that is not reserved is read with, as
/// htslib reads one (it says `Number=1`; [`Str`](crate::vcf::Str) reads either).
pub(crate) static UNDECLARED_INFO: InfoDef = InfoDef {
    number: Dot,
    typ: ValueType::String,
    description: SmolStr::new_static("Not declared in the header"),
    attributes: Vec::new(),
};

/// As [`UNDECLARED_INFO`], for FORMAT.
pub(crate) static UNDECLARED_FORMAT: FormatDef = FormatDef {
    number: Dot,
    typ: ValueType::String,
    description: SmolStr::new_static("Not declared in the header"),
    attributes: Vec::new(),
};

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// Each table names a key once, and the keys the reader depends on
    /// (END, SVLEN and LEN for the reference length; GT) are as VCF 4.5
    /// defines them.
    #[test]
    fn the_tables_hold_each_key_once() {
        for table in [
            INFO.iter().map(|(n, _)| *n).collect::<Vec<_>>(),
            FORMAT.iter().map(|(n, _)| *n).collect(),
        ] {
            assert_eq!(table.iter().collect::<HashSet<_>>().len(), table.len());
        }
        let info = |name| INFO.iter().find(|(n, _)| *n == name).map(|(_, d)| (d.number, d.typ));
        let format = |name| FORMAT.iter().find(|(n, _)| *n == name).map(|(_, d)| (d.number, d.typ));
        assert_eq!(info("END"), Some((ONE, ValueType::Integer)));
        assert_eq!(info("SVLEN"), Some((A, ValueType::Integer)));
        assert_eq!(info("SVTYPE"), Some((ONE, ValueType::String)));
        assert_eq!(format("LEN"), Some((ONE, ValueType::Integer)));
        assert_eq!(format("GT"), Some((ONE, ValueType::String)));
        assert_eq!(format("PL"), Some((G, ValueType::Integer)));
        assert_eq!(info("TRID"), None);
    }
}
