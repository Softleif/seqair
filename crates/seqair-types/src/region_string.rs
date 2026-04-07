use crate::{One, Pos, SmolStr};
use std::fmt;
use winnow::{
    Parser,
    ascii::dec_uint,
    combinator::{opt, preceded},
    error::{StrContext, StrContextValue},
    token::{literal, take_till},
};

/// A struct representing a genomic region string.
///
/// # Examples
///
/// ```rust
/// use std::str::FromStr;
/// use seqair_types::RegionString;
///
/// let entire_chr17 = RegionString::from_str("chr17").unwrap();
/// let chr17_from_100 = RegionString::from_str("chr17:100").unwrap();
/// let chr17_from_100_to_200 = RegionString::from_str("chr17:100-200").unwrap();
/// ```
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[must_use]
pub struct RegionString {
    /// The chromosome name, includes the "chr" prefix.
    pub chromosome: SmolStr,
    /// The start position of the region, inclusive (1-based).
    pub start: Option<Pos<One>>,
    /// The end position of the region, inclusive (1-based).
    pub end: Option<Pos<One>>,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for RegionString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.chromosome)?;
        if let Some(start) = self.start {
            write!(f, ":{}", start.get())?;
        }
        if let Some(end) = self.end {
            write!(f, "-{}", end.get())?;
        }
        Ok(())
    }
}

impl std::str::FromStr for RegionString {
    type Err = RegionStringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err(RegionStringError::EmptyInput);
        }
        if !s.is_ascii() {
            return Err(RegionStringError::InvalidAscii);
        }
        let res = parser.parse(s).map_err(|_| RegionStringError::Malformed)?;
        if let RegionString { start: Some(start), end: Some(end), .. } = res
            && start > end
        {
            return Err(RegionStringError::StartGreaterThanEnd);
        }

        Ok(res)
    }
}

#[allow(clippy::type_complexity)]
fn parser(input: &mut &str) -> winnow::Result<RegionString> {
    fn chromosome_name<'d>(input: &mut &'d str) -> winnow::Result<&'d str> {
        take_till(1.., |c| c == ':')
            .context(StrContext::Expected(StrContextValue::Description("chromosome name")))
            .parse_next(input)
    }

    fn index(input: &mut &str) -> winnow::Result<Pos<One>> {
        dec_uint::<&str, u32, _>
            .verify(|v| *v > 0 && *v != u32::MAX)
            .map(|v| Pos::<One>::new(v).expect("BUG: verified v > 0 && v != u32::MAX"))
            .context(StrContext::Expected(StrContextValue::Description("number greater than 0")))
            .parse_next(input)
    }

    (
        chromosome_name,
        opt((
            preceded(literal(":"), index.context(StrContext::Label("start")))
                .context(StrContext::Label("start2")),
            opt(preceded(literal("-"), index.context(StrContext::Label("end")))
                .context(StrContext::Label("end2"))),
        )
            .context(StrContext::Label("range"))),
    )
        .context(StrContext::Label("range"))
        .map(|(token, slice): (&str, Option<(Pos<One>, Option<Pos<One>>)>)| match slice {
            None => RegionString { chromosome: token.into(), start: None, end: None },
            Some((start, None)) => {
                RegionString { chromosome: token.into(), start: Some(start), end: None }
            }
            Some((start, end)) => {
                RegionString { chromosome: token.into(), start: Some(start), end }
            }
        })
        .parse_next(input)
}

/// Errors that can occur when parsing a [`RegionString`].
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RegionStringError {
    /// The input string was empty.
    #[error("Empty region string")]
    EmptyInput,
    /// The input string contained non-ASCII characters.
    #[error("Invalid ASCII string")]
    InvalidAscii,
    /// The input could not be parsed as a valid region string.
    #[error("Malformed region string")]
    Malformed,
    /// The start position was greater than the end position.
    #[error("Start position cannot be greater than end position")]
    StartGreaterThanEnd,
}

#[cfg(feature = "hts-compat")]
mod hts {
    use super::*;
    use rust_htslib::bam::FetchDefinition;

    impl<'a> From<&'a RegionString> for FetchDefinition<'a> {
        fn from(region: &'a RegionString) -> Self {
            let chromosome = region.chromosome.as_bytes();
            match (region.start, region.end) {
                (Some(start), Some(end)) => FetchDefinition::from((
                    chromosome,
                    start.to_zero_based().get() as i64,
                    end.get() as i64,
                )),
                (Some(start), None) => FetchDefinition::from((
                    chromosome,
                    start.to_zero_based().get() as i64,
                    i64::MAX,
                )),
                (None, None) => FetchDefinition::from(chromosome),
                (None, Some(_)) => {
                    unreachable!("End position cannot be specified without a start position")
                }
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[cfg(feature = "hts-compat")]
    use rust_htslib::bam::FetchDefinition;

    #[test]
    fn test_valid_region_strings() {
        // full chromosome
        let region = RegionString::from_str("chr1").unwrap();
        assert_eq!(region, RegionString { chromosome: "chr1".into(), start: None, end: None });

        // chromosome with start position
        let region = RegionString::from_str("chr2:100").unwrap();
        assert_eq!(
            region,
            RegionString { chromosome: "chr2".into(), start: Pos::<One>::new(100), end: None }
        );

        // chromosome with start and end positions
        let region = RegionString::from_str("chr3:100-200").unwrap();
        assert_eq!(
            region,
            RegionString {
                chromosome: "chr3".into(),
                start: Pos::<One>::new(100),
                end: Pos::<One>::new(200)
            }
        );
    }

    #[test]
    fn test_error_cases() {
        // empty input
        let err = RegionString::from_str("").unwrap_err();
        assert!(matches!(err, RegionStringError::EmptyInput));

        // invalid chromosome
        let err = RegionString::from_str(":100").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // invalid start position
        let err = RegionString::from_str("chr1:invalid").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // invalid end position
        let err = RegionString::from_str("chr1:100-invalid").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // start greater than end
        let err = RegionString::from_str("chr1:200-100").unwrap_err();
        assert!(matches!(err, RegionStringError::StartGreaterThanEnd));
    }

    #[test]
    fn test_edge_cases() {
        // with whitespace
        let region = RegionString::from_str(" chr4:150-250 ").unwrap();
        assert_eq!(
            region,
            RegionString {
                chromosome: "chr4".into(),
                start: Pos::<One>::new(150),
                end: Pos::<One>::new(250)
            }
        );

        // only whitespace in chromosome part
        let err = RegionString::from_str("  :100").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // empty range part
        let err = RegionString::from_str("chr1:  ").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // invalid characters in start
        let err = RegionString::from_str("chr1:xxx").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // invalid characters in end
        let err = RegionString::from_str("chr1:100-xxx").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));

        // non-ascii characters
        let err = RegionString::from_str("chrü1:100-200").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidAscii));

        // mad scientist name
        let region = RegionString::from_str("bacteria_1-strain-6#5").unwrap();
        assert_eq!(
            region,
            RegionString { chromosome: "bacteria_1-strain-6#5".into(), start: None, end: None }
        );
    }

    #[test]
    fn test_display() {
        // Test full chromosome
        let region = RegionString { chromosome: "chr1".into(), start: None, end: None };
        insta::assert_snapshot!(region, @"chr1");

        // Test chromosome with start position
        let region =
            RegionString { chromosome: "chr2".into(), start: Pos::<One>::new(100), end: None };
        insta::assert_snapshot!(region, @"chr2:100");

        // Test chromosome with start and end positions
        let region = RegionString {
            chromosome: "chr3".into(),
            start: Pos::<One>::new(100),
            end: Pos::<One>::new(200),
        };
        insta::assert_snapshot!(region, @"chr3:100-200");
    }

    proptest::proptest! {
        #[test]
        fn proptest_roundtrip_region_string(
            // Generate chromosome names with "chr" prefix and some alphanumeric characters
            chrom in "chr[0-9A-Za-z]{1,10}",
            // Ensure start is between 1 and 1,000,000
            start in 1u32..1_000_000u32,
            // Ensure end is >= start and within reasonable bounds
            end_offset in 0u32..1_000_000u32
        ) {
            let end = start + end_offset;

            // chromosome only
            let region_str = chrom.clone();
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, None);
            assert_eq!(parsed.end, None);
            #[cfg(feature = "hts-compat")]
            let _ = FetchDefinition::from(&parsed);

            // chromosome with start
            let region_str = format!("{chrom}:{start}");
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, Pos::<One>::new(start));
            assert_eq!(parsed.end, None);
            #[cfg(feature = "hts-compat")]
            let _ = FetchDefinition::from(&parsed);

            // chromosome with start and end
            let region_str = format!("{chrom}:{start}-{end}");
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, Pos::<One>::new(start));
            assert_eq!(parsed.end, Pos::<One>::new(end));
            #[cfg(feature = "hts-compat")]
            let _ = FetchDefinition::from(&parsed);
        }

        #[test]
        fn proptest_roundtrip_random_string(
            // Generate random strings with up to 100 printable characters
            random_str in r"\PC{0,100}"
        ) {
            let Ok(parsed) = RegionString::from_str(&random_str) else {
                // We're just checking that there is no panic, but errors are fine!
                return Ok(());
            };
            let display = parsed.to_string();
            if let Some(start) = parsed.start {
                proptest::prop_assert!(display.contains(':'), "display missing ':': {display}");
                let after_colon = display.split(':').nth(1).unwrap_or("");
                let start_str = after_colon.split('-').next().unwrap_or("");
                let parsed_start: u32 = start_str.parse().expect("start in display must be numeric");
                proptest::prop_assert_eq!(parsed_start, start.get());
            }
            if let Some(end) = parsed.end {
                proptest::prop_assert!(display.contains('-'), "display missing '-': {display}");
                let after_dash = display.split('-').nth(1).unwrap_or("");
                let parsed_end: u32 = after_dash.parse().expect("end in display must be numeric");
                proptest::prop_assert_eq!(parsed_end, end.get());
            }
            // can be used as a FetchDefinition for bam
            #[cfg(feature = "hts-compat")]
            let _ = FetchDefinition::from(&parsed);
        }
    }

    #[test]
    fn test_chromosome_name_with_underscore() {
        // Test chromosome name with underscore
        let region = RegionString::from_str("bacteriophage_lambda_CpG:1-48502").unwrap();
        assert_eq!(
            region,
            RegionString {
                chromosome: "bacteriophage_lambda_CpG".into(),
                start: Pos::<One>::new(1),
                end: Pos::<One>::new(48502)
            }
        );
    }

    // r[verify io.errors.typed_variants]
    #[test]
    fn test_malformed_is_unit_variant() {
        let err = RegionString::from_str(":100").unwrap_err();
        assert!(matches!(err, RegionStringError::Malformed));
        // Display output has no dynamic content
        assert_eq!(err.to_string(), "Malformed region string");
    }

    // r[verify io.minimal_public_api]
    #[test]
    fn test_error_accessible_from_crate_root() {
        // Verify the error type is accessible via the crate re-export
        let err: crate::RegionStringError = RegionStringError::Malformed;
        assert!(matches!(err, crate::RegionStringError::Malformed));
    }
}
