//! `VcfHeader::parse`: VCF header text into the model the writer builds.
//!
//! Follows htslib's `bcf_hdr_parse`: the dictionary is assigned as
//! `bcf_hdr_register_hrec` assigns it, and lines htslib only warns about
//! are accepted with a `warn!`.

use super::{
    Attribute, ContigDef, FilterDef, FormatDef, InfoDef, Line, LineKind, MAX_CONTIGS,
    MAX_DICT_ENTRIES, MAX_HEADER_TEXT, MAX_SAMPLES, Number, StringMap, ValueType, VcfHeader,
    VcfVersion,
};
use crate::vcf::error::{MetaLineError, VcfHeaderError};
use seqair_types::SmolStr;
use std::collections::HashSet;
use tracing::warn;

const FIXED_COLUMNS: [&str; 8] = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

impl VcfHeader {
    // r[impl vcf_header.parse]
    /// Parse VCF header text: the `##` meta-information lines and the
    /// `#CHROM` line, as at the start of a VCF file or in a BCF header's
    /// `l_text` block (trailing NULs are ignored).
    ///
    /// The header keeps the version its `##fileformat` line declares, its
    /// definitions in file order, the dictionary htslib would assign
    /// (`IDX=` honoured), and every line it does not interpret, verbatim.
    pub fn parse(text: &str) -> Result<Self, VcfHeaderError> {
        // r[impl vcf_header.parse.limits]
        if text.len() > MAX_HEADER_TEXT {
            return Err(VcfHeaderError::TextTooLarge { len: text.len(), max: MAX_HEADER_TEXT });
        }
        let mut lines = text
            .trim_end_matches('\0')
            .split('\n')
            .map(|line| line.strip_suffix('\r').unwrap_or(line))
            .enumerate()
            .map(|(i, line)| (i.saturating_add(1), line));

        let version = match lines.next() {
            Some((_, first)) => match first.strip_prefix("##fileformat=") {
                Some(value) => parse_version(value)?,
                None => return Err(VcfHeaderError::MissingFileFormat),
            },
            None => return Err(VcfHeaderError::MissingFileFormat),
        };

        let mut header = Self::empty(version);
        // r[impl vcf_header.line_order]
        header.canonical = false;
        let mut parser = Parser { header, pass_declared: false };
        while let Some((line_no, line)) = lines.next() {
            if line.starts_with("#CHROM") {
                parser.header.samples = parse_samples(line)?;
                if let Some((line, _)) = lines.find(|(_, rest)| !rest.is_empty()) {
                    return Err(VcfHeaderError::TrailingData { line });
                }
                return Ok(parser.header);
            }
            match line.strip_prefix("##") {
                Some(meta) => parser.meta_line(line_no, meta)?,
                None => {
                    warn!(line = line_no, "skipping a header line that is neither ## nor #CHROM")
                }
            }
        }
        Err(VcfHeaderError::MissingHeaderLine)
    }
}

// r[impl vcf_header.parse.version]
fn parse_version(value: &str) -> Result<VcfVersion, VcfHeaderError> {
    let malformed = || VcfHeaderError::MalformedVersion { value: SmolStr::from(value) };
    let (major, minor) =
        value.strip_prefix("VCFv").and_then(|v| v.split_once('.')).ok_or_else(malformed)?;
    let version = VcfVersion::new(
        parse_digits(major).ok_or_else(malformed)?,
        parse_digits(minor).ok_or_else(malformed)?,
    );
    if version.major != 4 {
        return Err(VcfHeaderError::UnsupportedVersion { version });
    }
    if version > VcfVersion::V4_5 {
        warn!(%version, "newer than VCFv4.5; reading it with the 4.5 rules");
    }
    Ok(version)
}

/// Decimal digits only: no sign, no space, not empty.
fn parse_digits<T: std::str::FromStr>(s: &str) -> Option<T> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse().ok()
}

struct Parser {
    header: VcfHeader,
    /// Whether the file declared PASS itself (the header starts with a default).
    pass_declared: bool,
}

impl Parser {
    /// One `##` line, without the `##`.
    fn meta_line(&mut self, line: usize, meta: &str) -> Result<(), VcfHeaderError> {
        let structured = meta.split_once('=').and_then(|(key, value)| {
            let body = value.strip_prefix('<')?;
            Some((key, body))
        });
        match structured {
            Some((key @ ("INFO" | "FORMAT" | "FILTER" | "contig"), body)) => {
                let attributes = parse_structured(body)
                    .map_err(|reason| VcfHeaderError::MalformedLine { line, reason })?;
                // r[impl vcf_header.line_form]
                let form = form(&attributes, |key| key != "IDX").into_boxed_str();
                let attributes = attributes.into_iter().map(|(attr, _)| attr).collect();
                self.definition(line, key, attributes, form)
            }
            // r[impl vcf_header.parse.preserve+2]
            // A second ##fileformat included: htslib keeps it as a line.
            _ => {
                if let Some(reason) = self.header.push_other_line(meta) {
                    warn!(line, %reason, "keeping a malformed structured header line verbatim");
                }
                Ok(())
            }
        }
    }

    fn definition(
        &mut self,
        line: usize,
        key: &str,
        mut attributes: Vec<Attribute>,
        form: Box<str>,
    ) -> Result<(), VcfHeaderError> {
        let id = take(&mut attributes, "ID")
            .ok_or_else(|| VcfHeaderError::MissingId { line, key: SmolStr::from(key) })?;
        let idx = take(&mut attributes, "IDX")
            .map(|value| {
                parse_digits::<u32>(&value)
                    .ok_or_else(|| VcfHeaderError::InvalidIdx { id: id.clone(), value })
            })
            .transpose()?;
        let header = &mut self.header;
        match key {
            "contig" => {
                if header.contigs.contains_key(&id) {
                    warn!(line, contig = %id, "duplicate contig definition; keeping the first");
                    return Ok(());
                }
                let length = take(&mut attributes, "length")
                    .map(|value| {
                        parse_digits::<u64>(&value).ok_or_else(|| {
                            VcfHeaderError::InvalidContigLength { name: id.clone(), value }
                        })
                    })
                    .transpose()?;
                assign(&mut header.contig_map, &id, idx, MAX_CONTIGS, || {
                    VcfHeaderError::TooManyContigs { max: MAX_CONTIGS }
                })?;
                let def_line =
                    Line::Def { kind: LineKind::Contig, id: id.clone(), form: Some(form) };
                header.lines.push(def_line);
                header.contigs.insert(id, ContigDef { length, attributes });
            }
            "FILTER" => {
                let is_pass = id == "PASS";
                let duplicate =
                    if is_pass { self.pass_declared } else { header.filters.contains_key(&id) };
                if duplicate {
                    warn!(line, filter = %id, "duplicate FILTER definition; keeping the first");
                    return Ok(());
                }
                let description = description(&mut attributes, line, &id);
                assign_field(&mut header.string_map, &id, idx)?;
                self.pass_declared |= is_pass;
                // r[impl vcf_header.line_order]
                // The file's PASS takes the default's place, after ##fileformat.
                match header.lines.iter_mut().find(|l| l.is_def(LineKind::Filter, "PASS")) {
                    Some(Line::Def { form: pass_form, .. }) if is_pass => *pass_form = Some(form),
                    _ => header.lines.push(Line::Def {
                        kind: LineKind::Filter,
                        id: id.clone(),
                        form: Some(form),
                    }),
                }
                header.filters.insert(id, FilterDef { description, attributes });
            }
            _ => {
                let is_info = key == "INFO";
                let declared = if is_info {
                    header.infos.contains_key(&id)
                } else {
                    header.formats.contains_key(&id)
                };
                if declared {
                    warn!(line, key, field = %id, "duplicate definition; keeping the first");
                    return Ok(());
                }
                let number = match take(&mut attributes, "Number") {
                    Some(value) => parse_number(&value)
                        .ok_or(VcfHeaderError::UnknownNumber { id: id.clone(), value })?,
                    None => {
                        warn!(line, key, field = %id, "no Number; assuming '.'");
                        Number::Unknown
                    }
                };
                let typ = match take(&mut attributes, "Type") {
                    Some(value) => parse_type(&value)
                        .ok_or(VcfHeaderError::UnknownType { id: id.clone(), value })?,
                    None => {
                        warn!(line, key, field = %id, "no Type; assuming String");
                        ValueType::String
                    }
                };
                // r[impl vcf_header.parse.defs]
                // Constraints the writer enforces, but which htslib only warns about.
                if typ == ValueType::Flag && number != Number::Count(0) {
                    warn!(line, key, field = %id, %number, "a Flag should have Number=0");
                }
                if typ == ValueType::Flag && !is_info {
                    warn!(line, field = %id, "a FORMAT field should not be a Flag");
                }
                let description = description(&mut attributes, line, &id);
                assign_field(&mut header.string_map, &id, idx)?;
                let kind = if is_info { LineKind::Info } else { LineKind::Format };
                header.lines.push(Line::Def { kind, id: id.clone(), form: Some(form) });
                if is_info {
                    header.infos.insert(id, InfoDef { number, typ, description, attributes });
                } else {
                    header.formats.insert(id, FormatDef { number, typ, description, attributes });
                }
            }
        }
        Ok(())
    }
}

fn description(attributes: &mut Vec<Attribute>, line: usize, id: &str) -> SmolStr {
    take(attributes, "Description").unwrap_or_else(|| {
        warn!(line, field = id, "no Description");
        SmolStr::default()
    })
}

/// Remove the first attribute named `key`, returning its value.
/// The first `key` attribute's value, every `key` attribute removed: a
/// repeat kept among the others would be written back before the one read,
/// and read in its place next time.
fn take(attributes: &mut Vec<Attribute>, key: &str) -> Option<SmolStr> {
    let at = attributes.iter().position(|attr| attr.key == key)?;
    let value = attributes.remove(at).value;
    attributes.retain(|attr| attr.key != key);
    Some(value)
}

// r[impl vcf_header.parse.dictionary]
fn assign_field(
    map: &mut StringMap,
    id: &SmolStr,
    idx: Option<u32>,
) -> Result<u32, VcfHeaderError> {
    assign(map, id, idx, MAX_DICT_ENTRIES, || VcfHeaderError::TooManyFields {
        max: MAX_DICT_ENTRIES,
    })
}

/// Give `id` its dictionary index as `bcf_hdr_register_hrec` does: an id seen
/// before (another category, or PASS) keeps its index; a new one takes its
/// `IDX`, or the next index past the highest.
fn assign(
    map: &mut StringMap,
    id: &SmolStr,
    idx: Option<u32>,
    max: usize,
    full: impl Fn() -> VcfHeaderError,
) -> Result<u32, VcfHeaderError> {
    match (map.get(id), idx) {
        (Some(first), Some(second)) if first != second as usize => {
            let first = u32::try_from(first).map_err(|_| full())?;
            Err(VcfHeaderError::InconsistentIdx { id: id.clone(), first, second })
        }
        (Some(existing), _) => u32::try_from(existing).map_err(|_| full()),
        (None, Some(idx)) => {
            if idx as usize >= max {
                return Err(full());
            }
            map.insert_at(id, idx).map_err(|first| VcfHeaderError::ConflictingIdx {
                idx,
                first,
                second: id.clone(),
            })?;
            Ok(idx)
        }
        (None, None) => map.push(id, max).ok_or_else(full),
    }
}

fn parse_number(value: &str) -> Option<Number> {
    Some(match value {
        "A" => Number::AlternateBases,
        "R" => Number::ReferenceAlternateBases,
        "G" => Number::Genotypes,
        "M" => Number::BaseModification,
        "P" => Number::Ploidy,
        "LA" => Number::LocalAlternateBases,
        "LR" => Number::LocalReferenceAlternateBases,
        "LG" => Number::LocalGenotypes,
        "." => Number::Unknown,
        n => Number::Count(parse_digits(n)?),
    })
}

fn parse_type(value: &str) -> Option<ValueType> {
    Some(match value {
        "Integer" => ValueType::Integer,
        "Float" => ValueType::Float,
        "Flag" => ValueType::Flag,
        "Character" => ValueType::Character,
        "String" => ValueType::String,
        _ => return None,
    })
}

/// The `#CHROM` line's sample names, checked.
fn parse_samples(line: &str) -> Result<Vec<SmolStr>, VcfHeaderError> {
    let mut columns = line.split('\t');
    for (column, expected) in FIXED_COLUMNS.iter().enumerate() {
        if columns.next() != Some(*expected) {
            return Err(VcfHeaderError::MalformedHeaderLine { column });
        }
    }
    match columns.next() {
        None => return Ok(Vec::new()),
        Some("FORMAT") => {}
        Some(_) => return Err(VcfHeaderError::MalformedHeaderLine { column: FIXED_COLUMNS.len() }),
    }
    // r[impl vcf_header.parse.limits]
    let n = columns.clone().count();
    if n > MAX_SAMPLES {
        return Err(VcfHeaderError::TooManySamples { max: MAX_SAMPLES });
    }
    let mut samples = Vec::with_capacity(n);
    let mut seen = HashSet::with_capacity(n);
    let first_sample_column = FIXED_COLUMNS.len().saturating_add(1);
    for (i, name) in columns.enumerate() {
        if name.is_empty() {
            return Err(VcfHeaderError::EmptySampleName {
                column: first_sample_column.saturating_add(i),
            });
        }
        if !seen.insert(name) {
            return Err(VcfHeaderError::DuplicateSample { name: SmolStr::from(name) });
        }
        samples.push(SmolStr::from(name));
    }
    Ok(samples)
}

/// The attributes of a structured line, `body` being everything after its
/// `<`. Grammar as htslib's `bcf_hdr_parse_line`: keys `[A-Za-z_][0-9A-Za-z_.]*`;
/// a value double-quoted (with `\"` and `\\` escapes), bracketed (`[a, b]`,
/// kept with its brackets) or bare up to the next `,` or `>` outside nested
/// `<…>`; spaces around keys and values ignored, and trailing ones inside
/// quotes. Each attribute comes with its value as written (escapes kept,
/// without the quotes).
pub(super) fn parse_structured(body: &str) -> Result<Vec<(Attribute, &str)>, MetaLineError> {
    let bytes = body.as_bytes();
    let at = |i: usize| bytes.get(i).copied();
    let skip_spaces = |mut i: usize| {
        while at(i) == Some(b' ') {
            i = i.saturating_add(1);
        }
        i
    };
    let slice = |from: usize, to: usize| body.get(from..to).unwrap_or_default();

    let mut attributes = Vec::new();
    let mut i = skip_spaces(0);
    loop {
        let key_start = i;
        if at(i).is_some_and(|b| b.is_ascii_alphabetic() || b == b'_') {
            i = i.saturating_add(1);
            while at(i).is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.') {
                i = i.saturating_add(1);
            }
        }
        if i == key_start {
            return Err(if at(i).is_none() {
                MetaLineError::Unterminated
            } else {
                MetaLineError::MissingKey
            });
        }
        let key = SmolStr::from(slice(key_start, i));
        i = skip_spaces(i);
        if at(i) != Some(b'=') {
            return Err(MetaLineError::MissingEquals);
        }
        i = skip_spaces(i.saturating_add(1));

        let (value, quoted, raw) = match at(i) {
            Some(b'"') => {
                let start = i.saturating_add(1);
                i = start;
                loop {
                    match at(i) {
                        None => return Err(MetaLineError::UnterminatedQuote),
                        // An escape: whatever follows cannot close the value.
                        Some(b'\\') => i = i.saturating_add(2),
                        Some(b'"') => break,
                        Some(_) => i = i.saturating_add(1),
                    }
                }
                // Trailing spaces dropped, as htslib drops them, unless that
                // would leave a backslash to escape the closing quote.
                let inner = slice(start, i);
                let trimmed = inner.trim_end_matches(' ');
                let inner = if trimmed.ends_with('\\') { inner } else { trimmed };
                i = i.saturating_add(1);
                (unescape(inner), true, inner)
            }
            Some(b'[') => {
                let start = i;
                let close = body
                    .get(start..)
                    .and_then(|rest| rest.find(']'))
                    .ok_or(MetaLineError::UnterminatedBracket)?;
                i = start.saturating_add(close).saturating_add(1);
                (SmolStr::from(slice(start, i)), false, slice(start, i))
            }
            _ => {
                let start = i;
                let mut depth = 0usize;
                loop {
                    match at(i) {
                        None => return Err(MetaLineError::Unterminated),
                        Some(b'<') => depth = depth.saturating_add(1),
                        Some(b'>' | b',') if depth == 0 => break,
                        Some(b'>') => depth = depth.saturating_sub(1),
                        Some(_) => {}
                    }
                    i = i.saturating_add(1);
                }
                let raw = slice(start, i).trim_end_matches(' ');
                (SmolStr::from(raw), false, raw)
            }
        };
        attributes.push((Attribute { key, value, quoted }, raw));

        i = skip_spaces(i);
        match at(i) {
            Some(b',') => i = skip_spaces(i.saturating_add(1)),
            Some(b'>') => {
                if !slice(i.saturating_add(1), bytes.len()).trim().is_empty() {
                    warn!("dropping trailing text after a structured header line's '>'");
                }
                return Ok(attributes);
            }
            None => return Err(MetaLineError::Unterminated),
            Some(_) => return Err(MetaLineError::JunkAfterValue),
        }
    }
}

// r[impl vcf_header.line_form]
/// The attributes `keep` keeps, as htslib's `bcf_hrec_format` writes them:
/// `key=value` with each value as written, joined by `,`.
pub(super) fn form(attributes: &[(Attribute, &str)], keep: impl Fn(&str) -> bool) -> String {
    let mut out = String::new();
    for (attr, raw) in attributes.iter().filter(|(attr, _)| keep(&attr.key)) {
        if !out.is_empty() {
            out.push(',');
        }
        out.push_str(&attr.key);
        out.push('=');
        if attr.quoted {
            out.push('"');
            out.push_str(raw);
            out.push('"');
        } else {
            out.push_str(raw);
        }
    }
    out
}

/// Undo `\"` and `\\`; any other backslash is kept as written.
fn unescape(raw: &str) -> SmolStr {
    if !raw.contains('\\') {
        return SmolStr::from(raw);
    }
    let mut out = String::with_capacity(raw.len());
    let mut chars = raw.chars();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.clone().next() {
                Some(next @ ('"' | '\\')) => {
                    out.push(next);
                    chars.next();
                }
                _ => out.push(c),
            }
        } else {
            out.push(c);
        }
    }
    SmolStr::from(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CHROM: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    fn parse(meta: &str) -> Result<VcfHeader, VcfHeaderError> {
        VcfHeader::parse(&format!("##fileformat=VCFv4.5\n{meta}{CHROM}\n"))
    }

    fn err(meta: &str) -> VcfHeaderError {
        parse(meta).expect_err("must not parse")
    }

    fn malformed(meta: &str) -> MetaLineError {
        match err(meta) {
            VcfHeaderError::MalformedLine { line: 2, reason } => reason,
            other => panic!("expected MalformedLine on line 2, got {other:?}"),
        }
    }

    // r[verify vcf_header.parse]
    #[test]
    fn quoted_values_hold_escapes_commas_and_equals() {
        let h = parse(
            "##INFO=<ID=X,Number=1,Type=String,Description=\"a \\\"b\\\", c=d \\\\ e\\x\",Source=\"s,1\">\n",
        )
        .unwrap();
        let def = &h.infos()["X"];
        assert_eq!(def.description, "a \"b\", c=d \\ e\\x");
        assert_eq!(
            def.attributes,
            [Attribute { key: "Source".into(), value: "s,1".into(), quoted: true }]
        );
    }

    // r[verify vcf_header.parse]
    // r[verify vcf_header.line_form]
    /// A structured other line is kept as htslib writes it: without the
    /// spaces around keys and values.
    #[test]
    fn spaces_brackets_and_nested_angles_are_read_as_htslib_reads_them() {
        let h = parse("##META=<ID=Assay, Type=String, Values=[WholeGenome, Exome]>\n##ALT=<ID=DEL,Description=\"<x>\">\n").unwrap();
        assert_eq!(
            h.other_lines(),
            [
                "META=<ID=Assay,Type=String,Values=[WholeGenome, Exome]>",
                "ALT=<ID=DEL,Description=\"<x>\">"
            ]
        );
        let attrs =
            parse_structured("ID=a , V = [x, y] ,N=<p,q>,D= \"a \\\"b\\x \",E=\"\\ \">").unwrap();
        assert_eq!(
            attrs.iter().map(|(attr, _)| attr.clone()).collect::<Vec<_>>(),
            [
                Attribute { key: "ID".into(), value: "a".into(), quoted: false },
                Attribute { key: "V".into(), value: "[x, y]".into(), quoted: false },
                Attribute { key: "N".into(), value: "<p,q>".into(), quoted: false },
                Attribute { key: "D".into(), value: "a \"b\\x".into(), quoted: true },
                Attribute { key: "E".into(), value: "\\ ".into(), quoted: true },
            ]
        );
        // r[verify vcf_header.line_form]
        // As htslib keeps them: escapes untouched, trailing spaces in quotes
        // dropped, but the one an escape would otherwise swallow the quote with.
        assert_eq!(
            form(&attrs, |_| true),
            "ID=a,V=[x, y],N=<p,q>,D=\"a \\\"b\\x\",E=\"\\ \"",
            "as written"
        );
    }

    // r[verify vcf_header.parse]
    #[test]
    fn each_malformed_structured_line_says_why() {
        assert_eq!(malformed("##INFO=<=1>\n"), MetaLineError::MissingKey);
        assert_eq!(malformed("##INFO=<ID>\n"), MetaLineError::MissingEquals);
        assert_eq!(
            malformed("##INFO=<ID=X,Description=\"open>\n"),
            MetaLineError::UnterminatedQuote
        );
        assert_eq!(malformed("##INFO=<ID=X,V=[a\n"), MetaLineError::UnterminatedBracket);
        assert_eq!(malformed("##INFO=<ID=X\n"), MetaLineError::Unterminated);
        assert_eq!(malformed("##INFO=<ID=\"X\"Y>\n"), MetaLineError::JunkAfterValue);
        // Not interpreted, so kept verbatim with a warning instead.
        assert_eq!(parse("##ALT=<ID=X\n").unwrap().other_lines(), ["ALT=<ID=X"]);
    }

    // r[verify vcf_header.parse.version]
    #[test]
    fn the_version_line_is_checked() {
        let bad = |text: &str| VcfHeader::parse(text).expect_err("must not parse");
        assert!(matches!(bad(""), VcfHeaderError::MissingFileFormat));
        assert!(matches!(bad(&format!("{CHROM}\n")), VcfHeaderError::MissingFileFormat));
        for value in ["VCF4.2", "VCFv4", "VCFv4.x", "VCFv+4.2", "VCFv4.256"] {
            let e = bad(&format!("##fileformat={value}\n{CHROM}\n"));
            assert!(
                matches!(e, VcfHeaderError::MalformedVersion { value: ref v } if v == value),
                "{value}: {e:?}"
            );
        }
        let e = bad(&format!("##fileformat=VCFv3.3\n{CHROM}\n"));
        assert!(
            matches!(e, VcfHeaderError::UnsupportedVersion { version } if version == VcfVersion::new(3, 3))
        );
        let newer = VcfHeader::parse(&format!("##fileformat=VCFv4.7\n{CHROM}\n")).unwrap();
        assert_eq!(newer.version(), VcfVersion::new(4, 7), "accepted, with a warning");
        assert!(newer.to_vcf_text().starts_with("##fileformat=VCFv4.7\n"));
    }

    // r[verify vcf_header.parse.defs]
    #[test]
    fn numbers_and_types_are_parsed_and_unknown_ones_named() {
        for (token, number) in [
            ("0", Number::Count(0)),
            ("12", Number::Count(12)),
            ("A", Number::AlternateBases),
            ("R", Number::ReferenceAlternateBases),
            ("G", Number::Genotypes),
            ("M", Number::BaseModification),
            ("P", Number::Ploidy),
            ("LA", Number::LocalAlternateBases),
            ("LR", Number::LocalReferenceAlternateBases),
            ("LG", Number::LocalGenotypes),
            (".", Number::Unknown),
        ] {
            let h =
                parse(&format!("##FORMAT=<ID=F,Number={token},Type=Integer,Description=\"\">\n"))
                    .unwrap();
            assert_eq!(h.formats()["F"].number, number, "{token}");
            assert_eq!(number.to_string(), token);
        }
        let e = err("##INFO=<ID=F,Number=-1,Type=Integer,Description=\"\">\n");
        assert!(
            matches!(e, VcfHeaderError::UnknownNumber { ref id, ref value } if id == "F" && value == "-1"),
            "{e:?}"
        );
        let e = err("##INFO=<ID=F,Number=1,Type=Int,Description=\"\">\n");
        assert!(
            matches!(e, VcfHeaderError::UnknownType { ref id, ref value } if id == "F" && value == "Int"),
            "{e:?}"
        );
    }

    // r[verify vcf_header.parse.defs]
    #[test]
    fn what_htslib_only_warns_about_is_accepted() {
        let h = parse(concat!(
            "##INFO=<ID=F,Number=1,Type=Flag,Description=\"flag with a count\">\n",
            "##FORMAT=<ID=G,Number=0,Type=Flag,Description=\"FORMAT flag\">\n",
            "##INFO=<ID=N>\n",
            "##INFO=<ID=F,Number=2,Type=Integer,Description=\"duplicate\">\n",
            "##FILTER=<ID=PASS,Description=\"Passed\">\n",
            "##FILTER=<ID=PASS,Description=\"again\">\n",
            "##contig=<ID=c,length=5>\n",
            "##contig=<ID=c,length=6>\n",
            "#comment\n",
            "##no equals sign\n",
            "##fileformat=VCFv4.1\n",
            "##fileformat=VCFv4.5\n",
        ))
        .unwrap();
        assert_eq!(
            h.infos()["F"],
            InfoDef::new(Number::Count(1), ValueType::Flag, "flag with a count")
        );
        assert_eq!(h.formats()["G"].typ, ValueType::Flag);
        assert_eq!(
            h.infos()["N"],
            InfoDef::new(Number::Unknown, ValueType::String, ""),
            "htslib's defaults"
        );
        assert_eq!(
            h.filters()["PASS"].description,
            "Passed",
            "the file's PASS replaces the default"
        );
        assert_eq!(h.contig_defs()["c"].length, Some(5), "the first definition wins");
        assert_eq!(
            h.other_lines(),
            ["no equals sign", "fileformat=VCFv4.1"],
            "a second ##fileformat is a line, as htslib keeps it; a repeat of the first is not"
        );
        assert_eq!(h.version(), VcfVersion::V4_5);
    }

    // r[verify vcf_header.line_form]
    /// A repeated other line is dropped as htslib's `bcf_hdr_add_hrec` drops
    /// it (`bcftools view -h` prints these lines once): the same `key=value`,
    /// or a structured line with the same key and `ID`. The builder drops
    /// them too.
    #[test]
    fn a_repeated_other_line_is_dropped_as_htslib_drops_it() {
        let h = parse(concat!(
            "##source=x\n##source=x\n##source=y\n",
            "##META=<ID=Assay, Type=String>\n##META=<ID=Assay,Type=Other>\n",
            "##ALT=<ID=Assay,Description=\"other key\">\n",
            "##META=<Type=no id>\n##META=<Type=no id>\n",
        ))
        .unwrap();
        let expected = [
            "source=x",
            "source=y",
            "META=<ID=Assay,Type=String>",
            "ALT=<ID=Assay,Description=\"other key\">",
            "META=<Type=no id>",
            "META=<Type=no id>",
        ];
        assert_eq!(h.other_lines(), expected);
        let mut builder = h.to_builder();
        builder.add_other_line("source=y");
        builder.add_other_line("META=<ID=Assay,Type=Third>");
        builder.add_other_line("source=z");
        let h = builder.build().unwrap();
        assert_eq!(h.other_lines().len(), expected.len() + 1);
        assert!(h.to_vcf_text().contains("##META=<Type=no id>\n##source=z\n#CHROM"));
    }

    // r[verify vcf_header.parse]
    #[test]
    fn a_definition_needs_an_id() {
        let e = err("##FILTER=<Description=\"x\">\n");
        assert!(
            matches!(e, VcfHeaderError::MissingId { line: 2, ref key } if key == "FILTER"),
            "{e:?}"
        );
    }

    #[test]
    fn contig_lengths_are_integers() {
        let h = parse("##contig=<ID=a>\n##contig=<ID=b,length=\"7\",md5=x>\n").unwrap();
        assert_eq!(h.contig_defs()["a"], ContigDef::new(None));
        assert_eq!(h.contig_defs()["b"].length, Some(7));
        let e = err("##contig=<ID=c,length=-1>\n");
        assert!(
            matches!(e, VcfHeaderError::InvalidContigLength { ref name, ref value } if name == "c" && value == "-1"),
            "{e:?}"
        );
    }

    // r[verify vcf_header.parse.dictionary]
    #[test]
    fn the_dictionary_follows_first_appearance_across_categories() {
        let h = parse(concat!(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n",
            "##FILTER=<ID=q10,Description=\"\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"\">\n",
            "##INFO=<ID=q10,Number=0,Type=Flag,Description=\"\">\n",
            "##contig=<ID=b>\n##contig=<ID=a>\n",
        ))
        .unwrap();
        let ids: Vec<_> = h.string_map().iter().collect();
        assert_eq!(ids, [(0, "PASS"), (1, "DP"), (2, "q10"), (3, "GT")]);
        assert_eq!(h.contig_map().iter().collect::<Vec<_>>(), [(0, "b"), (1, "a")]);
        let names: Vec<_> = h.contigs().map(|(c, _)| c.name().clone()).collect();
        assert_eq!(names, ["b", "a"], "dictionary order");
        // r[verify vcf_header.line_order]
        // Written back in the order read, the lines give the same dictionary
        // without IDX.
        let text = h.to_bcf_text();
        assert!(!text.contains("IDX="), "{text}");
        assert_eq!(text, h.to_vcf_text());
        let again = VcfHeader::parse(&text).unwrap();
        assert_eq!(again.string_map(), h.string_map());
        assert_eq!(again.contig_map(), h.contig_map());
    }

    // r[verify vcf_header.parse.preserve+2]
    /// An id read quoted that a bare value cannot hold is written quoted
    /// (fuzzer-found).
    #[test]
    fn an_id_a_bare_value_cannot_hold_is_written_quoted() {
        let h = parse("##FILTER=<ID=\"a, b>\",Description=\"\">\n##contig=<ID=\"[c\">\n").unwrap();
        let again = VcfHeader::parse(&h.to_vcf_text()).unwrap();
        assert_eq!(
            again.filters().keys().collect::<Vec<_>>(),
            h.filters().keys().collect::<Vec<_>>()
        );
        assert_eq!(again.contig_map(), h.contig_map());
    }

    // r[verify vcf_header.parse.dictionary]
    // r[verify vcf_header.parse.preserve+2]
    // r[verify vcf_header.line_form]
    /// A repeated `IDX` (or `ID`) is read once, the first, and not written
    /// back as an attribute that would win the next read (fuzzer-found).
    #[test]
    fn a_repeated_reserved_key_reads_and_writes_as_the_first() {
        let h = parse("##FILTER=<ID=f,Description=\"\",IDX=7,x=y,IDX=2,ID=g>\n").unwrap();
        assert_eq!(h.string_map().get("f"), Some(7));
        // As htslib writes it: every key as read, but IDX.
        let text = h.to_vcf_text();
        assert!(text.contains("##FILTER=<ID=f,Description=\"\",x=y,ID=g>\n"), "{text}");
        assert_eq!(VcfHeader::parse(&text).unwrap().filters().keys().nth(1).unwrap(), "f");
        let text = h.to_bcf_text();
        assert!(text.contains("##FILTER=<ID=f,Description=\"\",x=y,ID=g,IDX=7>"), "{text}");
        let again = VcfHeader::parse(&text).unwrap();
        assert_eq!(again.string_map(), h.string_map());
    }

    // r[verify vcf_header.parse.dictionary]
    #[test]
    fn idx_is_honoured_with_gaps_and_conflicts_are_errors() {
        let h = parse(concat!(
            "##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=5>\n",
            "##INFO=<ID=B,Number=1,Type=Integer,Description=\"\",IDX=2>\n",
            "##INFO=<ID=C,Number=1,Type=Integer,Description=\"\">\n",
            "##FORMAT=<ID=A,Number=1,Type=Integer,Description=\"\">\n",
            "##contig=<ID=c,IDX=3>\n",
        ))
        .unwrap();
        assert_eq!(
            h.string_map().iter().collect::<Vec<_>>(),
            [(0, "PASS"), (2, "B"), (5, "A"), (6, "C")]
        );
        assert_eq!(h.string_map().len(), 7);
        assert_eq!(h.string_map().name(1), None, "a hole");
        assert_eq!(h.contig_map().get("c"), Some(3));
        // r[verify vcf_header.contigs]
        let contigs: Vec<_> = h.contigs().map(|(c, _)| (c.index(), c.name().clone())).collect();
        assert_eq!(contigs, [(3, "c".into())], "by IDX, the holes skipped");
        // r[verify vcf_header.serialization+2]
        assert!(!h.to_vcf_text().contains("IDX="), "VCF text never carries IDX");
        let text = h.to_bcf_text();
        assert!(
            text.contains("##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>"),
            "{text}"
        );
        assert!(text.contains("##contig=<ID=c,IDX=3>"), "{text}");
        let again = VcfHeader::parse(&text).unwrap();
        assert_eq!((again.string_map(), again.contig_map()), (h.string_map(), h.contig_map()));

        let e = err(concat!(
            "##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=1>\n",
            "##INFO=<ID=B,Number=1,Type=Integer,Description=\"\",IDX=1>\n",
        ));
        assert!(
            matches!(e, VcfHeaderError::ConflictingIdx { idx: 1, ref first, ref second } if first == "A" && second == "B"),
            "{e:?}"
        );
        let e = err("##FILTER=<ID=PASS,Description=\"\",IDX=3>\n");
        assert!(
            matches!(e, VcfHeaderError::InconsistentIdx { ref id, first: 0, second: 3 } if id == "PASS"),
            "{e:?}"
        );
        let e = err(concat!(
            "##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=1>\n",
            "##FORMAT=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=2>\n",
        ));
        assert!(
            matches!(e, VcfHeaderError::InconsistentIdx { ref id, first: 1, second: 2 } if id == "A"),
            "{e:?}"
        );
        let e = err("##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=x>\n");
        assert!(
            matches!(e, VcfHeaderError::InvalidIdx { ref id, ref value } if id == "A" && value == "x"),
            "{e:?}"
        );
    }

    // r[verify vcf_header.parse.limits]
    #[test]
    fn limits_are_enforced() {
        let e = err("##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=65536>\n");
        assert!(matches!(e, VcfHeaderError::TooManyFields { max: MAX_DICT_ENTRIES }), "{e:?}");
        parse("##INFO=<ID=A,Number=1,Type=Integer,Description=\"\",IDX=65535>\n").unwrap();
        let e = err("##contig=<ID=c,IDX=1048576>\n");
        assert!(matches!(e, VcfHeaderError::TooManyContigs { max: MAX_CONTIGS }), "{e:?}");

        let mut text = String::from("##fileformat=VCFv4.5\n");
        text.extend(std::iter::repeat_n('#', MAX_HEADER_TEXT));
        let e = VcfHeader::parse(&text).unwrap_err();
        assert!(matches!(e, VcfHeaderError::TextTooLarge { max: MAX_HEADER_TEXT, .. }), "{e:?}");

        let mut line = format!("{CHROM}\tFORMAT");
        line.extend(std::iter::repeat_n("\ts", MAX_SAMPLES.saturating_add(1)));
        assert!(matches!(
            parse_samples(&line),
            Err(VcfHeaderError::TooManySamples { max: MAX_SAMPLES })
        ));
    }

    // r[verify vcf_header.parse]
    // r[verify vcf_header.sample_names]
    #[test]
    fn the_chrom_line_is_checked() {
        let header = |chrom: &str| VcfHeader::parse(&format!("##fileformat=VCFv4.5\n{chrom}\n"));
        assert!(header(CHROM).unwrap().samples().is_empty());
        assert!(header(&format!("{CHROM}\tFORMAT")).unwrap().samples().is_empty());
        assert_eq!(header(&format!("{CHROM}\tFORMAT\tA\tB")).unwrap().samples(), ["A", "B"]);
        assert_eq!(
            VcfHeader::parse(&format!("##fileformat=VCFv4.5\r\n{CHROM}\tFORMAT\tA\r\n\0\0"))
                .unwrap()
                .samples(),
            ["A"]
        );

        let bad = |chrom: &str| header(chrom).expect_err("must not parse");
        assert!(matches!(
            bad("#CHROM POS ID REF ALT QUAL FILTER INFO"),
            VcfHeaderError::MalformedHeaderLine { column: 0 }
        ));
        assert!(matches!(
            bad("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"),
            VcfHeaderError::MalformedHeaderLine { column: 7 }
        ));
        assert!(matches!(
            bad(&format!("{CHROM}\tA")),
            VcfHeaderError::MalformedHeaderLine { column: 8 }
        ));
        assert!(matches!(
            bad(&format!("{CHROM}\tFORMAT\tA\t")),
            VcfHeaderError::EmptySampleName { column: 10 }
        ));
        assert!(
            matches!(bad(&format!("{CHROM}\tFORMAT\tA\tA")), VcfHeaderError::DuplicateSample { ref name } if name == "A")
        );
        assert!(matches!(
            bad(&format!("{CHROM}\nchr1\t1")),
            VcfHeaderError::TrailingData { line: 3 }
        ));
        assert!(matches!(
            VcfHeader::parse("##fileformat=VCFv4.5\n").unwrap_err(),
            VcfHeaderError::MissingHeaderLine
        ));
    }
}
