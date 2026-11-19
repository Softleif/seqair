//! A VCF text record's reference length, as htslib's `get_rlen` computes it
//! from the record it encodes (`r[vcf_values.rlen]`).
//!
//! htslib reads INFO END, the SVLEN of `<DEL>`/`<DUP>`/`<INV>`/`<CNV>`
//! alleles (and their `:` subtypes) and FORMAT LEN of a gVCF `<*>` or
//! `<NON_REF>` allele whatever the file's version ("as version is not kept
//! properly updated in practice, it is ignored"), and so does this. Indexes
//! were built from that value (tabix's `tbx_parse1` is kept in sync with
//! it), so a reader that computes another end misses records.

use super::file::FileShared;
use super::keys::FieldTable;
use super::text_value::parse_i32;
use crate::io::text_scan::{for_each_byte, for_each_info_entry};
use crate::vcf::error::KeyKind;
use crate::vcf::record_encoder::ContigId;
use seqair_types::SmallVec;
use std::ops::ControlFlow;

/// Which of END, SVLEN and LEN a header declares so that htslib encodes
/// them as integers: only those take part in the length.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct RlenKeys {
    pub(crate) end: bool,
    pub(crate) svlen: bool,
    pub(crate) len: bool,
}

impl RlenKeys {
    // r[impl vcf_values.rlen+2]
    /// Those the field table gives `Type=Integer`: declared so or, in VCF
    /// text, not declared at all (VCF reserves them as integers).
    pub(crate) fn new(fields: &FieldTable) -> Self {
        use crate::vcf::ValueType::Integer;
        let info = |name| fields.declared_name(KeyKind::Info, name) == Some(Integer);
        let len = fields.declared_name(KeyKind::Format, "LEN") == Some(Integer);
        Self { end: info("END"), svlen: info("SVLEN"), len }
    }

    /// Those `header` declares `Type=Integer`: what a writer can write,
    /// since it writes declared fields only. For a declared field this is
    /// what [`Self::new`] gives the file it writes.
    pub(crate) fn from_header(header: &crate::vcf::VcfHeader) -> Self {
        use crate::vcf::ValueType::Integer;
        let info = |name| header.infos().get(name).is_some_and(|d| d.typ == Integer);
        let len = header.formats().get("LEN").is_some_and(|d| d.typ == Integer);
        Self { end: info("END"), svlen: info("SVLEN"), len }
    }
}

/// An ALT whose SVLEN is its reference length: `<CNV>`, `<DEL>`, `<DUP>`,
/// `<INV>` and their `:` subtypes (htslib's `svlen_on_ref_for_vcf_alt`).
fn svlen_on_ref(alt: &[u8]) -> bool {
    let kind = alt.get(..4);
    alt.len() >= 5
        && matches!(kind, Some(b"<CNV" | b"<DEL" | b"<DUP" | b"<INV"))
        && matches!(alt.get(4), Some(b'>' | b':'))
        && alt.last() == Some(&b'>')
}

/// An integer as htslib reads one for the length: `.`, anything that does
/// not parse, and anything outside `i32` or in BCF's reserved band below
/// `i32::MIN + 8` (htslib stores those as missing) count as 0.
fn int_or_zero(token: &[u8]) -> i64 {
    int(token).unwrap_or(0)
}

/// [`int_or_zero`], with `None` for what counts as 0 there.
fn int(token: &[u8]) -> Option<i64> {
    let reserved = i32::MIN.saturating_add(8);
    parse_i32(token).ok().flatten().filter(|&v| v >= reserved).map(i64::from)
}

/// The columns [`text_rlen`] reads.
pub(crate) struct RlenColumns<'a> {
    pub(crate) reference: &'a [u8],
    pub(crate) alt: &'a [u8],
    pub(crate) info: &'a [u8],
    /// The FORMAT column and the samples after it, TAB-separated.
    pub(crate) format: &'a [u8],
}

// r[impl vcf_values.rlen+2]
// r[impl vcf_values.rlen.end_before_pos]
/// The reference length htslib gives a VCF text record on `contig` at
/// 0-based `pos`; an END before POS is warned of, once per file.
pub(crate) fn text_rlen(
    cols: &RlenColumns<'_>,
    file: &FileShared,
    contig: &ContigId,
    pos: i64,
) -> i64 {
    let (rlen, end) = rlen_and_end(cols, file.rlen_keys, pos);
    if let Some(end) = end
        && end <= pos
    {
        file.warn_end_before_pos(contig, pos, end);
    }
    rlen
}

// r[impl vcf_values.rlen+2]
// r[impl record_encoder.rlen]
/// [`text_rlen`] for a header's `keys`: the writer's entry point, so that what
/// it indexes a record under is what a reader computes for it.
pub(crate) fn rlen(cols: &RlenColumns<'_>, keys: RlenKeys, pos: i64) -> i64 {
    rlen_and_end(cols, keys, pos).0
}

/// [`rlen`], and the first END as read (`None` if absent, missing or not
/// one integer).
fn rlen_and_end(cols: &RlenColumns<'_>, keys: RlenKeys, pos: i64) -> (i64, Option<i64>) {
    let len_ref = i64::try_from(cols.reference.len()).unwrap_or(i64::MAX);
    // Bit `a` for ALT `a` (REF is 0), as htslib's `svlenals`.
    let mut sv_alts: SmallVec<u64, 2> = SmallVec::new();
    let mut gvcf = false;
    let mut n_allele = 1usize;
    if cols.alt != b"." {
        for (i, alt) in cols.alt.split(|&b| b == b',').enumerate() {
            let allele = i.saturating_add(1);
            n_allele = allele.saturating_add(1);
            if alt.first() != Some(&b'<') {
                continue;
            }
            if svlen_on_ref(alt) {
                let word = allele / 64;
                if sv_alts.len() <= word {
                    sv_alts.resize(word.saturating_add(1), 0);
                }
                if let Some(w) = sv_alts.get_mut(word) {
                    *w |= 1 << (allele % 64);
                }
            } else if alt == b"<*>" || alt == b"<NON_REF>" {
                gvcf = true;
            }
        }
    }
    // htslib tests SVLEN element `i` against bit `(i + 1) & 7` of byte
    // `i >> 3`: ALT `i + 1`, except that every eighth element reads the bit
    // of the ALT seven before (ALT 0, REF, for the eighth). Matched exactly.
    let flagged = |i: usize| {
        let bit = (i & !7) | (i.saturating_add(1) & 7);
        sv_alts.get(bit / 64).is_some_and(|w| w & (1 << (bit % 64)) != 0)
    };
    let use_svlen = keys.svlen && sv_alts.iter().any(|&w| w != 0);
    let (mut end_read, mut svlen) = (None::<i64>, None::<i64>);
    if (keys.end || use_svlen) && cols.info != b"." {
        let (mut seen_end, mut seen_svlen) = (false, false);
        let _ = for_each_info_entry(cols.info, |e| -> ControlFlow<()> {
            let key = cols.info.get(e.key.clone()).unwrap_or_default();
            let value = e.value.and_then(|v| cols.info.get(v));
            if keys.end && !seen_end && key == b"END" {
                seen_end = true;
                // A vector, or no value, holds no single END.
                end_read = value.filter(|v| !v.contains(&b',')).and_then(int);
            } else if use_svlen && !seen_svlen && key == b"SVLEN" {
                seen_svlen = true;
                let mut longest = 0i64;
                let values = value.map(|v| v.split(|&b| b == b','));
                for (i, v) in values.into_iter().flatten().enumerate() {
                    if i.saturating_add(1) >= n_allele {
                        break;
                    }
                    if flagged(i) {
                        longest = longest.max(int_or_zero(v).saturating_abs());
                    }
                }
                svlen = Some(longest);
            }
            if (seen_end || !keys.end) && (seen_svlen || !use_svlen) {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        });
    }
    let end = end_read.unwrap_or(0);
    let after = |e: i64, d: i64| if e > pos { e.saturating_sub(pos).saturating_sub(d) } else { 0 };
    let len = match svlen {
        Some(l) if l != 0 => l,
        _ if end != 0 => after(end, 1),
        _ => 0,
    };
    let end_svlen = pos.saturating_add(len).saturating_add(1);
    let fmt_len = if gvcf && keys.len { format_len(cols.format) } else { None };
    let len = match fmt_len {
        Some(l) if l != 0 => l,
        _ if end != 0 => after(end, 0),
        _ => 0,
    };
    let end_fmtlen = pos.saturating_add(len);
    let hpos = end.max(end_svlen).max(end_fmtlen);
    (hpos.saturating_sub(pos).max(len_ref), end_read)
}

/// The largest FORMAT LEN value of any sample (0 if none is positive);
/// `None` if FORMAT has no LEN.
fn format_len(format: &[u8]) -> Option<i64> {
    let (keys, samples) = match format.iter().position(|&b| b == b'\t') {
        Some(tab) => (format.get(..tab)?, format.get(tab.saturating_add(1)..)?),
        None => (format, &[][..]),
    };
    let k = keys.split(|&b| b == b':').position(|key| key == b"LEN")?;
    let mut longest = 0i64;
    let mut start = 0usize;
    let mut column = |col: &[u8]| {
        if let Some(value) = col.split(|&b| b == b':').nth(k) {
            for v in value.split(|&b| b == b',') {
                longest = longest.max(int_or_zero(v));
            }
        }
    };
    let _ = for_each_byte(samples, b'\t', |tab| -> ControlFlow<()> {
        column(samples.get(start..tab).unwrap_or_default());
        start = tab.saturating_add(1);
        ControlFlow::Continue(())
    });
    column(samples.get(start..).unwrap_or_default());
    Some(longest)
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_values.rlen+2]
    /// The writer's keys (declared fields only) are the reader's for every
    /// header that declares the field; only a VCF text file whose header
    /// leaves END/SVLEN/LEN out reads them where a writer cannot write them.
    #[test]
    fn the_writers_keys_are_the_readers_for_declared_fields() {
        let keys = |k: RlenKeys| (k.end, k.svlen, k.len);
        let header = |lines: &str| {
            let text = format!(
                "##fileformat=VCFv4.3\n{lines}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            );
            crate::vcf::VcfHeader::parse(&text).unwrap()
        };
        let declared = header(concat!(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n",
            "##INFO=<ID=SVLEN,Number=A,Type=String,Description=\"s\">\n",
            "##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"l\">\n",
        ));
        for text in [true, false] {
            let read = RlenKeys::new(&FieldTable::new(&declared, text));
            assert_eq!(keys(read), keys(RlenKeys::from_header(&declared)));
            assert_eq!(keys(read), (true, false, true));
        }
        let silent = header("");
        assert_eq!(keys(RlenKeys::from_header(&silent)), (false, false, false));
        assert_eq!(keys(RlenKeys::new(&FieldTable::new(&silent, false))), (false, false, false));
        assert_eq!(keys(RlenKeys::new(&FieldTable::new(&silent, true))), (true, true, true));
    }

    /// Each `warn!`'s fields, `name=value` joined by spaces.
    #[derive(Default)]
    struct Warnings(std::sync::Mutex<Vec<String>>);

    impl tracing::Subscriber for Warnings {
        fn enabled(&self, _: &tracing::Metadata<'_>) -> bool {
            true
        }
        fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
            tracing::span::Id::from_u64(1)
        }
        fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}
        fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}
        fn event(&self, event: &tracing::Event<'_>) {
            struct Fields(String);
            impl tracing::field::Visit for Fields {
                fn record_debug(&mut self, f: &tracing::field::Field, v: &dyn std::fmt::Debug) {
                    use std::fmt::Write;
                    write!(self.0, "{}={v:?} ", f.name()).unwrap();
                }
            }
            if *event.metadata().level() == tracing::Level::WARN {
                let mut fields = Fields(String::new());
                event.record(&mut fields);
                self.0.lock().unwrap().push(fields.0);
            }
        }
        fn enter(&self, _: &tracing::span::Id) {}
        fn exit(&self, _: &tracing::span::Id) {}
    }

    // r[verify vcf_values.rlen.end_before_pos]
    /// An END before POS leaves the span at REF's length and is warned of
    /// once per file, with the path; an END at POS is a one-base span, and
    /// of two ENDs the first counts.
    #[test]
    fn an_end_before_pos_is_ignored_and_warned_of_once() {
        use crate::vcf::{Locus, VcfFile, VcfPos};
        use std::io::Write;
        let header = "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
            ##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let spans = |records: &[&str]| {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            write!(tmp, "{header}{}", records.join("")).unwrap();
            let warnings = std::sync::Arc::new(Warnings::default());
            let file = VcfFile::open(tmp.path()).unwrap();
            let spans: Vec<(VcfPos, VcfPos)> =
                tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
                    let reader = file.reader::<Locus>().open().unwrap();
                    reader.into_iter().map(|l| l.map(|l| (l.start(), l.end())).unwrap()).collect()
                });
            let at = |p: VcfPos| match p {
                VcfPos::At(p) => p.as_u32(),
                VcfPos::Telomere => 0,
            };
            let spans: Vec<(u32, u32)> = spans.into_iter().map(|(s, e)| (at(s), at(e))).collect();
            let warnings = warnings.0.lock().unwrap().clone();
            (spans, warnings, tmp)
        };
        let (got, warnings, tmp) = spans(&[
            "chr1\t20\t.\tA\t<DEL>\t.\t.\tEND=5\n",
            "chr1\t30\t.\tAC\t<DEL>\t.\t.\tEND=29\n",
            "chr1\t40\t.\tA\t<DEL>\t.\t.\tEND=40\n",
        ]);
        assert_eq!(got, [(20, 20), (30, 31), (40, 40)]);
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        let w = &warnings[0];
        assert!(w.contains("chrom=chr1 pos=20 end=5 "), "{w}");
        assert!(w.contains(&format!("{:?}", Some(tmp.path()))), "{w}");

        let (got, warnings, _) = spans(&[
            "chr1\t20\t.\tA\t<DEL>\t.\t.\tEND=50;END=5\n",
            "chr1\t30\t.\tA\t<DEL>\t.\t.\tEND=5;END=50\n",
        ]);
        assert_eq!(got, [(20, 50), (30, 30)]);
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert!(warnings[0].contains("pos=30 end=5 "), "{warnings:?}");
    }

    #[test]
    fn symbolic_alleles_whose_svlen_spans_the_reference() {
        for yes in [&b"<DEL>"[..], b"<DUP>", b"<INV>", b"<CNV>", b"<DEL:ME>", b"<CNV:TR>"] {
            assert!(svlen_on_ref(yes), "{yes:?}");
        }
        for no in [&b"<INS>"[..], b"<*>", b"<DELETION>", b"<DEL", b"<DE>", b"A", b"<NON_REF>"] {
            assert!(!svlen_on_ref(no), "{no:?}");
        }
    }
}
