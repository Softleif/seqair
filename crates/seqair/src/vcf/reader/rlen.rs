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
    let reserved = i32::MIN.saturating_add(8);
    parse_i32(token).ok().flatten().filter(|&v| v >= reserved).map_or(0, i64::from)
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
/// The reference length htslib gives a VCF text record at 0-based `pos`.
pub(crate) fn text_rlen(cols: &RlenColumns<'_>, file: &FileShared, pos: i64) -> i64 {
    let keys = file.rlen_keys;
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
    let (mut end, mut svlen) = (0i64, None::<i64>);
    if (keys.end || use_svlen) && cols.info != b"." {
        let (mut seen_end, mut seen_svlen) = (false, false);
        let _ = for_each_info_entry(cols.info, |e| -> ControlFlow<()> {
            let key = cols.info.get(e.key.clone()).unwrap_or_default();
            let value = e.value.and_then(|v| cols.info.get(v));
            if keys.end && !seen_end && key == b"END" {
                seen_end = true;
                // A vector, or no value, holds no single END.
                end = value.filter(|v| !v.contains(&b',')).map_or(0, int_or_zero);
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
    hpos.saturating_sub(pos).max(len_ref)
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
