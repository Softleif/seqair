//! Parse BAM headers and look up contigs. [`BamHeader`] maps between tid integers,
//! sequence names, and lengths; shared across forked readers via `Arc`.

use super::bgzf::{BgzfError, BgzfReader};
use seqair_types::SmolStr;
use std::{
    collections::HashMap,
    io::{Read, Seek},
    path::Path,
};
use tracing::instrument;

#[derive(Debug)]
pub struct BamHeader {
    header_text: String,
    targets: Vec<TargetInfo>,
    name_to_tid: HashMap<SmolStr, u32>,
}

#[derive(Debug, Clone)]
struct TargetInfo {
    name: SmolStr,
    length: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContigInfo {
    pub tid: u32,
    pub len: u64,
}

// r[impl bam_writer.header_add_pg]
/// Fields for a `@PG` header record.
#[derive(Debug, Clone)]
pub struct PgRecord {
    pub id: SmolStr,
    pub pn: Option<SmolStr>,
    pub vn: Option<SmolStr>,
    pub cl: Option<SmolStr>,
    pub ds: Option<SmolStr>,
    pub pp: Option<SmolStr>,
}

/// Access target name and length without exposing `TargetInfo` internals.
pub trait TargetInfoAccess {
    fn target_name(&self) -> &str;
    fn target_length(&self) -> u64;
}

impl TargetInfoAccess for TargetInfo {
    fn target_name(&self) -> &str {
        self.name.as_str()
    }
    fn target_length(&self) -> u64 {
        self.length
    }
}

impl AsRef<str> for TargetInfo {
    fn as_ref(&self) -> &str {
        self.name.as_str()
    }
}

// r[impl bam.header.errors]
#[derive(Debug, thiserror::Error)]
pub enum BamHeaderError {
    #[error("invalid BAM magic bytes (expected BAM\\x01)")]
    InvalidMagic,

    #[error("@SQ header line missing SN (sequence name) tag: {line}")]
    MissingSn { line: String },

    #[error("@SQ header line missing LN (sequence length) tag: {line}")]
    MissingLn { line: String },

    #[error("header text is not valid UTF-8 (BAM spec requires UTF-8)")]
    NonUtf8HeaderText { source: std::str::Utf8Error },

    #[error("reference name is not valid UTF-8 (BAM spec requires printable ASCII)")]
    NonUtf8ReferenceName { source: std::str::Utf8Error },

    #[error("no @SQ lines found in SAM header text")]
    NoSequences,

    #[error(
        "indexed region queries require coordinate-sorted input, but the header specifies \
         SO:{sort_order}. Re-sort the file with `samtools sort -o sorted.bam input.bam` \
         and re-index with `samtools index sorted.bam`."
    )]
    UnsupportedSortOrder { sort_order: SmolStr },

    #[error("contig `{name}` not found in BAM header")]
    ContigNotFound { name: SmolStr },

    #[error("contig `{name}` (tid={tid}) has no length in BAM header")]
    ContigNoLength { name: SmolStr, tid: u32 },

    #[error("BAM header field `{field}` has negative value {value}")]
    NegativeLength { field: &'static str, value: i32 },

    #[error("BAM header field `{field}` value {value} exceeds limit {limit} — data may be corrupt")]
    FieldTooLarge { field: &'static str, value: usize, limit: usize },

    #[error("BGZF error reading BAM header")]
    Bgzf {
        #[from]
        source: BgzfError,
    },
}

impl BamHeader {
    // r[impl bam.header.magic]
    // r[impl bam.header.text]
    // r[impl bam.header.references]
    /// Parse the BAM header from a BGZF reader positioned at the start of the file.
    #[instrument(level = "debug", skip(bgzf), err)]
    /// Parse the BAM binary header from a BGZF reader.
    #[cfg(feature = "fuzz")]
    #[instrument(level = "debug", skip(bgzf), err)]
    pub fn parse<R: Read + Seek>(bgzf: &mut BgzfReader<R>) -> Result<Self, BamHeaderError> {
        Self::parse_inner(bgzf)
    }

    #[cfg(not(feature = "fuzz"))]
    #[instrument(level = "debug", skip(bgzf), err)]
    pub(crate) fn parse<R: Read + Seek>(bgzf: &mut BgzfReader<R>) -> Result<Self, BamHeaderError> {
        Self::parse_inner(bgzf)
    }

    fn parse_inner<R: Read + Seek>(bgzf: &mut BgzfReader<R>) -> Result<Self, BamHeaderError> {
        // Magic: BAM\1
        let mut magic = [0u8; 4];
        bgzf.read_exact_into(&mut magic)?;
        if magic != *b"BAM\x01" {
            return Err(BamHeaderError::InvalidMagic);
        }

        // r[impl bam.header.non_negative_lengths]
        // Header text
        let l_text_raw = bgzf.read_i32()?;
        if l_text_raw < 0 {
            return Err(BamHeaderError::NegativeLength { field: "l_text", value: l_text_raw });
        }
        let l_text = l_text_raw as usize;
        if l_text > 256 * 1024 * 1024 {
            return Err(BamHeaderError::FieldTooLarge {
                field: "l_text",
                value: l_text,
                limit: 256 * 1024 * 1024,
            });
        }
        let mut text_buf = vec![0u8; l_text];
        bgzf.read_exact_into(&mut text_buf)?;
        let header_text = std::str::from_utf8(&text_buf)
            .map_err(|source| BamHeaderError::NonUtf8HeaderText { source })?
            .to_owned();

        // Reference sequences
        let n_ref_raw = bgzf.read_i32()?;
        if n_ref_raw < 0 {
            return Err(BamHeaderError::NegativeLength { field: "n_ref", value: n_ref_raw });
        }
        let n_ref = n_ref_raw as usize;
        if n_ref > 1_000_000 {
            return Err(BamHeaderError::FieldTooLarge {
                field: "n_ref",
                value: n_ref,
                limit: 1_000_000,
            });
        }
        let mut targets = Vec::with_capacity(n_ref);
        let mut name_to_tid = HashMap::with_capacity(n_ref);

        for tid in 0..n_ref {
            let l_name_raw = bgzf.read_i32()?;
            if l_name_raw < 0 {
                return Err(BamHeaderError::NegativeLength { field: "l_name", value: l_name_raw });
            }
            let l_name = l_name_raw as usize;
            if l_name > 256 * 1024 {
                return Err(BamHeaderError::FieldTooLarge {
                    field: "l_name",
                    value: l_name,
                    limit: 256 * 1024,
                });
            }
            let mut name_buf = vec![0u8; l_name];
            bgzf.read_exact_into(&mut name_buf)?;
            // Remove trailing NUL
            if name_buf.last() == Some(&0) {
                name_buf.pop();
            }
            let name = SmolStr::new(
                std::str::from_utf8(&name_buf)
                    .map_err(|source| BamHeaderError::NonUtf8ReferenceName { source })?,
            );

            let l_ref = bgzf.read_i32()? as u64;

            name_to_tid.insert(name.clone(), tid as u32);
            targets.push(TargetInfo { name, length: l_ref });
        }

        Ok(BamHeader { header_text, targets, name_to_tid })
    }

    /// Open a BAM file and parse just the header.
    #[instrument(level = "debug", fields(path = %path.display()), err)]
    pub fn from_bam_path(path: &Path) -> Result<Self, BamHeaderError> {
        let mut bgzf = BgzfReader::open(path)?;
        Self::parse(&mut bgzf)
    }

    // r[impl unified.header_from_sam_text]
    /// Parse a `BamHeader` from SAM header text (the `@HD`, `@SQ`, `@RG`, etc. lines).
    ///
    /// Extracts target names and lengths from `@SQ` lines (`SN` and `LN` tags).
    /// Returns an error if no `@SQ` lines are found.
    #[instrument(level = "debug", skip(text), fields(text_len = text.len()), err)]
    pub fn from_sam_text(text: &str) -> Result<Self, BamHeaderError> {
        let mut targets = Vec::new();
        let mut name_to_tid = HashMap::new();

        for line in text.lines() {
            if !line.starts_with("@SQ") {
                continue;
            }
            let mut name: Option<&str> = None;
            let mut length: Option<u64> = None;

            for field in line.split('\t').skip(1) {
                if let Some(val) = field.strip_prefix("SN:") {
                    name = Some(val);
                } else if let Some(val) = field.strip_prefix("LN:") {
                    length = val.parse().ok();
                }
            }

            match (name, length) {
                (Some(n), Some(l)) => {
                    let tid = targets.len() as u32;
                    name_to_tid.insert(SmolStr::new(n), tid);
                    targets.push(TargetInfo { name: SmolStr::new(n), length: l });
                }
                (None, _) => return Err(BamHeaderError::MissingSn { line: line.to_string() }),
                (Some(_), None) => {
                    return Err(BamHeaderError::MissingLn { line: line.to_string() });
                }
            }
        }

        if targets.is_empty() {
            return Err(BamHeaderError::NoSequences);
        }

        Ok(BamHeader { header_text: text.to_string(), targets, name_to_tid })
    }

    pub fn resolve_contig(&self, name: &str) -> Result<ContigInfo, BamHeaderError> {
        let tid = self
            .tid(name)
            .ok_or_else(|| BamHeaderError::ContigNotFound { name: SmolStr::new(name) })?;
        let len = self
            .target_len(tid)
            .ok_or_else(|| BamHeaderError::ContigNoLength { name: SmolStr::new(name), tid })?;
        Ok(ContigInfo { tid, len })
    }

    pub fn target_count(&self) -> usize {
        self.targets.len()
    }

    // r[impl bam.header.tid_lookup]
    pub fn tid(&self, name: &str) -> Option<u32> {
        self.name_to_tid.get(name).copied()
    }

    pub fn target_name(&self, tid: u32) -> Option<&str> {
        self.targets.get(tid as usize).map(|t| t.name.as_str())
    }

    // r[impl bam.header.target_len]
    pub fn target_len(&self, tid: u32) -> Option<u64> {
        self.targets.get(tid as usize).map(|t| t.length)
    }

    // r[impl bam.header.target_names]
    pub fn target_names(&self) -> impl Iterator<Item = &str> {
        self.targets.iter().map(|t| t.name.as_str())
    }

    // r[impl bam_writer.header_text_generation]
    pub fn header_text(&self) -> &str {
        &self.header_text
    }

    /// Access the full target list (name + length) for binary header serialization.
    pub fn targets(&self) -> &[impl AsRef<str> + TargetInfoAccess] {
        &self.targets
    }

    // r[impl bam_writer.header_from_template]
    /// Create a mutable copy from an existing header for use with BAM writing.
    pub fn from_template(other: &BamHeader) -> Self {
        Self {
            header_text: other.header_text.clone(),
            targets: other.targets.clone(),
            name_to_tid: other.name_to_tid.clone(),
        }
    }

    // r[impl bam_writer.header_add_pg]
    /// Append a `@PG` header record. If `pg.pp` is `None`, it is auto-set to the ID of the
    /// last existing `@PG` line (forming a provenance chain per [SAM1] §1.3).
    pub fn add_pg(&mut self, pg: PgRecord) {
        let pp = pg.pp.or_else(|| self.last_pg_id());

        let mut line = format!("@PG\tID:{}", pg.id);
        if let Some(pn) = &pg.pn {
            line.push_str(&format!("\tPN:{pn}"));
        }
        if let Some(vn) = &pg.vn {
            line.push_str(&format!("\tVN:{vn}"));
        }
        if let Some(cl) = &pg.cl {
            line.push_str(&format!("\tCL:{cl}"));
        }
        if let Some(ds) = &pg.ds {
            line.push_str(&format!("\tDS:{ds}"));
        }
        if let Some(pp_id) = &pp {
            line.push_str(&format!("\tPP:{pp_id}"));
        }
        line.push('\n');

        self.header_text.push_str(&line);
    }

    /// Find the ID of the last `@PG` line in the header text.
    fn last_pg_id(&self) -> Option<SmolStr> {
        self.header_text
            .lines()
            .rev()
            .find(|l| l.starts_with("@PG"))
            .and_then(|line| line.split('\t').find_map(|f| f.strip_prefix("ID:")).map(SmolStr::new))
    }

    // r[impl unified.sort_order]
    pub fn validate_sort_order(&self) -> Result<(), BamHeaderError> {
        let hd_line = match self.header_text.lines().find(|l| l.starts_with("@HD")) {
            Some(line) => line,
            None => return Ok(()),
        };

        let sort_order = match hd_line.split('\t').find_map(|f| f.strip_prefix("SO:")) {
            Some(so) => so,
            None => return Ok(()),
        };

        match sort_order {
            "unsorted" | "queryname" => {
                Err(BamHeaderError::UnsupportedSortOrder { sort_order: SmolStr::new(sort_order) })
            }
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on controlled values")]
mod tests {
    use super::*;

    fn make_bgzf_block(data: &[u8]) -> Vec<u8> {
        let mut compressor =
            libdeflater::Compressor::new(libdeflater::CompressionLvl::new(1).unwrap());
        let bound = compressor.deflate_compress_bound(data.len());
        let mut compressed = vec![0u8; bound];
        let compressed_len = compressor.deflate_compress(data, &mut compressed).unwrap();
        compressed.truncate(compressed_len);

        let mut crc = libdeflater::Crc::new();
        crc.update(data);

        let bsize = (18 + compressed_len + 8 - 1) as u16;
        let mut block = Vec::with_capacity(18 + compressed_len + 8);
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
        block.extend_from_slice(&[0; 4]);
        block.push(0);
        block.push(0xff);
        block.extend_from_slice(&6u16.to_le_bytes());
        block.extend_from_slice(&[b'B', b'C', 2, 0]);
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&compressed);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());
        block
    }

    fn make_bgzf_eof() -> Vec<u8> {
        vec![
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ]
    }

    // r[verify unified.header_from_sam_text]
    #[test]
    fn from_sam_text_parses_sq_lines() {
        let text = "@HD\tVN:1.6\tSO:coordinate\n\
                     @SQ\tSN:chr1\tLN:248956422\n\
                     @SQ\tSN:chr2\tLN:242193529\n\
                     @RG\tID:sample1\tSM:sample1\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        assert_eq!(header.target_count(), 2);
        assert_eq!(header.tid("chr1"), Some(0));
        assert_eq!(header.tid("chr2"), Some(1));
        assert_eq!(header.target_name(0), Some("chr1"));
        assert_eq!(header.target_len(0), Some(248_956_422));
        assert_eq!(header.target_len(1), Some(242_193_529));
        assert!(header.header_text().contains("@HD"));
    }

    #[test]
    fn from_sam_text_errors_on_no_sq() {
        let text = "@HD\tVN:1.6\n@RG\tID:x\tSM:x\n";
        let err = BamHeader::from_sam_text(text).unwrap_err();
        assert!(matches!(err, BamHeaderError::NoSequences));
    }

    #[test]
    fn from_sam_text_errors_on_malformed_sq() {
        let text = "@SQ\tSN:chr1\n"; // missing LN
        let err = BamHeader::from_sam_text(text).unwrap_err();
        match err {
            BamHeaderError::MissingLn { line } => {
                assert!(line.contains("SN:chr1"), "line field should contain the SQ line: {line}")
            }
            other => panic!("expected MissingLn, got {other:?}"),
        }
    }

    #[test]
    fn from_sam_text_errors_on_missing_sn() {
        let text = "@SQ\tLN:1000\n"; // missing SN
        let err = BamHeader::from_sam_text(text).unwrap_err();
        match err {
            BamHeaderError::MissingSn { line } => {
                assert!(line.contains("LN:1000"), "line field should contain the SQ line: {line}")
            }
            other => panic!("expected MissingSn, got {other:?}"),
        }
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_coordinate_accepted() {
        let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        header.validate_sort_order().unwrap();
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_absent_accepted() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        header.validate_sort_order().unwrap();
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_no_hd_accepted() {
        let text = "@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        header.validate_sort_order().unwrap();
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_unsorted_rejected() {
        let text = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        let err = header.validate_sort_order().unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("unsorted"), "error should mention the sort order: {msg}");
        assert!(msg.contains("samtools sort"), "error should suggest samtools sort: {msg}");
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_queryname_rejected() {
        let text = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        let err = header.validate_sort_order().unwrap_err();
        assert!(matches!(err, BamHeaderError::UnsupportedSortOrder { .. }));
        assert!(err.to_string().contains("queryname"));
    }

    // r[verify unified.sort_order]
    #[test]
    fn sort_order_unknown_accepted() {
        let text = "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();
        header.validate_sort_order().unwrap();
    }

    #[test]
    fn resolve_contig_returns_tid_and_len() {
        let text = "@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n";
        let header = BamHeader::from_sam_text(text).unwrap();

        let info = header.resolve_contig("chr1").unwrap();
        assert_eq!(info, ContigInfo { tid: 0, len: 248_956_422 });

        let info = header.resolve_contig("chr2").unwrap();
        assert_eq!(info, ContigInfo { tid: 1, len: 242_193_529 });
    }

    #[test]
    fn resolve_contig_not_found() {
        let text = "@SQ\tSN:chr1\tLN:1000\n";
        let header = BamHeader::from_sam_text(text).unwrap();

        let err = header.resolve_contig("chrX").unwrap_err();
        assert!(matches!(err, BamHeaderError::ContigNotFound { .. }));
        assert!(err.to_string().contains("chrX"));
    }

    // r[verify bam.header.non_negative_lengths]
    #[test]
    fn parse_rejects_negative_l_text() {
        use std::io::Write;

        // Build a minimal BGZF file with BAM magic followed by l_text = -1.
        let mut payload = Vec::new();
        payload.extend_from_slice(b"BAM\x01"); // magic
        payload.extend_from_slice(&(-1i32).to_le_bytes()); // l_text = -1

        let block = make_bgzf_block(&payload);
        let eof = make_bgzf_eof();

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("neg_ltext.bam");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(&block).unwrap();
        f.write_all(&eof).unwrap();
        drop(f);

        let mut bgzf = BgzfReader::open(&path).unwrap();
        let result = BamHeader::parse(&mut bgzf);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            matches!(err, BamHeaderError::NegativeLength { field: "l_text", value: -1 }),
            "expected NegativeLength for l_text, got {err:?}"
        );
    }

    // r[verify bam_writer.header_from_template]
    #[test]
    fn from_template_clones_all_fields() {
        let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n";
        let original = BamHeader::from_sam_text(text).unwrap();
        let copy = BamHeader::from_template(&original);

        assert_eq!(copy.target_count(), 2);
        assert_eq!(copy.tid("chr1"), Some(0));
        assert_eq!(copy.tid("chr2"), Some(1));
        assert_eq!(copy.target_len(0), Some(1000));
        assert_eq!(copy.target_len(1), Some(2000));
        assert_eq!(copy.header_text(), original.header_text());
    }

    // r[verify bam_writer.header_add_pg]
    #[test]
    fn add_pg_appends_line() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
        let mut header = BamHeader::from_sam_text(text).unwrap();
        header.add_pg(PgRecord {
            id: SmolStr::new("seqair"),
            pn: Some(SmolStr::new("seqair")),
            vn: Some(SmolStr::new("0.1.0")),
            cl: None,
            ds: None,
            pp: None,
        });

        let text = header.header_text();
        assert!(text.contains("@PG\tID:seqair\tPN:seqair\tVN:0.1.0\n"));
    }

    #[test]
    fn add_pg_auto_chains_pp() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:bwa\tPN:bwa\n";
        let mut header = BamHeader::from_sam_text(text).unwrap();
        header.add_pg(PgRecord {
            id: SmolStr::new("seqair"),
            pn: Some(SmolStr::new("seqair")),
            vn: None,
            cl: None,
            ds: None,
            pp: None, // should auto-set to "bwa"
        });

        let text = header.header_text();
        assert!(text.contains("PP:bwa"), "PP should be auto-set to last PG ID: {text}");
    }

    #[test]
    fn add_pg_explicit_pp_overrides_auto() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:bwa\tPN:bwa\n";
        let mut header = BamHeader::from_sam_text(text).unwrap();
        header.add_pg(PgRecord {
            id: SmolStr::new("seqair"),
            pn: None,
            vn: None,
            cl: None,
            ds: None,
            pp: Some(SmolStr::new("custom")),
        });

        let text = header.header_text();
        assert!(text.contains("PP:custom"), "explicit PP should override auto: {text}");
    }

    #[test]
    fn from_sam_text_matches_bam_header() {
        let bam_path =
            std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
        let bam_header = BamHeader::from_bam_path(bam_path).unwrap();
        let sam_header = BamHeader::from_sam_text(bam_header.header_text()).unwrap();

        assert_eq!(bam_header.target_count(), sam_header.target_count());
        for tid in 0..bam_header.target_count() as u32 {
            assert_eq!(bam_header.target_name(tid), sam_header.target_name(tid));
            assert_eq!(bam_header.target_len(tid), sam_header.target_len(tid));
        }
    }
}
