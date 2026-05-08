use crate::bam::bgzf::BgzfReader;

use super::ReaderError;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy)]
pub enum Format {
    Bam,
    Sam,
    Cram,
}

// r[impl unified.detect_format]
pub fn detect(path: &Path) -> Result<Format, ReaderError> {
    let mut magic = [0u8; 4];
    let file_len = std::fs::metadata(path)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?
        .len();

    if file_len < 4 {
        return Err(
            FormatDetectionError::FileTooShort { path: path.to_path_buf(), len: file_len }.into()
        );
    }

    let mut file = std::fs::File::open(path)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?;

    use std::io::Read;
    file.read_exact(&mut magic)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?;

    // CRAM: magic bytes "CRAM"
    if magic == *b"CRAM" {
        return Ok(Format::Cram);
    }

    // Uncompressed SAM: starts with '@'
    if magic[0] == b'@' {
        return Err(FormatDetectionError::UncompressedSam { path: path.to_path_buf() }.into());
    }

    // BGZF: gzip magic 1f 8b
    if magic[0] == 0x1f && magic[1] == 0x8b {
        // Decompress first block to check BAM magic vs SAM text
        let mut bgzf = BgzfReader::open(path)
            .map_err(|source| FormatDetectionError::NotBgzf { path: path.to_path_buf(), source })?;

        let mut first_bytes = [0u8; 4];
        bgzf.read_exact_into(&mut first_bytes)
            .map_err(|source| FormatDetectionError::NotBgzf { path: path.to_path_buf(), source })?;

        if first_bytes == *b"BAM\x01" {
            return Ok(Format::Bam);
        }

        if first_bytes[0] == b'@' {
            return Ok(Format::Sam);
        }

        return Err(FormatDetectionError::UnrecognizedBgzfContent {
            path: path.to_path_buf(),
            magic: first_bytes,
        }
        .into());
    }

    Err(FormatDetectionError::UnrecognizedFormat { path: path.to_path_buf(), magic }.into())
}

// r[impl io.non_exhaustive_enums]
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum FormatDetectionError {
    #[error("file too short to determine format ({len} bytes)")]
    FileTooShort { path: PathBuf, len: u64 },

    #[error(
        "uncompressed SAM detected. \
         Compress with `bgzip {path}` then index with `tabix -p sam {path}.gz`"
    )]
    UncompressedSam { path: PathBuf },

    #[error(
        "could not decompress first BGZF block from {path}. \
         If this is a plain gzip file, use `bgzip` instead of `gzip`"
    )]
    NotBgzf { path: PathBuf, source: crate::bam::bgzf::BgzfError },

    #[error(
        "BGZF-compressed file {path} does not contain BAM or SAM data \
         (first decompressed bytes: {magic:?}). \
         Supported formats: BAM (.bam), bgzf-compressed SAM (.sam.gz), CRAM (.cram)."
    )]
    UnrecognizedBgzfContent { path: PathBuf, magic: [u8; 4] },

    #[error(
        "unrecognized file format for {path} (magic bytes: {magic:02x?}). \
         Supported formats: BAM (.bam), bgzf-compressed SAM (.sam.gz), CRAM (.cram)."
    )]
    UnrecognizedFormat { path: PathBuf, magic: [u8; 4] },
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_bgzf(content: &[u8]) -> Vec<u8> {
        use bgzf::{CompressionLevel, Writer as BgzfWriter};
        let mut output = Vec::new();
        let mut writer =
            BgzfWriter::new(&mut output, CompressionLevel::new(1).expect("valid level"));
        writer.write_all(content).expect("write");
        writer.finish().expect("finish");
        output
    }

    // r[verify unified.detect_format]
    #[test]
    fn file_too_short() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"BA").expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::FileTooShort { len: 2, .. } }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn file_empty() {
        let f = NamedTempFile::new().expect("tempfile");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::FileTooShort { len: 0, .. } }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn uncompressed_sam() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"@HD\tVN:1.6\n").expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::UncompressedSam { .. } }
            ),
            "unexpected error: {err}"
        );
    }

    // TODO: `FormatDetectionError::NotBgzf` is only triggered when `BgzfReader::open`
    // fails with an I/O error. Because `detect_format` already opens the file once before
    // calling `BgzfReader::open`, this path is not reachable via content manipulation —
    // it would require a race condition (file deleted or permissions revoked between the
    // two opens). Portable unit testing of this variant is not feasible.
    //
    #[test]
    fn not_valid_bgzf_returns_not_bgzf() {
        // gzip magic but invalid BGZF header (byte 3 is 0x00, not 0x04 FEXTRA flag).
        // BgzfReader::open succeeds (pure file open); read_exact_into fails → NotBgzf.
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&[0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]).expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::NotBgzf { .. } }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn short_bgzf_content_returns_not_bgzf() {
        // A valid BGZF block with only 2 bytes of decompressed content — too short
        // to read 4 bytes, triggering NotBgzf.
        let compressed = write_bgzf(b"AB");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::NotBgzf { .. } }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn unrecognized_bgzf_content() {
        // A valid BGZF block whose decompressed content starts with neither
        // "BAM\x01" nor "@" — triggers UnrecognizedBgzfContent.
        let compressed = write_bgzf(b"\xDE\xAD\xBE\xEF and some more data");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format {
                    source: FormatDetectionError::UnrecognizedBgzfContent {
                        magic: [0xDE, 0xAD, 0xBE, 0xEF],
                        ..
                    }
                }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn bgzf_compressed_sam_detected() {
        // Sanity-check the Sam fast-path: "@"-prefixed BGZF content → Format::Sam
        let compressed = write_bgzf(b"@HD\tVN:1.6\n");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        // detect_format should succeed with Format::Sam (no error)
        let result = detect(f.path());
        assert!(matches!(result, Ok(Format::Sam)), "expected Sam, got {result:?}");
    }

    #[test]
    fn unrecognized_format() {
        // 4+ raw bytes that are not gzip magic, not "@", not "CRAM"
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&[0xDE, 0xAD, 0xBE, 0xEF]).expect("write");
        f.flush().expect("flush");
        let err = detect(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format {
                    source: FormatDetectionError::UnrecognizedFormat {
                        magic: [0xDE, 0xAD, 0xBE, 0xEF],
                        ..
                    }
                }
            ),
            "unexpected error: {err}"
        );
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;
    use std::io::Write as _;
    use tempfile::NamedTempFile;

    proptest! {
        #[test]
        fn fuzz_detect_format(data in proptest::collection::vec(any::<u8>(), 0..256 * 1024)) {
            // Just ensure detect_format doesn't panic on arbitrary input data.
            let mut f = NamedTempFile::new().expect("tempfile");
            f.write_all(data.as_slice()).expect("write");
            f.flush().expect("flush");
            let _ = detect(f.path());
        }
    }
}
