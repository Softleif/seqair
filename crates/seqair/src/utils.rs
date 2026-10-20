use tracing::Level;

pub trait TraceErr<T> {
    /// Like [`Option::and`], but traces a message when the value is `None`.
    fn trace_err(self, message: &'static str) -> Option<T>;
}

impl<T> TraceErr<T> for Option<T> {
    #[track_caller]
    #[inline(always)]
    fn trace_err(self, message: &'static str) -> Option<T> {
        if tracing::enabled!(Level::TRACE) && self.is_none() {
            tracing::trace!("{message}");
        }
        self
    }
}

pub trait TraceOk<T> {
    /// Like [`Result::ok`], but traces the error before discarding it.
    fn trace_ok(self, message: &'static str) -> Option<T>;
}

impl<T, E: std::fmt::Display> TraceOk<T> for Result<T, E> {
    #[track_caller]
    #[inline(always)]
    fn trace_ok(self, message: &'static str) -> Option<T> {
        match self {
            Ok(v) => Some(v),
            Err(error) => {
                if tracing::enabled!(Level::TRACE) {
                    tracing::trace!(%error, "{message}");
                }
                None
            }
        }
    }
}

// r[impl index.locate]
/// The index files `path` may have, in htslib's lookup order: for each of
/// `exts`, appended to the whole file name, then in place of its last
/// extension (`x.bam.csi`, `x.csi`, `x.bam.bai`, `x.bai`). Only the part after
/// the last dot is ever replaced, so a dotted stem (`sample.v2.bam`) keeps its
/// dots.
pub(crate) fn index_candidates<'a>(
    path: &'a std::path::Path,
    exts: &'a [&'a str],
) -> impl Iterator<Item = std::path::PathBuf> + 'a {
    let with = |base: &std::ffi::OsStr, ext: &str| {
        let mut name = base.to_os_string();
        name.push(".");
        name.push(ext);
        path.with_file_name(name)
    };
    let name = path.file_name();
    let stem = path.extension().and(path.file_stem());
    exts.iter().flat_map(move |ext| {
        let appended = name.map(|name| with(name, ext));
        let replaced = stem.map(|stem| with(stem, ext));
        appended.into_iter().chain(replaced)
    })
}

/// The first of [`index_candidates`] that exists.
pub(crate) fn find_index(path: &std::path::Path, exts: &[&str]) -> Option<std::path::PathBuf> {
    index_candidates(path, exts).find(|candidate| candidate.exists())
}

#[cfg(test)]
mod index_candidate_tests {
    use super::index_candidates;
    use std::path::{Path, PathBuf};

    fn names(path: &str, exts: &[&str]) -> Vec<PathBuf> {
        index_candidates(Path::new(path), exts).collect()
    }

    // r[verify index.locate]
    #[test]
    fn the_order_is_htslibs() {
        let expected: Vec<PathBuf> =
            ["d/x.bam.csi", "d/x.csi", "d/x.bam.bai", "d/x.bai"].map(PathBuf::from).into();
        assert_eq!(names("d/x.bam", &["csi", "bai"]), expected);
    }

    // r[verify index.locate]
    #[test]
    fn a_dotted_stem_keeps_its_dots() {
        let expected: Vec<PathBuf> =
            ["sample.v2.bam.crai", "sample.v2.crai"].map(PathBuf::from).into();
        assert_eq!(names("sample.v2.bam", &["crai"]), expected);
    }

    // r[verify index.locate]
    #[test]
    fn a_name_without_extension_has_only_appended_names() {
        assert_eq!(names("dir.v1/sample", &["bai"]), vec![PathBuf::from("dir.v1/sample.bai")]);
    }
}
