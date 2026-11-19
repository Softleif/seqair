//! A VCF/BCF file opened from bytes and every field of every record read
//! through `FullRecord`, then through narrower types, and every record
//! copied into a writer with the file's own header (`r[vcf_reader.fuzz+2]`):
//! nothing panics, however the file is malformed. BCF and VCF text alike:
//! the format is detected from the content.
//!
//! The first byte picks the framing: even, the rest is BGZF-compressed here
//! first (so mutations reach the records, not the deflate stream); odd, the
//! rest is the file as is.
#![no_main]

// r[verify vcf_reader.fuzz+2]

use libfuzzer_sys::fuzz_target;
use seqair::io::BgzfWriter;
use seqair::vcf::VcfFile;
use seqair_fuzz::vcf_fields::{copy_all, read_all, read_narrow};

fuzz_target!(|data: &[u8]| {
    let Some((&mode, rest)) = data.split_first() else { return };
    let bytes = if mode % 2 == 0 {
        let mut w = BgzfWriter::new(Vec::new());
        if w.write_all(rest).is_err() {
            return;
        }
        match w.finish() {
            Ok(b) => b,
            Err(_) => return,
        }
    } else {
        rest.to_vec()
    };
    if let Ok(file) = VcfFile::from_bytes(bytes, None) {
        std::hint::black_box(read_all(&file) + read_narrow(&file) + copy_all(&file));
    }
});
