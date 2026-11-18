//! VCF text records under a fixed header that declares every value kind,
//! END, SVLEN and LEN, read through every component and `RawRecord`'s views
//! (`r[vcf_reader.fuzz]`): nothing panics, however a line is malformed.
//!
//! The first byte picks the framing: even, BGZF (lines then span blocks
//! and are copied together); odd, plain text.
#![no_main]

// r[verify vcf_reader.fuzz]

use libfuzzer_sys::fuzz_target;
use seqair::io::BgzfWriter;
use seqair::vcf::VcfFile;
use seqair_fuzz::vcf_fields::{read_all, read_narrow};

const HEADER: &str = "##fileformat=VCFv4.5\n\
##FILTER=<ID=PASS,Description=\"p\">\n\
##FILTER=<ID=q10,Description=\"q\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"s\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"b\">\n\
##INFO=<ID=S,Number=.,Type=String,Description=\"s\">\n\
##INFO=<ID=C,Number=1,Type=Character,Description=\"c\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"l\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"g\">\n\
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"f\">\n\
##contig=<ID=chr1,length=1000000>\n\
##contig=<ID=chr2>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n";

fuzz_target!(|data: &[u8]| {
    let Some((&mode, lines)) = data.split_first() else { return };
    let mut text = HEADER.as_bytes().to_vec();
    text.extend_from_slice(lines);
    let bytes = if mode % 2 == 0 {
        let mut w = BgzfWriter::new(Vec::new());
        if w.write_all(&text).is_err() {
            return;
        }
        match w.finish() {
            Ok(b) => b,
            Err(_) => return,
        }
    } else {
        text
    };
    let file = VcfFile::from_bytes(bytes, None).expect("the fixed header parses");
    std::hint::black_box(read_all(&file) + read_narrow(&file));
});
