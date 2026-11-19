//! VCF header text parsed (`r[vcf_reader.fuzz]`): nothing panics, and a
//! header that parses writes BCF text that parses again, with the same
//! dictionaries and samples, and VCF text that writes itself again.
#![no_main]

// r[verify vcf_reader.fuzz]

use libfuzzer_sys::fuzz_target;
use seqair::vcf::VcfHeader;

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else { return };
    let Ok(header) = VcfHeader::parse(text) else { return };
    let again = VcfHeader::parse(&header.to_bcf_text()).expect("written header text parses");
    assert_eq!(again.string_map(), header.string_map());
    let contigs = |h: &VcfHeader| -> Vec<_> {
        h.contigs().map(|(c, def)| (c.name().clone(), c.index(), def.clone())).collect()
    };
    assert_eq!(contigs(&again), contigs(&header));
    assert_eq!(again.samples(), header.samples());
    let vcf = header.to_vcf_text();
    let again = VcfHeader::parse(&vcf).expect("written VCF header text parses");
    assert_eq!(again.to_vcf_text(), vcf, "VCF text is a fixed point");
});
