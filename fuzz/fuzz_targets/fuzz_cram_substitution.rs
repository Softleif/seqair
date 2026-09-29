#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::compression_header::SubstitutionMatrix;

#[derive(Arbitrary, Debug)]
struct SubstitutionInput {
    bytes: [u8; 5],
    ref_bases: Vec<u8>,
    codes: Vec<u8>,
}

fuzz_target!(|input: SubstitutionInput| {
    let matrix = SubstitutionMatrix::parse(&input.bytes);
    // Exercise substitute() with arbitrary ref bases and codes
    for (&ref_base, &code) in input.ref_bases.iter().zip(input.codes.iter()) {
        let _ = matrix.substitute(ref_base, code);
    }
});
