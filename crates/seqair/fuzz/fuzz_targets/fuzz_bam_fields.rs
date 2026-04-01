#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::record::compute_end_pos_from_raw;
use seqair::bam::{aux, cigar, seq};

#[derive(Arbitrary, Debug)]
struct SeqInput {
    data: Vec<u8>,
    len: u16,
}

#[derive(Arbitrary, Debug)]
struct AuxTagInput {
    data: Vec<u8>,
    tag: [u8; 2],
}

#[derive(Arbitrary, Debug)]
enum FuzzInput {
    Cigar(Vec<u8>),
    Seq(SeqInput),
    AuxFind(AuxTagInput),
    AuxIter(Vec<u8>),
    EndPos(Vec<u8>),
}

fuzz_target!(|input: FuzzInput| {
    match input {
        FuzzInput::Cigar(data) => {
            // CIGAR bytes must be 4-byte aligned for BAM packed CIGAR ops
            let aligned = &data[..data.len() - data.len() % 4];
            let _ = cigar::calc_matches_indels(aligned);
        }
        FuzzInput::Seq(s) => {
            let len = s.len as usize;
            let _ = seq::decode_seq(&s.data, len);
            let _ = seq::decode_bases(&s.data, len);
        }
        FuzzInput::AuxFind(a) => {
            let _ = aux::find_tag(&a.data, a.tag);
        }
        FuzzInput::AuxIter(data) => {
            for tag in aux::iter_tags(&data) {
                let _ = tag;
            }
        }
        FuzzInput::EndPos(data) => {
            let _ = compute_end_pos_from_raw(&data);
        }
    }
});
