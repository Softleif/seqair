#![no_main]

//! Differential fuzzing of the CRAM arithmetic coder (method 6): the
//! production decoder against the spec-transliterating reference decoder,
//! directly and inside tok3. Where both accept an input they must produce
//! the same bytes; their verdicts may differ (see the reference's docs).

use libfuzzer_sys::fuzz_target;
use seqair::cram::{arith, tok3};

fuzz_target!(|data: &[u8]| {
    // The size only matters to a stream without its own (NOSZ).
    let size = data.len();
    if let (Ok(ours), Some(theirs)) = (arith::decode(data, size), arith::reference::decode(data, size))
    {
        assert_eq!(ours, theirs, "production and reference arith decoders disagree");
    }

    // The same bytes as a tok3 block, its streams arith-coded.
    if data.get(8).is_some_and(|&use_arith| use_arith != 0)
        && let (Ok(ours), Ok(theirs)) = (tok3::decode(data), tok3::decode_with_arith_reference(data))
    {
        assert_eq!(ours, theirs, "tok3 over the two arith decoders disagrees");
    }
});
