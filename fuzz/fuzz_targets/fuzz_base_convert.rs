#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair_types::Base;

fuzz_target!(|data: Vec<u8>| {
    // Fuzz convert_ascii_in_place
    let mut in_place = data.clone();
    Base::convert_ascii_in_place(&mut in_place);

    // Fuzz from_ascii_vec (consumes the vec)
    let result = Base::from_ascii_vec(data.clone());

    // Both paths must agree: the in_place bytes must match from_ascii_vec discriminants
    assert_eq!(in_place.len(), result.len());
    for (byte, base) in in_place.iter().zip(result.iter()) {
        assert_eq!(*byte, *base as u8);
    }
});
