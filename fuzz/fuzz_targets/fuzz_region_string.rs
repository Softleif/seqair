#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair_types::RegionString;
use std::str::FromStr;

fuzz_target!(|data: &[u8]| {
    let Ok(s) = std::str::from_utf8(data) else {
        return;
    };
    let _ = RegionString::from_str(s);
});
