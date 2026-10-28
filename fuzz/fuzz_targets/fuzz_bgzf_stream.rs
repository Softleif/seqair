//! `BgzfStream` on arbitrary bytes: never panics or hangs, and inflating on
//! worker threads yields exactly what inflating on the caller does — the
//! same bytes, then the same error (`r[bgzf.stream.threads]`) — after which
//! every read is `Poisoned` (`r[bgzf.stream.poison]`).
//!
//! The first byte picks the thread count (1..=4) and EOF-marker handling,
//! the second the read size; the rest is the stream.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::BgzfError;
use seqair::bam::bgzf::{BgzfStream, EofMarker};

/// Every byte up to the end or the first error, read `chunk` bytes at a time,
/// and that error's description.
fn read_all(
    input: &[u8],
    threads: usize,
    eof: EofMarker,
    chunk: usize,
) -> (Vec<u8>, Option<String>) {
    let builder = BgzfStream::builder(input).eof_marker(eof);
    let mut stream = match u8::try_from(threads).ok().and_then(std::num::NonZeroU8::new) {
        Some(n) => builder.threads(n),
        None => builder,
    }
    .build();
    let mut out = Vec::new();
    let mut buf = vec![0u8; chunk];
    loop {
        match stream.read_up_to(&mut buf) {
            Ok(0) => {
                assert!(stream.is_finished(), "end of stream not finished");
                return (out, None);
            }
            Ok(n) => out.extend_from_slice(buf.get(..n).unwrap_or_default()),
            Err(e) => {
                assert!(
                    matches!(stream.read_up_to(&mut buf), Err(BgzfError::Poisoned)),
                    "not poisoned after {e:?}"
                );
                return (out, Some(format!("{e:?}")));
            }
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let [mode, size, input @ ..] = data else { return };
    let threads = usize::from(mode % 4) + 1;
    let eof = if mode & 0x80 == 0 { EofMarker::Required } else { EofMarker::Optional };
    let chunk = usize::from(*size).saturating_mul(257).max(1);
    let serial = read_all(input, 0, eof, chunk);
    let parallel = read_all(input, threads, eof, chunk);
    assert!(serial.0 == parallel.0, "data differs with {threads} threads");
    assert_eq!(serial.1, parallel.1, "error differs with {threads} threads");
});
