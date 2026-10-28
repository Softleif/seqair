//! `BamStreamReader` on arbitrary bytes as a BAM stream. Every record it
//! yields is a valid `RawBamRecord` (`r[bam.raw_record.validation]`): no
//! accessor panics, and writing the records back with `write_raw` and
//! streaming the result gives the same payloads (`r[bam_writer.write_raw]`).
//!
//! The first byte picks the thread count and EOF-marker handling; the rest
//! is the BGZF-compressed BAM.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::RawBamRecord;
use seqair::bam::bgzf::{BgzfStream, EofMarker};
use seqair::bam::stream::{BamStreamError, BamStreamReader};
use seqair::bam::writer::BamWriterBuilder;

/// Call every accessor; the values themselves only have to exist.
fn touch(record: &RawBamRecord) -> usize {
    let mut n = record.qname().len()
        + record.cigar().count()
        + record.cigar_len()
        + record.cigar_bytes().len()
        + record.packed_seq().len()
        + record.bases().len()
        + usize::from(record.bin())
        + usize::from(record.mapq().as_byte())
        + usize::from(record.flags().raw());
    let mut seq = Vec::new();
    record.decode_seq_into(&mut seq);
    assert_eq!(seq.len(), record.seq_len());
    n += record.qual().map_or(0, <[_]>::len);
    n += record.aux().iter_checked().take(1000).count();
    n += record.aux().iter_tags().take(1000).count();
    n += usize::from(record.has_cg_placeholder());
    n += record.end_pos().map_or(0, |p| p.as_u64() as usize & 1);
    n += record.tid().map_or(0, |t| t.as_usize() & 1);
    n += record.pos().map_or(0, |p| p.as_u64() as usize & 1);
    let mate = record.mate_fields();
    n += usize::from(mate.tid.is_some()) + usize::from(mate.pos.is_some());
    n + format!("{record:?}").len()
}

fuzz_target!(|data: &[u8]| {
    let [mode, input @ ..] = data else { return };
    let threads = usize::from(mode % 3);
    let eof = if mode & 0x80 == 0 { EofMarker::Required } else { EofMarker::Optional };
    let Ok(stream) = BgzfStream::builder(input).threads(threads).eof_marker(eof).build() else {
        return;
    };
    let Ok(mut reader) = BamStreamReader::from_bgzf(stream) else { return };

    let mut records = Vec::new();
    let mut record = RawBamRecord::default();
    for _ in 0..10_000 {
        match reader.read_record(&mut record) {
            Ok(true) => {
                std::hint::black_box(touch(&record));
                records.push(record.clone());
            }
            Ok(false) => break,
            Err(BamStreamError::InvalidRecord { .. }) => {
                assert_eq!(record, RawBamRecord::default(), "failed refill left a partial record");
            }
            Err(_) => {
                assert!(matches!(reader.read_record(&mut record), Err(_)), "not poisoned");
                break;
            }
        }
    }

    let header = reader.into_header();
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header)
        .build()
        .expect("a header seqair read must be one it can write");
    for record in &records {
        writer.write_raw(record).expect("write_raw of a validated record");
    }
    let (bytes, _) = writer.finish().expect("finish");
    let mut again = BamStreamReader::new(bytes.as_slice()).expect("re-read header");
    let mut back = RawBamRecord::default();
    for (i, record) in records.iter().enumerate() {
        assert!(again.read_record(&mut back).expect("re-read record"), "record {i} missing");
        assert!(back.as_payload() == record.as_payload(), "record {i} changed");
    }
    assert!(!again.read_record(&mut back).expect("re-read end"), "extra records");
});
