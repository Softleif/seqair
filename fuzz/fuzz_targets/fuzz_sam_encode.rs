#![no_main]
//! SAM text → BAM records: `SamStreamReader` over arbitrary input,
//! `SamRecordEncoder` on its first line against a fixed header, and
//! `parse_cigar` on the whole input. Never panics, and every record that
//! encodes passes `RawBamRecord` validation again unchanged.

use libfuzzer_sys::fuzz_target;
use seqair::bam::{BamHeader, CigarStr, RawBamRecord};
use seqair::sam::{SamRecordEncoder, SamStreamReader, parse_cigar};

fn revalidate(record: &RawBamRecord, header: &BamHeader) {
    let again = RawBamRecord::from_payload(record.as_payload().to_vec(), header)
        .expect("an encoded record re-validates");
    assert_eq!(again.as_payload(), record.as_payload());
    // Every accessor must work on it.
    let mut seq = Vec::new();
    record.decode_seq_into(&mut seq);
    assert_eq!(seq.len(), record.seq_len());
    std::hint::black_box((
        record.cigar().count(),
        record.end_pos(),
        record.qual(),
        record.aux().iter_checked().count(),
    ));
}

fuzz_target!(|data: &[u8]| {
    if let Ok(mut reader) = SamStreamReader::new(data) {
        let mut record = RawBamRecord::default();
        // Bounded: each call consumes at least one line or ends the stream.
        for _ in 0..10_000 {
            match reader.read_record(&mut record) {
                Ok(true) => revalidate(&record, reader.header()),
                Ok(false) => break,
                Err(_) => {
                    if reader.read_record(&mut record).is_err() {
                        break;
                    }
                }
            }
        }
    }

    let header = BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:2147483647\n@SQ\tSN:c2\tLN:100\n")
        .expect("fixed header");
    let line = data.split(|&b| b == b'\n').next().unwrap_or_default();
    let mut record = RawBamRecord::default();
    if SamRecordEncoder::new().encode(line, &header, &mut record).is_ok() {
        revalidate(&record, &header);
    }

    let mut ops = Vec::new();
    if parse_cigar(data, &mut ops).is_ok()
        && ops.iter().all(|op| op.op_code() <= 8)
        && !ops.is_empty()
    {
        let text = CigarStr(&ops).to_string();
        let mut again = Vec::new();
        parse_cigar(text.as_bytes(), &mut again).expect("CigarStr output parses");
        assert_eq!(again, ops);
    }
});
