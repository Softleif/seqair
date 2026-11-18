//! Arbitrary bytes as a BAM stream; the records it yields are written as
//! CRAM against a reference made up for the header's targets, read back,
//! compared with what was written (`r[cram_writer.fidelity]`) and formatted
//! as SAM.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::stream::BamStreamReader;
use seqair::cram::writer::{CramProfile, CramVersion};
use seqair_fuzz::cram_roundtrip;

fuzz_target!(|data: &[u8]| {
    let Ok(mut reader) = BamStreamReader::new(data) else { return };
    let records: Vec<_> = reader.records().map_while(Result::ok).take(200).collect();
    let header = reader.into_header();
    // A reference for every target but every third, whose records embed
    // one built from them, if they are small enough to write.
    let targets: Vec<_> = header.targets().collect();
    if targets.iter().map(|t| t.length).sum::<u64>() > 200_000 {
        return;
    }
    let contigs: Vec<(String, Vec<u8>)> = targets
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 3 != 2)
        .map(|(_, t)| {
            let bases = (0..t.length).map(|i| b"ACGTTGCAN"[(i % 9) as usize]).collect();
            (t.name.to_string(), bases)
        })
        .collect();
    if contigs.iter().any(|(name, _)| name.contains(|c: char| c.is_whitespace())) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let fasta = cram_roundtrip::write_fasta(dir.path(), &contigs);
    cram_roundtrip::check(
        &header,
        &fasta,
        &records,
        0,
        10_000,
        CramVersion::V3_1(CramProfile::Normal),
    );
});
