//! Types for the `fuzz_reader_indexed` target — shared between the fuzz target
//! and the seed generator so seeds always match the deserialization format.
//!
//! Uses a simple hand-rolled binary format instead of Arbitrary derive, because
//! Arbitrary's length encoding wastes most of the seed data on bad Vec splits.

/// Header layout (16 bytes):
///   [0]    format: 0=BAM, 1=SAM(BGZF), 2=CRAM, 3=PlainSAM
///   [1..5] data1_len: u32 LE (BAM/SAM/CRAM file bytes)
///   [5..9] data2_len: u32 LE (BAI/TBI/CRAI index bytes; unused for PlainSAM)
///   [9..11] fai_len: u16 LE (FAI text bytes)
///   [11..13] gzi_len: u16 LE (GZI binary bytes)
///   [13..16] reserved (0)
///
/// Body: [data1][data2][fasta_gz][fai][gzi]
///   fasta_gz gets whatever bytes remain after the other fields.
pub const HEADER_SIZE: usize = 16;

#[derive(Debug)]
pub enum Format {
    Bam,
    /// BGZF-compressed SAM with tabix index.
    Sam,
    Cram,
    /// Plain (uncompressed) SAM — no BGZF, no index, linear scan.
    PlainSam,
}

#[derive(Debug)]
pub struct Input<'a> {
    pub format: Format,
    pub data1: &'a [u8], // alignment file (BAM/SAM/CRAM)
    pub data2: &'a [u8], // index file (BAI/TBI/CRAI)
    pub fasta_gz: &'a [u8],
    pub fai: &'a [u8],
    pub gzi: &'a [u8],
}

impl<'a> Input<'a> {
    /// Parse from raw bytes. Returns None if too short or sizes don't fit.
    pub fn parse(data: &'a [u8]) -> Option<Self> {
        if data.len() < HEADER_SIZE {
            return None;
        }

        let format = match data[0] % 4 {
            0 => Format::Bam,
            1 => Format::Sam,
            2 => Format::Cram,
            _ => Format::PlainSam,
        };

        let data1_len = u32::from_le_bytes([data[1], data[2], data[3], data[4]]) as usize;
        let data2_len = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
        let fai_len = u16::from_le_bytes([data[9], data[10]]) as usize;
        let gzi_len = u16::from_le_bytes([data[11], data[12]]) as usize;

        let body = &data[HEADER_SIZE..];
        let fixed_total = data1_len + data2_len + fai_len + gzi_len;
        if fixed_total > body.len() {
            return None;
        }

        let mut off = 0;
        let data1 = body.get(off..off + data1_len)?;
        off += data1_len;
        let data2 = body.get(off..off + data2_len)?;
        off += data2_len;

        // fasta_gz gets everything between data2 and fai
        let fasta_end = body.len().checked_sub(fai_len + gzi_len)?;
        let fasta_gz = body.get(off..fasta_end)?;

        let fai = body.get(fasta_end..fasta_end + fai_len)?;
        let gzi = body.get(fasta_end + fai_len..)?;

        Some(Input { format, data1, data2, fasta_gz, fai, gzi })
    }

    /// Serialize an Input into bytes matching the parse format.
    pub fn encode(
        format: u8,
        data1: &[u8],
        data2: &[u8],
        fasta_gz: &[u8],
        fai: &[u8],
        gzi: &[u8],
    ) -> Vec<u8> {
        let mut buf = Vec::with_capacity(
            HEADER_SIZE + data1.len() + data2.len() + fasta_gz.len() + fai.len() + gzi.len(),
        );

        // Header
        buf.push(format);
        buf.extend_from_slice(&(data1.len() as u32).to_le_bytes());
        buf.extend_from_slice(&(data2.len() as u32).to_le_bytes());
        buf.extend_from_slice(&(fai.len() as u16).to_le_bytes());
        buf.extend_from_slice(&(gzi.len() as u16).to_le_bytes());
        buf.extend_from_slice(&[0, 0, 0]); // reserved

        // Body
        buf.extend_from_slice(data1);
        buf.extend_from_slice(data2);
        buf.extend_from_slice(fasta_gz);
        buf.extend_from_slice(fai);
        buf.extend_from_slice(gzi);

        buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let data1 = b"hello BAM";
        let data2 = b"hello BAI";
        let fasta = b"hello FASTA";
        let fai = b"chr1\t100\t5\t70\t71\n";
        let gzi = b"\x01\x00\x00\x00\x00\x00\x00\x00";

        let encoded = Input::encode(0, data1, data2, fasta, fai, gzi);
        let parsed = Input::parse(&encoded).unwrap();

        assert!(matches!(parsed.format, Format::Bam));
        assert_eq!(parsed.data1, data1);
        assert_eq!(parsed.data2, data2);
        assert_eq!(parsed.fasta_gz, fasta);
        assert_eq!(parsed.fai, fai);
        assert_eq!(parsed.gzi, gzi);
    }
}
