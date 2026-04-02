//! Parse the CRAM compression header from a container. [`CompressionHeader`] holds the
//! encoding descriptor for each data series, used by the slice decoder to reconstruct records.

use super::{
    encoding::{ByteArrayEncoding, ByteEncoding, IntEncoding},
    reader::CramError,
    varint,
};
use rustc_hash::FxHashMap;

/// Parsed compression header from a CRAM container.
#[derive(Debug)]
pub struct CompressionHeader {
    pub preservation: PreservationMap,
    pub data_series: DataSeriesEncodings,
    pub tag_encodings: FxHashMap<i32, ByteArrayEncoding>,
}

/// Preservation map from the compression header.
#[derive(Debug)]
pub struct PreservationMap {
    pub read_names_included: bool,
    pub ap_delta: bool,
    pub reference_required: bool,
    pub substitution_matrix: SubstitutionMatrix,
    pub tag_dictionary: Vec<Vec<TagDictEntry>>,
}

/// An entry in the tag dictionary: 2-byte tag name + BAM type character.
#[derive(Debug, Clone)]
pub struct TagDictEntry {
    pub tag: [u8; 2],
    pub bam_type: u8,
}

/// Substitution matrix: for each reference base, the 4 possible substitution
/// alternatives ordered by their 2-bit code.
#[derive(Debug, Clone)]
pub struct SubstitutionMatrix {
    /// `matrix[ref_base_index]` gives 4 alternatives.
    /// ref_base_index: 0=A, 1=C, 2=G, 3=T, 4=N
    pub matrix: [[u8; 4]; 5],
}

impl SubstitutionMatrix {
    /// Possible read bases for each reference base (excluding the ref itself).
    const READ_BASES: [[u8; 4]; 5] = [
        [b'C', b'G', b'T', b'N'], // ref A
        [b'A', b'G', b'T', b'N'], // ref C
        [b'A', b'C', b'T', b'N'], // ref G
        [b'A', b'C', b'G', b'N'], // ref T
        [b'A', b'C', b'G', b'T'], // ref N
    ];

    // r[impl cram.compression.substitution_matrix]
    #[allow(
        clippy::indexing_slicing,
        reason = "ref_idx ≤ 4 (5-element zip), inner index = 2-bit value (0–3), matrix is [[u8;4];5]"
    )]
    pub fn parse(bytes: &[u8; 5]) -> Self {
        // Each byte encodes the ordering of substitution alternatives.
        // The 2-bit codes in the byte are positions into the substitutions
        // array. `read_bases[i]` are placed at position `code` in the output.
        let mut matrix = [[0u8; 4]; 5];
        for (ref_idx, (&codes_byte, read_bases)) in
            bytes.iter().zip(Self::READ_BASES.iter()).enumerate()
        {
            matrix[ref_idx][(codes_byte >> 6 & 0x03) as usize] = read_bases[0];
            matrix[ref_idx][(codes_byte >> 4 & 0x03) as usize] = read_bases[1];
            matrix[ref_idx][(codes_byte >> 2 & 0x03) as usize] = read_bases[2];
            matrix[ref_idx][(codes_byte & 0x03) as usize] = read_bases[3];
        }
        Self { matrix }
    }

    /// Look up the substitution base for a given reference base and code.
    pub fn substitute(&self, ref_base: u8, code: u8) -> u8 {
        let ref_idx = match ref_base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        let code_idx = (code & 0x03) as usize;
        self.matrix.get(ref_idx).and_then(|row| row.get(code_idx)).copied().unwrap_or(b'N')
    }
}

/// Encodings for all CRAM data series.
#[derive(Debug)]
pub struct DataSeriesEncodings {
    pub bam_flags: IntEncoding,           // BF
    pub cram_flags: IntEncoding,          // CF
    pub ref_id: IntEncoding,              // RI
    pub read_length: IntEncoding,         // RL
    pub alignment_pos: IntEncoding,       // AP
    pub read_group: IntEncoding,          // RG
    pub read_name: ByteArrayEncoding,     // RN
    pub mate_flags: IntEncoding,          // MF
    pub next_segment_ref: IntEncoding,    // NS
    pub next_mate_pos: IntEncoding,       // NP
    pub template_size: IntEncoding,       // TS
    pub next_fragment: IntEncoding,       // NF
    pub tag_line: IntEncoding,            // TL
    pub feature_count: IntEncoding,       // FN
    pub feature_code: ByteEncoding,       // FC
    pub feature_pos: IntEncoding,         // FP
    pub deletion_length: IntEncoding,     // DL
    pub mapping_quality: IntEncoding,     // MQ
    pub base: ByteEncoding,               // BA
    pub quality_score: ByteEncoding,      // QS
    pub base_sub: ByteEncoding,           // BS
    pub insertion: ByteArrayEncoding,     // IN
    pub soft_clip: ByteArrayEncoding,     // SC
    pub ref_skip: IntEncoding,            // RS
    pub padding: IntEncoding,             // PD
    pub hard_clip: IntEncoding,           // HC
    pub bases_block: ByteArrayEncoding,   // BB
    pub quality_block: ByteArrayEncoding, // QQ
}

impl CompressionHeader {
    /// Parse a compression header from block data.
    pub fn parse(data: &[u8]) -> Result<Self, CramError> {
        let mut cursor: &[u8] = data;

        let preservation = parse_preservation_map(&mut cursor)?;
        let data_series = parse_data_series_map(&mut cursor)?;
        let tag_encodings = parse_tag_encoding_map(&mut cursor)?;

        Ok(CompressionHeader { preservation, data_series, tag_encodings })
    }
}

// r[impl cram.compression.preservation]
fn parse_preservation_map(cursor: &mut &[u8]) -> Result<PreservationMap, CramError> {
    let map_size = varint::read_itf8_from(cursor)
        .ok_or(CramError::Truncated { context: "preservation map size" })?
        as usize;

    let map_data =
        cursor.get(..map_size).ok_or(CramError::Truncated { context: "preservation map data" })?;
    let mut mcur: &[u8] = map_data;

    let num_entries = varint::read_itf8_from(&mut mcur)
        .ok_or(CramError::Truncated { context: "preservation map entry count" })?;

    let mut read_names_included = true;
    let mut ap_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = SubstitutionMatrix { matrix: [[b'N'; 4]; 5] };
    let mut tag_dictionary = Vec::new();

    for _ in 0..num_entries {
        let key0 = *mcur.first().ok_or(CramError::Truncated { context: "preservation map key" })?;
        let key1 = *mcur.get(1).ok_or(CramError::Truncated { context: "preservation map key" })?;
        mcur = mcur
            .get(2..)
            .ok_or(CramError::Truncated { context: "preservation map key advance" })?;

        match (key0, key1) {
            (b'R', b'N') => {
                read_names_included =
                    *mcur.first().ok_or(CramError::Truncated { context: "RN value" })? != 0;
                mcur = mcur.get(1..).ok_or(CramError::Truncated { context: "RN advance" })?;
            }
            (b'A', b'P') => {
                ap_delta = *mcur.first().ok_or(CramError::Truncated { context: "AP value" })? != 0;
                mcur = mcur.get(1..).ok_or(CramError::Truncated { context: "AP advance" })?;
            }
            (b'R', b'R') => {
                reference_required =
                    *mcur.first().ok_or(CramError::Truncated { context: "RR value" })? != 0;
                mcur = mcur.get(1..).ok_or(CramError::Truncated { context: "RR advance" })?;
            }
            (b'S', b'M') => {
                let sm_bytes: &[u8] =
                    mcur.get(..5).ok_or(CramError::Truncated { context: "SM bytes" })?;
                substitution_matrix = SubstitutionMatrix::parse(
                    sm_bytes
                        .try_into()
                        .map_err(|_| CramError::Truncated { context: "SM bytes conversion" })?,
                );
                mcur = mcur.get(5..).ok_or(CramError::Truncated { context: "SM advance" })?;
            }
            (b'T', b'D') => {
                let td_len = varint::read_itf8_from(&mut mcur)
                    .ok_or(CramError::Truncated { context: "TD length" })?
                    as usize;
                let td_data =
                    mcur.get(..td_len).ok_or(CramError::Truncated { context: "TD data" })?;
                tag_dictionary = parse_tag_dictionary(td_data);
                mcur = mcur.get(td_len..).ok_or(CramError::Truncated { context: "TD advance" })?;
            }
            _ => {
                // Unknown key — skip 1 byte value (best effort)
                mcur = mcur.get(1..).unwrap_or(&[]);
            }
        }
    }

    *cursor = cursor
        .get(map_size..)
        .ok_or(CramError::Truncated { context: "preservation map advance" })?;

    Ok(PreservationMap {
        read_names_included,
        ap_delta,
        reference_required,
        substitution_matrix,
        tag_dictionary,
    })
}

fn parse_tag_dictionary(data: &[u8]) -> Vec<Vec<TagDictEntry>> {
    let mut result = Vec::new();
    let mut current_set = Vec::new();

    let mut i = 0;
    while i < data.len() {
        if data.get(i).copied() == Some(0) {
            // Null separator: end of this tag set
            result.push(std::mem::take(&mut current_set));
            i += 1;
            continue;
        }
        if let Some(&[a, b, t]) = data.get(i..i + 3) {
            current_set.push(TagDictEntry { tag: [a, b], bam_type: t });
            i += 3;
        } else {
            break;
        }
    }
    if !current_set.is_empty() {
        result.push(current_set);
    }
    result
}

// r[impl cram.compression.ds_encodings]
fn parse_data_series_map(cursor: &mut &[u8]) -> Result<DataSeriesEncodings, CramError> {
    let map_size = varint::read_itf8_from(cursor)
        .ok_or(CramError::Truncated { context: "data series map size" })?
        as usize;

    let map_data =
        cursor.get(..map_size).ok_or(CramError::Truncated { context: "data series map data" })?;
    let mut mcur: &[u8] = map_data;

    let num_entries = varint::read_itf8_from(&mut mcur)
        .ok_or(CramError::Truncated { context: "data series entry count" })?;

    let mut encodings: FxHashMap<[u8; 2], RawEncoding> = FxHashMap::default();

    for _ in 0..num_entries {
        let key0 = *mcur.first().ok_or(CramError::Truncated { context: "ds key" })?;
        let key1 = *mcur.get(1).ok_or(CramError::Truncated { context: "ds key" })?;
        mcur = mcur.get(2..).ok_or(CramError::Truncated { context: "ds key advance" })?;

        let key = [key0, key1];

        // Determine type from key and parse accordingly
        match &key {
            // ByteArray encodings
            b"RN" | b"IN" | b"SC" | b"BB" | b"QQ" => {
                let enc = ByteArrayEncoding::parse(&mut mcur)?;
                encodings.insert(key, RawEncoding::ByteArray(enc));
            }
            // Byte encodings
            b"FC" | b"BA" | b"QS" | b"BS" => {
                let enc = ByteEncoding::parse(&mut mcur)?;
                encodings.insert(key, RawEncoding::Byte(enc));
            }
            // Everything else is Int
            _ => {
                let enc = IntEncoding::parse(&mut mcur)?;
                encodings.insert(key, RawEncoding::Int(enc));
            }
        }
    }

    *cursor = cursor
        .get(map_size..)
        .ok_or(CramError::Truncated { context: "data series map advance" })?;

    fn take_int(map: &mut FxHashMap<[u8; 2], RawEncoding>, key: &[u8; 2]) -> IntEncoding {
        match map.remove(key) {
            Some(RawEncoding::Int(e)) => e,
            _ => IntEncoding::Null,
        }
    }
    fn take_byte(map: &mut FxHashMap<[u8; 2], RawEncoding>, key: &[u8; 2]) -> ByteEncoding {
        match map.remove(key) {
            Some(RawEncoding::Byte(e)) => e,
            _ => ByteEncoding::Null,
        }
    }
    fn take_byte_array(
        map: &mut FxHashMap<[u8; 2], RawEncoding>,
        key: &[u8; 2],
    ) -> ByteArrayEncoding {
        match map.remove(key) {
            Some(RawEncoding::ByteArray(e)) => e,
            _ => ByteArrayEncoding::Null,
        }
    }

    Ok(DataSeriesEncodings {
        bam_flags: take_int(&mut encodings, b"BF"),
        cram_flags: take_int(&mut encodings, b"CF"),
        ref_id: take_int(&mut encodings, b"RI"),
        read_length: take_int(&mut encodings, b"RL"),
        alignment_pos: take_int(&mut encodings, b"AP"),
        read_group: take_int(&mut encodings, b"RG"),
        read_name: take_byte_array(&mut encodings, b"RN"),
        mate_flags: take_int(&mut encodings, b"MF"),
        next_segment_ref: take_int(&mut encodings, b"NS"),
        next_mate_pos: take_int(&mut encodings, b"NP"),
        template_size: take_int(&mut encodings, b"TS"),
        next_fragment: take_int(&mut encodings, b"NF"),
        tag_line: take_int(&mut encodings, b"TL"),
        feature_count: take_int(&mut encodings, b"FN"),
        feature_code: take_byte(&mut encodings, b"FC"),
        feature_pos: take_int(&mut encodings, b"FP"),
        deletion_length: take_int(&mut encodings, b"DL"),
        mapping_quality: take_int(&mut encodings, b"MQ"),
        base: take_byte(&mut encodings, b"BA"),
        quality_score: take_byte(&mut encodings, b"QS"),
        base_sub: take_byte(&mut encodings, b"BS"),
        insertion: take_byte_array(&mut encodings, b"IN"),
        soft_clip: take_byte_array(&mut encodings, b"SC"),
        ref_skip: take_int(&mut encodings, b"RS"),
        padding: take_int(&mut encodings, b"PD"),
        hard_clip: take_int(&mut encodings, b"HC"),
        bases_block: take_byte_array(&mut encodings, b"BB"),
        quality_block: take_byte_array(&mut encodings, b"QQ"),
    })
}

#[derive(Debug)]
enum RawEncoding {
    Int(IntEncoding),
    Byte(ByteEncoding),
    ByteArray(ByteArrayEncoding),
}

// r[impl cram.compression.tag_encodings]
fn parse_tag_encoding_map(
    cursor: &mut &[u8],
) -> Result<FxHashMap<i32, ByteArrayEncoding>, CramError> {
    let map_size = varint::read_itf8_from(cursor)
        .ok_or(CramError::Truncated { context: "tag encoding map size" })?
        as usize;

    let map_data =
        cursor.get(..map_size).ok_or(CramError::Truncated { context: "tag encoding map data" })?;
    let mut mcur: &[u8] = map_data;

    let num_entries = varint::read_itf8_from(&mut mcur)
        .ok_or(CramError::Truncated { context: "tag encoding entry count" })?;

    let mut map = FxHashMap::default();
    for _ in 0..num_entries {
        let key = varint::read_itf8_from(&mut mcur)
            .ok_or(CramError::Truncated { context: "tag encoding key" })? as i32;
        let enc = ByteArrayEncoding::parse(&mut mcur)?;
        map.insert(key, enc);
    }

    *cursor = cursor
        .get(map_size..)
        .ok_or(CramError::Truncated { context: "tag encoding map advance" })?;

    Ok(map)
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "test-only arithmetic on known-valid CRAM file data"
)]
mod tests {
    use super::*;

    fn load_first_data_container_compression_header() -> CompressionHeader {
        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        use super::super::{block, container::ContainerHeader};

        // Skip file def + header container
        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];
        let header_container = ContainerHeader::parse(after_file_def).unwrap();
        let data_start = header_container.header_size + header_container.length as usize;

        // Parse first data container
        #[allow(clippy::indexing_slicing)]
        let data_container_bytes = &after_file_def[data_start..];
        let data_container = ContainerHeader::parse(data_container_bytes).unwrap();

        // First block in data container is the compression header
        #[allow(clippy::indexing_slicing)]
        let block_data = &data_container_bytes[data_container.header_size..];
        let (comp_block, _) = block::parse_block(block_data).unwrap();
        assert_eq!(comp_block.content_type, block::ContentType::CompressionHeader);

        CompressionHeader::parse(&comp_block.data).unwrap()
    }

    #[test]
    fn parse_compression_header_from_real_cram() {
        let ch = load_first_data_container_compression_header();

        // Verify preservation map has sensible defaults
        assert!(ch.preservation.ap_delta, "AP should be delta-coded by default");
        assert!(ch.preservation.reference_required, "RR should be true by default");

        // Substitution matrix should have non-N entries
        let sm = &ch.preservation.substitution_matrix;
        // For ref=A (index 0), first substitution should be C, G, or T (not A)
        assert_ne!(sm.matrix[0][0], b'A', "substitution for A should not be A");
    }

    #[test]
    fn data_series_encodings_present() {
        let ch = load_first_data_container_compression_header();
        let ds = &ch.data_series;

        // Essential data series should not be Null
        assert!(!matches!(ds.bam_flags, IntEncoding::Null), "BF should be present");
        assert!(!matches!(ds.cram_flags, IntEncoding::Null), "CF should be present");
        assert!(!matches!(ds.read_length, IntEncoding::Null), "RL should be present");
        assert!(!matches!(ds.alignment_pos, IntEncoding::Null), "AP should be present");
        assert!(!matches!(ds.feature_count, IntEncoding::Null), "FN should be present");
        assert!(!matches!(ds.mapping_quality, IntEncoding::Null), "MQ should be present");
    }

    #[test]
    fn tag_dictionary_parsed() {
        let ch = load_first_data_container_compression_header();
        // There should be at least one tag set in the dictionary
        assert!(!ch.preservation.tag_dictionary.is_empty(), "tag dictionary should have entries");
    }

    #[test]
    fn tag_encodings_present() {
        let ch = load_first_data_container_compression_header();
        // There should be tag encodings for the tags used in the file
        assert!(!ch.tag_encodings.is_empty(), "should have tag encodings");
    }

    #[test]
    fn substitution_matrix_lookup() {
        // Default samtools matrix: for ref A, subs are C(0) G(1) T(2) N(3)
        let sm = SubstitutionMatrix {
            matrix: [
                [b'C', b'G', b'T', b'N'], // ref A
                [b'A', b'G', b'T', b'N'], // ref C
                [b'A', b'C', b'T', b'N'], // ref G
                [b'A', b'C', b'G', b'N'], // ref T
                [b'A', b'C', b'G', b'T'], // ref N
            ],
        };
        assert_eq!(sm.substitute(b'A', 0), b'C');
        assert_eq!(sm.substitute(b'A', 1), b'G');
        assert_eq!(sm.substitute(b'C', 0), b'A');
        assert_eq!(sm.substitute(b'N', 3), b'T');
        assert_eq!(sm.substitute(b'a', 0), b'C'); // case insensitive
    }

    #[test]
    fn substitution_matrix_matches_noodles_test_vector() {
        // Test vector from noodles: § 10.6.4 "Mapped reads: Substitution Matrix Format"
        let bytes = [0x93, 0x1b, 0x6c, 0xb1, 0xc6];
        let sm = SubstitutionMatrix::parse(&bytes);

        // Expected from noodles test:
        // ref A: [T, G, C, N]
        // ref C: [A, G, T, N]
        // ref G: [N, A, C, T]
        // ref T: [G, N, A, C]
        // ref N: [C, G, T, A]
        assert_eq!(sm.substitute(b'A', 0), b'T');
        assert_eq!(sm.substitute(b'A', 1), b'G');
        assert_eq!(sm.substitute(b'A', 2), b'C');
        assert_eq!(sm.substitute(b'A', 3), b'N');
        assert_eq!(sm.substitute(b'C', 0), b'A');
        assert_eq!(sm.substitute(b'C', 1), b'G');
        assert_eq!(sm.substitute(b'C', 2), b'T');
        assert_eq!(sm.substitute(b'C', 3), b'N');
        assert_eq!(sm.substitute(b'G', 0), b'N');
        assert_eq!(sm.substitute(b'G', 1), b'A');
        assert_eq!(sm.substitute(b'G', 2), b'C');
        assert_eq!(sm.substitute(b'G', 3), b'T');
        assert_eq!(sm.substitute(b'T', 0), b'G');
        assert_eq!(sm.substitute(b'T', 1), b'N');
        assert_eq!(sm.substitute(b'T', 2), b'A');
        assert_eq!(sm.substitute(b'T', 3), b'C');
        assert_eq!(sm.substitute(b'N', 0), b'C');
        assert_eq!(sm.substitute(b'N', 1), b'G');
        assert_eq!(sm.substitute(b'N', 2), b'T');
        assert_eq!(sm.substitute(b'N', 3), b'A');
    }

    proptest::proptest! {
        #[test]
        fn substitution_matrix_code_in_bounds(ref_idx in 0u8..5, code in 0u8..4) {
            let sm = SubstitutionMatrix {
                matrix: [
                    [b'C', b'G', b'T', b'N'],
                    [b'A', b'G', b'T', b'N'],
                    [b'A', b'C', b'T', b'N'],
                    [b'A', b'C', b'G', b'N'],
                    [b'A', b'C', b'G', b'T'],
                ],
            };
            let ref_base = match ref_idx {
                0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => b'N',
            };
            let result = sm.substitute(ref_base, code);
            proptest::prop_assert!(
                matches!(result, b'A' | b'C' | b'G' | b'T' | b'N'),
                "substitute(ref={}, code={}) = 0x{:02x} '{}', must be in {{A,C,G,T,N}}",
                ref_base as char, code, result, result as char
            );
            // A substitution must differ from the reference base.
            // For ref A/C/G/T, READ_BASES excludes the ref itself.
            // For ref N, READ_BASES is [A,C,G,T], so N is never returned.
            proptest::prop_assert_ne!(
                result, ref_base,
                "substitute(ref={}, code={}) returned ref base itself",
                ref_base as char, code
            );
        }
    }
}
