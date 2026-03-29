//! Tests for RecordStore: slab-based BAM record storage.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use seqair::bam::record_store::RecordStore;
use seqair_types::Base;

// r[verify record_store.push_raw+2]
// r[verify record_store.field_access]
#[test]
fn decode_record_into_slabs() {
    // Minimal BAM record: tid=0, pos=100, mapq=60, 1 CIGAR op (4M), seq=ACGT, qual=[30,30,30,30]
    let raw =
        make_test_record(0, 100, 0x63, 60, b"read1", &[(4 << 4)], &[0x12, 0x48], &[30; 4], &[]);

    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).expect("decode");

    assert_eq!(store.len(), 1);
    let rec = store.record(idx);
    assert_eq!(rec.pos, 100);
    assert_eq!(rec.mapq, 60);
    assert_eq!(rec.flags, 0x63);
    assert_eq!(rec.seq_len, 4);

    assert_eq!(store.qname(idx), b"read1");
    assert_eq!(store.seq_at(idx, 0), Base::A);
    assert_eq!(store.seq_at(idx, 1), Base::C);
    assert_eq!(store.seq_at(idx, 2), Base::G);
    assert_eq!(store.seq_at(idx, 3), Base::T);
    assert_eq!(store.qual(idx).len(), 4);
    assert!(store.qual(idx).iter().all(|&q| q == 30));
}

// r[verify record_store.push_raw+2]
#[test]
fn multiple_records_share_slabs() {
    let raw1 =
        make_test_record(0, 100, 0x63, 60, b"read1", &[(4 << 4)], &[0x12, 0x48], &[30; 4], &[]);
    let raw2 =
        make_test_record(0, 200, 0x63, 50, b"read2", &[(4 << 4)], &[0x12, 0x48], &[25; 4], &[]);

    let mut store = RecordStore::new();
    let idx1 = store.push_raw(&raw1).unwrap();
    let idx2 = store.push_raw(&raw2).unwrap();

    assert_eq!(store.len(), 2);
    assert_eq!(store.record(idx1).pos, 100);
    assert_eq!(store.record(idx2).pos, 200);
    assert_eq!(store.qname(idx1), b"read1");
    assert_eq!(store.qname(idx2), b"read2");
    assert_eq!(store.qual(idx1)[0], 30);
    assert_eq!(store.qual(idx2)[0], 25);
}

// r[verify record_store.clear+2]
#[test]
fn clear_retains_capacity() {
    let raw =
        make_test_record(0, 100, 0x63, 60, b"read1", &[(4 << 4)], &[0x12, 0x48], &[30; 4], &[]);

    let mut store = RecordStore::new();
    for _ in 0..100 {
        store.push_raw(&raw).unwrap();
    }

    let records_cap = store.records_capacity();
    let names_cap = store.names_capacity();
    let data_cap = store.data_capacity();

    store.clear();
    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() >= records_cap);
    assert!(store.names_capacity() >= names_cap);
    assert!(store.data_capacity() >= data_cap);
}

// r[verify record_store.capacity]
#[test]
fn with_capacity_hint_preallocates() {
    let store = RecordStore::with_byte_hint(100_000);
    assert!(store.records_capacity() > 0);
    assert!(store.names_capacity() > 0);
    assert!(store.data_capacity() > 0);
}

// r[verify record_store.field_access]
#[test]
fn aux_tag_accessible_from_store() {
    // Record with an RG:Z:group1 aux tag
    let aux = b"RGZgroup1\0";
    let raw =
        make_test_record(0, 100, 0x63, 60, b"read1", &[(4 << 4)], &[0x12, 0x48], &[30; 4], aux);

    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).unwrap();

    let aux_bytes = store.aux(idx);
    assert!(!aux_bytes.is_empty());
    assert!(aux_bytes.starts_with(b"RG"));
}

// r[verify record_store.no_rc]
// r[verify record_store.region_scoped]
#[test]
fn integration_with_real_bam() {
    use seqair::bam::reader::IndexedBamReader;
    use std::path::Path;

    let bam_path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = IndexedBamReader::open(bam_path).expect("open BAM");
    let tid = reader.header().tid("chr19").expect("tid");

    let mut store = RecordStore::new();
    let count = reader.fetch_into(tid, 6_105_700, 6_105_800, &mut store).expect("fetch");

    assert!(count > 0);
    for i in 0..store.len() as u32 {
        let rec = store.record(i);
        assert!(rec.pos >= 0);
        assert!(rec.seq_len > 0);
        assert_eq!(store.qual(i).len(), rec.seq_len as usize);
    }
}

// r[verify record_store.push_fields]
// r[verify unified.push_fields_equivalence]
// r[verify unified.record_store_push]
#[test]
fn push_fields_matches_push_raw() -> Result<(), Box<dyn std::error::Error>> {
    let cigar_ops: &[u32] = &[(4 << 4)]; // 4M
    let packed_seq = &[0x12, 0x48]; // ACGT in 4-bit packed
    let qual = &[30u8; 4];
    let aux = b"RGZgroup1\0";

    let raw = make_test_record(0, 100, 0x63, 60, b"read1", cigar_ops, packed_seq, qual, aux);

    let mut store_raw = RecordStore::new();
    let idx_raw = store_raw.push_raw(&raw)?;

    // Build the same record from pre-parsed fields
    let cigar_packed: Vec<u8> = cigar_ops.iter().flat_map(|op| op.to_le_bytes()).collect();
    let bases = [Base::A, Base::C, Base::G, Base::T];

    let mut store_fields = RecordStore::new();
    let idx_fields = store_fields.push_fields(
        100, // pos
        103, // end_pos (pos + 4M - 1)
        0x63,
        60,
        store_raw.record(idx_raw).matching_bases,
        store_raw.record(idx_raw).indel_bases,
        b"read1",
        &cigar_packed,
        &bases,
        qual,
        aux,
    )?;

    // Compare fixed fields
    let r = store_raw.record(idx_raw);
    let f = store_fields.record(idx_fields);
    assert_eq!(r.pos, f.pos);
    assert_eq!(r.end_pos, f.end_pos);
    assert_eq!(r.flags, f.flags);
    assert_eq!(r.mapq, f.mapq);
    assert_eq!(r.seq_len, f.seq_len);
    assert_eq!(r.n_cigar_ops, f.n_cigar_ops);
    assert_eq!(r.matching_bases, f.matching_bases);
    assert_eq!(r.indel_bases, f.indel_bases);

    // Compare slab contents
    assert_eq!(store_raw.qname(idx_raw), store_fields.qname(idx_fields));
    assert_eq!(store_raw.seq(idx_raw), store_fields.seq(idx_fields));
    assert_eq!(store_raw.qual(idx_raw), store_fields.qual(idx_fields));
    assert_eq!(store_raw.cigar(idx_raw), store_fields.cigar(idx_fields));
    assert_eq!(store_raw.aux(idx_raw), store_fields.aux(idx_fields));
    Ok(())
}

// r[verify unified.push_fields_equivalence]
#[test]
fn push_fields_with_real_bam_records() -> Result<(), Box<dyn std::error::Error>> {
    use seqair::bam::reader::IndexedBamReader;
    use std::path::Path;

    let bam_path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = IndexedBamReader::open(bam_path)?;
    let tid = reader.header().tid("chr19").ok_or("chr19 not found")?;

    let mut store = RecordStore::new();
    reader.fetch_into(tid, 6_105_700, 6_105_800, &mut store)?;
    assert!(!store.is_empty());

    // Re-push every record via push_fields and compare
    let mut store2 = RecordStore::new();
    for i in 0..store.len() as u32 {
        let rec = store.record(i);
        store2.push_fields(
            rec.pos,
            rec.end_pos,
            rec.flags,
            rec.mapq,
            rec.matching_bases,
            rec.indel_bases,
            store.qname(i),
            store.cigar(i),
            store.seq(i),
            store.qual(i),
            store.aux(i),
        )?;
    }

    assert_eq!(store.len(), store2.len());
    for i in 0..store.len() as u32 {
        assert_eq!(store.record(i).pos, store2.record(i).pos, "rec {i}: pos");
        assert_eq!(store.record(i).flags, store2.record(i).flags, "rec {i}: flags");
        assert_eq!(store.record(i).seq_len, store2.record(i).seq_len, "rec {i}: seq_len");
        assert_eq!(store.qname(i), store2.qname(i), "rec {i}: qname");
        assert_eq!(store.seq(i), store2.seq(i), "rec {i}: seq");
        assert_eq!(store.qual(i), store2.qual(i), "rec {i}: qual");
        assert_eq!(store.cigar(i), store2.cigar(i), "rec {i}: cigar");
        assert_eq!(store.aux(i), store2.aux(i), "rec {i}: aux");
    }
    Ok(())
}

/// Build a minimal raw BAM record for testing.
#[allow(clippy::too_many_arguments)]
fn make_test_record(
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    qname: &[u8],
    cigar_ops: &[u32],
    packed_seq: &[u8],
    qual: &[u8],
    aux: &[u8],
) -> Vec<u8> {
    let name_len = qname.len() + 1; // +1 for NUL
    let n_cigar_ops = cigar_ops.len();
    let cigar_bytes = n_cigar_ops * 4;
    let seq_len = qual.len() as u32;
    let seq_bytes = packed_seq.len();

    let total = 32 + name_len + cigar_bytes + seq_bytes + qual.len() + aux.len();
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[12..14].copy_from_slice(&(n_cigar_ops as u16).to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());

    let mut off = 32;
    raw[off..off + qname.len()].copy_from_slice(qname);
    raw[off + qname.len()] = 0; // NUL
    off += name_len;

    for &op in cigar_ops {
        raw[off..off + 4].copy_from_slice(&op.to_le_bytes());
        off += 4;
    }

    raw[off..off + seq_bytes].copy_from_slice(packed_seq);
    off += seq_bytes;

    raw[off..off + qual.len()].copy_from_slice(qual);
    off += qual.len();

    raw[off..off + aux.len()].copy_from_slice(aux);

    raw
}
