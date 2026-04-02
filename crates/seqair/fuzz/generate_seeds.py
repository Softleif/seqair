#!/usr/bin/env python3
"""Generate fuzz seed corpora from test data files.

Run from the repo root:
    python3 crates/seqair/fuzz/generate_seeds.py
"""

import struct
import os
import sys
import zlib

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(REPO_ROOT, "tests", "data")
SEED_DIR = os.path.join(REPO_ROOT, "crates", "seqair", "fuzz", "seeds")


def write_seed(target: str, name: str, data: bytes) -> None:
    d = os.path.join(SEED_DIR, target)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name)
    with open(path, "wb") as f:
        f.write(data)
    print(f"  {target}/{name}: {len(data)} bytes")


def decompress_bgzf_blocks(bam: bytes, max_raw: int = 100_000):
    """Yield (block_bytes, decompressed_bytes) from a BGZF stream."""
    pos = 0
    total_raw = 0
    while pos < len(bam) and total_raw < max_raw:
        if bam[pos : pos + 2] != b"\x1f\x8b":
            break
        xlen = struct.unpack_from("<H", bam, pos + 10)[0]
        extra = bam[pos + 12 : pos + 12 + xlen]
        ep = 0
        bsize = None
        while ep + 4 <= len(extra):
            si1, si2 = extra[ep], extra[ep + 1]
            slen = struct.unpack_from("<H", extra, ep + 2)[0]
            if si1 == ord("B") and si2 == ord("C") and slen == 2:
                bsize = struct.unpack_from("<H", extra, ep + 4)[0]
                break
            ep += 4 + slen
        if bsize is None:
            break
        total = bsize + 1
        block_bytes = bam[pos : pos + total]
        deflate_start = 12 + xlen
        footer_start = total - 8
        try:
            decompressed = zlib.decompress(block_bytes[deflate_start:footer_start], -15)
        except Exception:
            decompressed = b""
        yield block_bytes, decompressed
        total_raw += len(decompressed)
        pos += total


def extract_bam_records(raw_bam: bytearray, max_records: int = 5):
    """Extract individual BAM records from decompressed BAM data."""
    if raw_bam[:4] != b"BAM\x01":
        return
    header_len = struct.unpack_from("<i", raw_bam, 4)[0]
    p = 8 + header_len
    n_ref = struct.unpack_from("<i", raw_bam, p)[0]
    p += 4
    for _ in range(n_ref):
        name_len = struct.unpack_from("<i", raw_bam, p)[0]
        p += 4 + name_len + 4
    for i in range(max_records):
        if p + 4 >= len(raw_bam):
            break
        block_size = struct.unpack_from("<i", raw_bam, p)[0]
        if block_size <= 0 or p + 4 + block_size > len(raw_bam):
            break
        yield i, bytes(raw_bam[p + 4 : p + 4 + block_size])
        p += 4 + block_size


def extract_sam_header(raw_bam: bytearray):
    """Extract SAM header text from decompressed BAM data."""
    if raw_bam[:4] != b"BAM\x01":
        return None
    header_len = struct.unpack_from("<i", raw_bam, 4)[0]
    return bytes(raw_bam[8 : 8 + header_len])


def main():
    if not os.path.isdir(DATA_DIR):
        print(f"Error: test data directory not found: {DATA_DIR}", file=sys.stderr)
        sys.exit(1)

    print("Generating fuzz seeds from test data...\n")

    # --- BAM-derived seeds ---
    bam_path = os.path.join(DATA_DIR, "test.bam")
    bam = open(bam_path, "rb").read()

    raw_bam = bytearray()
    bgzf_blocks = []
    for block_bytes, decompressed in decompress_bgzf_blocks(bam):
        bgzf_blocks.append(block_bytes)
        raw_bam.extend(decompressed)

    # fuzz_bam_record: individual records
    for i, record in extract_bam_records(raw_bam, max_records=5):
        write_seed("fuzz_bam_record", f"record_{i}", record)

    # fuzz_bgzf: individual BGZF blocks (first 3)
    for i, block in enumerate(bgzf_blocks[:3]):
        write_seed("fuzz_bgzf", f"block_{i}", block)

    # fuzz_sam_header: header text + individual lines
    header_text = extract_sam_header(raw_bam)
    if header_text:
        write_seed("fuzz_sam_header", "full_header", header_text)
        for i, line in enumerate(header_text.split(b"\n")[:5]):
            if line:
                write_seed("fuzz_sam_header", f"line_{i}", line + b"\n")

    # --- BAI seed (raw file is the exact format) ---
    bai = open(os.path.join(DATA_DIR, "test.bam.bai"), "rb").read()
    write_seed("fuzz_bam_index", "test_bai", bai)

    # --- FAI seeds ---
    for name in ["test.fasta.fai", "test.fasta.gz.fai"]:
        path = os.path.join(DATA_DIR, name)
        if os.path.exists(path):
            write_seed("fuzz_fasta_index", name.replace(".", "_"), open(path, "rb").read())

    # --- GZI seed (trimmed to keep small) ---
    gzi = open(os.path.join(DATA_DIR, "test.fasta.gz.gzi"), "rb").read()
    count = struct.unpack_from("<Q", gzi, 0)[0]
    trimmed_count = min(count, 16)
    write_seed("fuzz_gzi", "trimmed", struct.pack("<Q", trimmed_count) + gzi[8 : 8 + trimmed_count * 16])
    # Also the full file (small enough)
    write_seed("fuzz_gzi", "full", gzi)

    # --- CRAM-derived seeds ---
    for cram_name in ["test.cram", "test_v30.cram", "test_gzip.cram"]:
        cram_path = os.path.join(DATA_DIR, cram_name)
        if not os.path.exists(cram_path):
            continue
        cram = open(cram_path, "rb").read()
        # Skip 26-byte file definition
        after_filedef = cram[26:]
        tag = cram_name.replace(".", "_")
        write_seed("fuzz_cram_container", f"{tag}_container", after_filedef[:128])
        write_seed("fuzz_cram_block", f"{tag}_blocks", after_filedef[:1024])
        write_seed("fuzz_cram_header", f"{tag}_header", after_filedef[:1024])

    # --- CRAM varint seeds ---
    for val, encoded in [
        (0, b"\x00"),
        (127, b"\x7f"),
        (128, b"\x80\x80"),
        (16383, b"\xbf\xff"),
        (16384, b"\xc0\x40\x00"),
        (2097151, b"\xdf\xff\xff"),
    ]:
        write_seed("fuzz_cram_varint", f"itf8_{val}", encoded)

    # --- Region string seeds ---
    for s in ["chr1:1-1000", "chrX", "chr1:100", "MT:1-16569", "chr1:1,000-2,000"]:
        write_seed("fuzz_region_string", s.replace(":", "_").replace(",", ""), s.encode())

    # --- rANS codec seed (noodles test vector) ---
    rans = bytes(
        [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49,
            0x65, 0x00, 0x82, 0x49, 0x6C, 0x82, 0x49, 0x6E, 0x82, 0x49, 0x6F, 0x00,
            0x84, 0x92, 0x73, 0x82, 0x49, 0x00, 0xE2, 0x06, 0x83, 0x18, 0x74, 0x7B,
            0x41, 0x0C, 0x2B, 0xA9, 0x41, 0x0C, 0x25, 0x31, 0x80, 0x03,
        ]
    )
    # Wrap as Rans4x8 Arbitrary enum variant (discriminant 0)
    write_seed("fuzz_cram_codecs", "rans4x8_order0", bytes([0]) + rans)

    total = sum(1 for _ in _walk_seeds())
    print(f"\nDone: {total} seed files in {SEED_DIR}")


def _walk_seeds():
    for root, _, files in os.walk(SEED_DIR):
        yield from files


if __name__ == "__main__":
    main()
