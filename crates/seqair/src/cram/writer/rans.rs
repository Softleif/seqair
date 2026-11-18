//! rANS 4×8 encoder, order 0 and order 1: the inverse of
//! [`crate::cram::rans`], laid out as htscodecs' `rANS_static.c` writes it.

use crate::cram::varint::write_itf8;

// r[impl cram_writer.rans4x8]

#[derive(Debug, Clone, Copy)]
pub(super) enum Order {
    Zero,
    One,
}

/// htscodecs normalizes to one less than 4096 (`fsum++`); its decoder
/// accepts 4095 or 4096.
const TOTAL: u32 = 4095;
const LOWER_BOUND: u32 = 1 << 23;

/// A symbol's slot range in a frequency table.
#[derive(Debug, Clone, Copy, Default)]
struct Symbol {
    start: u32,
    freq: u32,
}

type Table = [Symbol; 256];

/// `data` compressed, or `None` when it is empty, too large for the
/// format's 32-bit sizes, or too small for order 1.
pub(super) fn encode(data: &[u8], order: Order) -> Option<Vec<u8>> {
    // An empty table would read as symbol 0 followed by the states.
    if data.is_empty() {
        return None;
    }
    let raw_len = u32::try_from(data.len()).ok()?;
    let mut table = Vec::new();
    let mut rev = Vec::with_capacity(data.len());
    let states = match order {
        Order::Zero => encode_order_0(data, &mut table, &mut rev),
        Order::One => encode_order_1(data, &mut table, &mut rev)?,
    };
    let payload_len = table.len().saturating_add(16).saturating_add(rev.len());
    let mut out = Vec::with_capacity(payload_len.saturating_add(9));
    out.push(match order {
        Order::Zero => 0,
        Order::One => 1,
    });
    out.extend_from_slice(&u32::try_from(payload_len).ok()?.to_le_bytes());
    out.extend_from_slice(&raw_len.to_le_bytes());
    out.extend_from_slice(&table);
    for state in states {
        out.extend_from_slice(&state.to_le_bytes());
    }
    out.extend(rev.iter().rev());
    Some(out)
}

fn encode_order_0(data: &[u8], table: &mut Vec<u8>, rev: &mut Vec<u8>) -> [u32; 4] {
    let mut counts = [0u32; 256];
    for &b in data {
        let c = counts.get_mut(usize::from(b));
        if let Some(c) = c {
            *c = c.saturating_add(1);
        }
    }
    let symbols = normalize(&counts);
    write_frequencies(table, &symbols);

    let mut states = [LOWER_BOUND; 4];
    // The decoder reads byte `i` with state `i % 4`, renormalizing after
    // each; the encoder runs it backwards.
    for (i, &b) in data.iter().enumerate().rev() {
        if let (Some(state), Some(&sym)) = (states.get_mut(i % 4), symbols.get(usize::from(b))) {
            put(state, rev, sym);
        }
    }
    states
}

fn encode_order_1(data: &[u8], table: &mut Vec<u8>, rev: &mut Vec<u8>) -> Option<[u32; 4]> {
    // htscodecs reads the byte before the remainder from the last quarter.
    if data.len() < 4 {
        return None;
    }
    let quarter = data.len() / 4;
    // Each byte's context: the byte before it in its quarter, 0 at a
    // quarter's start; the remainder continues the last quarter.
    let context = |i: usize| -> u8 {
        if i < quarter.saturating_mul(4) && i.is_multiple_of(quarter) {
            0
        } else {
            data.get(i.wrapping_sub(1)).copied().unwrap_or(0)
        }
    };

    let mut counts = vec![[0u32; 256]; 256];
    for (i, &b) in data.iter().enumerate() {
        if let Some(c) =
            counts.get_mut(usize::from(context(i))).and_then(|row| row.get_mut(usize::from(b)))
        {
            *c = c.saturating_add(1);
        }
    }
    let mut tables = vec![[Symbol::default(); 256]; 256];
    let mut active = [false; 256];
    for ((row, symbols), on) in counts.iter().zip(tables.iter_mut()).zip(active.iter_mut()) {
        if row.iter().any(|&c| c > 0) {
            *symbols = normalize(row);
            *on = true;
        }
    }
    write_runs(table, &active, |table, ctx| {
        if let Some(symbols) = tables.get(ctx) {
            write_frequencies(table, symbols);
        }
    });

    let mut states = [LOWER_BOUND; 4];
    let mut put_at = |state: usize, i: usize| {
        let b = data.get(i).copied().unwrap_or(0);
        let sym = tables
            .get(usize::from(context(i)))
            .and_then(|row| row.get(usize::from(b)))
            .copied()
            .unwrap_or_default();
        if let Some(state) = states.get_mut(state) {
            put(state, rev, sym);
        }
    };
    // Decode order: position by position across the four quarters, then
    // the remainder with state 3. Encode it backwards.
    for i in (quarter.saturating_mul(4)..data.len()).rev() {
        put_at(3, i);
    }
    for pos in (0..quarter).rev() {
        for q in (0..4).rev() {
            put_at(q, q.saturating_mul(quarter).saturating_add(pos));
        }
    }
    Some(states)
}

/// `RansEncPutSymbol`: renormalize (bytes written backwards into `rev`),
/// then encode.
#[allow(clippy::arithmetic_side_effects, reason = "freq ≥ 1, so the division cannot fail")]
fn put(x: &mut u32, rev: &mut Vec<u8>, sym: Symbol) {
    debug_assert!(sym.freq > 0, "encoding a symbol without frequency");
    let freq = sym.freq.max(1);
    let x_max = ((LOWER_BOUND >> 12) << 8).saturating_mul(freq);
    while *x >= x_max {
        #[allow(clippy::cast_possible_truncation, reason = "the low byte")]
        rev.push(*x as u8);
        *x >>= 8;
    }
    *x = ((*x / freq) << 12).wrapping_add(*x % freq).wrapping_add(sym.start);
}

/// Scale `counts` to frequencies summing to [`TOTAL`], each present symbol
/// at least 1, and give each its slot range.
fn normalize(counts: &[u32; 256]) -> Table {
    let total: u64 = counts.iter().map(|&c| u64::from(c)).sum();
    let mut freq = [0u32; 256];
    for (f, &c) in freq.iter_mut().zip(counts) {
        if c > 0 {
            let scaled =
                u64::from(c).saturating_mul(u64::from(TOTAL)).checked_div(total).unwrap_or(0);
            *f = u32::try_from(scaled).unwrap_or(TOTAL).max(1);
        }
    }
    let mut sum: u32 = freq.iter().sum();
    // At most 256 symbols at 1 each, so the target is always reachable by
    // taking from the largest.
    while sum != TOTAL {
        let Some(max) = freq.iter_mut().max_by_key(|f| **f) else { break };
        if sum < TOTAL {
            *max = max.saturating_add(TOTAL.saturating_sub(sum));
            sum = TOTAL;
        } else {
            let take = sum.saturating_sub(TOTAL).min(max.saturating_sub(1));
            *max = max.saturating_sub(take);
            sum = sum.saturating_sub(take);
            if take == 0 {
                break;
            }
        }
    }
    let mut table = [Symbol::default(); 256];
    let mut start = 0u32;
    for (sym, &f) in table.iter_mut().zip(&freq) {
        *sym = Symbol { start, freq: f };
        start = start.saturating_add(f);
    }
    table
}

/// An order-0 frequency table: the present symbols, run-length coded as
/// htscodecs writes them, each with its frequency as ITF8; then 0.
fn write_frequencies(out: &mut Vec<u8>, symbols: &Table) {
    let mut present = [false; 256];
    for (p, s) in present.iter_mut().zip(symbols) {
        *p = s.freq > 0;
    }
    write_runs(out, &present, |out, sym| {
        let freq = symbols.get(sym).map_or(0, |s| s.freq);
        write_itf8(out, i32::try_from(freq).unwrap_or(i32::MAX));
    });
}

/// The alphabet coding both table levels share: a present symbol is written
/// unless a run covers it; one following a present symbol starts a run,
/// the count of present symbols right after it. `value` writes each one's
/// payload. The list ends with 0.
#[allow(clippy::cast_possible_truncation, reason = "symbols and run lengths are below 256")]
fn write_runs(
    out: &mut Vec<u8>,
    present: &[bool; 256],
    mut value: impl FnMut(&mut Vec<u8>, usize),
) {
    let mut run = 0usize;
    for j in 0..256usize {
        if !present.get(j).copied().unwrap_or(false) {
            continue;
        }
        if run > 0 {
            run = run.saturating_sub(1);
        } else {
            out.push(j as u8);
            if j > 0 && present.get(j.wrapping_sub(1)).copied().unwrap_or(false) {
                run = present.iter().skip(j.saturating_add(1)).take_while(|&&p| p).count();
                out.push(run as u8);
            }
        }
        value(out, j);
    }
    out.push(0);
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;
    use crate::cram::rans::decode;
    use hegel::prelude::*;

    // r[verify cram_writer.rans4x8]
    #[hegel::test]
    fn every_input_round_trips_through_the_decoder(tc: TestCase) {
        // Small alphabets and long runs as well as uniform bytes: the
        // frequency tables' runs and the renormalization both get exercised.
        let alphabet: u8 = tc.draw(gs::integers::<u8>().min_value(1));
        let len: usize = tc.draw(gs::integers::<usize>().max_value(5000));
        let seed: u64 = tc.draw(gs::integers::<u64>());
        let mut x = seed | 1;
        let data: Vec<u8> = (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                #[allow(clippy::cast_possible_truncation, reason = "below alphabet ≤ 255")]
                let b = (x % u64::from(alphabet)) as u8;
                b
            })
            .collect();
        for order in [Order::Zero, Order::One] {
            if let Some(encoded) = encode(&data, order) {
                assert_eq!(decode(&encoded).unwrap(), data, "{order:?}");
            }
        }
    }

    #[test]
    fn edge_inputs_round_trip() {
        let inputs: [&[u8]; 6] =
            [b"", b"a", b"abcd", b"abcde", &[0; 1000], &[255, 0, 255, 0, 1, 254, 254, 254, 254]];
        for data in inputs {
            for order in [Order::Zero, Order::One] {
                if let Some(encoded) = encode(data, order) {
                    assert_eq!(decode(&encoded).unwrap(), data, "{order:?} {data:?}");
                }
            }
        }
        let every_byte: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        for order in [Order::Zero, Order::One] {
            let encoded = encode(&every_byte, order).unwrap();
            assert_eq!(decode(&encoded).unwrap(), every_byte);
        }
    }

    #[test]
    fn skewed_input_compresses() {
        let data: Vec<u8> = (0..10_000u32).map(|i| if i % 50 == 0 { b'B' } else { b'A' }).collect();
        let encoded = encode(&data, Order::Zero).unwrap();
        assert!(encoded.len() < data.len() / 4, "{} bytes", encoded.len());
    }
}
