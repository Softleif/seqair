//! The rANS 4×8 encoder's reference counting and coding loops: a byte at a
//! time, each byte's context worked out on its own, dividing by each
//! symbol's reciprocal. The fast path in the parent module must write
//! these bytes for every input.

use super::{LOWER_BOUND, Symbol, Table, normalize, write_frequencies, write_runs};

pub(super) fn encode_order_0(data: &[u8], table: &mut Vec<u8>, rev: &mut Vec<u8>) -> [u32; 4] {
    let mut counts = [0u32; 256];
    for &b in data {
        let c = counts.get_mut(usize::from(b));
        if let Some(c) = c {
            *c = c.saturating_add(1);
        }
    }
    let symbols = symbols(&normalize(&counts));
    write_frequencies(table, &symbols.map(|s| s.freq));

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

pub(super) fn encode_order_1(
    data: &[u8],
    table: &mut Vec<u8>,
    rev: &mut Vec<u8>,
) -> Option<[u32; 4]> {
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
            *symbols = self::symbols(&normalize(row));
            *on = true;
        }
    }
    write_runs(table, &active, |table, ctx| {
        if let Some(symbols) = tables.get(ctx) {
            write_frequencies(table, &symbols.map(|s| s.freq));
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

/// Each symbol's slot range and reciprocal.
fn symbols(freq: &[u32; 256]) -> Table {
    let mut table = [Symbol::default(); 256];
    let mut start = 0u32;
    for (sym, &f) in table.iter_mut().zip(freq) {
        *sym = Symbol::new(start, f);
        start = start.saturating_add(f);
    }
    table
}

/// `RansEncPutSymbol`: renormalize (bytes written backwards into `rev`),
/// then encode, dividing by the symbol's reciprocal.
#[inline]
fn put(x: &mut u32, rev: &mut Vec<u8>, sym: Symbol) {
    debug_assert!(sym.freq > 0, "encoding a symbol without frequency");
    let x_max = ((LOWER_BOUND >> 12) << 8).saturating_mul(sym.freq.max(1));
    while *x >= x_max {
        #[allow(clippy::cast_possible_truncation, reason = "the low byte")]
        rev.push(*x as u8);
        *x >>= 8;
    }
    // x < 2^31 after renormalization: (q << 12) + x % f + start, as
    // x + start + q·(4096 − f).
    let q = sym.quotient(*x);
    *x = x
        .wrapping_add(sym.start)
        .wrapping_add(q.wrapping_mul((1u32 << 12).saturating_sub(sym.freq)));
}
