//! A TBI or CSI parsed into bins keyed by id, so two indexes compare equal
//! whatever order their bins were written in (htslib writes its hash
//! table's), and a comparison that names the first bin that differs.
#![allow(
    dead_code,
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::cast_sign_loss,
    reason = "test support"
)]

/// A bin's CSI `loffset` (TBI: `None`) and its chunks.
type Bin = (Option<u64>, Vec<(u64, u64)>);

/// One reference's index: bins (id → loffset for CSI, chunks), linear index.
#[derive(Debug, PartialEq, Eq)]
struct RefIndex {
    bins: std::collections::BTreeMap<u32, Bin>,
    linear: Vec<u64>,
}

/// A parsed TBI or CSI: the references' indexes and the `n_no_coor` trailer.
#[derive(Debug, PartialEq, Eq)]
pub struct ParsedIndex {
    refs: Vec<RefIndex>,
    n_no_coor: Option<u64>,
}

struct Cursor<'a>(&'a [u8]);
impl Cursor<'_> {
    fn take(&mut self, n: usize) -> &[u8] {
        let (head, tail) = self.0.split_at(n);
        self.0 = tail;
        head
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn u32(&mut self) -> u32 {
        u32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn u64(&mut self) -> u64 {
        u64::from_le_bytes(self.take(8).try_into().unwrap())
    }
}

/// Parse a decompressed TBI (`csi == false`) or CSI. Bins are keyed by id, so
/// htslib's hash-table bin order does not matter.
pub fn parse_index(data: &[u8], csi: bool) -> ParsedIndex {
    let mut c = Cursor(data);
    let n_ref = if csi {
        assert_eq!(c.take(4), b"CSI\x01");
        let _min_shift = c.i32();
        let _depth = c.i32();
        let l_aux = c.i32() as usize;
        c.take(l_aux);
        c.i32()
    } else {
        assert_eq!(c.take(4), b"TBI\x01");
        let n_ref = c.i32();
        c.take(24); // format, col_seq, col_beg, col_end, meta, skip
        let l_nm = c.i32() as usize;
        c.take(l_nm);
        n_ref
    };
    let refs = (0..n_ref)
        .map(|_| {
            let n_bin = c.i32();
            let bins = (0..n_bin)
                .map(|_| {
                    let id = c.u32();
                    let loffset = csi.then(|| c.u64());
                    let n_chunk = c.i32();
                    let chunks = (0..n_chunk).map(|_| (c.u64(), c.u64())).collect();
                    (id, (loffset, chunks))
                })
                .collect();
            let linear = if csi {
                Vec::new()
            } else {
                let n_intv = c.i32();
                (0..n_intv).map(|_| c.u64()).collect()
            };
            RefIndex { bins, linear }
        })
        .collect();
    let n_no_coor = (c.0.len() == 8).then(|| c.u64());
    assert!(c.0.is_empty(), "{} bytes left over", c.0.len());
    ParsedIndex { refs, n_no_coor }
}

/// Equal, or a panic naming the first bin (or linear-index entry) that differs.
pub fn assert_same_index(ours: &ParsedIndex, theirs: &ParsedIndex, what: &str) {
    assert_eq!(ours.refs.len(), theirs.refs.len(), "{what}: n_ref");
    for (i, (a, b)) in ours.refs.iter().zip(&theirs.refs).enumerate() {
        let ids: std::collections::BTreeSet<_> = a.bins.keys().chain(b.bins.keys()).collect();
        for id in ids {
            assert_eq!(a.bins.get(id), b.bins.get(id), "{what}: ref {i} bin {id} (ours, theirs)");
        }
        assert_eq!(a.linear, b.linear, "{what}: ref {i} linear index");
    }
    assert_eq!(ours.n_no_coor, theirs.n_no_coor, "{what}: n_no_coor");
}
