//! htscodecs' decoders, from the static htslib that the rust-htslib
//! dev-dependency links, as the baseline for the CRAM codec benchmarks.
//!
//! Included via `#[path]` by the `cram` and `fqzcomp` bench targets. Each
//! decode returns a buffer htscodecs `malloc`ed, freed on drop, so a timed
//! call allocates and releases its output just as seqair's `Vec`-returning
//! decoders do (and as htslib does per CRAM block).

use std::ffi::{c_char, c_uint};
use std::ptr::NonNull;

// Pulls in hts-sys' static libhts, which contains htscodecs.
use rust_htslib as _;

/// `arith_dynamic.h`, `fqzcomp_qual.h` and `tokenise_name3.h`.
mod ffi {
    use std::ffi::{c_char, c_int, c_uint};

    unsafe extern "C" {
        pub(super) fn arith_uncompress_to(
            input: *mut u8,
            in_size: c_uint,
            out: *mut u8,
            out_size: *mut c_uint,
        ) -> *mut u8;

        pub(super) fn fqz_decompress(
            input: *mut c_char,
            in_size: usize,
            out_size: *mut usize,
            lengths: *mut c_int,
            nlengths: c_int,
        ) -> *mut c_char;

        pub(super) fn tok3_decode_names(
            input: *mut u8,
            sz: c_uint,
            out_len: *mut c_uint,
        ) -> *mut u8;
    }
}

/// `len` bytes htscodecs `malloc`ed; `free`d on drop.
pub struct MallocBuf {
    ptr: NonNull<u8>,
    len: usize,
}

impl MallocBuf {
    /// # Safety
    /// A non-null `ptr` must point at `len` initialised bytes from `malloc`
    /// that nothing else frees.
    unsafe fn from_raw(ptr: *mut u8, len: usize) -> Option<Self> {
        NonNull::new(ptr).map(|ptr| Self { ptr, len })
    }

    pub fn as_slice(&self) -> &[u8] {
        // SAFETY: `from_raw`'s contract; the buffer lives until `self` drops.
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for MallocBuf {
    fn drop(&mut self) {
        // SAFETY: `from_raw`'s contract: malloc'd and owned only by `self`.
        unsafe { libc::free(self.ptr.as_ptr().cast()) };
    }
}

/// `arith_uncompress_to` with a null `out`, so htscodecs sizes and allocates
/// the output; `None` if it rejects the stream.
pub fn arith_uncompress(input: &[u8]) -> Option<MallocBuf> {
    let in_size = c_uint::try_from(input.len()).ok()?;
    let mut out_size: c_uint = 0;
    // SAFETY: `input` is live for the call and `in_size` bytes long;
    // htscodecs only reads it despite the non-const signature. With a null
    // `out` it mallocs the output and stores its length in `out_size`.
    let ptr = unsafe {
        ffi::arith_uncompress_to(
            input.as_ptr().cast_mut(),
            in_size,
            std::ptr::null_mut(),
            &raw mut out_size,
        )
    };
    let len = usize::try_from(out_size).expect("a u32 fits a usize");
    // SAFETY: a non-null return is the malloc'd output of `out_size` bytes.
    unsafe { MallocBuf::from_raw(ptr, len) }
}

/// `fqz_decompress` without record lengths; `None` if it rejects the stream.
pub fn fqz_decompress(input: &[u8]) -> Option<MallocBuf> {
    let mut out_size = 0usize;
    // SAFETY: `input` is live for the call and `input.len()` bytes long;
    // htscodecs only reads it despite the non-const signature. `lengths` may
    // be null when `nlengths` is 0.
    let ptr = unsafe {
        ffi::fqz_decompress(
            input.as_ptr().cast_mut().cast::<c_char>(),
            input.len(),
            &raw mut out_size,
            std::ptr::null_mut(),
            0,
        )
    };
    // SAFETY: a non-null return is the malloc'd output of `out_size` bytes.
    unsafe { MallocBuf::from_raw(ptr.cast(), out_size) }
}

/// `tok3_decode_names`; `None` if it rejects the block.
pub fn tok3_decode_names(input: &[u8]) -> Option<MallocBuf> {
    let size = c_uint::try_from(input.len()).ok()?;
    let mut out_len: c_uint = 0;
    // SAFETY: `input` is live for the call and `size` bytes long; htscodecs
    // only reads it despite the non-const signature.
    let ptr = unsafe { ffi::tok3_decode_names(input.as_ptr().cast_mut(), size, &raw mut out_len) };
    let len = usize::try_from(out_len).expect("a u32 fits a usize");
    // SAFETY: a non-null return is the malloc'd output of `out_len` bytes.
    unsafe { MallocBuf::from_raw(ptr, len) }
}
