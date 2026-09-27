//! Latency of dependent instruction chains, on the machine it runs on.
//!
//! `compair_latency NAME N` runs `N` iterations of a ten-deep chain of the
//! named sequence, each link consuming the previous one's result, and prints
//! nanoseconds per link. To get cycles:
//!
//! - x86-64: `perf stat -e cycles:u` and divide by `10 * N`;
//!   `taskset -c` one idle core.
//! - AArch64 (no user cycle counter on macOS): divide by the `iadd` chain's
//!   time, which is one cycle per link. Take the minimum of several runs;
//!   the first one after idle catches the clock ramping.
//!
//! The chains are the ones on compair's carried paths (notes section 14):
//! the flush (`cmp` + `andn`), the lane shift in its three x86 spellings,
//! mul + add, and the FP/integer-domain bypasses between them.
//!
//! ```bash
//! cargo build --release --manifest-path tools/compair_latency/Cargo.toml
//! perf stat -x, -e cycles:u tools/compair_latency/target/release/compair_latency vpermps 20000000
//! ```
use std::arch::asm;
use std::hint::black_box;
use std::time::Instant;

macro_rules! rep10 { ($s:literal) => { concat!($s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n", $s, "\n") }; }

#[cfg(target_arch = "x86_64")]
fn run(name: &str, n: u64) -> bool {
    unsafe {
        macro_rules! chain { ($body:literal) => {{
            asm!("vxorps ymm0, ymm0, ymm0", "vxorps ymm1, ymm1, ymm1", "vxorps ymm2, ymm2, ymm2",
                 "vpcmpeqd ymm3, ymm3, ymm3", "vxorps ymm4, ymm4, ymm4",
                 "2:", rep10!($body), "dec {n}", "jnz 2b",
                 n = inout(reg) black_box(n) => _, out("ymm0") _, out("ymm1") _, out("ymm2") _, out("ymm3") _, out("ymm4") _, out("ymm5") _, options(nostack));
        }}}
        match name {
            "vaddps" => chain!("vaddps ymm0, ymm0, ymm1"),
            "iadd" => { let mut x: u64 = 0; asm!("2:", rep10!("add {x}, 1"), "dec {n}", "jnz 2b", x = inout(reg) x, n = inout(reg) black_box(n) => _, options(nostack)); black_box(x); }
            "vmulps" => chain!("vmulps ymm0, ymm0, ymm1"),
            "vandps" => chain!("vandps ymm0, ymm0, ymm3"),
            "vmaxps" => chain!("vmaxps ymm0, ymm0, ymm1"),
            "vcmpltps" => chain!("vcmpltps ymm0, ymm0, ymm1"),
            // flush: x = andn(x < c, x)
            "flush" => chain!("vcmpltps ymm2, ymm0, ymm1\nvandnps ymm0, ymm2, ymm0"),
            // integer-domain flush: x = andn(c > x (as int), x)
            "flush_int" => chain!("vpcmpgtd ymm2, ymm1, ymm0\nvpandn ymm0, ymm2, ymm0"),
            "vpcmpgtd" => chain!("vpcmpgtd ymm0, ymm0, ymm1"),
            "vblendvps" => chain!("vblendvps ymm0, ymm0, ymm1, ymm3"),
            "vblendps" => chain!("vblendps ymm0, ymm0, ymm1, 1"),
            "vpermps" => chain!("vpermps ymm0, ymm4, ymm0"),
            "shift_perm" => chain!("vpermps ymm0, ymm4, ymm0\nvblendps ymm0, ymm0, ymm1, 1"),
            // vinsertf128 + vpalignr: ymm5 = [bcast | lo(x)], x = alignr(x, ymm5, 12)
            "shift_alignr" => chain!("vinsertf128 ymm5, ymm1, xmm0, 1\nvpalignr ymm0, ymm0, ymm5, 12"),
            "shift_perm2f" => chain!("vperm2f128 ymm5, ymm0, ymm1, 0x02\nvpalignr ymm0, ymm0, ymm5, 12"),
            // vinsertf128 + two vshufps, all FP domain: [t3,t3,x0,x0] then [t3,x0,x1,x2]
            "shift_shufps" => chain!("vinsertf128 ymm5, ymm1, xmm0, 1\nvshufps ymm5, ymm5, ymm0, 0x0f\nvshufps ymm0, ymm5, ymm0, 0x98"),
            // as the kernel sees it: FP op in, shift, FP op out
            "and_alignr_mul" => chain!("vandps ymm0, ymm0, ymm3\nvinsertf128 ymm5, ymm1, xmm0, 1\nvpalignr ymm0, ymm0, ymm5, 12\nvmulps ymm0, ymm0, ymm2"),
            "and_shufps_mul" => chain!("vandps ymm0, ymm0, ymm3\nvinsertf128 ymm5, ymm1, xmm0, 1\nvshufps ymm5, ymm5, ymm0, 0x0f\nvshufps ymm0, ymm5, ymm0, 0x98\nvmulps ymm0, ymm0, ymm2"),
            "and_perm_mul" => chain!("vandps ymm0, ymm0, ymm3\nvpermps ymm0, ymm4, ymm0\nvblendps ymm0, ymm0, ymm1, 1\nvmulps ymm0, ymm0, ymm2"),
            "and_mul" => chain!("vandps ymm0, ymm0, ymm3\nvmulps ymm0, ymm0, ymm2"),
            "vpalignr" => chain!("vpalignr ymm0, ymm0, ymm1, 12"),
            "vinsertf128" => chain!("vinsertf128 ymm0, ymm1, xmm0, 1"),
            "vperm2f128" => chain!("vperm2f128 ymm0, ymm0, ymm1, 0x02"),
            // mul then add, the recurrence's a*b + c
            "muladd" => chain!("vmulps ymm0, ymm0, ymm1\nvaddps ymm0, ymm0, ymm2"),
            // fp result into integer op and back: bypass
            "add_pandn" => chain!("vaddps ymm0, ymm0, ymm1\nvpandn ymm0, ymm3, ymm0"),
            "add_andnps" => chain!("vaddps ymm0, ymm0, ymm1\nvandnps ymm0, ymm2, ymm0"),
            // running.last() round trip
            "last_rt" => chain!("vextractf128 xmm5, ymm0, 1\nvshufps xmm5, xmm5, xmm5, 0xff\nvbroadcastss ymm5, xmm5\nvblendps ymm0, ymm0, ymm5, 0x80"),
            "vextract_store_reload" => chain!("vextractf128 xmm5, ymm0, 1\nvaddps ymm0, ymm0, ymm5"),
            _ => return false,
        }
    }
    true
}

#[cfg(target_arch = "aarch64")]
fn run(name: &str, n: u64) -> bool {
    unsafe {
        macro_rules! chain { ($body:literal) => {{
            asm!("movi v0.16b, #0", "movi v1.16b, #0", "movi v2.16b, #0", "movi v3.16b, #0xff", "movi v4.16b, #0",
                 "2:", rep10!($body), "subs {n}, {n}, #1", "b.ne 2b",
                 n = inout(reg) black_box(n) => _, out("v0") _, out("v1") _, out("v2") _, out("v3") _, out("v4") _, out("v5") _, options(nostack));
        }}}
        match name {
            "fadd" => chain!("fadd v0.4s, v0.4s, v1.4s"),
            "iadd" => { let mut x: u64 = 0; asm!("2:", rep10!("add {x}, {x}, #1"), "subs {n}, {n}, #1", "b.ne 2b", x = inout(reg) x, n = inout(reg) black_box(n) => _, options(nostack)); black_box(x); }
            "fcmgt_bic_fmul" => chain!("fmul v0.4s, v0.4s, v1.4s\nfcmgt v2.4s, v1.4s, v0.4s\nbic v0.16b, v0.16b, v2.16b"),
            "fadd_and" => chain!("fadd v0.4s, v0.4s, v1.4s\nand v0.16b, v0.16b, v3.16b"),
            "fmul_fadd_indep" => chain!("fmul v0.4s, v0.4s, v1.4s\nfadd v0.4s, v0.4s, v2.4s\nfmul v0.4s, v0.4s, v1.4s"),
            "fmul" => chain!("fmul v0.4s, v0.4s, v1.4s"),
            "fmla" => chain!("fmla v0.4s, v1.4s, v2.4s"),
            "and" => chain!("and v0.16b, v0.16b, v3.16b"),
            "fmax" => chain!("fmax v0.4s, v0.4s, v1.4s"),
            "fcmgt" => chain!("fcmgt v0.4s, v1.4s, v0.4s"),
            "flush" => chain!("fcmgt v2.4s, v1.4s, v0.4s\nbic v0.16b, v0.16b, v2.16b"),
            "flush_int" => chain!("cmgt v2.4s, v1.4s, v0.4s\nbic v0.16b, v0.16b, v2.16b"),
            "bsl" => chain!("bsl v0.16b, v1.16b, v2.16b"),
            "ext" => chain!("ext v0.16b, v1.16b, v0.16b, #12"),
            // f32x8 as two q regs: shift_in = ext hi from lo, ext lo from scalar
            "ins" => chain!("mov v0.s[0], v1.s[0]"),
            "muladd" => chain!("fmul v0.4s, v0.4s, v1.4s\nfadd v0.4s, v0.4s, v2.4s"),
            "dup_lane" => chain!("dup v0.4s, v0.s[3]"),
            _ => return false,
        }
    }
    true
}

fn main() {
    let mut args = std::env::args().skip(1);
    let name = args.next().unwrap_or_default();
    let n: u64 = args.next().and_then(|s| s.parse().ok()).unwrap_or(100_000_000);
    let t = Instant::now();
    if !run(&name, n) {
        eprintln!("unknown test {name}");
        std::process::exit(2);
    }
    let ns = t.elapsed().as_secs_f64() * 1e9 / (10.0 * n as f64);
    println!("{name} {ns:.4} ns/op");
}
