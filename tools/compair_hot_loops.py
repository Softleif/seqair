#!/usr/bin/env python3
"""Static size of compair's hot loops, from `objdump -d` of the real build.

Do not use cargo-show-asm for this: it builds with -C codegen-units=1, which
changes what LLVM does to these kernels -- it hid ten reloaded Vec headers
per batch step that the default 16-unit build, the one benches link, had.
This reads the rlib the benches actually link.

    cargo build --release -p compair --lib --target x86_64-unknown-linux-gnu
    objdump -d --demangle --x86-asm-syntax=intel --no-show-raw-insn \
        target/x86_64-unknown-linux-gnu/release/deps/libcompair-*.rlib > dis
    python3 tools/compair_hot_loops.py dis [--dump NAME-SUBSTRING]

Strip kernels: the unmasked middle phase, unrolled two steps (so /2 per
step). Batch and pairs: the column loop. Diagonal: the full-chunk loop.
"""
import re
import sys
from collections import Counter

path = sys.argv[1]
dump = sys.argv[sys.argv.index("--dump") + 1] if "--dump" in sys.argv else None

text = open(path).read().splitlines()
# `#` starts a comment in x86 dumps and an immediate in AArch64 ones, whose
# comments start with `;` (Mach-O) or `//` (ELF).
x86 = any(re.search(r"\b[xy]mm\d+\b", line) for line in text)


def uncomment(ins):
    if x86:
        return ins.split("#")[0].strip()
    return re.split(r"\s;|//", ins)[0].strip()


functions = {}
current = None
for line in text:
    m = re.match(r"^[0-9a-f]+ <(.*)>:$", line)
    if m:
        current = m.group(1)
        functions[current] = []
        continue
    m = re.match(r"^\s*([0-9a-f]+):\s+(.*)$", line)
    if m and current is not None:
        functions[current].append((int(m.group(1), 16), uncomment(m.group(2))))

WANTED = {
    "strip-avx2": "vectorize_avx2::<compair::simd::strip_kernel_at",
    "batch-avx2": "vectorize_avx2::<compair::simd::batch_kernel_at",
    "pairs-avx2": "vectorize_avx2::<compair::simd::pairs_kernel_at",
    "banded-avx2": "vectorize_avx2::<compair::simd::banded_kernel_at",
    "strip-avx512": "vectorize_avx512::<compair::simd::strip_kernel_at",
    "batch-avx512": "vectorize_avx512::<compair::simd::batch_kernel_at",
    "pairs-avx512": "vectorize_avx512::<compair::simd::pairs_kernel_at",
    "banded-avx512": "vectorize_avx512::<compair::simd::banded_kernel_at",
    "strip-sse4.2": "vectorize_sse4_2::<compair::simd::strip_kernel_at",
    "batch-sse4.2": "vectorize_sse4_2::<compair::simd::batch_kernel_at",
    # NEON on aarch64, SSE2 on x86-64: the baseline level is inlined into the
    # dispatching function itself.
    "strip-base": "compair::simd::strip_kernel_at",
    "batch-base": "compair::simd::batch_kernel_at",
    "pairs-base": "compair::simd::pairs_kernel_at",
    "banded-base": "compair::simd::banded_kernel_at",
    # the `compair` branch head, for comparison
    "hand-strip": "compair::intrinsics::strip_kernel_avx2",
    "hand-batch": "compair::intrinsics::batch_kernel_avx2",
    "wide-strip": "compair::strips::strip_kernel_wide",
    "wide-batch": "compair::batch::batch_kernel_wide",
}


def loops(body):
    index = {a: i for i, (a, _) in enumerate(body)}
    found = {}
    for i, (a, ins) in enumerate(body):
        m = re.match(r"^(j\w+|b(\.\w+)?|cbn?z|tbn?z)\s+(?:[\w#]+,\s*)*0x([0-9a-f]+)", ins)
        if m:
            t = int(m.group(3), 16)
            if t <= a and t in index:
                s = index[t]
                if s not in found or i < found[s]:
                    found[s] = i
    spans = sorted(found.items())
    return [(a, b) for a, b in spans if not any((c, d) != (a, b) and a <= c and d <= b for c, d in spans)]


if not any("vectorize_avx2" in n for n in functions):
    # aarch64: the NEON instance is inlined into the dispatching function.
    WANTED = {k: v for k, v in WANTED.items() if "avx" not in v and "sse" not in v}
    WANTED["hand-strip"] = "compair::intrinsics::strip_kernel_intrinsics"
    WANTED["hand-batch"] = "compair::intrinsics::batch_kernel_intrinsics"

for label, pattern in WANTED.items():
    names = [n for n in functions if n == pattern or ("::<" in pattern and pattern in n)]
    if not names:
        continue
    body = functions[names[0]]
    kind = "strip" if "strip" in label else ("banded" if "banded" in label else "batch")
    best = None
    for a, b in loops(body):
        ins = [x for _, x in body[a : b + 1]]
        ops = [x.split()[0] for x in ins]
        if not ({"vmulps", "fmul.4s", "mulps"} & set(ops)) or {"call", "bl", "blr"} & set(ops):
            continue
        if kind == "strip" and (len(ins) < 100 or {"vcvtsi2ss", "ucvtf", "scvtf", "cvtsi2ss"} & set(ops)):
            continue
        if kind in ("batch", "banded") and len(ins) < 50:
            continue
        if best is None or len(ins) < len(best):
            best = ins
    if best is None:
        print(f"{label}: no loop")
        continue
    per = len(best) / (2 if kind == "strip" else 1)
    stack = sum(1 for x in best if re.search(r"\[rsp", x))
    mem_gpr = sum(1 for x in best if re.match(r"^(mov|cmp)\s+\w+, qword ptr \[(?!rsp)", x) or re.match(r"^cmp\s+\w+, qword ptr", x))
    calls = sum(1 for _, x in body if x.startswith("call"))
    print(f"{label:<16} {len(best):>4} instr, {per:>5.1f}/step, rsp refs {stack:>3}, gpr loads/cmp-mem {mem_gpr:>3}, calls in fn {calls}")
    if dump and dump in label:
        print("\n".join("    " + x for x in best))
