#!/usr/bin/env python3
"""A hot loop from `compair_hot_loops.py --dump`, as llvm-mca input.

    python3 tools/compair_hot_loops.py dis --dump avx2 > loops
    python3 tools/compair_mca.py loops strip-avx2 > strip.s
    llvm-mca -mtriple=x86_64-unknown-linux-gnu -mcpu=znver2 -iterations=300 \\
        -timeline -bottleneck-analysis strip.s

AArch64 (a Mach-O rlib disassembles in Apple's `fmul.4s v0, v1, v2[0]`
syntax, which llvm-mca does not parse; it is rewritten to the generic
`fmul v0.4s, v1.4s, v2.s[0]`):

    python3 tools/compair_mca.py loops batch-base > batch.s
    llvm-mca -mtriple=aarch64 -mcpu=apple-m4 -iterations=300 batch.s

llvm-mca is in Homebrew's `llvm` (and `llvm@21`) keg, not in rustup's
llvm-tools. Branch targets become one local label; the loop is the whole
region, taken branches and all, so a loop whose hot path skips a block
(the batch kernel's `summing`) is over-counted by that block.

Treat the models with care (see the notes, section 14): the znver2 model
dispatches four wide where Zen 2 dispatches six, and the Apple models
predate the M4's latencies. Use them for the dependency graph and the
per-resource split, and measure cycles.
"""
import re
import sys

EL = {"4s": "s", "2s": "s", "2d": "d", "1d": "d", "16b": "b", "8b": "b", "8h": "h", "4h": "h"}


def generic(line):
    """Apple NEON syntax to generic: `op.4s v0, v1[2]` -> `op v0.4s, v1.s[2]`."""
    m = re.match(r"^(\w+)\.(\w+)(\s+)(.*)$", line)
    if not m:
        return line
    op, arrangement, space, args = m.groups()
    if arrangement in EL:
        el = EL[arrangement]
        args = re.sub(r"\bv(\d+)\[(\d+)\]", lambda x: f"v{x.group(1)}.{el}[{x.group(2)}]", args)
        args = re.sub(r"\bv(\d+)\b(?![.\[])", lambda x: f"v{x.group(1)}.{arrangement}", args)
    elif arrangement in ("s", "d", "b", "h"):
        # lane moves and single-lane loads/stores: `mov.s v0[1], v1[0]`,
        # `st1.s { v3 }[3], [x5]`
        args = re.sub(r"\bv(\d+)\[(\d+)\]", lambda x: f"v{x.group(1)}.{arrangement}[{x.group(2)}]", args)
        args = re.sub(r"\{ v(\d+) \}", lambda x: f"{{ v{x.group(1)}.{arrangement} }}", args)
    else:
        return line
    return f"{op}{space}{args}"


def main():
    path, label = sys.argv[1], sys.argv[2]
    lines, grab = [], False
    for line in open(path):
        if re.match(r"^\S", line):
            grab = line.startswith(label + " ")
            continue
        if grab and line.startswith("    "):
            lines.append(line.strip())
    if not lines:
        sys.exit(f"no loop labelled {label} in {path}")
    x86 = any(re.search(r"\b[xy]mm\d+\b", l) for l in lines)
    comment = "#" if x86 else "//"
    out = [".intel_syntax noprefix"] if x86 else []
    out += [f"{comment} LLVM-MCA-BEGIN {label}", ".Ltop:"]
    for l in lines:
        l = re.sub(r"\s*<.*$", "", l)
        l = re.sub(r"^((?:j\w+|b(?:\.\w+)?|cbn?z|tbn?z)\s+(?:[\w#]+,\s*)*)0x[0-9a-f]+$", r"\1.Ltop", l)
        out.append(l if x86 else generic(l))
    out.append(f"{comment} LLVM-MCA-END")
    print("\n".join(out))


if __name__ == "__main__":
    main()
