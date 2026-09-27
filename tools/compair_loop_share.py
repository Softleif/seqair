#!/usr/bin/env python3
"""Measured cycles per loop: each loop's share of the samples in a function.

A whole-run `perf stat` gives cycles per vector step averaged over
everything the run does (plan fill, renormalisation, prologues). This splits
it by loop: every backward branch in the function's disassembly defines a
loop span, and the samples whose instruction pointer falls in the span are
that loop's share. Multiply by the run's cycles (or wall time) and divide by
the loop's iteration count to get cycles per iteration.

Linux (perf):

    perf record -e cycles:u -c 20011 -o p.data -- BINARY ARGS
    perf script -F ip,sym,symoff -i p.data > p.ips
    objdump -d --no-show-raw-insn -M intel -C BINARY > p.dis
    python3 tools/compair_loop_share.py perf p.dis p.ips batch_kernel_at [CYCLES]

The binary is position-independent, so the sample addresses are rebased
through the first sample that `perf script` resolved into the function.

macOS (samply, raw addresses relative to the binary's __TEXT):

    samply record --save-only -r 4000 -o p.json -- BINARY ARGS
    llvm-objdump -d --no-show-raw-insn --demangle BINARY > p.dis
    python3 tools/compair_loop_share.py samply p.dis p.json strip_kernel [SECONDS]

Only loops above 1 % of the samples are printed. Nested spans all print; the
innermost one with the hot share is the loop.
"""
import json
import re
import sys
from collections import Counter

MACHO_TEXT = 0x100000000


def functions(dis):
    found, current = {}, None
    for line in open(dis):
        m = re.match(r"^([0-9a-f]+) <(.*)>:$", line)
        if m:
            current = m.group(2)
            found[current] = []
            continue
        m = re.match(r"^\s*([0-9a-f]+):\s+(.*)$", line)
        if m and current:
            found[current].append((int(m.group(1), 16), m.group(2)))
    return found


def perf_samples(ips, funcs, want):
    raw = []
    for line in open(ips):
        parts = line.split()
        if len(parts) >= 2:
            m = re.search(r"\+0x([0-9a-f]+)$", parts[-1])
            raw.append((int(parts[0], 16), " ".join(parts[1:]), int(m.group(1), 16) if m else None))
    starts = {name: body[0][0] for name, body in funcs.items() if body}
    base = 0
    for ip, sym, off in raw:
        if off is not None and want in sym:
            matches = [start for name, start in starts.items() if want in name and name in sym]
            matches = matches or [start for name, start in starts.items() if want in name]
            if matches:
                base = ip - off - matches[0]
                break
    return Counter(ip - base for ip, _, _ in raw)


def samply_samples(profile):
    data = json.load(open(profile))
    leaf = Counter()
    for thread in data["threads"]:
        stacks, frames = thread["stackTable"], thread["frameTable"]
        for stack in thread["samples"]["stack"]:
            if stack is not None:
                leaf[frames["address"][stacks["frame"][stack]] + MACHO_TEXT] += 1
    return leaf


def main():
    mode, dis, samples_path, want = sys.argv[1:5]
    scale = float(sys.argv[5]) if len(sys.argv) > 5 else None
    funcs = functions(dis)
    samples = perf_samples(samples_path, funcs, want) if mode == "perf" else samply_samples(samples_path)
    total = sum(samples.values())
    branch = re.compile(r"^(j\w+|b(?:\.\w+)?|cbn?z|tbn?z)\s+(?:[\w#]+,\s*)*(?:0x)?([0-9a-f]+)\b")
    for name, body in funcs.items():
        if want not in name or not body:
            continue
        lo, hi = body[0][0], body[-1][0]
        in_function = sum(c for a, c in samples.items() if lo <= a <= hi)
        if not in_function:
            continue
        print(f"{name[:120]}\n  function: {in_function / total:.1%}")
        spans = set()
        for address, ins in body:
            m = branch.match(ins)
            if m:
                target = int(m.group(2), 16)
                if lo <= target <= address:
                    spans.add((target, address))
        for start, end in sorted(spans):
            share = sum(c for a, c in samples.items() if start <= a <= end) / total
            if share > 0.01:
                count = sum(1 for a, _ in body if start <= a <= end)
                scaled = f" = {share * scale:.4g}" if scale else ""
                print(f"  loop {start:x}-{end:x} ({count} instr): {share:.1%}{scaled}")


if __name__ == "__main__":
    main()
