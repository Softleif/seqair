#!/usr/bin/env python3
"""A/B two compair trees: interleaved profile_150x48 runs, then criterion.

    python3 tools/compair_ab.py BASE_TREE NEW_TREE [race] [criterion]

Both trees are seqair checkouts (a worktree, or an rsync copy on the bench
box). Each is built with `-C llvm-args=-align-loops=64` into its own target
directory: on Zen 2 an unchanged kernel moved 9% in cycles when a change
elsewhere shifted its loop by 16 bytes, and aligning loops removes that layout
noise from the comparison. Set ALIGN_LOOPS=0 to build without it.

Environment:
    CARGO   the cargo command (default `cargo`; on the 3950X box
            CARGO="$HOME/.cargo/bin/cargo +1.98.1")
    LOCK    a command prefix every timed run is wrapped in, e.g.
            LOCK="flock $HOME/bench.lock" or LOCK="lockf -k /path/bench.lock"
    REPS    interleaved repetitions (default 5)
    ROUNDS  profile_150x48 rounds per run (default 20000, x 8 haplotypes)
    ARMS    the race's arms, `label=base_kernel:new_kernel` separated by
            spaces; the default pairs the hand-written lane on `compair`
            with the fearless_simd lane after it

The race reports min and median wall time, and cycles and instructions under
`perf stat` where `perf` exists (Linux); min is the number to quote.
"""

import os
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

CARGO = shlex.split(os.environ.get("CARGO", "cargo"))
LOCK = shlex.split(os.environ.get("LOCK", ""))
REPS = int(os.environ.get("REPS", "5"))
ROUNDS = os.environ.get("ROUNDS", "20000")
ALIGN = os.environ.get("ALIGN_LOOPS", "64")
PERF = shutil.which("perf") is not None
DEFAULT_ARMS = (
    "strips=strips-intrinsics-taps-workspace:strips-simd-taps-workspace "
    "batch=batch-taps-workspace:batch-taps-workspace "
    "candidates=candidates-taps-workspace:candidates-taps-workspace "
    "banded=simd-taps-workspace:simd-taps-workspace "
    "strips-wide=strips-simd-taps-workspace:strips-simd-taps-workspace"
)


def env():
    flags = f"-C llvm-args=-align-loops={ALIGN}" if ALIGN != "0" else ""
    return {**os.environ, "RUSTFLAGS": flags, "CARGO_TARGET_DIR": f"target/ab-align{ALIGN}"}


def build(tree):
    subprocess.run(
        LOCK + CARGO + ["build", "--release", "-p", "compair", "--example", "profile_150x48"],
        cwd=tree, env=env(), check=True,
    )
    return Path(tree) / f"target/ab-align{ALIGN}/release/examples/profile_150x48"


def measure(binary, kernel):
    cmd = [str(binary), kernel, ROUNDS]
    if PERF:
        cmd = ["perf", "stat", "-x,", "-e", "cycles,instructions", "--"] + cmd
    start = time.perf_counter()
    result = subprocess.run(LOCK + cmd, capture_output=True, text=True, check=True)
    wall = time.perf_counter() - start
    counters = {}
    for line in result.stderr.splitlines():
        parts = line.split(",")
        if len(parts) > 2 and parts[2].split(":")[0] in ("cycles", "instructions"):
            try:
                counters[parts[2].split(":")[0]] = float(parts[0])
            except ValueError:
                pass
    return wall, counters.get("cycles"), counters.get("instructions"), result.stdout.split()[-1]


def race(base, new):
    binaries = {"base": build(base), "new": build(new)}
    arms = []
    for arm in os.environ.get("ARMS", DEFAULT_ARMS).split():
        label, kernels = arm.split("=")
        base_kernel, new_kernel = kernels.split(":")
        arms += [(label, "base", base_kernel), (label, "new", new_kernel)]
    samples = {(label, side): [] for label, side, _ in arms}
    for rep in range(REPS):
        for label, side, kernel in arms:
            samples[(label, side)].append(measure(binaries[side], kernel))
        print(f"rep {rep + 1}/{REPS}", file=sys.stderr, flush=True)
    print(f"\n## interleaved, {ROUNDS} rounds x 8 haplotypes, {REPS} reps, align-loops={ALIGN}")
    print(f"{'arm':<12} {'side':<5} {'wall min':>9} {'wall med':>9} {'cycles min':>12} {'instr min':>12}  checksum")
    for (label, side), rows in samples.items():
        walls = [r[0] for r in rows]
        cycles = [r[1] for r in rows if r[1]]
        instrs = [r[2] for r in rows if r[2]]
        counters = f"{min(cycles):>12.4e} {min(instrs):>12.4e}" if cycles and instrs else f"{'-':>12} {'-':>12}"
        print(f"{label:<12} {side:<5} {min(walls):>9.4f} {statistics.median(walls):>9.4f} {counters}  {rows[0][3]}")


def criterion(base, new):
    print("\n## criterion (one run per tree, the lock held for the whole bench)")
    for side, tree in (("base", base), ("new", new)):
        for bench, filt in (("tenspeed", "compair/"), ("align", "workspace")):
            out = subprocess.run(
                LOCK + CARGO + ["bench", "-p", "compair", "--bench", bench, "--", filt],
                cwd=tree, env=env(), capture_output=True, text=True, check=True,
            ).stdout.splitlines()
            for i, line in enumerate(out):
                if line.startswith(("10s/", "align/")):
                    timing = next((t.strip() for t in out[i:i + 3] if "time:" in t), "")
                    print(f"{side:<5} {line.split()[0]:<58} {timing}")


def main():
    base, new = sys.argv[1], sys.argv[2]
    steps = sys.argv[3:] or ["race", "criterion"]
    print(f"# compair A/B: base {base}, new {new}, perf {PERF}, {time.ctime()}")
    if "race" in steps:
        race(base, new)
    if "criterion" in steps:
        criterion(base, new)


main()
