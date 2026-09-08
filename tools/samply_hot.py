#!/usr/bin/env python3
"""Symbolicate and summarise a samply profile recorded with `--save-only`.

A `--save-only` profile carries raw addresses, so the Firefox-profiler UI (and
any tool that just reads `funcTable.name`) shows `0x115b8` instead of a
function. Record with the sidecar as well and this joins the two:

    samply record --save-only --unstable-presymbolicate -o p.json -- ./prog
    python3 tools/samply_hot.py p.json                 # hot self time
    python3 tools/samply_hot.py p.json --callers foo   # who reaches foo
    python3 tools/samply_hot.py p.json --tree          # inclusive tree

Self time is the leaf frame of each sample; inclusive time counts a function
once per sample whose stack contains it (recursion is not double counted).
"""

import argparse
import bisect
import json
import re
from collections import defaultdict
from pathlib import Path


class Symbols:
    """Address → function name, from a `--unstable-presymbolicate` sidecar."""

    def __init__(self, path: Path):
        self.by_debug_id: dict[str, tuple[list[int], list[tuple[int, str]]]] = {}
        if not path.exists():
            return
        blob = json.loads(path.read_text())
        strings = blob["string_table"]
        for lib in blob["data"]:
            table = sorted(lib["symbol_table"], key=lambda e: e["rva"])
            starts = [e["rva"] for e in table]
            entries = [(e["rva"] + e.get("size", 0), strings[e["symbol"]]) for e in table]
            self.by_debug_id[lib["debug_id"].replace("-", "").upper()] = (starts, entries)

    def lookup(self, debug_id: str, address: int) -> str | None:
        found = self.by_debug_id.get(debug_id)
        if found is None:
            return None
        starts, entries = found
        i = bisect.bisect_right(starts, address) - 1
        if i < 0:
            return None
        end, name = entries[i]
        # A sample past the symbol's extent belongs to no symbol at all;
        # naming it after the previous one would invent a hot function.
        return name if address < end else None


def short(name: str) -> str:
    """Drop generic arguments and hash suffixes that make names unreadable."""
    name = re.sub(r"::h[0-9a-f]{16}$", "", name)
    out, depth = [], 0
    for ch in name:
        if ch == "<":
            depth += 1
            if depth == 1:
                out.append("<…>")
        elif ch == ">":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def resolve(profile: dict, syms: Symbols) -> list[dict]:
    """Return one dict per thread: name, and a symbolicated frame-name list."""
    lib_debug_id = [lib["breakpadId"][:32].upper() for lib in profile["libs"]]
    lib_name = [lib["name"] for lib in profile["libs"]]
    threads = []
    for thread in profile["threads"]:
        strings = thread["stringArray"]
        funcs = thread["funcTable"]
        resources = thread["resourceTable"]
        frames = thread["frameTable"]
        names, libs = [], []
        for i in range(frames["length"]):
            func = frames["func"][i]
            res = resources["lib"][funcs["resource"][func]] if funcs["resource"][func] >= 0 else None
            address = frames["address"][i]
            name = None
            if res is not None and address is not None and address >= 0:
                name = syms.lookup(lib_debug_id[res], address)
            names.append(short(name) if name else strings[funcs["name"][func]])
            libs.append(lib_name[res] if res is not None else "?")
        threads.append({"name": thread["name"], "frames": names, "libs": libs, "thread": thread})
    return threads


def walk(thread: dict):
    """Yield (leaf_frame, [frames root→leaf]) per sample."""
    stacks = thread["thread"]["stackTable"]
    prefix, frame_of = stacks["prefix"], stacks["frame"]
    chains: list[list[int]] = []
    for i in range(stacks["length"]):
        parent = prefix[i]
        chains.append((chains[parent] if parent is not None else []) + [frame_of[i]])
    for stack in thread["thread"]["samples"]["stack"]:
        if stack is None:
            continue
        chain = chains[stack]
        if chain:
            yield chain[-1], chain


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("profile", type=Path)
    ap.add_argument("--syms", type=Path, help="sidecar (default: <profile>.syms.json)")
    ap.add_argument("-n", "--top", type=int, default=25)
    ap.add_argument("-t", "--thread", help="only threads whose name contains this")
    ap.add_argument("-l", "--lib", help="only frames from libraries matching this")
    ap.add_argument("-c", "--callers", help="show what reaches functions matching this")
    ap.add_argument("--inclusive", action="store_true", help="rank by inclusive time")
    args = ap.parse_args()

    profile = json.loads(args.profile.read_text())
    syms_path = args.syms or args.profile.with_suffix("").with_suffix(".syms.json")
    if not syms_path.exists():
        syms_path = Path(str(args.profile).replace(".json", ".syms.json"))
    syms = Symbols(syms_path)
    threads = resolve(profile, syms)
    if args.thread:
        threads = [t for t in threads if args.thread in t["name"]]

    self_time: dict[str, int] = defaultdict(int)
    incl_time: dict[str, int] = defaultdict(int)
    callers: dict[str, int] = defaultdict(int)
    total = 0
    per_lib: dict[str, int] = defaultdict(int)

    for thread in threads:
        names, libs = thread["frames"], thread["libs"]
        for leaf, chain in walk(thread):
            total += 1
            per_lib[libs[leaf]] += 1
            if args.lib is None or args.lib in libs[leaf]:
                self_time[names[leaf]] += 1
            for name in {names[f] for f in chain}:
                incl_time[name] += 1
            if args.callers:
                for depth, f in enumerate(chain):
                    if args.callers in names[f] and depth:
                        callers[names[chain[depth - 1]]] += 1
                        break

    if not total:
        print("no samples matched")
        return

    print(f"threads: {', '.join(t['name'] for t in threads)}")
    print(f"samples: {total}\n")
    print(f"{'samples':>9} {'share':>7}  library")
    for lib, n in sorted(per_lib.items(), key=lambda kv: -kv[1])[:8]:
        print(f"{n:>9} {n / total:>6.1%}  {lib}")

    if args.callers:
        print(f"\ncallers of *{args.callers}*")
        for name, n in sorted(callers.items(), key=lambda kv: -kv[1])[: args.top]:
            print(f"{n:>9} {n / total:>6.1%}  {name}")
        return

    ranked = incl_time if args.inclusive else self_time
    label = "inclusive" if args.inclusive else "self"
    print(f"\n{'samples':>9} {'share':>7}  function ({label})")
    for name, n in sorted(ranked.items(), key=lambda kv: -kv[1])[: args.top]:
        print(f"{n:>9} {n / total:>6.1%}  {name}")


if __name__ == "__main__":
    main()
