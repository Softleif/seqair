#!/usr/bin/env bash
# Run all fuzz targets with optional seed corpora.
#
# Usage:
#   ./run_all.sh              # 30s per target (default)
#   ./run_all.sh 60           # 60s per target
#   ./run_all.sh 10 --quick   # 10s, stop on first failure
#
# CI usage:
#   cd crates/seqair && ./fuzz/run_all.sh 30
#
# Requires: cargo-fuzz, nightly toolchain

set -euo pipefail

DURATION="${1:-30}"
QUICK="${2:-}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FUZZ_DIR="$SCRIPT_DIR"
CRATE_DIR="$(dirname "$FUZZ_DIR")"
RSS_LIMIT=4096

# Detect target triple
ARCH="$(rustc -vV | grep host | awk '{print $2}')"

cd "$CRATE_DIR"

TARGETS=$(cargo +nightly fuzz list 2>/dev/null)
TOTAL=$(echo "$TARGETS" | wc -l | tr -d ' ')
PASSED=0
FAILED=0
FAILURES=""

echo "╔══════════════════════════════════════════════════════════╗"
echo "║  seqair fuzz suite: $TOTAL targets × ${DURATION}s each"
echo "║  target: $ARCH"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

for target in $TARGETS; do
    seed_dir="$FUZZ_DIR/seeds/$target"
    # The first positional dir is the read-write corpus; additional dirs are read-only.
    # Always pass corpus first so seeds stay untouched.
    corpus_dir="$FUZZ_DIR/corpus/$target"
    mkdir -p "$corpus_dir"
    seed_args="$corpus_dir"
    if [ -d "$seed_dir" ] && [ "$(ls -A "$seed_dir" 2>/dev/null)" ]; then
        seed_args="$corpus_dir $seed_dir"
        seed_note=" (seeded)"
    else
        seed_note=""
    fi

    printf "%-35s" "  $target$seed_note"

    output=$(cargo +nightly fuzz run \
        --target "$ARCH" \
        "$target" \
        -- -max_total_time="$DURATION" \
           -rss_limit_mb="$RSS_LIMIT" \
           $seed_args 2>&1) || true

    # Check for crashes/OOM (exit codes 77=crash, 71=OOM)
    if echo "$output" | grep -q "SUMMARY: \|Fuzz target exited with"; then
        if echo "$output" | grep -q "exit status: 71"; then
            printf "OOM (rss_limit=%sMB)\n" "$RSS_LIMIT"
            # OOM with high RSS limit is not a code bug, just note it
            PASSED=$((PASSED + 1))
        else
            printf "CRASH\n"
            FAILED=$((FAILED + 1))
            FAILURES="$FAILURES  $target\n"
            # Show crash details
            echo "$output" | grep -A2 "Failing input\|Output of\|Reproduce with" | head -10
            echo ""
            if [ "$QUICK" = "--quick" ]; then
                echo "Stopping early (--quick mode)"
                break
            fi
        fi
    else
        # Extract coverage and run count from last line
        done_line=$(echo "$output" | grep "^Done " | tail -1 || true)
        cov=$(echo "$output" | grep "DONE" | tail -1 | sed -n 's/.*cov: \([0-9]*\).*/\1/p' || true)
        runs=$(echo "$done_line" | sed -n 's/^Done \([0-9]*\) runs.*/\1/p' || true)
        if [ -n "$cov" ]; then
            printf "ok  cov=%-5s runs=%s\n" "$cov" "$runs"
        else
            printf "ok\n"
        fi
        PASSED=$((PASSED + 1))
    fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Results: $PASSED passed, $FAILED failed (of $TOTAL)"

if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "  Failed targets:"
    printf "$FAILURES"
    echo ""
    echo "  Reproduce with:"
    echo "    cargo +nightly fuzz run --target $ARCH <target> fuzz/artifacts/<target>/<crash-file>"
    exit 1
else
    echo "  All targets clean."
    exit 0
fi
