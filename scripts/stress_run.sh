#!/usr/bin/env bash
set -euo pipefail

# Simple stability stress test: run the same netlist repeatedly.
# Usage:
#   scripts/stress_run.sh [netlist] [repeats]

NETLIST="${1:-tests/buffer.sp}"
REPEATS="${2:-200}"
EXE="${EXE:-build/CircuitSimulator}"

if [[ ! -x "$EXE" ]]; then
  echo "ERROR: executable not found: $EXE"
  exit 2
fi

if [[ ! -f "$NETLIST" ]]; then
  echo "ERROR: netlist not found: $NETLIST"
  exit 3
fi

echo "Netlist: $NETLIST"
echo "Repeats: $REPEATS"
echo "Exe    : $EXE"

fail=0
for i in $(seq 1 "$REPEATS"); do
  echo "== Run $i/$REPEATS =="
  if ! "$EXE" "$NETLIST" >/dev/null; then
    echo "Run $i failed"
    fail=1
    break
  fi
done

if [[ "$fail" -ne 0 ]]; then
  echo "Stress run failed"
  exit 1
fi

echo "Stress run passed"
