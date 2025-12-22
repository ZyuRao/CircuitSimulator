#!/usr/bin/env bash
set -euo pipefail

# Repeat profiling runs on the same build/netlist and report medians of key stages.
#
# Usage:
#   chmod +x scripts/repeat_profile.sh
#   scripts/repeat_profile.sh [netlist] [repeats]
#
# Env (optional):
#   EXE=build/CircuitSimulator
#   OUT_DIR=out
#   BENCH_DIR=out/bench
#   LU_BLOCK=24
#   LU_THRESHOLD=64
#   TIMING_GLOB="*timing*.csv"   # if your timing file name differs
#
# Notes:
# - Does NOT change plot behavior. If netlist has .PROBE, hb_plot will appear.
# - Assumes you already built with -DCSIM_ENABLE_PROF=ON.

NETLIST="${1:-tests/buffer.sp}"
REPEATS="${2:-5}"

EXE="${EXE:-build/CircuitSimulator}"
OUT_DIR="${OUT_DIR:-out}"
BENCH_DIR="${BENCH_DIR:-$OUT_DIR/bench}"
LU_BLOCK="${LU_BLOCK:-16}"
LU_THRESHOLD="${LU_THRESHOLD:-64}"
TIMING_GLOB="${TIMING_GLOB:-*timing*.csv}"

mkdir -p "$BENCH_DIR"

if [[ ! -x "$EXE" ]]; then
  echo "ERROR: executable not found: $EXE"
  echo "Build first with: cmake -DCSIM_ENABLE_PROF=ON -S . -B build && cmake --build build -j"
  exit 2
fi

if [[ ! -f "$NETLIST" ]]; then
  echo "ERROR: netlist not found: $NETLIST"
  exit 3
fi

echo "Netlist     : $NETLIST"
echo "Repeats     : $REPEATS"
echo "Executable  : $EXE"
echo "LU params   : CSIM_LU_BLOCK=$LU_BLOCK CSIM_LU_THRESHOLD=$LU_THRESHOLD"
echo "Bench dir   : $BENCH_DIR"
echo ""

# Find newest timing csv under OUT_DIR
find_latest_timing_csv() {
  local latest
  latest=$(find "$OUT_DIR" -maxdepth 3 -type f -name "$TIMING_GLOB" -printf '%T@ %p\n' 2>/dev/null \
          | sort -nr | head -n 1 | awk '{print $2}')
  echo "${latest:-}"
}

RUNS_CSV="$BENCH_DIR/repeat_runs.csv"
echo "run,netlist,lu_block,lu_threshold,timing_csv,total_wall_ms,hb_run_ms,J_add_ms,assemble_ms,factorize_ms,hb_plot_ms" > "$RUNS_CSV"

parse_one_timing() {
  local csv="$1"
  python3 - <<'PY' "$csv"
import csv, sys
path = sys.argv[1]
want = {
  "total_wall": None,
  "hb_run": None,
  "HB.J_add_nonlinear": None,
  "HB.assemble_J": None,
  "HB.solve_factorize": None,
  "hb_plot": None,
}
with open(path, newline="") as f:
  r = csv.DictReader(f)
  for row in r:
    st = row.get("stage","").strip()
    if st in want:
      want[st] = float(row["total_ms"])
def g(k):
  v = want[k]
  return "nan" if v is None else f"{v:.6f}"
print(",".join([
  g("total_wall"),
  g("hb_run"),
  g("HB.J_add_nonlinear"),
  g("HB.assemble_J"),
  g("HB.solve_factorize"),
  g("hb_plot"),
]))
PY
}

for i in $(seq 1 "$REPEATS"); do
  echo "== Run $i/$REPEATS =="
  CSIM_LU_BLOCK="$LU_BLOCK" CSIM_LU_THRESHOLD="$LU_THRESHOLD" "$EXE" "$NETLIST" >/dev/null

  timing="$(find_latest_timing_csv)"
  if [[ -z "$timing" ]]; then
    echo "ERROR: timing csv not found under $OUT_DIR (glob=$TIMING_GLOB)"
    exit 4
  fi

  base="$(basename "$NETLIST" .sp)"
  dst="$BENCH_DIR/timing_${base}_b${LU_BLOCK}_t${LU_THRESHOLD}_r${i}.csv"
  cp "$timing" "$dst"

  metrics="$(parse_one_timing "$dst")"
  IFS=',' read -r total_wall hb_run jadd assemble fac hb_plot <<<"$metrics"

  echo "  total_wall=${total_wall} ms | hb_run=${hb_run} ms | J_add=${jadd} ms | assemble=${assemble} ms | factorize=${fac} ms | hb_plot=${hb_plot} ms"
  echo "$i,$NETLIST,$LU_BLOCK,$LU_THRESHOLD,$dst,$total_wall,$hb_run,$jadd,$assemble,$fac,$hb_plot" >> "$RUNS_CSV"
done

echo ""
echo "Wrote run log: $RUNS_CSV"

# Compute medians
python3 - <<'PY' "$RUNS_CSV"
import csv, statistics, sys, math
path = sys.argv[1]
cols = ["total_wall_ms","hb_run_ms","J_add_ms","assemble_ms","factorize_ms","hb_plot_ms"]
vals = {c: [] for c in cols}

with open(path, newline="") as f:
  r = csv.DictReader(f)
  for row in r:
    for c in cols:
      try:
        v = float(row[c])
        if not math.isnan(v):
          vals[c].append(v)
      except:
        pass

def med(c):
  v = vals[c]
  return statistics.median(v) if v else float("nan")

print("\n=== Median summary ===")
for c in cols:
  print(f"{c:14s}: {med(c):10.3f} ms")
PY
