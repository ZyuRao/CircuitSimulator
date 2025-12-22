#!/usr/bin/env bash
set -euo pipefail

# Thread sweep for LU trailing parallel experiment (CSIM_LU_TRAIL_PAR).
#
# Usage:
#   chmod +x scripts/threadsweep_trailpar.sh
#   scripts/threadsweep_trailpar.sh [netlist] [repeats] [threads_list...]
#
# Examples:
#   scripts/threadsweep_trailpar.sh tests/buffer.sp 5 1 2 4 8 16
#   scripts/threadsweep_trailpar.sh tests/buffer.sp 7 1 2 4 8
#
# Env (optional):
#   EXE=build/CircuitSimulator
#   OUT_DIR=out
#   BENCH_DIR=out/thread_sweep_trailpar
#   LU_BLOCK=16
#   LU_THRESHOLD=64
#   CSIM_PARALLEL=1
#   TIMING_GLOB="*timing*.csv"

NETLIST="${1:-tests/buffer.sp}"
REPEATS="${2:-5}"
shift $(( $#>=2 ? 2 : $# )) || true

THREADS=("$@")
if [[ ${#THREADS[@]} -eq 0 ]]; then
  THREADS=(1 2 4 8 16)
fi

EXE="${EXE:-build/CircuitSimulator}"
OUT_DIR="${OUT_DIR:-out}"
BENCH_DIR="${BENCH_DIR:-$OUT_DIR/thread_sweep_trailpar}"
LU_BLOCK="${LU_BLOCK:-16}"
LU_THRESHOLD="${LU_THRESHOLD:-64}"
CSIM_PARALLEL="${CSIM_PARALLEL:-1}"
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

echo "Netlist      : $NETLIST"
echo "Repeats      : $REPEATS"
echo "Threads list : ${THREADS[*]}"
echo "Executable   : $EXE"
echo "LU params    : CSIM_LU_BLOCK=$LU_BLOCK CSIM_LU_THRESHOLD=$LU_THRESHOLD"
echo "HB parallel  : CSIM_PARALLEL=$CSIM_PARALLEL"
echo "Bench dir    : $BENCH_DIR"
echo ""

find_latest_timing_csv() {
  local latest
  latest=$(find "$OUT_DIR" -maxdepth 3 -type f -name "$TIMING_GLOB" -printf '%T@ %p\n' 2>/dev/null \
          | sort -nr | head -n 1 | awk '{print $2}')
  echo "${latest:-}"
}

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
  "LU.trailing_update": None,
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
  g("HB.solve_factorize"),
  g("LU.trailing_update"),
  g("HB.J_add_nonlinear"),
  g("HB.assemble_J"),
  g("hb_plot"),
]))
PY
}

RUNS_CSV="$BENCH_DIR/runs.csv"
echo "trail_par,threads,run,netlist,lu_block,lu_threshold,hb_parallel,timing_csv,total_wall_ms,hb_run_ms,factorize_ms,trail_update_ms,J_add_ms,assemble_ms,hb_plot_ms" > "$RUNS_CSV"

BASE="$(basename "$NETLIST" .sp)"

for par in 0 1; do
  echo "============================="
  echo "CSIM_LU_TRAIL_PAR=$par"
  echo "============================="
  for t in "${THREADS[@]}"; do
    echo "----- threads=$t -----"
    for r in $(seq 1 "$REPEATS"); do
      echo "  - run $r/$REPEATS"
      CSIM_LU_TRAIL_PAR="$par" \
      CSIM_LU_BLOCK="$LU_BLOCK" \
      CSIM_LU_THRESHOLD="$LU_THRESHOLD" \
      CSIM_PARALLEL="$CSIM_PARALLEL" \
      CSIM_THREADS="$t" \
        "$EXE" "$NETLIST" >/dev/null

      timing="$(find_latest_timing_csv)"
      if [[ -z "$timing" ]]; then
        echo "ERROR: timing csv not found under $OUT_DIR (glob=$TIMING_GLOB)"
        exit 4
      fi

      dst="$BENCH_DIR/timing_${BASE}_trail${par}_thr${t}_b${LU_BLOCK}_t${LU_THRESHOLD}_r${r}.csv"
      cp "$timing" "$dst"

      metrics="$(parse_one_timing "$dst")"
      IFS=',' read -r total_wall hb_run fac trail jadd assemble hb_plot <<<"$metrics"

      echo "    hb_run=${hb_run} ms | fac=${fac} ms | trail=${trail} ms | J_add=${jadd} ms | wall=${total_wall} ms"
      echo "$par,$t,$r,$NETLIST,$LU_BLOCK,$LU_THRESHOLD,$CSIM_PARALLEL,$dst,$total_wall,$hb_run,$fac,$trail,$jadd,$assemble,$hb_plot" >> "$RUNS_CSV"
    done
  done
done

echo ""
echo "Wrote run log: $RUNS_CSV"

python3 - <<'PY' "$RUNS_CSV" "$BENCH_DIR"
import csv, statistics, sys, math, os
path, outdir = sys.argv[1], sys.argv[2]

def fnum(x):
  try: return float(x)
  except: return float("nan")

rows=[]
with open(path, newline="") as f:
  r=csv.DictReader(f)
  rows=list(r)

groups={}
for row in rows:
  key=(int(row["trail_par"]), int(row["threads"]))
  groups.setdefault(key, []).append(row)

def median(vals):
  vals=[v for v in vals if not math.isnan(v)]
  return statistics.median(vals) if vals else float("nan")

out=[]
for (par,t),rs in sorted(groups.items()):
  out.append({
    "trail_par": par,
    "threads": t,
    "n": len(rs),
    "hb_run_med_ms": median([fnum(r["hb_run_ms"]) for r in rs]),
    "factorize_med_ms": median([fnum(r["factorize_ms"]) for r in rs]),
    "trail_update_med_ms": median([fnum(r["trail_update_ms"]) for r in rs]),
    "J_add_med_ms": median([fnum(r["J_add_ms"]) for r in rs]),
    "assemble_med_ms": median([fnum(r["assemble_ms"]) for r in rs]),
    "hb_plot_med_ms": median([fnum(r["hb_plot_ms"]) for r in rs]),
    "total_wall_med_ms": median([fnum(r["total_wall_ms"]) for r in rs]),
  })

# Pretty print: sort each mode by hb_run
print("\n=== Median summary: LU trailing parallel OFF (trail_par=0) ===")
off=[d for d in out if d["trail_par"]==0]
off.sort(key=lambda d:d["hb_run_med_ms"])
print("thr  n   hb_run_med  factorize_med  trail_update_med  J_add_med  total_wall_med")
for d in off:
  print(f"{d['threads']:>3} {d['n']:>2} {d['hb_run_med_ms']:>11.3f} {d['factorize_med_ms']:>13.3f} "
        f"{d['trail_update_med_ms']:>16.3f} {d['J_add_med_ms']:>9.3f} {d['total_wall_med_ms']:>14.3f}")

print("\n=== Median summary: LU trailing parallel ON (trail_par=1) ===")
on=[d for d in out if d["trail_par"]==1]
on.sort(key=lambda d:d["hb_run_med_ms"])
print("thr  n   hb_run_med  factorize_med  trail_update_med  J_add_med  total_wall_med")
for d in on:
  print(f"{d['threads']:>3} {d['n']:>2} {d['hb_run_med_ms']:>11.3f} {d['factorize_med_ms']:>13.3f} "
        f"{d['trail_update_med_ms']:>16.3f} {d['J_add_med_ms']:>9.3f} {d['total_wall_med_ms']:>14.3f}")

# Write CSV
out_csv=os.path.join(outdir,"medians.csv")
with open(out_csv,"w",newline="") as f:
  w=csv.DictWriter(f, fieldnames=list(out[0].keys()))
  w.writeheader(); w.writerows(out)

print(f"\nWrote: {out_csv}")


echo ""
echo "Done. All timing CSVs are under: $BENCH_DIR"
echo "Median summary CSV: $BENCH_DIR/medians.csv"

