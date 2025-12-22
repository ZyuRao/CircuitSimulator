#!/usr/bin/env bash
set -euo pipefail

# ========= user-configurable defaults =========
BUILD_DIR="${BUILD_DIR:-build}"
NETLIST="${1:-tests/buffer.sp}"
REPEATS="${REPEATS:-3}"          # 每组参数跑几次，取中位数
EXE="${EXE:-$BUILD_DIR/CircuitSimulator}"
OUT_DIR="${OUT_DIR:-out}"
BENCH_DIR="${BENCH_DIR:-$OUT_DIR/bench}"
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"

# 默认 sweep 组合（你也可以改成两数组分别笛卡尔积，当前按“成对组合”跑）
CONFIGS=(
  "16 64"
  "24 64"
  "32 64"
  "48 96"
  "64 128"
)

# timing csv 文件名不确定时，用这个 glob 找最新的
TIMING_GLOB="${TIMING_GLOB:-*timing*.csv}"

usage() {
  echo "Usage: $0 <netlist.sp>"
  echo "Env vars:"
  echo "  REPEATS=3 BUILD_DIR=build EXE=build/CircuitSimulator OUT_DIR=out"
  echo "  TIMING_GLOB='*timing*.csv' (if your timing filename differs)"
  exit 1
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; fi

mkdir -p "$BENCH_DIR"

# ========= helper: build with profiling =========
echo "[1/4] Configure & build (profiling ON) ..."
cmake -S . -B "$BUILD_DIR" -DCSIM_ENABLE_PROF=ON -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE"
cmake --build "$BUILD_DIR" -j

if [[ ! -x "$EXE" ]]; then
  echo "ERROR: executable not found/executable: $EXE"
  echo "Set EXE=... or check your build output."
  exit 2
fi

if [[ ! -f "$NETLIST" ]]; then
  echo "ERROR: netlist not found: $NETLIST"
  exit 3
fi

# ========= warn about .PROBE =========
if grep -Eiq '^\s*\.probe\b' "$NETLIST"; then
  echo "WARNING: Netlist contains .PROBE. Your run may include plotting time (hb_plot)."
  echo "         For clean LU benchmarking, consider using a netlist without .PROBE."
fi

# ========= helper: find newest timing csv =========
find_latest_timing_csv() {
  # search OUT_DIR for newest matching csv after each run
  # GNU find: print mtime epoch + path
  local latest
  latest=$(find "$OUT_DIR" -maxdepth 3 -type f -name "$TIMING_GLOB" -printf '%T@ %p\n' 2>/dev/null \
          | sort -nr | head -n 1 | awk '{print $2}')
  echo "${latest:-}"
}

# ========= helper: parse timing csv -> key metrics (no pandas) =========
# outputs: total_wall,hb_run,solve_factorize,assemble_J,hb_plot
parse_timing_csv() {
  local csv="$1"
  python3 - <<'PY' "$csv"
import csv, sys, math
path = sys.argv[1]
want = {
  "total_wall": None,
  "hb_run": None,
  "HB.solve_factorize": None,
  "HB.assemble_J": None,
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
  g("HB.assemble_J"),
  g("hb_plot"),
]))
PY
}

# ========= main run loop =========
SUMMARY_CSV="$BENCH_DIR/summary_runs.csv"
echo "block,threshold,run,netlist,total_wall_ms,hb_run_ms,solve_factorize_ms,assemble_J_ms,hb_plot_ms,timing_csv" > "$SUMMARY_CSV"

echo "[2/4] Running sweep on netlist: $NETLIST"
echo "      repeats per config: $REPEATS"
echo "      results -> $BENCH_DIR"

# mark time to filter out old timing files (optional safety)
START_TS=$(date +%s)

for cfg in "${CONFIGS[@]}"; do
  block=$(echo "$cfg" | awk '{print $1}')
  thr=$(echo "$cfg" | awk '{print $2}')
  echo ""
  echo "== Config: CSIM_LU_BLOCK=$block CSIM_LU_THRESHOLD=$thr =="

  for r in $(seq 1 "$REPEATS"); do
    echo "  - Run $r/$REPEATS ..."
    # run
    CSIM_LU_BLOCK="$block" CSIM_LU_THRESHOLD="$thr" "$EXE" "$NETLIST" >/dev/null

    # find newest timing csv
    timing_csv=$(find_latest_timing_csv)
    if [[ -z "$timing_csv" ]]; then
      echo "ERROR: could not locate timing csv under '$OUT_DIR' (glob=$TIMING_GLOB)"
      echo "       Set TIMING_GLOB to match your filename, e.g. TIMING_GLOB='timing*.csv'"
      exit 4
    fi

    # copy it to bench dir with unique name
    base="$(basename "$NETLIST" .sp)"
    dst="$BENCH_DIR/timing_${base}_b${block}_t${thr}_r${r}.csv"
    cp "$timing_csv" "$dst"

    # parse metrics
    metrics=$(parse_timing_csv "$dst")
    IFS=',' read -r total_wall hb_run solve_fac assemble_J hb_plot <<<"$metrics"

    echo "    total_wall=${total_wall} ms | hb_run=${hb_run} ms | factorize=${solve_fac} ms | assemble_J=${assemble_J} ms | hb_plot=${hb_plot} ms"
    echo "$block,$thr,$r,$NETLIST,$total_wall,$hb_run,$solve_fac,$assemble_J,$hb_plot,$dst" >> "$SUMMARY_CSV"
  done
done

# ========= summarize medians =========
echo ""
echo "[3/4] Computing medians (no pandas) ..."
python3 - <<'PY' "$SUMMARY_CSV"
import csv, statistics, sys
path = sys.argv[1]

rows = []
with open(path, newline="") as f:
  r = csv.DictReader(f)
  for row in r:
    rows.append(row)

def fnum(x):
  try: return float(x)
  except: return float("nan")

groups = {}
for row in rows:
  key = (row["block"], row["threshold"])
  groups.setdefault(key, []).append(row)

# compute medians per group
out = []
for (b,t), rs in groups.items():
  def med(col):
    vals = [fnum(r[col]) for r in rs]
    vals = [v for v in vals if v==v]  # drop NaN
    return statistics.median(vals) if vals else float("nan")
  out.append({
    "block": b,
    "threshold": t,
    "n": len(rs),
    "total_wall_ms_med": med("total_wall_ms"),
    "hb_run_ms_med": med("hb_run_ms"),
    "solve_factorize_ms_med": med("solve_factorize_ms"),
    "assemble_J_ms_med": med("assemble_J_ms"),
    "hb_plot_ms_med": med("hb_plot_ms"),
  })

# sort by hb_run median first, then factorize
out.sort(key=lambda d: (d["hb_run_ms_med"], d["solve_factorize_ms_med"]))

print("\n=== Median summary (sorted by hb_run_ms_med) ===")
print("block thr  n   hb_run_med(ms)  factorize_med(ms)  assembleJ_med(ms)  hb_plot_med(ms)  total_wall_med(ms)")
for d in out:
  print(f"{d['block']:>5} {d['threshold']:>3} {d['n']:>2} "
        f"{d['hb_run_ms_med']:>14.3f} {d['solve_factorize_ms_med']:>17.3f} "
        f"{d['assemble_J_ms_med']:>16.3f} {d['hb_plot_ms_med']:>14.3f} "
        f"{d['total_wall_ms_med']:>16.3f}")

# write summary csv
out_path = path.replace("summary_runs.csv", "summary_medians.csv")
with open(out_path, "w", newline="") as f:
  w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
  w.writeheader()
  w.writerows(out)
print(f"\nWrote: {out_path}")
PY

echo ""
echo "[4/4] Done."
echo "Run-level log : $SUMMARY_CSV"
echo "Timing copies  : $BENCH_DIR/timing_*.csv"
echo "Median summary : $BENCH_DIR/summary_medians.csv"