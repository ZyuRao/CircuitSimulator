#!/usr/bin/env bash
set -euo pipefail

# =======================
# CircuitSimulator Bench
# - Default: Release build + CSIM_FAST_MATH=ON
# - Measures wall time via /usr/bin/time
# =======================

BUILD_DIR="build"
EXE=""                  # default: <build_dir>/CircuitSimulator
RUNS=5
WARMUP=1
OUT_CSV="bench_summary.csv"
OUT_DIR="bench_raw"
CLEAN=0
NO_BUILD=0
FAST_MATH=1             # default ON (auto enable fast-math)

usage() {
  cat <<EOF
Usage:
  ./bench.sh [options] [netlists...]

Options:
  -b <build_dir>       Build dir (default: build)
  -e <exe_path>        Executable path (default: <build_dir>/CircuitSimulator)
  -r <runs>            Measured runs per netlist (default: 5)
  -w <warmup>          Warmup runs per netlist (default: 1)
  -o <csv_path>        Output summary CSV (default: bench_summary.csv)
  -d <out_dir>         Output raw times dir (default: bench_raw)
  -c                   Clean build dir before build
  --no-build            Do not build; only run benchmark
  --no-fast-math        Disable fast-math (default is ON)
  -h, --help            Show help

CPU pinning to reduce jitter (recommended for stable results):
  TASKSET="taskset -c 0-7" ./bench.sh -c -w 3 -r 20 tests/dbmixer.sp
EOF
}

# Parse args
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) BUILD_DIR="$2"; shift 2;;
    -e) EXE="$2"; shift 2;;
    -r) RUNS="$2"; shift 2;;
    -w) WARMUP="$2"; shift 2;;
    -o) OUT_CSV="$2"; shift 2;;
    -d) OUT_DIR="$2"; shift 2;;
    -c) CLEAN=1; shift 1;;
    --no-build) NO_BUILD=1; shift 1;;
    --no-fast-math) FAST_MATH=0; shift 1;;
    -h|--help) usage; exit 0;;
    *) ARGS+=("$1"); shift 1;;
  esac
done

# Default executable path
if [[ -z "${EXE}" ]]; then
  EXE="${BUILD_DIR}/CircuitSimulator"
fi

# Decide netlists
NETLISTS=("${ARGS[@]}")
if [[ ${#NETLISTS[@]} -eq 0 ]]; then
  if compgen -G "tests/*.sp" > /dev/null; then
    # shellcheck disable=SC2207
    NETLISTS=($(ls -1 tests/*.sp))
  else
    echo "[ERR] No netlists provided and tests/*.sp not found."
    echo "      Example: ./bench.sh tests/buffer.sp"
    exit 1
  fi
fi

# Optional CPU pinning cmd injection
# Example:
#   TASKSET="taskset -c 0-7" ./bench.sh ...
TASKSET_CMD="${TASKSET:-}"

# Clean if requested
if [[ "$CLEAN" -eq 1 ]]; then
  echo "[INFO] Cleaning build dir: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

# Build (Release) unless --no-build
if [[ "$NO_BUILD" -eq 0 ]]; then
  echo "[INFO] Configuring (Release) ..."
  CMAKE_EXTRA=()
  if [[ "$FAST_MATH" -eq 1 ]]; then
    CMAKE_EXTRA+=("-DCSIM_FAST_MATH=ON")
  else
    CMAKE_EXTRA+=("-DCSIM_FAST_MATH=OFF")
  fi

  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release "${CMAKE_EXTRA[@]}"

  echo "[INFO] Building ..."
  cmake --build "${BUILD_DIR}" -j "$(nproc)"
else
  echo "[INFO] --no-build enabled: skip build."
fi

# Check executable
if [[ ! -x "${EXE}" ]]; then
  echo "[ERR] Executable not found or not executable: ${EXE}"
  exit 1
fi

mkdir -p "${OUT_DIR}"

# Environment info (for reproducibility)
META_FILE="${OUT_DIR}/meta.txt"
{
  echo "timestamp: $(date -Is)"
  echo "exe: ${EXE}"
  echo "runs: ${RUNS}"
  echo "warmup: ${WARMUP}"
  echo "build_dir: ${BUILD_DIR}"
  echo "fast_math: $([[ "$FAST_MATH" -eq 1 ]] && echo ON || echo OFF)"
  echo "git_commit: $(git rev-parse HEAD 2>/dev/null || echo N/A)"
  echo "compiler: $(${CXX:-g++} --version 2>/dev/null | head -n 1 || echo N/A)"
  echo "uname: $(uname -a)"
  echo "cpu: $(lscpu 2>/dev/null | grep -E 'Model name|CPU\\(s\\)|Thread|Core|Socket|MHz' || true)"
} > "${META_FILE}"

# Write CSV header
echo "timestamp,netlist,runs,warmup,min_s,mean_s,max_s,std_s" > "${OUT_CSV}"

stats_from_file() {
  awk '
    BEGIN { sum=0; sumsq=0; }
    {
      x=$1;
      if (NR==1) { min=x; max=x; }
      if (x<min) min=x;
      if (x>max) max=x;
      sum += x;
      sumsq += x*x;
    }
    END {
      n=NR;
      mean = sum/n;
      var = (sumsq/n) - (mean*mean);
      if (var < 0) var = 0;
      std = sqrt(var);
      printf "%.6f %.6f %.6f %.6f\n", min, mean, max, std;
    }
  ' "$1"
}

run_one() {
  local net="$1"
  local timefile="$2"
  local errfile="$3"
  /usr/bin/time -f "%e" -o "${timefile}" ${TASKSET_CMD} "${EXE}" "${net}" > /dev/null 2> "${errfile}"
}

echo "[INFO] Benchmark start."
echo "[INFO] Summary CSV: ${OUT_CSV}"
echo "[INFO] Raw times dir: ${OUT_DIR}"
echo "[INFO] Meta saved: ${META_FILE}"
echo

for net in "${NETLISTS[@]}"; do
  if [[ ! -f "$net" ]]; then
    echo "[WARN] Netlist not found: $net (skip)"
    continue
  fi

  base="$(basename "$net")"
  raw="${OUT_DIR}/${base}.times.txt"
  err="${OUT_DIR}/${base}.stderr.txt"

  : > "${raw}"
  : > "${err}"

  echo "[BENCH] ${net}"

  # Warmup
  for ((i=1; i<=WARMUP; i++)); do
    tmp_time="$(mktemp)"
    tmp_err="$(mktemp)"
    run_one "$net" "$tmp_time" "$tmp_err" || true
    rm -f "$tmp_time" "$tmp_err"
    echo "  warmup ${i}/${WARMUP} done"
  done

  # Measured runs
  for ((i=1; i<=RUNS; i++)); do
    tmp_time="$(mktemp)"
    tmp_err="$(mktemp)"
    if run_one "$net" "$tmp_time" "$tmp_err"; then
      t="$(tail -n 1 "$tmp_time")"
      echo "$t" >> "${raw}"
      echo "  run ${i}/${RUNS}: ${t}s"
    else
      echo "[WARN] run failed on ${net} (see stderr log)." >> "${err}"
      cat "$tmp_err" >> "${err}"
      echo "  run ${i}/${RUNS}: FAILED"
    fi
    rm -f "$tmp_time" "$tmp_err"
  done

  if [[ ! -s "${raw}" ]]; then
    echo "  [WARN] No valid timing data for ${net}, skip stats."
    continue
  fi

  read -r min mean max std < <(stats_from_file "${raw}")
  ts="$(date -Is)"

  echo "  stats: min=${min}s mean=${mean}s max=${max}s std=${std}s"
  echo "${ts},${net},${RUNS},${WARMUP},${min},${mean},${max},${std}" >> "${OUT_CSV}"
  echo
done

echo "[INFO] Done."
echo "[INFO] Open summary: ${OUT_CSV}"
