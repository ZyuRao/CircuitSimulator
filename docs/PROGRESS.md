# Performance/Profiling Progress

## Step 0 - Baseline audit (no code changes)
- Changed files: none
- Purpose: Record which optimizations already exist and what is missing before changes.
- Expected timing impact: none.
- Verify: n/a

## Step 1 - Route TimingRegistry into transient paths
- Changed files: include/analysis.hpp, src/tanalisis.cpp, src/runner.cpp, src/hbanalysis.cpp
- Purpose: Allow transient runs (including HB init transient) to receive the same optional TimingRegistry used by DC/HB.
- Expected timing impact: none when profiling is off; enables future TRAN.* timers without additional plumbing.
- Verify: Build with `-DCSIM_ENABLE_PROF=ON`, run `./build/CircuitSimulator tests/buffer.sp`, confirm no functional output change and timing CSV still generated when enabled.

## Step 2 - Cache transient element lists and reuse per-iteration buffers
- Changed files: src/tanalisis.cpp, include/element.hpp
- Purpose: Remove per-iteration dynamic_pointer_cast in transient loops, reuse G/I buffers, and keep branch-current output ordering stable.
- Expected timing impact: reduces TR/BE per-iteration overhead; no change to Newton steps or residuals beyond floating-point reordering.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp` and compare `out/buffer_tran_raw.csv` columns/order and basic waveform sanity; compare iteration warnings and final output file presence.

## Step 3 - HB LU solve split timers + in-place factor/solve
- Changed files: include/solver.hpp, src/hbanalysis.cpp
- Purpose: Split `HB.solve_linear` into `HB.solve_factorize` and `HB.solve_backsolve`, log system size once, and route HB to an in-place LU path with explicit pivot vector.
- Expected timing impact: clearer attribution for factor/backsolve; removes an extra LU copy in HB without changing Newton flow or residual criteria.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp`, compare HB iteration count/residual log sequence and `out/buffer_hb_raw.csv` (or filtered CSV) for consistency.

## Step 4 - LU inner-loop memory access optimization
- Changed files: include/solver.hpp
- Purpose: Reorder elimination loops to favor column-major access in Eigen’s default storage, reducing cache misses in factorization without altering pivoting or arithmetic.
- Expected timing impact: `HB.solve_factorize` (and other LU callers) should drop modestly for mid/large n.
- Verify: Re-run `./build/CircuitSimulator tests/buffer.sp` and compare HB iteration/residual logs and output CSVs to prior step.

## Step 5 - Blocked LU factorization (panel + trailing update)
- Changed files: include/solver.hpp
- Purpose: Add a blocked LU path (`A22 -= L21*U12`) with the same partial pivoting, improving cache behavior for larger systems.
- Expected timing impact: `HB.solve_factorize` should drop more noticeably for n >= ~64; backsolve unchanged.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp` and compare HB iteration/residual logs and output CSVs; optionally run a larger netlist to see factorize timing change.

## Step 6 - Plot gating during profiling
- Changed files: include/runtime.hpp, src/runtime.cpp, src/runner.cpp
- Purpose: Disable plotting by default when profiling is enabled; allow runtime override via `CSIM_PLOT=1`/`0`.
- Expected timing impact: reduces `hb_plot`/`tran_plot` wall time and noise in profiling runs.
- Verify: Run with and without `CSIM_PLOT=1`; ensure no plot PNGs are created by default when profiling is on and that outputs/iterations remain unchanged.

## Step 7 - LU block-size/threshold tuning + noalias trailing update
- Changed files: include/solver.hpp, include/runtime.hpp, src/runtime.cpp
- Purpose: Make block size/threshold configurable via `CSIM_LU_BLOCK`/`CSIM_LU_THRESHOLD` and use `noalias` for `A22 -= L21*U12` to avoid temporaries.
- Expected timing impact: tune `HB.solve_factorize` via block size/threshold; reduce overhead in trailing update.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp` with varying `CSIM_LU_BLOCK`/`CSIM_LU_THRESHOLD` values (e.g., 16/24/32/48/64 and n>=64/96/128), compare HB iteration/residual sequence and output CSVs.

## Step 8 - HB.assemble_J index precompute
- Changed files: src/hbanalysis.cpp
- Purpose: Precompute and cache per-variable harmonic/node/real-imag indices and row-base offsets to reduce integer work inside `HB.assemble_J`.
- Expected timing impact: reduce `HB.assemble_J` CPU time without changing residual/Jacobian math.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp` and compare HB iteration/residual sequence and `out/*_hb_raw.csv` column order/values.

## Step 9 - Plot gating via .PROBE only
- Changed files: include/sim.hpp, src/parser.cpp, src/runner.cpp, include/runtime.hpp, src/runtime.cpp
- Purpose: Default plotting off; enable plotting only when `.PROBE` is present in the netlist; profiling no longer affects plotting.
- Expected timing impact: remove `hb_plot`/`tran_plot` from runs without `.PROBE`; allow timing clarity without behavior changes from profiling.
- Verify: Run a netlist without `.PROBE` (no plot PNGs, no `hb_plot` stage); run a netlist with `.PROBE` (plot PNGs created and `hb_plot` appears).

## Step 10 - HB assemble/Residual sub-timers
- Changed files: src/hbanalysis.cpp
- Purpose: Split HB Jacobian/residual assembly into sub-stages for diagnosis (`HB.J_init`, `HB.F_init`, `HB.J_add_nonlinear`, `HB.F_add_nonlinear`, `HB.gmin_diag`).
- Expected timing impact: instrumentation only; no numerical change beyond FP ordering from separated accumulation.
- Verify: Run the same netlist, compare HB iteration count/residual sequence and `out/*_hb_raw.csv` column order/values; confirm new timing stages appear.

## Step 11 - HB.J_add_nonlinear active-index packing
- Changed files: src/hbanalysis.cpp
- Purpose: Prepack active q indices and per-q dv table pointers to reduce inner-loop branching in `HB.J_add_nonlinear`.
- Expected timing impact: reduce `HB.J_add_nonlinear` time without changing Jacobian math.
- Verify: Run the same netlist and compare HB iteration count/residual sequence and `out/*_hb_raw.csv` values; compare timing CSV for `HB.J_add_nonlinear`.

## Step 12 - HB.J_add_nonlinear column-pointer writes
- Changed files: src/hbanalysis.cpp
- Purpose: Use direct column pointers for Jacobian updates in `HB.J_add_nonlinear` to reduce per-entry indexing overhead.
- Expected timing impact: small reduction in `HB.J_add_nonlinear` for large n; no numerical change.
- Verify: Run the same netlist and compare HB iteration count/residual sequence and `out/*_hb_raw.csv` values; compare timing CSV for `HB.J_add_nonlinear`.

## Step 13 - HB.J_add_nonlinear precomputed active metadata
- Changed files: src/hbanalysis.cpp
- Purpose: Precompute active-q node indices, dv pointers, column pointers, and row-base imag offsets to eliminate inner-loop mapping/branching in `HB.J_add_nonlinear`.
- Expected timing impact: reduce `HB.J_add_nonlinear` by cutting index/branch overhead.
- Verify: Run `./build/CircuitSimulator tests/buffer.sp` (with `.PROBE`) and compare HB iteration/residual sequence and output CSV/plot files; compare timing CSV for `HB.J_add_nonlinear`.

## Step 14 - HB.J_add_nonlinear cached Gnl column pointers
- Changed files: src/hbanalysis.cpp
- Purpose: Precompute `Gnl_t_vec[n_t].col(iq).data()` pointers for every time sample and active q to remove inner-loop indexing/lookup during `HB.J_add_nonlinear`.
- Expected timing impact: reduce `HB.J_add_nonlinear` by avoiding matrix indexing and branch checks.
- Verify: `scripts/repeat_profile.sh tests/buffer.sp 5` (with `.PROBE`, plots must appear) and compare HB iteration/residual sequence and outputs.
- Median comparison (tests/buffer.sp, 5 runs, LU_BLOCK=24, LU_THRESHOLD=64):
  - Before (Step 13 baseline): J_add=3457.041ms, assemble=3458.799ms, factorize=1422.543ms, hb_run=5034.231ms, total_wall=5820.257ms
  - After (Step 14): J_add=3752.182ms, assemble=3754.101ms, factorize=1627.531ms, hb_run=5571.841ms, total_wall=6389.803ms

## Step 14b - Roll back Gnl column pointer table
- Changed files: src/hbanalysis.cpp
- Purpose: Remove the per-(n_t, q) `Gnl` column pointer table (regression vs Step 13) and restore baseline access pattern.
- Expected timing impact: undo the Step 14 regression in `HB.J_add_nonlinear`.
- Verify: `scripts/repeat_profile.sh tests/buffer.sp 5` (with `.PROBE`) and confirm timings return near Step 13 baseline.

## Step 15 - Direct Gnl column pointer arithmetic (no pointer table)
- Changed files: src/hbanalysis.cpp
- Purpose: Access `Gnl` columns via `Gnl.data() + colIdx * N` to avoid `Gnl.col()` and 2D pointer tables inside `HB.J_add_nonlinear`.
- Expected timing impact: reduce inner-loop index/lookup cost in `HB.J_add_nonlinear` without changing math.
- Verify: `scripts/repeat_profile.sh tests/buffer.sp 10` (with `.PROBE`) and compare HB iteration/residual sequence and outputs.
- Median comparison (tests/buffer.sp, 10 runs, LU_BLOCK=24, LU_THRESHOLD=64):
  - Before (Step 13 baseline): J_add=3457.041ms, assemble=3458.799ms, factorize=1422.543ms, hb_run=5034.231ms, total_wall=5820.257ms
  - After (Step 15): J_add=3607.273ms, assemble=3609.122ms, factorize=1516.867ms, hb_run=5296.917ms, total_wall=6081.455ms, hb_plot=666.175ms

## Step 16 - HB.J_add_nonlinear per-column delta buffer (serial)
- Changed files: src/hbanalysis.cpp
- Purpose: Accumulate nonlinear Jacobian contributions into a contiguous `deltaCol` and apply a single sequential column add, reducing scattered writes.
- Expected timing impact: reduce `HB.J_add_nonlinear` by improving write locality; no change to Newton logic.
- Verify: `scripts/repeat_profile.sh tests/buffer.sp 5` (with `.PROBE`) and compare HB iteration/residual sequence and outputs.
- Median comparison (tests/buffer.sp, 5 runs, LU_BLOCK=24, LU_THRESHOLD=64):
  - Before (Step 13 baseline): J_add=3457.041ms, assemble=3458.799ms, factorize=1422.543ms, hb_run=5034.231ms, total_wall=5820.257ms
  - After (Step 16): J_add=3530.492ms, assemble=3532.545ms, factorize=1528.626ms, hb_run=5242.390ms, total_wall=6040.718ms, hb_plot=639.140ms

## Step 17 - HB.J_add_nonlinear parallel by column (per-q)
- Changed files: src/hbanalysis.cpp, src/runtime.cpp
- Purpose: Parallelize per-column nonlinear Jacobian assembly with per-thread scratch buffers; add `CSIM_PARALLEL`/`CSIM_THREADS` runtime controls.
- Expected timing impact: reduce `HB.J_add_nonlinear` and `HB.assemble_J` when threads >1; no change to Newton flow or plot behavior.
- Verify (serial, 1 thread): `CSIM_THREADS=1 scripts/repeat_profile.sh tests/buffer.sp 5`
- Verify (parallel, 8 threads): `CSIM_THREADS=8 scripts/repeat_profile.sh tests/buffer.sp 5`
- Verify (iteration/residual sequence): `CSIM_THREADS=1 ./build/CircuitSimulator tests/buffer.sp > out/hb_serial.log` and `CSIM_THREADS=8 ./build/CircuitSimulator tests/buffer.sp > out/hb_parallel.log`, then `rg "^\[HB\] iter" ...` and `diff -u` (no diff).
- Median comparison (tests/buffer.sp, 5 runs, LU_BLOCK=24, LU_THRESHOLD=64):
  - Serial (CSIM_THREADS=1): J_add=3566.198ms, assemble=3567.966ms, factorize=1524.616ms, hb_run=5258.115ms, total_wall=6066.452ms, hb_plot=646.005ms
  - Parallel (CSIM_THREADS=8): J_add=664.454ms, assemble=666.235ms, factorize=1703.118ms, hb_run=2540.757ms, total_wall=3320.688ms, hb_plot=663.591ms

## Step 18 - LU factorize sub-timers (pivot/panel/trailing)
- Changed files: include/solver.hpp, src/hbanalysis.cpp
- Purpose: Split `HB.solve_factorize` into LU sub-stages (`LU.pivot_search`, `LU.pivot_apply`, `LU.panel_factor`, `LU.trailing_update`) for bottleneck diagnosis without changing math.
- Expected timing impact: profiling-only overhead; no change to Newton flow or plot behavior.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - HB.solve_factorize: 1637.276ms
  - LU.pivot_search: 41.669ms (2.5%)
  - LU.pivot_apply: 100.405ms (6.1%)
  - LU.panel_factor: 206.754ms (12.6%)
  - LU.trailing_update: 1280.953ms (78.2%)

## Step 19 - LU trailing_update block view hoist + noalias audit
- Changed files: include/solver.hpp
- Purpose: Hoist `A22/L21/U12` block views (single construction per block) and enforce `A22.noalias() -= L21 * U12` in the trailing update.
- Expected timing impact: reduce `LU.trailing_update` overhead slightly; no change to Newton flow or plot behavior.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - HB.solve_factorize: 1599.582ms
  - LU.pivot_search: 40.873ms (2.6%)
  - LU.pivot_apply: 99.397ms (6.2%)
  - LU.panel_factor: 191.688ms (12.0%)
  - LU.trailing_update: 1259.350ms (78.7%)

## Step 20 - LU trailing_update GEMM blocking (CSIM_LU_GEMM_BLOCK)
- Changed files: include/runtime.hpp, src/runtime.cpp, include/solver.hpp
- Purpose: Add `CSIM_LU_GEMM_BLOCK` (default = `CSIM_LU_BLOCK`) and use it to tile the trailing update `A22 -= L21*U12`.
- Expected timing impact: tune `LU.trailing_update` without changing math; no change to Newton flow or plot behavior.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 CSIM_LU_GEMM_BLOCK=<val> scripts/repeat_profile.sh tests/buffer.sp 5`
- Verify (iteration/residual sequence): compare `[HB] iter` logs vs baseline (`CSIM_LU_GEMM_BLOCK=16`); no diff observed for 8/16/24/32/48.
- Median sweep (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - gemm_block=8:  HB.solve_factorize=3435.665ms, LU.trailing_update=3066.585ms
  - gemm_block=16: HB.solve_factorize=2419.981ms, LU.trailing_update=2032.634ms
  - gemm_block=24: HB.solve_factorize=2070.106ms, LU.trailing_update=1682.861ms
  - gemm_block=32: HB.solve_factorize=1948.800ms, LU.trailing_update=1572.798ms
  - gemm_block=48: HB.solve_factorize=1838.082ms, LU.trailing_update=1456.881ms

## Step 21 - Disable GEMM tiling by default (experimental-only)
- Changed files: include/runtime.hpp, src/runtime.cpp, include/solver.hpp
- Purpose: Keep `A22.noalias() -= L21*U12` as the default trailing update; only enable tiling/micro-kernel in explicit experimental mode.
- Expected timing impact: restore baseline performance when experiment flag is off.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - HB.solve_factorize: 1621.108ms
  - HB.assemble_J: 434.958ms
  - HB.J_add_nonlinear: 433.037ms
  - hb_run: 2219.807ms
  - total_wall: 3024.610ms

## Step 22 - Experimental nb=16 trailing_update micro-kernel (disabled by default)
- Changed files: include/solver.hpp
- Purpose: Evaluate a manual `nb==16` trailing-update micro-kernel under `CSIM_LU_GEMM_EXPERIMENT=1` (no tiling when `CSIM_LU_GEMM_BLOCK=0`).
- Expected timing impact: reduce `LU.trailing_update` if the micro-kernel is beneficial; no change to Newton flow or plot behavior when experiment is off.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 CSIM_LU_GEMM_EXPERIMENT=1 CSIM_LU_GEMM_BLOCK=0 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64, CSIM_LU_GEMM_EXPERIMENT=1):
  - HB.solve_factorize: 4285.067ms (regression)
  - HB.assemble_J: 440.318ms
  - HB.J_add_nonlinear: 438.346ms
  - hb_run: 4898.867ms
  - total_wall: 5697.600ms
- Residual/iteration check: baseline vs experimental log compare shows mismatched residual sequence; micro-kernel is not acceptable and remains experimental-only.

## Step 23 - Remove micro-kernel path + stabilize defaults
- Changed files: include/solver.hpp, include/runtime.hpp, src/runtime.cpp, README.md
- Purpose: Remove the nb=16 micro-kernel path (performance regression + residual mismatch) and ensure default runs always use `A22.noalias() -= L21 * U12`; set recommended default runtime parameters.
- Expected timing impact: restore baseline LU behavior; improve out-of-box performance with tuned defaults.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - HB.solve_factorize: 1570.257ms
  - HB.assemble_J: 429.889ms
  - HB.J_add_nonlinear: 428.320ms
  - hb_run: 2159.515ms
  - total_wall: 2962.007ms

## Step 24 - Experimental LU.trailing_update parallel by columns
- Changed files: CMakeLists.txt, include/runtime.hpp, src/runtime.cpp, include/solver.hpp, README.md
- Purpose: Add `CSIM_LU_TRAIL_PAR` experimental switch (default off) to parallelize trailing update by disjoint column blocks; keep serial baseline by default.
- Expected timing impact: reduce `LU.trailing_update` when enabled; no change to Newton flow or plot behavior when disabled.
- Verify: `CSIM_THREADS=16 CSIM_LU_BLOCK=16 CSIM_LU_THRESHOLD=64 CSIM_LU_TRAIL_PAR=0|1 scripts/repeat_profile.sh tests/buffer.sp 5`
- Median summary (tests/buffer.sp, 5 runs, CSIM_THREADS=16, CSIM_LU_BLOCK=16, CSIM_LU_THRESHOLD=64):
  - Serial (CSIM_LU_TRAIL_PAR=0): HB.solve_factorize=1631.730ms, hb_run=2234.613ms, total_wall=3041.481ms
  - Parallel (CSIM_LU_TRAIL_PAR=1): HB.solve_factorize=880.439ms, hb_run=1608.997ms, total_wall=2401.166ms
- Iter/residual sequence: serial vs parallel logs match exactly (no diff).

## Step 25 - Update recommended benchmark params + add stability tooling
- Changed files: README.md, scripts/stress_run.sh, docs/PROGRESS.md
- Purpose: Update recommended benchmark env (trail parallel + threads=8) and add a 200x stress-run script + sanitizer build instructions.
- Expected timing impact: none (docs/tools only).
- Verify: `scripts/stress_run.sh tests/buffer.sp 200` (optional); sanitizer build per README.
