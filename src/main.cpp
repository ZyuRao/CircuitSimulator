// main.cpp
// Entry point for CircuitSimulator using Runner (auto/manual + timing handled in runner.cpp)
//
// Usage (recommended):
//   sim <netlist.sp> [op|tran|hb|all] [options]
//
// Examples:
//   sim buffer.sp
//   sim buffer.sp tran --tran-method=be --dc-solver=lu
//   sim buffer.sp --run=op,hb --dc-solver=gs
//
// Output directory:
//   Default: "out"
//   You may override via environment variable: SIM_OUTDIR=out_test
//
// Note:
//   - If runner.hpp provides Runner::run(argc, argv), CLI options work.
//   - If not, this main will fall back to Runner::run(netlistPath) (you can still use env vars supported by runner.cpp).

#include "runner.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

// Compile-time friendly dispatch: prefer run(argc, argv) if available; otherwise fallback to run(netlistPath).
template <typename T>
auto runRunner(T& r, int argc, char** argv, int) -> decltype(r.run(argc, argv)) {
    return r.run(argc, argv);
}
template <typename T>
int runRunner(T& r, int argc, char** argv, long) {
    // fallback: only netlistPath; extra CLI options are ignored
    return r.run(std::string(argv[1]));
}

int main(int argc, char** argv) {
    if (argc < 2 || argv == nullptr) {
        std::cerr << "Usage:\n"
                  << "  sim <netlist.sp> [op|tran|hb|all] [options]\n"
                  << "Try: sim --help\n";
        return 1;
    }

    // Optional: output directory via env var (doesn't affect Runner::run argument parsing)
    std::string outDir = "out";
    if (const char* v = std::getenv("SIM_OUTDIR")) {
        if (*v != '\0') outDir = v;
    }

    Runner runner(outDir);
    return runRunner(runner, argc, argv, 0);
}
