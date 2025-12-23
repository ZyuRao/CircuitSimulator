#include <chrono>
#include <iostream>
#include <string>

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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [--verbose] <netlist.sp>\n";
        return 1;
    }

    bool verbose = false;
    std::string netlistPath;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") {
            verbose = true;
            continue;
        }
        if (!netlistPath.empty()) {
            std::cerr << "Unexpected extra argument: " << arg << "\n";
            return 1;
        }
        netlistPath = arg;
    }

    if (netlistPath.empty()) {
        std::cerr << "Usage: " << argv[0] << " [--verbose] <netlist.sp>\n";
        return 1;
    }

    Runner runner(std::string("out"), verbose);
    auto start = std::chrono::steady_clock::now();
    int rc = runner.run(netlistPath);
    auto end = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "[TIME] total_wall_ms=" << ms << "\n";
    return rc;
}
