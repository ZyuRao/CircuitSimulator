#include <iostream>
#include <string>

#include "runner.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <netlist.sp>\n";
        return 1;
    }

    std::string netlistPath = argv[1];
    Runner runner(std::string("out"));
    return runner.run(netlistPath);
}
