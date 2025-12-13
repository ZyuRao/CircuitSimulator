// hb_main.cpp
#include <iostream>
#include <iomanip>
#include <string>

#include "parser.hpp"
#include "circuit.hpp"
#include "analysis.hpp"
#include "sim.hpp"
#include "element.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: CircuitSimulator <netlist.sp> [hb_time.csv]\n";
        return 1;
    }

    const std::string netlistFile = argv[1];
    const std::string outCsv      = (argc >= 3) ? argv[2] : "hb_time.csv";

    Circuit ckt;
    SimulationConfig sim;

    std::cout << "Reading netlist: " << netlistFile << "\n";
    if (!parseNetlist(netlistFile, ckt, sim)) {
        std::cerr << "ERROR: parseNetlist() failed.\n";
        return 2;
    }

    // 分配方程号（节点电压 + 电压源/电感支路电流）
    ckt.assignEquationIndices();

    std::cout << "\n==== Circuit summary ====\n";
    std::cout << "Nodes        : " << ckt.nodes.size()    << "\n";
    std::cout << "Elements     : " << ckt.elements.size() << "\n";
    std::cout << "Unknowns     : " << ckt.numUnknowns()
              << "  (nodeEq="   << ckt.numNodeEquations()
              << ", branchEq=" << ckt.numVoltageBranches() << ")\n";

    // ---------- DC operating point ----------
    std::cout << "\nRunning DC operating point...\n";
    Eigen::VectorXd xdc;
    try {
        DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel); // 你也可以改成 LU
        xdc = dc.run();
    } catch (const std::exception& e) {
        std::cerr << "DC solve failed: " << e.what() << "\n";
        return 3;
    }

    if (xdc.size() != ckt.numUnknowns()) {
        std::cerr << "DC solution size mismatch: xdc.size()=" << xdc.size()
                  << ", numUnknowns=" << ckt.numUnknowns() << "\n";
        return 4;
    }

    // 可选：打印 DC 节点电压
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n==== DC node voltages ====\n";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            std::cout << "V(" << node.name << ") = " << xdc(node.eqIndex)
                      << " V   [eqIndex=" << node.eqIndex << "]\n";
        } else {
            std::cout << "V(" << node.name << ") = 0.000000 V   [GND]\n";
        }
    }

    // ---------- HB ----------
    if (!sim.hb.enabled) {
        std::cerr << "\nNo .hb card found. Skip HB.\n";
        return 0;
    }

    std::cout << "\nRunning HB: f0=" << sim.hb.f0
              << " Hz, nHarm=" << sim.hb.nHarm
              << ", out=" << outCsv << "\n";

    Eigen::VectorXd xhb;
    try {
        HbAnalysis hb(ckt, sim, xdc);
        bool ok = hb.run(xhb, outCsv);
        if (!ok) {
            std::cerr << "HB failed to converge.\n";
            return 5;
        }
    } catch (const std::exception& e) {
        std::cerr << "HB exception: " << e.what() << "\n";
        return 6;
    }

    std::cout << "HB done. CSV written: " << outCsv << "\n";
    std::cout << "Plot: python3 plot_tran.py " << outCsv << "\n";
    return 0;
}
