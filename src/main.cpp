// main.cpp -- DC 工作点 + 瞬态仿真（后向欧拉）

#include <iostream>
#include <iomanip>
#include <string>
#include <memory>

#include "parser.hpp"
#include "circuit.hpp"
#include "dcanalysis.hpp"
#include "sim.hpp"
#include "element.hpp"
#include "tanalisis.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: mysim.exe <netlist.sp> [tran_out.csv]\n";
        return 1;
    }

    std::string netlistFile = argv[1];
    std::string tranOutFile = (argc >= 3) ? argv[2] : "tran_out.csv";

    Circuit ckt;
    SimulationConfig sim;

    std::cout << "Reading netlist: " << netlistFile << "\n";

    if (!parseNetlist(netlistFile, ckt, sim)) {
        std::cerr << "parseNetlist() failed.\n";
        return 1;
    }

    ckt.assignEquationIndices();

    std::cout << "\n==== Circuit summary ====\n";
    std::cout << "Node count   : " << ckt.nodes.size()    << "\n";
    std::cout << "Element count: " << ckt.elements.size() << "\n";
    std::cout << "Unknowns     : " << ckt.numUnknowns()
              << "  (nodeEq=" << ckt.numNodeEquations()
              << ", branchEq=" << ckt.numVoltageBranches() << ")\n";

    // ===== DC 工作点 =====
    std::cout << "\nRunning DC operating point...\n";

    Eigen::VectorXd xdc;
    try {
        xdc = computeDcOperatingPoint(ckt);
    } catch (const std::exception& e) {
        std::cerr << "DC solve failed: " << e.what() << "\n";
        return 1;
    }

    if (xdc.size() != ckt.numUnknowns()) {
        std::cerr << "DC solution size mismatch.\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n==== DC node voltages ====\n";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            double v = xdc(node.eqIndex);
            std::cout << "V(" << node.name << ") = " << v << " V"
                      << "   [eqIndex=" << node.eqIndex << "]\n";
        } else {
            std::cout << "V(" << node.name << ") = 0.000000 V   [GND]\n";
        }
    }

    std::cout << "\n==== DC branch currents (voltage sources / inductors) ====\n";
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            int k = vs->getBranchEqIndex();
            double I = (k >= 0 && k < xdc.size()) ? xdc(k) : 0.0;
            std::cout << "I(" << vs->getName()
                      << ", +" << ckt.nodes[vs->getNodeIds()[0]].name
                      << " -> -" << ckt.nodes[vs->getNodeIds()[1]].name
                      << ") = " << I << " A"
                      << "   [branchEq=" << k << "]\n";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            int k = L->getBranchEqIndex();
            double I = (k >= 0 && k < xdc.size()) ? xdc(k) : 0.0;
            std::cout << "I(" << L->getName()
                      << ", " << ckt.nodes[L->getNodeIds()[0]].name
                      << " -> "  << ckt.nodes[L->getNodeIds()[1]].name
                      << ") = " << I << " A"
                      << "   [branchEq=" << k << "]\n";
        }
    }

    std::cout << "\nDC analysis finished.\n";

    // ===== 瞬态仿真（如果 .TRAN 启用） =====
    if (sim.tran.enabled) {
        std::cout << "\nRunning transient analysis (Backward Euler)...\n";
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  .TRAN: tstep=" << sim.tran.tstep
                  << ", tstop=" << sim.tran.tstop
                  << ", tstart=" << sim.tran.tstart << "\n";
        std::cout << "  output file: " << tranOutFile << "\n";

        try {
            runTransientAnalysisBackwardEuler(ckt, sim, tranOutFile);
        } catch (const std::exception& e) {
            std::cerr << "Transient failed: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cout << "\nNo .TRAN card; transient analysis skipped.\n";
    }

    return 0;
}
