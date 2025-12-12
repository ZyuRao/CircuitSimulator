#include <iostream>
#include <iomanip>
#include <string>
#include <memory>

#include "analysis.hpp"
#include "parser.hpp"
#include "circuit.hpp"
#include "sim.hpp"
#include "element.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: dc_test <netlist.sp>\n";
        return 1;
    }

    std::string netlistFile = argv[1];

    Circuit ckt;
    SimulationConfig sim;

    std::cout << "Reading netlist: " << netlistFile << "\n";

    if (!parseNetlist(netlistFile, ckt, sim)) {
        std::cerr << "ERROR: parseNetlist() failed.\n";
        return 1;
    }

    // 分配方程号
    ckt.assignEquationIndices();

    std::cout << "\n==== Circuit summary ====\n";
    std::cout << "Nodes        : " << ckt.nodes.size()    << "\n";
    std::cout << "Elements     : " << ckt.elements.size() << "\n";
    std::cout << "Unknowns     : " << ckt.numUnknowns()
              << "  (nodeEq=" << ckt.numNodeEquations()
              << ", branchEq=" << ckt.numVoltageBranches() << ")\n";

    // ========== DC 工作点 ==========
    std::cout << "\nRunning DC operating point (Gauss-Seidel / Newton)...\n";

    Eigen::VectorXd xdc;
    try {
        DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
        xdc = dc.run();
    } catch (const std::exception& e) {
        std::cerr << "DC solve failed: " << e.what() << "\n";
        return 1;
    }

    if (xdc.size() != ckt.numUnknowns()) {
        std::cerr << "DC solution size mismatch: xdc.size() = "
                  << xdc.size()
                  << ", numUnknowns = " << ckt.numUnknowns() << "\n";
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
   for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            double v = xdc(node.eqIndex);
            std::cout << "V(" << node.name << ") = " << v << " V"
                      << "   [eqIndex=" << node.eqIndex << "]\n";
        } else {
            // 地节点
            std::cout << "V(" << node.name << ") = 0.000000 V   [GND]\n";
        }
    }

    // 顺便把电压源 / 电感的支路电流也打印出来，方便 debug
    std::cout << "\n==== DC branch currents (voltage sources / inductors) ====\n";

    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            int k = vs->getBranchEqIndex();
            double I = (k >= 0 && k < xdc.size()) ? xdc(k) : 0.0;
            // 约定：电流方向是从正端节点流向负端（和 stamp 的定义一致）
            std::cout << "I(" << vs->getName()
                      << ", from +" << ckt.nodes[vs->getNodeIds()[0]].name
                      << " to -"   << ckt.nodes[vs->getNodeIds()[1]].name
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
    return 0;
}
