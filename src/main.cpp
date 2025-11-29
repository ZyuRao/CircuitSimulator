// 咱就简单测试试试
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "parser.hpp"
#include "circuit.hpp"
#include "sim.hpp"
#include "dcanalysis.hpp"

int main(int argc, char** argv) {
    // 默认用 buffer.sp，也可以在命令行传其他 netlist
    std::string netlistFile = "buffer.sp";
    if (argc >= 2) {
        netlistFile = argv[1];
    }

    Circuit ckt;
    SimulationConfig sim;

    // 解析 netlist
    if (!parseNetlist(netlistFile, ckt, sim)) {
        std::cerr << "Failed to parse netlist: " << netlistFile << std::endl;
        return 1;
    }

    std::cout << "Parsed netlist: " << netlistFile << std::endl;
    std::cout << "Number of unknowns: " << ckt.numUnknowns() << std::endl;

    // 做一次 DC 工作点求解（默认用你改好的 Gauss-Seidel 版本）
    Eigen::VectorXd x = dcSolve(ckt);

    // 打印每个节点的电压
    std::cout << "\nDC operating point (node voltages):\n";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex < 0) {
            // eqIndex < 0 代表接地节点
            std::cout << "Node " << node.name
                      << " (id=" << node.id << "): 0 (ground)\n";
        } else {
            double v = x(node.eqIndex);
            std::cout << "Node " << node.name
                      << " (id=" << node.id
                      << ", eq=" << node.eqIndex
                      << "): " << v << " V\n";
        }
    }

    return 0;
}
