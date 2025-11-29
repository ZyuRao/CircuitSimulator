// main.cpp
#include <iostream>
#include <string>

#include "parser.hpp"
#include "circuit.hpp"
#include "dcanalysis.hpp"

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " netlist.sp\n";
        return 1;
    }

    std::string netlistFile = argv[1];

    Circuit ckt;
    SimulationConfig sim;

    // 1) 解析网表，填充 Circuit 和 SimulationConfig
    if (!parseNetlist(netlistFile, ckt, sim)) {  // parser.hpp 里的 inline 函数
        std::cerr << "Parse netlist failed.\n";
        return 1;
    }

    // 2) 分配 MNA 方程编号：
    //    - 每个非 GND 节点一个电压未知数
    //    - 每个电压源 / 电感一个支路电流未知数
    ckt.assignEquationIndices();
    int N = ckt.numUnknowns();

    std::cout << "DC analysis on netlist: " << netlistFile << "\n";
    std::cout << "Total unknowns (node voltages + branch currents): " << N << "\n";

    // 3) 做一次 DC 工作点求解
    //    - 有 MOS: 调用 dcSolveNewton()
    //    - 无 MOS: 调用 dcSolveDirect()
    Eigen::VectorXd x = dcSolve(ckt);

    if (x.size() == 0) {
        std::cerr << "No unknowns, nothing to solve.\n";
        return 0;
    }

    // 4) 打印节点电压
    std::cout << "\n=== Node voltages (DC operating point) ===\n";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex < 0) {
            // 约定：eqIndex == -1 表示地节点
            std::cout << "V(" << node.name << ") = 0.0 V\n";
        } else if (node.eqIndex < x.size()) {
            std::cout << "V(" << node.name << ") = " << x(node.eqIndex) << " V\n";
        }
    }

    // 5) 打印电压源 / 电感的支路电流
    //    MNA 里这些元件各自有一个额外未知量，对应 eqIndex = 节点数 + k
    std::cout << "\n=== Branch currents of voltage sources / inductors ===\n";
    for (const auto& e : ckt.elements) {
        // 电压源
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            int k = vs->getBranchEqIndex();
            if (k >= 0 && k < x.size()) {
                std::cout << "I(" << vs->getName()
                          << ") = " << x(k) << " A  (current from + to - terminal)\n";
            }
        }

        // 电感：DC 中等价为 0V 电压源 + 支路电流未知数
        if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            int k = L->getBranchEqIndex();
            if (k >= 0 && k < x.size()) {
                std::cout << "I(" << L->getName()
                          << ") = " << x(k) << " A  (DC current through inductor)\n";
            }
        }
    }

    return 0;
}
