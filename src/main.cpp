// 以下测试所用，证明了网表读取无误
// 但是DC分析仍不收敛

#include <iostream>
#include <string>
#include <memory>

#include "parser.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "sim.hpp"

// 小工具：把枚举转成字符串，方便打印
static const char* waveformTypeName(WaveformType t) {
    switch (t) {
        case WaveformType::NONE:  return "NONE";
        case WaveformType::PULSE: return "PULSE";
        case WaveformType::SIN:   return "SIN";
        case WaveformType::PWL:   return "PWL";
    }
    return "UNKNOWN";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: mysim.exe <netlist.sp>\n";
        return 1;
    }

    std::string netlistFile = argv[1];

    Circuit ckt;
    SimulationConfig sim;

    std::cout << "Reading netlist: " << netlistFile << "\n";

    if (!parseNetlist(netlistFile, ckt, sim)) {
        std::cerr << "parseNetlist() returned false\n";
        return 1;
    }

    std::cout << "\n==== Circuit summary ====\n";
    std::cout << "Node count   : " << ckt.nodes.size()    << "\n";
    std::cout << "Element count: " << ckt.elements.size() << "\n";
    std::cout << "MOS model cnt: " << ckt.mosModels.size() << "\n";

    // 给节点 & 电压源/电感分配方程号（便于后面打印 eqIndex）
    ckt.assignEquationIndices();
    std::cout << "numUnknowns  : " << ckt.numUnknowns()
              << "  (nodeEq=" << ckt.numNodeEquations()
              << ", branchEq=" << ckt.numVoltageBranches() << ")\n";

    // 打印 MOS 模型表
    if (!ckt.mosModels.empty()) {
        std::cout << "\n==== MOS models ====\n";
        for (const auto& kv : ckt.mosModels) {
            const auto& m = kv.second;
            std::cout << "ModelId=" << kv.first
                      << "  name=" << m.name
                      << "  type=" << (m.isP ? "PMOS" : "NMOS")
                      << "  VT=" << m.VT
                      << "  MU=" << m.MU
                      << "  COX=" << m.COX
                      << "  LAMBDA=" << m.LAMBDA
                      << "  CJO=" << m.CJO
                      << "\n";
        }
    }

    // 打印节点列表和连接的元件
    std::cout << "\n==== Nodes ====\n";
    for (const auto& node : ckt.nodes) {
        std::cout << "Node id=" << node.id
                  << "  name=" << node.name
                  << "  eqIndex=" << node.eqIndex
                  << "  attached=[";
        for (std::size_t k = 0; k < node.attachedElements.size(); ++k) {
            int ei = node.attachedElements[k];
            if (ei >= 0 && ei < (int)ckt.elements.size()) {
                std::cout << ckt.elements[ei]->getName();
            } else {
                std::cout << "?" << ei;
            }
            if (k + 1 < node.attachedElements.size()) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    // 打印元件详细信息
    std::cout << "\n==== Elements ====\n";
    for (std::size_t i = 0; i < ckt.elements.size(); ++i) {
        auto e = ckt.elements[i];
        std::cout << "[" << i << "] " << e->getName() << "  type=";

        // 判断类型
        if (std::dynamic_pointer_cast<Resistor>(e)) {
            std::cout << "R";
        } else if (std::dynamic_pointer_cast<CapacitorElement>(e)) {
            std::cout << "C";
        } else if (std::dynamic_pointer_cast<Inductor>(e)) {
            std::cout << "L";
        } else if (std::dynamic_pointer_cast<VoltageSource>(e)) {
            std::cout << "V";
        } else if (std::dynamic_pointer_cast<CurrentSource>(e)) {
            std::cout << "I";
        } else if (std::dynamic_pointer_cast<NMosElement>(e)) {
            std::cout << "NMOS";
        } else if (std::dynamic_pointer_cast<PMosElement>(e)) {
            std::cout << "PMOS";
        } else {
            std::cout << "Unknown";
        }

        // 连接的节点名
        std::cout << "  nodes=[";
        const auto& nodeIds = e->getNodeIds();
        for (std::size_t k = 0; k < nodeIds.size(); ++k) {
            int nid = nodeIds[k];
            if (nid >= 0 && nid < (int)ckt.nodes.size()) {
                std::cout << ckt.nodes[nid].name;
            } else {
                std::cout << "?" << nid;
            }
            if (k + 1 < nodeIds.size()) std::cout << ", ";
        }
        std::cout << "]";

        // 如果是电压源，把 SourceSpec 打印出来
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            const SourceSpec& spec = vs->getSpec();
            std::cout << "  [Vsrc: dc=" << spec.dcValue
                      << ", acMag=" << spec.acMag
                      << ", acPhase=" << spec.acPhaseDeg
                      << "deg, tranType=" << waveformTypeName(spec.tran.type);

            if (spec.tran.type == WaveformType::SIN) {
                const SinSpec& s = spec.tran.sine;
                std::cout << ", SIN(v0=" << s.v0
                          << ", va=" << s.va
                          << ", freq=" << s.freq
                          << ", td=" << s.td
                          << ", phi=" << s.phi << ")";
            }
            std::cout << "]";
        }

        // 如果是电流源，同样打印 SourceSpec
        if (auto isrc = std::dynamic_pointer_cast<CurrentSource>(e)) {
            const SourceSpec& spec = isrc->setSpec(); // 注意：接口名叫 setSpec
            std::cout << "  [Isrc: dc=" << spec.dcValue
                      << ", acMag=" << spec.acMag
                      << ", acPhase=" << spec.acPhaseDeg
                      << "deg, tranType=" << waveformTypeName(spec.tran.type)
                      << "]";
        }

        std::cout << "\n";
    }

    // 打印仿真配置（.OP / .DC / .TRAN / .AC）
    std::cout << "\n==== SimulationConfig ====\n";
    std::cout << "doOp      : " << (sim.doOp ? "true" : "false") << "\n";
    std::cout << "DC sweeps : " << sim.dcSweeps.size() << "\n";
    for (std::size_t i = 0; i < sim.dcSweeps.size(); ++i) {
        const auto& dc = sim.dcSweeps[i];
        std::cout << "  [" << i << "] source=" << dc.sourceName
                  << "  start=" << dc.start
                  << "  stop=" << dc.stop
                  << "  step=" << dc.step << "\n";
    }

    std::cout << "TRAN      : enabled=" << (sim.tran.enabled ? "true" : "false")
              << "  tstart=" << sim.tran.tstart
              << "  tstop=" << sim.tran.tstop
              << "  tstep=" << sim.tran.tstep << "\n";

    std::cout << "AC        : enabled=" << (sim.ac.enabled ? "true" : "false")
              << "  nPoints=" << sim.ac.nPoints
              << "  fstart=" << sim.ac.fstart
              << "  fstop=" << sim.ac.fstop << "\n";

    std::cout << "\nDone.\n";
    return 0;
}
