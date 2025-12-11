#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <cmath>

#include "circuit.hpp"
#include "parser.hpp"
#include "analysis.hpp"
#include "sim.hpp"
#include "utils.hpp"

using std::cerr;
using std::cout;
using std::endl;

// ====== 工具函数：把 HB 解 xHB 转成时域波形并写 CSV ======
static void exportHbSolutionToCsv(const Circuit&          ckt,
                                  const SimulationConfig& sim,
                                  const Eigen::VectorXd&  xHB,
                                  const std::string&      outFile)
{
    const int N = ckt.numUnknowns();
    if (N <= 0) {
        cerr << "[HB] export: circuit has no unknowns.\n";
        return;
    }

    if (!sim.hb.enabled) {
        cerr << "[HB] export: HB not enabled in SimulationConfig.\n";
        return;
    }

    const double f0 = sim.hb.f0;
    if (f0 <= 0.0) {
        cerr << "[HB] export: invalid f0 in sim.hb.\n";
        return;
    }

    const int K = sim.hb.nHarm;
    const int H = K + 1;
    const int expectedSize = 2 * H * N;
    if (xHB.size() != expectedSize) {
        cerr << "[HB] export: xHB size mismatch, expect " << expectedSize
             << ", got " << xHB.size() << "\n";
        return;
    }

    using Cplx = std::complex<double>;
    const double T      = 1.0 / f0;
    const double omega0 = 2.0 * M_PI * f0;

    // 布局：p = h * (2N) + offset
    // offset < N      -> 实部
    // offset in [N,2N)-> 虚部
    std::vector< std::vector<Cplx> > Vk(H, std::vector<Cplx>(N, Cplx(0.0, 0.0)));

    const int block = 2 * N;
    for (int p = 0; p < xHB.size(); ++p) {
        int h      = p / block;
        int offset = p % block;
        if (h < 0 || h >= H) continue;

        if (offset < N) {
            int idx = offset;
            Vk[h][idx].real(xHB[p]);
        } else {
            int idx = offset - N;
            Vk[h][idx].imag(xHB[p]);
        }
    }

    // 采样点数：和 HbAnalysis 里类似，足够细
    int minSamples   = std::max(64, 8 * (2 * K + 1));
    int nTimeSamples = 1;
    while (nTimeSamples < minSamples) nTimeSamples <<= 1;

    std::ofstream ofs(outFile);
    if (!ofs) {
        cerr << "[HB] export: cannot open output file '" << outFile << "'\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    // CSV 头：time, V(node1), V(node2), ...
    ofs << "time";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    ofs << "\n";

    // 逐个采样点计算 v_i(t) = Re( sum_{h=0..K} Vk[h][i] * e^{j h ω0 t} )
    for (int n = 0; n < nTimeSamples; ++n) {
        double t = (static_cast<double>(n) / nTimeSamples) * T;
        ofs << t;

        for (const auto& node : ckt.nodes) {
            if (node.eqIndex < 0) continue;
            int idx = node.eqIndex;
            if (idx < 0 || idx >= N) {
                ofs << ",0.0";
                continue;
            }

            Cplx v(0.0, 0.0);
            for (int h = 0; h < H; ++h) {
                double ang = h * omega0 * t;
                Cplx ej( std::cos(ang), std::sin(ang) );
                v += Vk[h][idx] * ej;
            }
            ofs << "," << v.real();
        }

        ofs << "\n";
    }

    cout << "[HB] Time-domain waveform written to '" << outFile << "'\n";
}

// ====== 运行瞬态 (后向欧拉)，输出 tran_out.csv ======
static void runTransient(const Circuit& ckt, const SimulationConfig& sim)
{
    if (!sim.tran.enabled) {
        cerr << "[TRAN] .TRAN not enabled in netlist, skip transient.\n";
        return;
    }
    cout << "[TRAN] Running transient (Backward Euler) -> tran_out.csv\n";
    TransientAnalysis tran(ckt, sim, "tran_out.csv");
    tran.runBackwardEuler();
}

// ====== 运行 HB + 导出时域波形 hb_out.csv ======
static void runHb(const Circuit& ckt, const SimulationConfig& sim)
{
    if (!sim.hb.enabled) {
        cerr << "[HB] .hb not enabled in netlist, skip HB.\n";
        return;
    }

    cout << "[HB] Solving DC operating point for HB initial guess...\n";
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
    Eigen::VectorXd xdc = dc.run();

    cout << "[HB] Running harmonic balance solver...\n";
    HbAnalysis hb(ckt, sim, xdc);
    Eigen::VectorXd xHB;
    bool ok = hb.run(xHB);
    if (!ok) {
        cerr << "[HB] Harmonic balance did not converge, no CSV output.\n";
        return;
    }

    exportHbSolutionToCsv(ckt, sim, xHB, "hb_out.csv");
}

// ====== 运行 Shooting PSS，输出 shoot_out.csv ======
static void runShooting(const Circuit& ckt, const SimulationConfig& sim)
{
    if (!sim.hb.enabled) {
        cerr << "[PSS] .hb card required to specify f0 for Shooting (use: .hb f0 nHarm).\n";
        return;
    }

    ShootingPssConfig cfg;
    cfg.periodT = 1.0 / sim.hb.f0;      // 周期由 .hb 的 f0 决定
    cfg.tstep   = 0.0;                   // 让 runPssShootingAnalysis 内部决定 dt
    cfg.maxIters = 50;
    cfg.tol      = 1e-6;
    cfg.relax    = 0.5;

    cout << "[PSS] Running Shooting method PSS -> shoot_out.csv\n";
    runPssShootingAnalysis(ckt, sim, cfg, "shoot_out.csv");
}

// ====== main 函数 ======
int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage:\n"
             << "  " << argv[0] << " <netlist> [mode]\n"
             << "    mode = auto  (default, 根据 .TRAN / .hb 自动判断)\n"
             << "         = tran  (只跑瞬态，输出 tran_out.csv)\n"
             << "         = hb    (只跑 HB，输出 hb_out.csv)\n"
             << "         = shoot (只跑 Shooting PSS，输出 shoot_out.csv)\n";
        return 1;
    }

    std::string netlistFile = argv[1];
    std::string mode        = (argc >= 3) ? argv[2] : "auto";

    Circuit ckt;
    SimulationConfig sim;

    NetlistParser parser(ckt, sim);
    if (!parser.parseFile(netlistFile)) {
        cerr << "Failed to parse netlist: " << netlistFile << "\n";
        return 1;
    }

    ckt.assignEquationIndices();

    if (mode == "tran") {
        runTransient(ckt, sim);
    } else if (mode == "hb") {
        runHb(ckt, sim);
    } else if (mode == "shoot") {
        runShooting(ckt, sim);
    } else if (mode == "auto") {
        if (sim.tran.enabled) {
            runTransient(ckt, sim);
        }
        if (sim.hb.enabled) {
            runShooting(ckt, sim);  // 这里默认优先测试 Shooting
            // 如果你也想自动跑 HB，可以再加一行：
            // runHb(ckt, sim);
        }
    } else {
        cerr << "Unknown mode: " << mode << "\n";
        return 1;
    }

    return 0;
}
