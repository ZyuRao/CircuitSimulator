#include "analysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "sim.hpp"
#include "solver.hpp"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <functional>
#include <limits>
#include <algorithm>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ================== PSS 专用的 BE 伴随模型（目前主流程不再用，仅为兼容接口） ==================

// 电容后向欧拉伴随模型
// C 接在 eq1 与 eq2 之间，vPrev = V(eq1)^n - V(eq2)^n
void pssStampCapBE(int eq1, int eq2, double C, double dt, double vPrev,
                   MatrixXd& G, VectorXd& I)
{
    if (C <= 0.0 || dt <= 0.0) return;

    double Gc = C / dt;

    if (eq1 >= 0) G(eq1, eq1) += Gc;
    if (eq2 >= 0) G(eq2, eq2) += Gc;
    if (eq1 >= 0 && eq2 >= 0) {
        G(eq1, eq2) -= Gc;
        G(eq2, eq1) -= Gc;
    }

    // 历史电流源：I_hist 从 eq1 -> eq2
    double I_hist = -Gc * vPrev;
    if (eq1 >= 0) I(eq1) -= I_hist;
    if (eq2 >= 0) I(eq2) += I_hist;
}

// 电感后向欧拉 Thévenin 伴随模型
// eqP, eqM 为电感两端节点方程，k 为支路电流方程索引
void pssStampIndBE(int eqP, int eqM, int k, double Lval, double dt, double iPrev,
                   MatrixXd& G, VectorXd& I)
{
    if (Lval <= 0.0 || dt <= 0.0) return;
    if (k < 0 || k >= G.rows()) return;

    // V = L di/dt ≈ L (i_{n+1} - i_n)/dt
    // 等效为：R_eq = L/dt，历史电压源 V_hist = -R_eq * iPrev
    double R_eq   = Lval / dt;
    double V_hist = -R_eq * iPrev;

    // KCL 方程：I_L 从 p -> m
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;

    // 支路方程：V(p) - V(m) - R_eq * I_L = V_hist
    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;
    G(k, k) += -R_eq;
    I(k)    += V_hist;
}

// 全局 gmin 导纳矩阵（非模板版本，供 DC / TRAN / SHOOT 使用）
void stampGlobalGmin(const Circuit& ckt, MatrixXd& G, double gmin)
{
    if (gmin <= 0.0) return;
    int N = static_cast<int>(G.rows());
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < N) {
            G(eq, eq) += gmin;
        }
    }
}

// ================== 调用瞬态里的 integrateOnePeriodBE ==================

VectorXd integrateOnePeriodPssBE(const Circuit& ckt,
                                 const SimulationConfig& sim,
                                 double dt,
                                 double T,
                                 const VectorXd& x0,
                                 const std::function<void(double, const VectorXd&)>& dumpRow)
{
    // 这里的 outFile 对 integrateOnePeriodBE 没有实际作用，只是构造函数需要
    TransientAnalysis tran(ckt, sim, "shoot_dummy.csv");
    return tran.integrateOnePeriodBE(x0, 0.0, dt, T, dumpRow);
}

// ================== Shooting 外层迭代 ==================

void runPssShootingAnalysis(const Circuit& ckt,
                            const SimulationConfig& sim,
                            const ShootingPssConfig& cfg,
                            const std::string& outFile)
{
    if (cfg.periodT <= 0.0) {
        std::cerr << "PSS-SHOOT: periodT must be > 0\n";
        return;
    }

    // dt 优先用 cfg.tstep，默认 T/200
    double dt = cfg.tstep;
    if (dt <= 0.0) {
        dt = cfg.periodT / 200.0;
    }

    const int N = ckt.numUnknowns();
    if (N <= 0) {
        std::cerr << "PSS-SHOOT: circuit has no unknowns\n";
        return;
    }

    std::cout << "[PSS-SHOOT] T = " << cfg.periodT
              << ", dt = " << dt
              << ", maxIters = " << cfg.maxIters
              << ", tol = " << cfg.tol
              << ", relax = " << cfg.relax << "\n";

    // 1. 先求 DC 工作点，作为 Shooting 的初始猜测 x0
    std::cout << "[PSS-SHOOT] Solving DC operating point as initial guess...\n";
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
    VectorXd x0 = dc.run();
    if (x0.size() != N) {
        x0.conservativeResize(N);
    }

    // 2. Shooting 外层迭代：H(x0) = x(T; x0) - x0 = 0
    VectorXd xInit = x0;
    for (int k = 0; k < cfg.maxIters; ++k) {
        std::cout << "[PSS-SHOOT] Iteration " << (k + 1)
                  << " / " << cfg.maxIters << "...\n";

        // 用瞬态 BE 逻辑从 xInit 穷举一个周期
        VectorXd xT = integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, nullptr);

        VectorXd residual = xT - xInit;
        double error = residual.norm();

        std::cout << "             ||x(T) - x(0)|| = " << error << "\n";

        // 收敛判断
        if (error < cfg.tol) {
            std::cout << "[PSS-SHOOT] Converged after "
                      << (k + 1) << " iterations.\n";
            xInit = xT; // 收敛时把末态作为最终周期解的“起点”
            break;
        }

        // 简单松弛迭代：x0_{k+1} = x0_k + relax * (xT - x0_k)
        xInit += cfg.relax * residual;

        if (k == cfg.maxIters - 1) {
            std::cerr << "[PSS-SHOOT] WARNING: not converged within maxIters; "
                      << "using last iterate for waveform output.\n";
        }
    }

    // 3. 用收敛（或最后一次）xInit 再积分一个周期，写出稳态波形到 CSV
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "PSS-SHOOT: cannot open output file '" << outFile << "'\n";
        return;
    }

    ofs << "time";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    ofs << "\n";

    auto dumpRow = [&](double t, const VectorXd& x) {
        ofs << t;
        for (const auto& node : ckt.nodes) {
            int eq = node.eqIndex;
            if (eq >= 0 && eq < x.size()) {
                ofs << "," << x(eq);
            }
        }
        ofs << "\n";
    };

    integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, dumpRow);

    std::cout << "[PSS-SHOOT] Steady-state waveform written to '"
              << outFile << "'.\n";
}
