#include "tanalisis.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <cmath>

#include "element.hpp"
#include "solver.hpp"
#include "dcanalysis.hpp"   // dcSolve
#include "utils.hpp"        // 如果需要一些工具函数的话

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ======= 小工具：从解向量中取节点电压 =======
static double getNodeVoltage(const Circuit& ckt,
                             const VectorXd& x,
                             int nodeId)
{
    int eq = ckt.nodes[nodeId].eqIndex;
    if (eq >= 0 && eq < x.size()) return x(eq);
    return 0.0;
}

// ======= 小工具: 增加微小导纳 =======
static void stampGlobalGmin(const Circuit& ckt,
                            MatrixXd& G,
                            double gmin)
{
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < G.rows()) {
            G(eq, eq) += gmin;
        }
    }
}

// ======= 对外：DC 工作点 =======
Eigen::VectorXd computeDcOperatingPoint(const Circuit& ckt)
{
    // 直接复用 DC 分析
    return dcSolve(ckt);
}

// ======= 主函数：后向欧拉瞬态 + 输出 CSV =======
void runTransientAnalysisBackwardEuler(const Circuit& ckt,
                                       const SimulationConfig& sim,
                                       const std::string& outFile)
{
    const TranConfig& cfg = sim.tran;

    if (!cfg.enabled) {
        std::cerr << "Transient analysis is not enabled in SimulationConfig.\n";
        return;
    }

    if (cfg.tstep <= 0.0 || cfg.tstop <= 0.0) {
        std::cerr << "Invalid .TRAN: tstep and tstop must be > 0.\n";
        return;
    }

    const double dt     = cfg.tstep;
    const double tstop  = cfg.tstop;
    const double tstart = cfg.tstart;   // 只影响输出起始时间

    const int N = ckt.numUnknowns();
    if (N <= 0) {
        std::cerr << "Transient analysis: circuit has no unknowns.\n";
        return;
    }

    // ===== 0. DC 工作点，作为 t=0 初值 =====
    VectorXd xdc;
    try {
        xdc = computeDcOperatingPoint(ckt);
    } catch (const std::exception& e) {
        std::cerr << "computeDcOperatingPoint failed: "
                  << e.what() << "\n";
        return;
    }

    if (xdc.size() != N) {
        std::cerr << "Transient: DC solution size mismatch: xdc.size()="
                  << xdc.size() << ", expected N=" << N << "\n";
        return;
    }

    // ===== 1. 收集所有 C / L 元件，建立状态 =====
    std::vector<std::shared_ptr<CapacitorElement>> caps;
    std::vector<std::shared_ptr<Inductor>>        inds;

    for (const auto& e : ckt.elements) {
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
            caps.push_back(c);
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            inds.push_back(L);
        }
    }

    // C: 记录 v_C^n；L: 记录 i_L^n
    std::unordered_map<const CapacitorElement*, double> capVprev;
    std::unordered_map<const Inductor*, double>         indIprev;

    // 用 DC 解初始化 t=0 状态
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, xdc, n1);
        double v2 = getNodeVoltage(ckt, xdc, n2);
        capVprev[c.get()] = v1 - v2;
    }
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        double iL0 = 0.0;
        if (k >= 0 && k < xdc.size()) {
            iL0 = xdc(k);
        }
        indIprev[L.get()] = iL0;
    }

    // ===== 2. 打开输出文件，写表头 =====
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Failed to open transient output file '"
                  << outFile << "'\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    ofs << "time";

    // 输出所有有方程号的节点电压（除了 GND）
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }

    // 输出所有电压源 / 电感的支路电流
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            ofs << ",I(" << vs->getName() << ")";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            ofs << ",I(" << L->getName() << ")";
        }
    }
    ofs << "\n";

    // 输出一行数据的 lambda（只输出 t >= tstart）
    auto dumpRow = [&](double t, const VectorXd& x) {
        if (t < tstart) return;

        ofs << t;

        // 节点电压
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) {
                ofs << "," << x(node.eqIndex);
            }
        }

        // 电压源 / 电感支路电流
        for (const auto& e : ckt.elements) {
            if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
                int k = vs->getBranchEqIndex();
                double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << Ibr;
            } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
                int k = L->getBranchEqIndex();
                double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << Ibr;
            }
        }

        ofs << "\n";
    };

    // ===== 3. 时间步循环：后向欧拉 + Newton =====

    // 打印一下配置，方便 debug
    std::cout << "[TRAN] tstep="  << std::scientific << cfg.tstep
              << ", tstop="       << cfg.tstop
              << ", tstart="      << cfg.tstart << "\n";

    // t=0 输出一行（如果 tstart <= 0）
    dumpRow(0.0, xdc);

    const int    maxNewtonIters = 100;
    const double tol            = 1e-6;
    const double gmin           = 1e-8;   // 比 DC 稍小一点就够

    int nSteps = static_cast<int>(std::floor(tstop / dt + 1e-9));
    std::cout << "[TRAN] total steps = " << nSteps << "\n";

    VectorXd x = xdc;     // 当前步解
    MatrixXd G(N, N);
    VectorXd I(N);

    for (int step = 0; step < nSteps; ++step) {
        double tNow = (step + 1) * dt;

        // 以上一个时间点的解作为初值（课件推荐）
        // x 此时已经是上一步的解

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 非 C / 非 L 元件 stamp（电阻/MOS/电压源/电流源等）
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e))        continue;
                e->stamp(G, I, ckt, x, ctx);
            }

            // 2) C 的后向欧拉 Norton 伴随模型
            // i_C^{n+1} = C (v^{n+1} - v^n)/dt
            // => 等效电导 g = C/dt，历史电流源 I_hist = -g * v^n
            for (auto& c : caps) {
                double Cval = c->getC();
                if (Cval <= 0.0) continue;

                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];

                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;

                double g = Cval / dt;

                // 电导部分
                if (eq1 >= 0) G(eq1, eq1) += g;
                if (eq2 >= 0) G(eq2, eq2) += g;
                if (eq1 >= 0 && eq2 >= 0) {
                    G(eq1, eq2) -= g;
                    G(eq2, eq1) -= g;
                }

                // 历史电流源 I_hist
                double vPrev = capVprev[c.get()];
                double I_hist = -g * vPrev;   // 从 n1 -> n2 的电流源

                if (eq1 >= 0) I(eq1) -= I_hist;
                if (eq2 >= 0) I(eq2) += I_hist;
            }

            // 3) L 的后向欧拉 Thevenin 伴随模型
            // v^{n+1} = L (i^{n+1} - i^n)/dt
            // => v^{n+1} - R_eq i^{n+1} = -R_eq i^n,  R_eq = L/dt
            for (auto& L : inds) {
                double Lval = L->getL();
                if (Lval <= 0.0) continue;

                const auto& nodes = L->getNodeIds();
                int np = nodes[0];
                int nm = nodes[1];

                int eqP = ckt.nodes[np].eqIndex;
                int eqM = ckt.nodes[nm].eqIndex;
                int k   = L->getBranchEqIndex();

                if (k < 0 || k >= N) {
                    std::cerr << "Transient: invalid branchEqIndex for inductor "
                              << L->getName() << "\n";
                    continue;
                }

                double R_eq  = Lval / dt;
                double iPrev = indIprev[L.get()];
                double V_hist = -R_eq * iPrev;

                // KCL：节点方程，对应 0V 源那两行（电流从 p 流向 m 为正）
                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;

                // 支路方程：Vp - Vm - R_eq * I_L = V_hist
                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
                G(k, k) += -R_eq;
                I(k)    += V_hist;
            }

            // 4) gmin：和 DC 一样给对角线加一点
            stampGlobalGmin(ckt, G, gmin);

            // 5) 解 G x_new = I
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);

            if (!xNew.allFinite()) {
                throw std::runtime_error(
                    "Transient Newton: linear solve produced NaN/Inf.");
            }

            // 6) 阻尼更新（仿照 dcSolveNewtonLU）
            const double alpha = 0.45;

            xNew = x + alpha * (xNew - x);
            double err = (xNew - x).norm();
            x = xNew;

            if (err < tol) {
                break;  // 本时间点收敛
            }

            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: transient Newton did not converge at t="
                          << std::scientific << tNow
                          << " (err=" << err << ")\n";
            }
        }

        // ===== 4) 步进成功：更新 C/L 状态，并输出一行 =====
        for (auto& c : caps) {
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            capVprev[c.get()] = v1 - v2;
        }

        for (auto& L : inds) {
            int k = L->getBranchEqIndex();
            double iL = 0.0;
            if (k >= 0 && k < x.size()) {
                iL = x(k);
            }
            indIprev[L.get()] = iL;
        }

        dumpRow(tNow, x);
    }

    std::cout << "Transient analysis (Backward Euler) finished. Results written to '"
              << outFile << "'.\n";
}
