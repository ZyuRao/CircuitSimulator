
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>

#include "circuit.hpp"
#include "element.hpp"
#include "sim.hpp"
#include "solver.hpp"
#include "dcanalysis.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ========= 小工具：取节点电压 =========
static double getNodeVoltage(const Circuit& ckt,
                             const VectorXd& x,
                             int nodeId)
{
    int eq = ckt.nodes[nodeId].eqIndex;
    if (eq >= 0 && eq < x.size()) return x(eq);
    return 0.0;
}

// ========= 小工具：全局 gmin-to-ground =========
static void stampGlobalGmin(const Circuit& ckt,
                            MatrixXd& G,
                            double gmin)
{
    int N = G.rows();
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < N) {
            G(eq, eq) += gmin;
        }
    }
}

// ========= 对外：求 DC 工作点 =========
VectorXd computeDcOperatingPoint(const Circuit& ckt)
{
    return dcSolve(ckt);
}

// MOS 寄生电容状态：上一时刻的电压差
struct MosCapState {
    double vgsPrev = 0.0;
    double vgdPrev = 0.0;
    double vsbPrev = 0.0;
    double vdbPrev = 0.0;
};

// 在 G / I 中 stamp 一个“电容 + 历史电流源”
// C 接在 eq1 与 eq2 之间，vPrev = V(eq1)^n - V(eq2)^n
static void stampCapBE(int eq1, int eq2,
                       double C, double dt,
                       double vPrev,
                       MatrixXd& G,
                       VectorXd& I)
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

// ========= 主函数：后向欧拉瞬态 + MOS 寄生电容 =========
void runTransientAnalysisBackwardEuler(const Circuit& ckt,
                                       const SimulationConfig& sim,
                                       const std::string& outFile)
{
    const TranConfig& cfg = sim.tran;

    if (!cfg.enabled) {
        std::cerr << "Transient analysis is not enabled (.TRAN missing).\n";
        return;
    }

    if (cfg.tstep <= 0.0 || cfg.tstop <= 0.0) {
        std::cerr << "Invalid .TRAN card: tstep and tstop must be > 0.\n";
        return;
    }

    const double dt     = cfg.tstep;
    const double tstop  = cfg.tstop;
    const double tstart = cfg.tstart;

    const int N = ckt.numUnknowns();
    if (N <= 0) {
        std::cerr << "Transient: circuit has no unknowns.\n";
        return;
    }

    // ===== 0. 先求 DC 工作点，作为 t=0 初值 =====
    VectorXd xdc;
    try {
        xdc = computeDcOperatingPoint(ckt);
    } catch (const std::exception& e) {
        std::cerr << "DC operating point failed: " << e.what() << "\n";
        return;
    }

    if (xdc.size() != N) {
        std::cerr << "Transient: DC solution size mismatch.\n";
        return;
    }

    // ===== 1. 收集 C / L / MOS 元件，并建立初始状态 =====
    std::vector<std::shared_ptr<CapacitorElement>> caps;
    std::vector<std::shared_ptr<Inductor>>        inds;
    std::vector<std::shared_ptr<MosfetBase>>      mosfets;

    for (const auto& e : ckt.elements) {
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
            caps.push_back(c);
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            inds.push_back(L);
        } else if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) {
            mosfets.push_back(m);
        }
    }

    // 显式电容：v^{0} = V(n1) - V(n2)
    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, xdc, n1);
        double v2 = getNodeVoltage(ckt, xdc, n2);
        capVprev[c.get()] = v1 - v2;
    }

    // 电感：i_L^{0} = DC 解中的支路电流
    std::unordered_map<const Inductor*, double> indIprev;
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        double i0 = 0.0;
        if (k >= 0 && k < xdc.size()) {
            i0 = xdc(k);
        }
        indIprev[L.get()] = i0;
    }

    // MOS 寄生电容状态：从 DC 工作点初始化
    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    for (auto& m : mosfets) {
        const auto& nodes = m->getNodeIds();
        int nD = nodes[0];
        int nG = nodes[1];
        int nS = nodes[2];
        int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

        double vD = getNodeVoltage(ckt, xdc, nD);
        double vG = getNodeVoltage(ckt, xdc, nG);
        double vS = getNodeVoltage(ckt, xdc, nS);
        double vB = getNodeVoltage(ckt, xdc, nB);

        MosCapState st;
        st.vgsPrev = vG - vS;
        st.vgdPrev = vG - vD;
        st.vsbPrev = vS - vB;
        st.vdbPrev = vD - vB;
        mosPrev[m.get()] = st;
    }

    // ===== 2. 打开输出文件，写表头 =====
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Cannot open transient output file '" << outFile << "'.\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    ofs << "time";
    // 节点电压
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    // 电压源与电感电流
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            ofs << ",I(" << vs->getName() << ")";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            ofs << ",I(" << L->getName() << ")";
        }
    }
    ofs << "\n";

    auto dumpRow = [&](double t, const VectorXd& x) {
        if (t < tstart) return;

        ofs << t;
        // 节点电压
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) {
                ofs << "," << x(node.eqIndex);
            }
        }
        // 电压源与电感电流
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

    // ===== 3. 时间步循环：Newton + LU + 后向欧拉 =====
    std::cout << "[TRAN] tstep="  << std::scientific << cfg.tstep
              << ", tstop="       << cfg.tstop
              << ", tstart="      << cfg.tstart << "\n";

    int nSteps = static_cast<int>(std::floor(tstop / dt + 1e-12));
    std::cout << "[TRAN] total steps = " << nSteps << "\n";

    const int    maxNewtonIters = 50;
    const double tol            = 1e-6;
    const double gmin           = 1e-6;
    const double alpha          = 0.45;   // 与 DC 保持一致

    // 当前解，从 DC 解开始
    VectorXd x = xdc;

    // t = 0 的一行
    dumpRow(0.0, xdc);

    MatrixXd G(N, N);
    VectorXd I(N);

    for (int step = 0; step < nSteps; ++step) {
        double tNow = (step + 1) * dt;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 先 stamp 所有“非 C / 非 L / 非 MOS”的部分
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e))        continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e))      continue;
                e->stamp(G, I, ckt, x, ctx);
            }

            // 2) 线性 MOS 导电部分（Ids、gm、gds 等）
            for (const auto& m : mosfets) {
                m->stamp(G, I, ckt, x, ctx);
            }

            // 3) 显式电容的后向欧拉伴随模型
            for (const auto& c : caps) {
                double Cval = c->getC();
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c.get()]; // V(n1) - V(n2) at previous step
                stampCapBE(eq1, eq2, Cval, dt, vPrev, G, I);
            }

            // 4) 电感的后向欧拉 Thévenin 伴随模型
            for (const auto& L : inds) {
                double Lval = L->getL();
                if (Lval <= 0.0) continue;

                const auto& nodes = L->getNodeIds();
                int np = nodes[0];
                int nm = nodes[1];
                int eqP = ckt.nodes[np].eqIndex;
                int eqM = ckt.nodes[nm].eqIndex;
                int k   = L->getBranchEqIndex();
                if (k < 0 || k >= N) continue;

                double R_eq  = Lval / dt;
                double iPrev = indIprev[L.get()];
                double V_hist = -R_eq * iPrev;

                // KCL：节点电流 = I_L
                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;

                // 支路方程：Vp - Vm - R_eq * I_L = V_hist
                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
                G(k, k) += -R_eq;
                I(k)    += V_hist;
            }

            // 5) MOS 寄生电容：Cgs, Cgd, Cs, Cd
            for (const auto& m : mosfets) {
                const auto& nodes = m->getNodeIds();
                int nD = nodes[0];
                int nG = nodes[1];
                int nS = nodes[2];
                int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

                int eqD = ckt.nodes[nD].eqIndex;
                int eqG = ckt.nodes[nG].eqIndex;
                int eqS = ckt.nodes[nS].eqIndex;
                int eqB = ckt.nodes[nB].eqIndex;

                // 目前简单实现：用 Cj0 近似构造 Cgs/Cgd/Cs/Cd
                // 你可以在 MosfetBase 里加 Cox、W、L 等，把这里改成课件中的：
                // Cgs = 0.5 * Cox * W * L, Cgd 同理，Cs=Cd=Cj0
                double Cj0 = m->getCj0();
                double Cgs = 0.5 * Cj0;
                double Cgd = 0.5 * Cj0;
                double CsJ = Cj0;
                double CdJ = Cj0;

                const MosCapState& stPrev = mosPrev[m.get()];

                // Gate-source 电容
                stampCapBE(eqG, eqS, Cgs, dt, stPrev.vgsPrev, G, I);
                // Gate-drain 电容
                stampCapBE(eqG, eqD, Cgd, dt, stPrev.vgdPrev, G, I);
                // Source-bulk 结电容
                stampCapBE(eqS, eqB, CsJ, dt, stPrev.vsbPrev, G, I);
                // Drain-bulk 结电容
                stampCapBE(eqD, eqB, CdJ, dt, stPrev.vdbPrev, G, I);
            }

            // 6) gmin 到地，和 DC 保持一致
            stampGlobalGmin(ckt, G, gmin);

            // 7) 解线性化方程 G x_new = I
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);
            if (!xNew.allFinite()) {
                throw std::runtime_error("Transient: LU produced NaN/Inf.");
            }

            // 8) Newton 阻尼更新
            xNew = x + alpha * (xNew - x);
            double err = (xNew - x).norm();
            x = xNew;

            if (err < tol) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: transient Newton did not converge at t="
                          << std::scientific << tNow
                          << " (err=" << err << ")\n";
            }
        }

        // ===== 4. 步进成功：更新所有状态，并输出 =====
        // 显式电容
        for (const auto& c : caps) {
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            capVprev[c.get()] = v1 - v2;
        }
        // 电感
        for (const auto& L : inds) {
            int k = L->getBranchEqIndex();
            double iL = 0.0;
            if (k >= 0 && k < x.size()) {
                iL = x(k);
            }
            indIprev[L.get()] = iL;
        }
        // MOS 寄生电容状态
        for (const auto& m : mosfets) {
            const auto& nodes = m->getNodeIds();
            int nD = nodes[0];
            int nG = nodes[1];
            int nS = nodes[2];
            int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

            double vD = getNodeVoltage(ckt, x, nD);
            double vG = getNodeVoltage(ckt, x, nG);
            double vS = getNodeVoltage(ckt, x, nS);
            double vB = getNodeVoltage(ckt, x, nB);

            MosCapState st;
            st.vgsPrev = vG - vS;
            st.vgdPrev = vG - vD;
            st.vsbPrev = vS - vB;
            st.vdbPrev = vD - vB;
            mosPrev[m.get()] = st;
        }

        dumpRow(tNow, x);
    }

    std::cout << "Transient analysis (Backward Euler) finished. "
              << "Results written to '" << outFile << "'.\n";
}
