#define _USE_MATH_DEFINES
#include "analysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "sim.hpp"
#include "solver.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <functional>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ========== 电容的 BE 伴随模型 stamp（PSS / Shooting 用） ==========
void pssStampCapBE(int eq1, int eq2,
                   double C, double dt, double vPrev,
                   MatrixXd& G, VectorXd& I)
{
    if (C <= 0.0 || dt <= 0.0) return;

    double Gc = C / dt;  // 等效导纳

    // KCL: C dv/dt = I_c，相当于节点到地的导纳 + 历史电流源
    if (eq1 >= 0) G(eq1, eq1) += Gc;
    if (eq2 >= 0) G(eq2, eq2) += Gc;
    if (eq1 >= 0 && eq2 >= 0) {
        G(eq1, eq2) -= Gc;
        G(eq2, eq1) -= Gc;
    }

    // 历史源：I_hist = -C/dt * v_prev
    double I_hist = -Gc * vPrev;   // 从 eq1 流向 eq2

    if (eq1 >= 0) I(eq1) -= I_hist;
    if (eq2 >= 0) I(eq2) += I_hist;
}

// ========== 电感的 BE 伴随模型 stamp（PSS / Shooting 用） ==========
// 方程形式：
//   KCL at P:  +I_L + ... = 0
//   KCL at M:  -I_L + ... = 0
//   branch eq: Vp - Vm - R_eq * I_L = V_hist
// 其中 R_eq = L/dt,  V_hist = -R_eq * I_prev
void pssStampIndBE(int eqP, int eqM, int k,
                   double Lval, double dt, double iPrev,
                   MatrixXd& G, VectorXd& I)
{
    if (Lval <= 0.0 || dt <= 0.0) return;
    if (k < 0 || k >= G.rows()) {
        std::cerr << "PSS-SHOOT: invalid branchEqIndex for inductor (k="
                  << k << ")\n";
        return;
    }

    double R_eq   = Lval / dt;
    double V_hist = -R_eq * iPrev;

    // KCL at nodes P / M
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;

    // branch equation row
    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;
    G(k, k) += -R_eq;
    I(k)    += V_hist;
}

// ========== 全局 gmin，专门给 MatrixXd 用 ==========
void stampGlobalGmin(const Circuit& ckt, MatrixXd& G, double gmin)
{
    if (gmin <= 0.0) return;
    int N = (int)G.rows();
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < N) {
            G(eq, eq) += gmin;
        }
    }
}

// ========== MOS 寄生电容的历史状态（和瞬态里的结构一致） ==========
struct MosCapState {
    double vgsPrev = 0.0;
    double vgdPrev = 0.0;
    double vsbPrev = 0.0;
    double vdbPrev = 0.0;
};

// ========== 在一个周期内，用后向欧拉积分状态方程（真正和瞬态一致） ==========
VectorXd integrateOnePeriodPssBE(const Circuit& ckt,
                                 const SimulationConfig& sim,
                                 double dt, double T,
                                 const VectorXd& x0,
                                 const std::function<void(double,
                                                          const VectorXd&)>& dumpRow)
{
    const int N = ckt.numUnknowns();
    if (N <= 0) {
        throw std::runtime_error("PSS-SHOOT: circuit has no unknowns");
    }
    if (x0.size() != N) {
        throw std::runtime_error("PSS-SHOOT: x0 size mismatch");
    }
    if (dt <= 0.0 || T <= 0.0) {
        throw std::runtime_error("PSS-SHOOT: dt and T must be > 0");
    }

    double stepsDouble = T / dt;
    if (!std::isfinite(stepsDouble) || stepsDouble <= 0.0) {
        throw std::runtime_error("PSS-SHOOT: T/dt <= 0");
    }
    if (stepsDouble > 1e7) {
        throw std::runtime_error("PSS-SHOOT: T/dt too large (>1e7)");
    }

    int nSteps = (int)std::floor(stepsDouble + 1e-12);
    if (nSteps <= 0) {
        throw std::runtime_error("PSS-SHOOT: T/dt gives zero steps");
    }

    // === 收集动态元件 ===
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

    // === 根据 x0 初始化历史量 ===
    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, x0, n1);
        double v2 = getNodeVoltage(ckt, x0, n2);
        capVprev[c.get()] = v1 - v2;
    }

    std::unordered_map<const Inductor*, double> indIprev;
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        double i0 = 0.0;
        if (k >= 0 && k < x0.size()) {
            i0 = x0(k);
        }
        indIprev[L.get()] = i0;
    }

    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    for (auto& m : mosfets) {
        const auto& nodes = m->getNodeIds();
        int nD = nodes[0];
        int nG = nodes[1];
        int nS = nodes[2];
        int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

        double vD = getNodeVoltage(ckt, x0, nD);
        double vG = getNodeVoltage(ckt, x0, nG);
        double vS = getNodeVoltage(ckt, x0, nS);
        double vB = getNodeVoltage(ckt, x0, nB);

        MosCapState st;
        st.vgsPrev = vG - vS;
        st.vgdPrev = vG - vD;
        st.vsbPrev = vS - vB;
        st.vdbPrev = vD - vB;
        mosPrev[m.get()] = st;
    }

    VectorXd x = x0;
    MatrixXd G = MatrixXd::Zero(N, N);
    VectorXd I = VectorXd::Zero(N);

    const int    maxNewtonIters = 50;
    const double tol            = 1e-6;
    double       gmin           = 1e-6;
    const double alpha          = 0.45;   // 与 DC / TR 保持一致的阻尼

    // t = 0 的一行（如果需要输出）
    if (dumpRow) {
        dumpRow(0.0, x);
    }

    // === 时间步进：每个 dt 用 BE + Newton ===
    for (int step = 0; step < nSteps; ++step) {
        double tNow = (step + 1) * dt;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;  // 重用瞬态的 stamp 行为
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 除 C/L/MOS 以外的所有元件
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e))        continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e))      continue;
                e->stamp(G, I, ckt, x, ctx);
            }

            // 2) MOS 导电部分
            for (const auto& m : mosfets) {
                m->stamp(G, I, ckt, x, ctx);
            }

            // 3) 显式电容 C 的 BE 伴随模型
            for (const auto& c : caps) {
                double Cval = c->getC();
                if (Cval <= 0.0) continue;
                const auto& nodes = c->getNodeIds();
                int n1  = nodes[0];
                int n2  = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c.get()];
                pssStampCapBE(eq1, eq2, Cval, dt, vPrev, G, I);
            }

            // 4) 电感 L 的 BE 伴随模型
            for (const auto& L : inds) {
                double Lval = L->getL();
                if (Lval <= 0.0) continue;

                const auto& nodes = L->getNodeIds();
                int nP  = nodes[0];
                int nM  = nodes[1];
                int eqP = ckt.nodes[nP].eqIndex;
                int eqM = ckt.nodes[nM].eqIndex;
                int k   = L->getBranchEqIndex();

                double iPrev = indIprev[L.get()];
                pssStampIndBE(eqP, eqM, k, Lval, dt, iPrev, G, I);
            }

            // 5) MOS 寄生电容：用 Cj0 简易近似成 Cgs、Cgd、Csj、Cdj
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

                double Cj0 = m->getCj0();
                double Cg0 = m->getCg0();

                double Cgs = Cg0;
                double Cgd = Cg0;
                double CsJ = Cj0;
                double CdJ = Cj0;

                const MosCapState& stPrev = mosPrev[m.get()];

                pssStampCapBE(eqG, eqS, Cgs, dt, stPrev.vgsPrev, G, I);
                pssStampCapBE(eqG, eqD, Cgd, dt, stPrev.vgdPrev, G, I);
                pssStampCapBE(eqS, eqB, CsJ, dt, stPrev.vsbPrev, G, I);
                pssStampCapBE(eqD, eqB, CdJ, dt, stPrev.vdbPrev, G, I);
            }

            // 6) 全局 gmin
            stampGlobalGmin(ckt, G, gmin);

            // 7) 解 G x_new = I
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);
            if (!xNew.allFinite()) {
                gmin = std::min(gmin * 10.0, 1e-4);
                continue;
            }

            // 8) Newton 阻尼
            xNew = x + alpha * (xNew - x);
            double err = (xNew - x).norm();
            x = xNew;

            if (err < tol) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: BE-period Newton did not converge at t="
                          << std::scientific << tNow
                          << " (err=" << err << ")\n";
            }
        }

        // === 时步收敛：更新历史量 ===
        for (const auto& c : caps) {
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            capVprev[c.get()] = v1 - v2;
        }

        for (const auto& L : inds) {
            int k = L->getBranchEqIndex();
            double iL = 0.0;
            if (k >= 0 && k < x.size()) {
                iL = x(k);
            }
            indIprev[L.get()] = iL;
        }

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

        if (dumpRow) {
            dumpRow(tNow, x);
        }
    }

    return x;
}

// ========== Shooting-Newton 外层：解 F(x0) = x(T;x0) - x0 = 0 ==========
void runPssShootingAnalysis(const Circuit& ckt,
                            const SimulationConfig& sim,
                            const ShootingPssConfig& cfg,
                            const std::string& outFile)
{
    if (cfg.periodT <= 0.0) {
        std::cerr << "PSS-SHOOT: periodT must be > 0.\n";
        return;
    }

    double dt = cfg.tstep;
    if (dt <= 0.0) {
        dt = cfg.periodT / 200.0;   // 默认一个周期 200 个点
    }

    const int N = ckt.numUnknowns();
    if (N <= 0) {
        std::cerr << "PSS-SHOOT: circuit has no unknowns.\n";
        return;
    }

    // 1) DC 工作点作为初始 shooting 向量
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
    VectorXd x0;
    try {
        x0 = dc.run();
    } catch (const std::exception& e) {
        std::cerr << "PSS-SHOOT: DC operating point failed: " << e.what() << "\n";
        return;
    }
    if (x0.size() != N) {
        std::cerr << "PSS-SHOOT: DC solution size mismatch.\n";
        return;
    }

    std::cout << "PSS-SHOOT: Shooting Method Analysis\n";
    std::cout << "  Period T = " << std::scientific << cfg.periodT << "\n";
    std::cout << "  dt       = " << dt << "\n";
    std::cout << "  maxIters = " << cfg.maxIters
              << ", tol = " << cfg.tol
              << ", relax = " << cfg.relax << "\n";

    VectorXd xInit = x0;

    // 2) Shooting-Newton 外层迭代
    const double relPert = 1e-3;
    const double absPert = 1e-6;

    for (int it = 0; it < cfg.maxIters; ++it) {
        // 一次射击：xInit -> xT
        VectorXd xT;
        try {
            xT = integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, nullptr);
        } catch (const std::exception& e) {
            std::cerr << "PSS-SHOOT: integrateOnePeriodPssBE failed at iter "
                      << it << ": " << e.what() << "\n";
            return;
        }

        // 残差 F = xT - xInit
        VectorXd F = xT - xInit;
        double err = F.norm();
        std::cout << "[PSS-SHOOT] Iter " << it
                  << "  ||x(T) - x(0)|| = " << std::scientific << err << "\n";

        if (!std::isfinite(err)) {
            std::cerr << "PSS-SHOOT: non-finite residual, abort.\n";
            return;
        }

        if (err < cfg.tol) {
            std::cout << "PSS-SHOOT: converged after " << it << " iterations.\n";
            break;
        }

        // 构造 Jacobian: J = dF/dx0 ≈ (Phi(x0+he_j) - Phi(x0))/h - I
        MatrixXd J(N, N);
        J.setZero();

        for (int j = 0; j < N; ++j) {
            double scale = std::max(std::abs(xInit(j)), 1.0);
            double h     = std::max(absPert, relPert * scale);

            VectorXd xPert = xInit;
            xPert(j) += h;

            VectorXd xT_pert;
            try {
                xT_pert = integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT,
                                                  xPert, nullptr);
            } catch (const std::exception& e) {
                std::cerr << "PSS-SHOOT: integrateOnePeriodPssBE failed "
                             "when building Jacobian at col "
                          << j << ": " << e.what() << "\n";
                return;
            }

            VectorXd dPhi = (xT_pert - xT) / h;
            J.col(j) = dPhi;
        }

        // F(x) = Phi(x) - x，所以 J = dPhi/dx - I
        for (int i = 0; i < N; ++i) {
            J(i, i) -= 1.0;
        }

        // 解 J dx = -F
        VectorXd rhs = -F;
        VectorXd dx  = Solver::solveLinearSystemLU(J, rhs);
        if (!dx.allFinite()) {
            std::cerr << "PSS-SHOOT: Newton linear solve produced NaN/Inf, "
                         "fallback to simple relaxation.\n";
            xInit += cfg.relax * F;
        } else {
            double alphaShoot = cfg.relax;   // 外层阻尼
            xInit += alphaShoot * dx;
        }

        if (it == cfg.maxIters - 1) {
            std::cerr << "PSS-SHOOT: not converged within maxIters.\n";
        }
    }

    // 3) 用最终 xInit 再射一次，输出一个周期的波形到 CSV
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "PSS-SHOOT: cannot open output file '"
                  << outFile << "'.\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    // 表头：time, V(node), I(Vsrc), I(L)
    ofs << "time";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            ofs << ",I(" << vs->getName() << ")";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            ofs << ",I(" << L->getName() << ")";
        }
    }
    ofs << "\n";

    auto dumpRow = [&](double t, const VectorXd& x) {
        ofs << t;
        // 节点电压
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) {
                ofs << "," << x(node.eqIndex);
            }
        }
        // 电压源 / 电感电流
        for (const auto& e : ckt.elements) {
            if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
                int k  = vs->getBranchEqIndex();
                double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << Ibr;
            } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
                int k  = L->getBranchEqIndex();
                double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << Ibr;
            }
        }
        ofs << "\n";
    };

    // t = 0 行
    dumpRow(0.0, xInit);

    try {
        integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, dumpRow);
    } catch (const std::exception& e) {
        std::cerr << "PSS-SHOOT: final integrateOnePeriodPssBE failed: "
                  << e.what() << "\n";
        return;
    }

    std::cout << "PSS-SHOOT: results written to " << outFile << "\n";
}
