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
#include <stdexcept>
#include <string>
#include <algorithm>


using Eigen::MatrixXd;
using Eigen::VectorXd;

// ==========================================================
// 辅助函数：Stamp 动态元件的 "质量矩阵" B = C/dt
// 用于灵敏度分析 RHS: J * S_new = B * S_old
// ==========================================================
void stampDynamicMatrix(const Circuit& ckt, MatrixXd& B, double dt) {
    if (dt <= 0.0) return;
    double invDt = 1.0 / dt;

    // 1. 独立电容
    for (const auto& e : ckt.elements) {
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
            double val = c->getC() * invDt;
            const auto& nodes = c->getNodeIds();
            int eq1 = ckt.nodes[nodes[0]].eqIndex;
            int eq2 = ckt.nodes[nodes[1]].eqIndex;
            
            if (eq1 >= 0) B(eq1, eq1) += val;
            if (eq2 >= 0) B(eq2, eq2) += val;
            if (eq1 >= 0 && eq2 >= 0) {
                B(eq1, eq2) -= val;
                B(eq2, eq1) -= val;
            }
        }
    }

    // 2. 电感 (Branch equation contribution)
    // 电感方程: vP - vM - (L/dt)*iL = ...
    // 对 i_prev 的导数出现在 RHS，系数为 -(L/dt)
    // 注意方程符号： Row k: vP - vM - R_eq * iL = ...
    // 对应 RHS 项是 -R_eq * i_prev。
    // 所以 B(k, k) 应该填入 -L/dt
    for (const auto& e : ckt.elements) {
        if (auto ind = std::dynamic_pointer_cast<Inductor>(e)) {
            double val = - (ind->getL() * invDt); 
            int k = ind->getBranchEqIndex();
            if (k >= 0 && k < B.rows()) {
                B(k, k) += val;
            }
        }
    }

    // 3. MOSFET 寄生电容 (近似为常数 Cg0, Cj0，与 pssStampCapBE 一致)
    for (const auto& e : ckt.elements) {
        if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) {
            const auto& nodes = m->getNodeIds();
            int eqD = ckt.nodes[nodes[0]].eqIndex;
            int eqG = ckt.nodes[nodes[1]].eqIndex;
            int eqS = ckt.nodes[nodes[2]].eqIndex;
            int eqB = ckt.nodes[(nodes.size() > 3) ? nodes[3] : nodes[2]].eqIndex;

            double Cj0 = m->getCj0();
            double Cg0 = m->getCg0();
            
            auto stampC = [&](int i, int j, double cVal) {
                if (cVal <= 0) return;
                double g = cVal * invDt;
                if (i >= 0) B(i, i) += g;
                if (j >= 0) B(j, j) += g;
                if (i >= 0 && j >= 0) {
                    B(i, j) -= g;
                    B(j, i) -= g;
                }
            };

            stampC(eqG, eqS, Cg0); // Cgs
            stampC(eqG, eqD, Cg0); // Cgd
            stampC(eqS, eqB, Cj0); // Csb
            stampC(eqD, eqB, Cj0); // Cdb
        }
    }
}

// ========== 复用原有的 BE Stamp 函数 ==========
void pssStampCapBE(int eq1, int eq2, double C, double dt, double vPrev, MatrixXd& G, VectorXd& I) {
    if (C <= 0.0 || dt <= 0.0) return;
    double Gc = C / dt;
    if (eq1 >= 0) G(eq1, eq1) += Gc;
    if (eq2 >= 0) G(eq2, eq2) += Gc;
    if (eq1 >= 0 && eq2 >= 0) { G(eq1, eq2) -= Gc; G(eq2, eq1) -= Gc; }
    double I_hist = -Gc * vPrev;
    if (eq1 >= 0) I(eq1) -= I_hist;
    if (eq2 >= 0) I(eq2) += I_hist;
}

void pssStampIndBE(int eqP, int eqM, int k, double Lval, double dt, double iPrev, MatrixXd& G, VectorXd& I) {
    if (Lval <= 0.0 || dt <= 0.0) return;
    double R_eq = Lval / dt;
    double V_hist = -R_eq * iPrev;
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;
    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;
    G(k, k) += -R_eq;
    I(k) += V_hist;
}

void stampGlobalGmin(const Circuit& ckt, MatrixXd& G, double gmin) {
    if (gmin <= 0.0) return;
    int N = (int)G.rows();
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < N) G(eq, eq) += gmin;
    }
}

struct MosCapState {
    double vgsPrev = 0.0;
    double vgdPrev = 0.0;
    double vsbPrev = 0.0;
    double vdbPrev = 0.0;
};

static bool solveOneStepBE(const Circuit& ckt,
                           const SimulationConfig& sim,
                           double tNow, double dt,
                           VectorXd& x, // in: x_n (as initial guess), out: x_{n+1}
                           const std::unordered_map<const CapacitorElement*, double>& capVprev,
                           const std::unordered_map<const Inductor*, double>& indIprev,
                           const std::unordered_map<const MosfetBase*, MosCapState>& mosPrev,
                           MatrixXd& outA)
{
    (void)sim;
    const int N = ckt.numUnknowns();
    const int maxNewtonIters = 50;
    const double tol = 1e-7;
    const double damp = 1.0; // 不稳就改成 0.5~0.9

    MatrixXd G(N, N);
    VectorXd I(N);

    // gmin continuation inside one step
    const int maxGminTries = 6;
    double gmin = 1e-12;

    for (int gtry = 0; gtry < maxGminTries; ++gtry) {
        VectorXd xk = x;

        bool converged = false;
        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type = AnalysisType::TRAN;
            ctx.time = tNow; // backward-euler at t_{n+1}

            // 1) static elements: R, independent sources, etc.
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e)) continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e)) continue;
                e->stamp(G, I, ckt, xk, ctx);
            }

            // 2) MOS conduction (Newton linearization)
            for (const auto& e : ckt.elements) {
                auto m = std::dynamic_pointer_cast<MosfetBase>(e);
                if (!m) continue;
                m->stamp(G, I, ckt, xk, ctx);
            }

            // 3) BE dynamic stamps (caps / inductors / MOS caps) using HISTORY from x_n
            for (const auto& e : ckt.elements) {
                if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
                    const auto& nodes = c->getNodeIds();
                    const int eq1 = ckt.nodes[nodes[0]].eqIndex;
                    const int eq2 = ckt.nodes[nodes[1]].eqIndex;
                    const double vPrev = capVprev.at(c.get());
                    pssStampCapBE(eq1, eq2, c->getC(), dt, vPrev, G, I);
                } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
                    const auto& nodes = L->getNodeIds();
                    const int eqP = ckt.nodes[nodes[0]].eqIndex;
                    const int eqM = ckt.nodes[nodes[1]].eqIndex;
                    const int k = L->getBranchEqIndex();
                    if (k < 0 || k >= N) continue;
                    const double iPrev = indIprev.at(L.get());
                    pssStampIndBE(eqP, eqM, k, L->getL(), dt, iPrev, G, I);
                } else if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) {
                    const auto& nodes = m->getNodeIds();
                    const int eqD = ckt.nodes[nodes[0]].eqIndex;
                    const int eqG = ckt.nodes[nodes[1]].eqIndex;
                    const int eqS = ckt.nodes[nodes[2]].eqIndex;
                    const int eqB = ckt.nodes[(nodes.size() > 3) ? nodes[3] : nodes[2]].eqIndex;

                    const MosCapState& st = mosPrev.at(m.get());
                    const double Cj0 = m->getCj0();
                    const double Cg0 = m->getCg0();

                    pssStampCapBE(eqG, eqS, Cg0, dt, st.vgsPrev, G, I); // Cgs
                    pssStampCapBE(eqG, eqD, Cg0, dt, st.vgdPrev, G, I); // Cgd
                    pssStampCapBE(eqS, eqB, Cj0, dt, st.vsbPrev, G, I); // Csb
                    pssStampCapBE(eqD, eqB, Cj0, dt, st.vdbPrev, G, I); // Cdb
                }
            }

            stampGlobalGmin(ckt, G, gmin);

            // Solve the linearized system
            VectorXd xsol = Solver::solveLinearSystemLU(G, I);
            if (!xsol.allFinite()) {
                converged = false;
                break;
            }

            VectorXd xnew = xk + damp * (xsol - xk);
            double stepNorm = (xnew - xk).norm();
            xk.swap(xnew);

            if (stepNorm < tol * (1.0 + xk.norm())) {
                converged = true;
                break;
            }
        }

        if (converged) {
            // Re-stamp one more time at the final x to get accurate Jacobian A
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type = AnalysisType::TRAN;
            ctx.time = tNow;

            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e)) continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e)) continue;
                e->stamp(G, I, ckt, xk, ctx);
            }
            for (const auto& e : ckt.elements) {
                auto m = std::dynamic_pointer_cast<MosfetBase>(e);
                if (!m) continue;
                m->stamp(G, I, ckt, xk, ctx);
            }
            for (const auto& e : ckt.elements) {
                if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
                    const auto& nodes = c->getNodeIds();
                    const int eq1 = ckt.nodes[nodes[0]].eqIndex;
                    const int eq2 = ckt.nodes[nodes[1]].eqIndex;
                    const double vPrev = capVprev.at(c.get());
                    pssStampCapBE(eq1, eq2, c->getC(), dt, vPrev, G, I);
                } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
                    const auto& nodes = L->getNodeIds();
                    const int eqP = ckt.nodes[nodes[0]].eqIndex;
                    const int eqM = ckt.nodes[nodes[1]].eqIndex;
                    const int k = L->getBranchEqIndex();
                    if (k < 0 || k >= N) continue;
                    const double iPrev = indIprev.at(L.get());
                    pssStampIndBE(eqP, eqM, k, L->getL(), dt, iPrev, G, I);
                } else if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) {
                    const auto& nodes = m->getNodeIds();
                    const int eqD = ckt.nodes[nodes[0]].eqIndex;
                    const int eqG = ckt.nodes[nodes[1]].eqIndex;
                    const int eqS = ckt.nodes[nodes[2]].eqIndex;
                    const int eqB = ckt.nodes[(nodes.size() > 3) ? nodes[3] : nodes[2]].eqIndex;

                    const MosCapState& st = mosPrev.at(m.get());
                    const double Cj0 = m->getCj0();
                    const double Cg0 = m->getCg0();

                    pssStampCapBE(eqG, eqS, Cg0, dt, st.vgsPrev, G, I);
                    pssStampCapBE(eqG, eqD, Cg0, dt, st.vgdPrev, G, I);
                    pssStampCapBE(eqS, eqB, Cj0, dt, st.vsbPrev, G, I);
                    pssStampCapBE(eqD, eqB, Cj0, dt, st.vdbPrev, G, I);
                }
            }

            stampGlobalGmin(ckt, G, gmin);

            x = xk;
            outA = G;
            return true;
        }

        gmin *= 10.0;
    }

    return false;
}


// =======================================================================
// 积分一个周期，并可选地计算 Monodromy 矩阵 (Sensitivity Analysis)
// =======================================================================
VectorXd integrateOnePeriodPssBE(const Circuit& ckt,
                                 const SimulationConfig& sim,
                                 double dt, double T,
                                 const VectorXd& x0,
                                 MatrixXd* outMonodromy,
                                 const std::function<void(double, const VectorXd&)>& dumpRow)
{
    const int N = ckt.numUnknowns();
    VectorXd x = x0;

    // Collect element lists and init history from x0
    std::vector<std::shared_ptr<CapacitorElement>> caps;
    std::vector<std::shared_ptr<Inductor>> inds;
    std::vector<std::shared_ptr<MosfetBase>> mosfets;

    for (const auto& e : ckt.elements) {
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) caps.push_back(c);
        else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) inds.push_back(L);
        else if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) mosfets.push_back(m);
    }

    std::unordered_map<const CapacitorElement*, double> capVprev;
    capVprev.reserve(caps.size());
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        capVprev[c.get()] = getNodeVoltage(ckt, x0, nodes[0]) - getNodeVoltage(ckt, x0, nodes[1]);
    }

    std::unordered_map<const Inductor*, double> indIprev;
    indIprev.reserve(inds.size());
    for (auto& L : inds) {
        const int k = L->getBranchEqIndex();
        indIprev[L.get()] = (k >= 0 && k < x0.size()) ? x0(k) : 0.0;
    }

    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    mosPrev.reserve(mosfets.size());
    for (auto& m : mosfets) {
        const auto& nodes = m->getNodeIds();
        const double vD = getNodeVoltage(ckt, x0, nodes[0]);
        const double vG = getNodeVoltage(ckt, x0, nodes[1]);
        const double vS = getNodeVoltage(ckt, x0, nodes[2]);
        const double vB = getNodeVoltage(ckt, x0, (nodes.size() > 3) ? nodes[3] : nodes[2]);
        MosCapState st;
        st.vgsPrev = vG - vS;
        st.vgdPrev = vG - vD;
        st.vsbPrev = vS - vB;
        st.vdbPrev = vD - vB;
        mosPrev[m.get()] = st;
    }

    MatrixXd S;
    if (outMonodromy) S = MatrixXd::Identity(N, N);

    if (dumpRow) dumpRow(0.0, x);

    double t = 0.0;
    const double eps = 1e-15;

    while (t < T - eps) {
        const double dtTry = std::min(dt, T - t);
        const double tNext = t + dtTry;

        // 1) Solve x_{n+1} and get step Jacobian A
        MatrixXd A(N, N);
        VectorXd xGuess = x;
        bool ok = solveOneStepBE(ckt, sim, tNext, dtTry, xGuess, capVprev, indIprev, mosPrev, A);
        if (!ok) {
            throw std::runtime_error("PSS shooting: BE step Newton failed (t=" + std::to_string(tNext) + ")");
        }
        x = xGuess;

        // 2) Sensitivity recursion:  A * S_{n+1} = B * S_n
        if (outMonodromy) {
            MatrixXd B = MatrixXd::Zero(N, N);
            stampDynamicMatrix(ckt, B, dtTry);

            MatrixXd RHS = B * S;

            // Factorize A once and solve for all RHS columns
            Eigen::PartialPivLU<MatrixXd> lu(A);
            S = lu.solve(RHS);
        }

        // 3) Update history using x (as x_n for next step)
        for (auto& c : caps) {
            const auto& nodes = c->getNodeIds();
            capVprev[c.get()] = getNodeVoltage(ckt, x, nodes[0]) - getNodeVoltage(ckt, x, nodes[1]);
        }
        for (auto& L : inds) {
            const int k = L->getBranchEqIndex();
            indIprev[L.get()] = (k >= 0 && k < x.size()) ? x(k) : 0.0;
        }
        for (auto& m : mosfets) {
            const auto& nodes = m->getNodeIds();
            const double vD = getNodeVoltage(ckt, x, nodes[0]);
            const double vG = getNodeVoltage(ckt, x, nodes[1]);
            const double vS = getNodeVoltage(ckt, x, nodes[2]);
            const double vB = getNodeVoltage(ckt, x, (nodes.size() > 3) ? nodes[3] : nodes[2]);
            MosCapState st;
            st.vgsPrev = vG - vS;
            st.vgdPrev = vG - vD;
            st.vsbPrev = vS - vB;
            st.vdbPrev = vD - vB;
            mosPrev[m.get()] = st;
        }

        if (dumpRow) dumpRow(tNext, x);
        t = tNext;
    }

    if (outMonodromy) *outMonodromy = S;
    return x;
}

// ==========================================================
// PSS Shooting 主函数 (Sensitivty + Broyden)
// ==========================================================
void runPssShootingAnalysis(const Circuit& ckt,
                            const SimulationConfig& sim,
                            const ShootingPssConfig& cfg,
                            const std::string& outFile)
{
    if (cfg.periodT <= 0.0) return;

    const double T = cfg.periodT;
    const double dt = (cfg.tstep > 0.0) ? cfg.tstep : (T / 400.0);

    const int N = ckt.numUnknowns();
    if (N <= 0) return;

    // 1) DC init
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
    VectorXd xInit = dc.run();
    if (sim.verbose) {
    std::cout << "PSS-SHOOT (Sensitivity Matrix / Monodromy)\n"
              << "  Period T = " << T << "\n"
              << "  dt       = " << dt << "\n"
              << "  N        = " << N << "\n";
    }
    const int maxIters = std::max(1, cfg.maxIters);
    const double tol = std::max(1e-14, cfg.tol);

    // Outer Newton on H(x0) = x(T;x0) - x0
    for (int it = 0; it < maxIters; ++it) {
        MatrixXd M(N, N);
        VectorXd xT = integrateOnePeriodPssBE(ckt, sim, dt, T, xInit, &M, nullptr);

        VectorXd H = xT - xInit;
        double hnorm = H.norm();
        if (sim.verbose) {
            std::cout << "  Iter " << it << "  ||H|| = " << std::scientific << hnorm << "\n";
            if (hnorm < tol) {
                std::cout << "  Converged.\n";
                break;
            }
        }
        MatrixXd J = M - MatrixXd::Identity(N, N);

        // Solve J * dx = -H
        VectorXd dx = Solver::solveLinearSystemLU(J, -H);
        if (!dx.allFinite()) {
            throw std::runtime_error("PSS shooting: linear solve failed (Jacobian singular/NaN)");
        }

        // Damping + backtracking (robust)
        double alpha = (cfg.relax > 0.0 && cfg.relax <= 1.0) ? cfg.relax : 1.0;
        bool accepted = false;

        for (int ls = 0; ls < 8; ++ls) {
            VectorXd xTrial = xInit + alpha * dx;
            VectorXd xT_trial = integrateOnePeriodPssBE(ckt, sim, dt, T, xTrial, nullptr, nullptr);
            VectorXd H_trial = xT_trial - xTrial;
            double n_trial = H_trial.norm();

            if (std::isfinite(n_trial) && n_trial < hnorm) {
                xInit.swap(xTrial);
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if (!accepted) {
            xInit += 1e-2 * dx; // fallback
        }
    }

    // 2) Output one-period waveform
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Cannot open output file: " << outFile << "\n";
        return;
    }
    ofs << std::setprecision(12);

    // Header: time, node voltages, then branch currents (V sources + L)
    ofs << "time";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) ofs << ",V(" << node.name << ")";
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
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) ofs << "," << x(node.eqIndex);
        }
        for (const auto& e : ckt.elements) {
            if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
                int k = vs->getBranchEqIndex();
                double val = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << val;
            } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
                int k = L->getBranchEqIndex();
                double val = (k >= 0 && k < x.size()) ? x(k) : 0.0;
                ofs << "," << val;
            }
        }
        ofs << "\n";
    };

    dumpRow(0.0, xInit);
    integrateOnePeriodPssBE(ckt, sim, dt, T, xInit, nullptr, dumpRow);

    std::cout << "PSS waveform written to: " << outFile << "\n";
}
