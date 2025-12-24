
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>

#include "analysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "sim.hpp"
#include "solver.hpp"
#include "utils.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;


TransientAnalysis::TransientAnalysis(const Circuit& ckt_,
                                     const SimulationConfig& sim_,
                                     const std::string& outFile_)
    : ckt(ckt_), sim(sim_), outFile(outFile_) {}

VectorXd TransientAnalysis::computeDcOperatingPoint() const {
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel); // 或LU
    return dc.run();
}


// 在 G / I 中 stamp 一个“电容 + 历史电流源”
// C 接在 eq1 与 eq2 之间，vPrev = V(eq1)^n - V(eq2)^n
void TransientAnalysis::stampCapBE(int eq1, int eq2,
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

// 梯形法版本的电容伴随模型
// C 接在 eq1 与 eq2 之间，vPrev = V(eq1)^n - V(eq2)^n，iPrev 为上一步从 eq1->eq2 的电流
void TransientAnalysis::stampCapTR(int eq1, int eq2,
                       double C, double dt,
                       double vPrev, double iPrev,
                       MatrixXd& G, VectorXd& I)
{
    if (C <= 0.0 || dt <= 0.0) return;

    // 梯形法等效电导
    double Gc = 2.0 * C / dt;

    if (eq1 >= 0) G(eq1, eq1) += Gc;
    if (eq2 >= 0) G(eq2, eq2) += Gc;
    if (eq1 >= 0 && eq2 >= 0) {
        G(eq1, eq2) -= Gc;
        G(eq2, eq1) -= Gc;
    }

    // 梯形法的历史电流源：
    //   i^{n+1} = Gc * v^{n+1} + I_hist
    // 推导得：I_hist = -Gc * vPrev - iPrev
    double I_hist = -Gc * vPrev - iPrev;
    if (eq1 >= 0) I(eq1) -= I_hist;
    if (eq2 >= 0) I(eq2) += I_hist;
}


// ========= 主函数：后向欧拉瞬态 + MOS 寄生电容 =========
void TransientAnalysis::runBackwardEuler() {
    const TranConfig& cfg = sim.tran;

    if (!cfg.enabled) {
        std::cerr << "Transient analysis is not enabled (.TRAN missing).\n";
        return;
    }

    if (cfg.tstep <= 0.0 || cfg.tstop <= 0.0) {
        std::cerr << "Invalid .TRAN card: tstep and tstop must be > 0.\n";
        return;
    }

    const double dtInit = cfg.tstep;          // 初始 dt
    double       dt     = dtInit;             // 当前 dt
    const double dtMin  = dtInit / 64.0;      // 最小 dt
    const double dtMax  = std::min(dtInit * 16.0, cfg.tstop / 400);     // 最大 dt
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
        xdc = computeDcOperatingPoint();
    } catch (const std::exception& e) {
        std::cerr << "DC operating point failed: " << e.what() << "\n";
        return;
    }

    if (xdc.size() != N) {
        std::cerr << "Transient: DC solution size mismatch.\n";
        return;
    }

    // 计算截断误差阈值：遍历电压源，取最大电压值的0.2%
    double maxVoltage = 0.0;
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            double Vac = 0.0;
            if (vs->getSpec().tran.type == WaveformType::SIN) {
                Vac = vs->getSpec().tran.sine.va;
            }
            double v = std::max(std::abs(vs->getSpec().dcValue), Vac);
            if (v > maxVoltage) maxVoltage = v;
        }
    }
    printf("bound for truncation error: %.6e V\n", maxVoltage * 0.2);
    double bound = maxVoltage * 0.2;

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

    // ===== 3. 自适应步长：收敛 + 阈值判定 =====
    std::cout << "[TRAN] tstep="  << std::scientific << cfg.tstep
              << ", tstop="       << cfg.tstop
              << ", tstart="      << cfg.tstart << "\n";

    int nStepsInit = static_cast<int>(std::floor(tstop / dtInit + 1e-12));
    std::cout << "[TRAN] total steps (based on initial dt) = " << nStepsInit << "\n";

    const int    maxNewtonIters = 50;
    const double tol            = 1e-6;
    const double gmin           = 1e-6;
    const double alpha          = 0.45;   // 与原版保持一致


    // 当前接受解（从 DC 开始）
    VectorXd x = xdc;
    VectorXd xPrevStep = xdc; // 上一个接受解（用于预测）
    double   dtPrev = dtInit;

    dumpRow(0.0, xdc);

    MatrixXd G(N, N);
    VectorXd I(N);

    // “一步 BE 求解”封装成 lambda：只负责 Newton/LU + 阻尼，不改状态表
    auto solveOneStepBE = [&](double tNow, double dtTry, VectorXd& xSolve) -> bool {
        double lastErr = 1e300;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 非 C / 非 L / 非 MOS
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e))        continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e))      continue;
                e->stamp(G, I, ckt, xSolve, ctx);
            }

            // 2) MOS 导电
            for (const auto& m : mosfets) {
                m->stamp(G, I, ckt, xSolve, ctx);
            }

            // 3) 显式电容 BE
            for (const auto& c : caps) {
                double Cval = c->getC();
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c.get()];
                stampCapBE(eq1, eq2, Cval, dtTry, vPrev, G, I);
            }

            // 4) 电感 BE
            for (const auto& L : inds) {
                double Lval = L->getL();
                if (Lval <= 0.0) continue;

                const auto& nodes = L->getNodeIds();
                int np  = nodes[0];
                int nm  = nodes[1];
                int eqP = ckt.nodes[np].eqIndex;
                int eqM = ckt.nodes[nm].eqIndex;
                int k   = L->getBranchEqIndex();
                if (k < 0 || k >= N) continue;

                double iPrev = indIprev[L.get()];
                double G_eq = dt / Lval;
                double I_hist = iPrev; 
                if (eqP >= 0) {
                    G(eqP, eqP) += G_eq;
                    I(eqP) += -I_hist;
                }
                if (eqM >= 0) {
                    G(eqM, eqM) += G_eq; 
                    I(eqM) += I_hist;
                }
                if (eqP >= 0 && eqM >= 0) {
                    G(eqP, eqM) -= G_eq;
                    G(eqM, eqP) -= G_eq;
                }
                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;
                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
            }


            // 5) MOS 寄生电容 BE
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

                double Cg0 = m->getCg0();
                double Cj0 = m->getCj0();
                double Cgs = Cg0;
                double Cgd = Cg0;
                double CsJ = Cj0;
                double CdJ = Cj0;

                const MosCapState& stPrev = mosPrev[m.get()];

                if (Cgs > 0.0) stampCapBE(eqG, eqS, Cgs, dtTry, stPrev.vgsPrev, G, I);
                if (Cgd > 0.0) stampCapBE(eqG, eqD, Cgd, dtTry, stPrev.vgdPrev, G, I);
                if (CsJ > 0.0) stampCapBE(eqS, eqB, CsJ, dtTry, stPrev.vsbPrev, G, I);
                if (CdJ > 0.0) stampCapBE(eqD, eqB, CdJ, dtTry, stPrev.vdbPrev, G, I);
            }

            // 6) gmin
            stampGlobalGmin(ckt, G, gmin);

            // 7) 解线性化方程
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);
            if (!xNew.allFinite()) {
                throw std::runtime_error("Transient: LU produced NaN/Inf.");
            }

            // 8) Newton 阻尼（与原版一致）
            xNew = xSolve + alpha * (xNew - xSolve);
            double err = (xNew - xSolve).norm();
            xSolve = xNew;

            lastErr = err;
            if (err < tol) return true;
        }

        return false;
    };

    auto passTruncation = [&](const VectorXd& xPred, const VectorXd& xSol) -> bool {
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex < 0) continue;
            int k = node.eqIndex;
            if (k < 0 || k >= xSol.size() || k >= xPred.size()) continue;

            double v  = xSol(k);
            double vp = xPred(k);
            double av = std::abs(v);
            double ap = std::abs(vp);
            double vmax = (av > ap) ? av : ap;

            if (std::abs(v - vp) > bound) return false;
        }
        return true;
    };

    double tNow = 0.0;
    int    step = 0;

    while (tNow + 1e-15 < tstop) {
        bool accepted = false;
        double dtTry  = dt;

        // 允许最后一个“余量步”小于 dtMin（只要不是通过不断减半逼到的）
        const double dtRemain0 = tstop - tNow;
        const bool allowFinalBelowDtMin = (dtRemain0 < dtMin * (1.0 - 1e-15));

        while (!accepted) {
            if (!allowFinalBelowDtMin && dtTry < dtMin * (1.0 - 1e-15)) {
                throw std::runtime_error("Transient(BE): dt reduced below dtMin, abort.");
            }

            double dtRemain = tstop - tNow;
            if (dtTry > dtRemain) dtTry = dtRemain;

            // 保存上一次“接受的解”
            VectorXd xBefore = x;

            // 预测：x_pred = x + (dtTry/dtPrev) * (x - xPrevStep)
            VectorXd xPred = xBefore;
            if (step > 0) {
                double ratio = (dtPrev > 0.0) ? (dtTry / dtPrev) : 1.0;
                xPred = xBefore + ratio * (xBefore - xPrevStep);
            }

            // 用预测作为初值
            x = xPred;

            // 先尝试求解
            bool converged = solveOneStepBE(tNow + dtTry, dtTry, x);
            if (!converged) {
                if (!allowFinalBelowDtMin && dtTry <= dtMin * (1.0 + 1e-15)) {
                    throw std::runtime_error("Transient(BE): reached dtMin but still not converged, abort.");
                }
                dtTry *= 0.5;
                x = xBefore;
                continue;
            }

            // 第一次：只要收敛就接受；从第二次开始做阈值判断
            if (step > 0) {
                if (!passTruncation(xPred, x)) {
                    if (!allowFinalBelowDtMin && dtTry <= dtMin * (1.0 + 1e-15)) {
                        throw std::runtime_error("Transient(BE): reached dtMin but truncation test still fails, abort.");
                    }
                    dtTry *= 0.5;
                    x = xBefore;
                    continue;
                }
            }

            // ===== 接受该步：更新状态并输出 =====
            tNow += dtTry;

            // 显式电容 vPrev
            for (const auto& c : caps) {
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                double v1 = getNodeVoltage(ckt, x, n1);
                double v2 = getNodeVoltage(ckt, x, n2);
                capVprev[c.get()] = v1 - v2;
            }
            // 电感电流
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

            // 维护预测历史
            xPrevStep = xBefore;
            dtPrev    = dtTry;

            // 当前 dt翻倍
            dt = dtTry >= 0.5 * dtMax ? dt : dtTry * 2.0;

            ++step;
            accepted = true;
        }
    }
    if (sim.verbose) {
    std::cout << "Transient analysis (Backward Euler) finished. "
              << "Results written to '" << outFile << "'.\n";
    }
}

// ========= 梯形法瞬态（Trapezoidal Rule） =========
void TransientAnalysis::runTrapezoidal(){
    const TranConfig& cfg = sim.tran;

    if (!cfg.enabled) {
        std::cerr << "Transient analysis (TR) is not enabled (.TRAN missing).\n";
        return;
    }

    if (cfg.tstep <= 0.0 || cfg.tstop <= 0.0) {
        std::cerr << "Invalid .TRAN card: tstep and tstop must be > 0.\n";
        return;
    }

    const double dtInit = cfg.tstep;      // 初始 dt（来自网表）
    double       dt     = dtInit;         // 当前 dt（只减不增）
    const double dtMin  = dtInit / 64.0;  // 最小 dt
    const double dtMax  = std::min(dtInit * 16.0, cfg.tstop / 400);     // 最大 dt
    const double tstop  = cfg.tstop;
    const double tstart = cfg.tstart;

    const int N = ckt.numUnknowns();
    if (N <= 0) {
        std::cerr << "Transient (TR): circuit has no unknowns.\n";
        return;
    }

    // ===== 0. DC 工作点作为 t=0 初值 =====
    VectorXd xdc;
    try {
        xdc = computeDcOperatingPoint();
    } catch (const std::exception& e) {
        std::cerr << "DC operating point failed: " << e.what() << "\n";
        return;
    }
    if (xdc.size() != N) {
        std::cerr << "Transient (TR): DC solution size mismatch.\n";
        return;
    }

    // 计算截断误差阈值：遍历电压源，取最大电压值的0.2%
    double maxVoltage = 0.0;
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            double Vac = 0.0;
            if (vs->getSpec().tran.type == WaveformType::SIN) {
                Vac = vs->getSpec().tran.sine.va;
            }
            double v = std::max(std::abs(vs->getSpec().dcValue), Vac);
            if (v > maxVoltage) maxVoltage = v;
        }
    }
    printf("bound for truncation error: %.6e V\n", maxVoltage * 0.2);
    double bound = maxVoltage * 0.2;

    // ===== 1. 收集 C / L / MOS，并建立 TR 状态 =====
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

    // 电容：上一时刻电压/电流（TR 需要记 i^n）
    struct CapTrapState {
        double vPrev = 0.0;  // V(n1)-V(n2) at t^n
        double iPrev = 0.0;  // i 从 n1->n2, at t^n
    };
    std::unordered_map<const CapacitorElement*, CapTrapState> capState;
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, xdc, n1);
        double v2 = getNodeVoltage(ckt, xdc, n2);
        CapTrapState st;
        st.vPrev = v1 - v2;
        st.iPrev = 0.0;
        capState[c.get()] = st;
    }

    // 电感：上一时刻电压/电流
    struct IndTrapState {
        double vPrev = 0.0;
        double iPrev = 0.0;
    };
    std::unordered_map<const Inductor*, IndTrapState> indState;
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        if (k < 0 || k >= N) continue;
        const auto& nodes = L->getNodeIds();
        int np = nodes[0];
        int nm = nodes[1];
        double vP = getNodeVoltage(ckt, xdc, np);
        double vM = getNodeVoltage(ckt, xdc, nm);
        IndTrapState st;
        st.vPrev = vP - vM;
        st.iPrev = xdc(k);
        indState[L.get()] = st;
    }

    // MOS 寄生电容：4 个电压/电流
    struct MosCapTrapState {
        double vgsPrev = 0.0, igsPrev = 0.0;
        double vgdPrev = 0.0, igdPrev = 0.0;
        double vsbPrev = 0.0, isbPrev = 0.0;
        double vdbPrev = 0.0, idbPrev = 0.0;
    };
    std::unordered_map<const MosfetBase*, MosCapTrapState> mosPrev;
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

        MosCapTrapState st;
        st.vgsPrev = vG - vS;
        st.vgdPrev = vG - vD;
        st.vsbPrev = vS - vB;
        st.vdbPrev = vD - vB;
        st.igsPrev = st.igdPrev = st.isbPrev = st.idbPrev = 0.0;
        mosPrev[m.get()] = st;
    }

    // ===== 2. 打开输出文件，写表头 =====
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Cannot open transient (TR) output file '" << outFile << "'.\n";
        return;
    }
    ofs << std::scientific << std::setprecision(9);

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
        if (t < tstart) return;

        ofs << t;
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) {
                ofs << "," << x(node.eqIndex);
            }
        }
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

    // ===== 3. 主时间步循环（第一步 BE，之后 TR；自适应 dt） =====
    std::cout << "[TRAN-TR] tstep=" << std::scientific << dtInit
              << ", tstop="        << tstop
              << ", tstart="       << tstart << "\n";

    int totalStepsInit = static_cast<int>(std::floor(tstop / dtInit + 1e-12));
    if (sim.verbose) {
        std::cout << "[TRAN-TR] total steps (based on initial dt) = " << totalStepsInit << "\n";
    }
    VectorXd x         = xdc;   // 当前接受解
    VectorXd xPrevStep = xdc;   // 上一个接受解（用于预测）
    double   dtPrev    = dtInit;

    dumpRow(0.0, xdc);

    ConvController ctrl;
    const double tol            = 1e-6;
    const int    maxNewtonIters = 60;

    // “一步求解”（BE/TR）封装：只负责 Newton + ConvController，不改状态表
    auto solveOneStep = [&](double tNow, double dtTry, bool useBE, VectorXd& xSolve) -> bool {
        double gmin    = ctrl.baseGmin(1.0);
        double alpha   = ctrl.initialAlphaLU();
        double prevErr = std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 非 C / 非 L / 非 MOS
            for (const auto& e : ckt.elements) {
                if (std::dynamic_pointer_cast<CapacitorElement>(e)) continue;
                if (std::dynamic_pointer_cast<Inductor>(e))        continue;
                if (std::dynamic_pointer_cast<MosfetBase>(e))      continue;
                e->stamp(G, I, ckt, xSolve, ctx);
            }

            // 2) MOS 导电
            for (const auto& m : mosfets) {
                m->stamp(G, I, ckt, xSolve, ctx);
            }

            if (useBE) {
                // ===== BE：电容/电感/MOS 电容 =====
                for (const auto& c : caps) {
                    double Cval = c->getC();
                    if (Cval <= 0.0) continue;

                    const auto& nodes = c->getNodeIds();
                    int n1 = nodes[0];
                    int n2 = nodes[1];
                    int eq1 = ckt.nodes[n1].eqIndex;
                    int eq2 = ckt.nodes[n2].eqIndex;

                    const CapTrapState& st = capState[c.get()];
                    stampCapBE(eq1, eq2, Cval, dtTry, st.vPrev, G, I);
                }

                // i_n+1 = dt/2L * v_n+1 + i_n + dt/2L * v_n
                for (const auto& L : inds) {
                    double Lval = L->getL();
                    if (Lval <= 0.0) continue;

                    const auto& nodes = L->getNodeIds();
                    int np  = nodes[0];
                    int nm  = nodes[1];
                    int eqP = ckt.nodes[np].eqIndex;
                    int eqM = ckt.nodes[nm].eqIndex;
                    int k   = L->getBranchEqIndex();
                    if (k < 0 || k >= N) continue;

                    const IndTrapState& st = indState[L.get()];

                    // TR: alpha = dt/(2L)
                    double alpha = dtTry / (2.0 * Lval);

                    // RHS = i_n + alpha * v_n
                    double rhs = st.iPrev + alpha * st.vPrev;

                    // 节点 KCL：保持你原来的方向约定（i 从 p -> m）
                    if (eqP >= 0) G(eqP, k) += 1.0;
                    if (eqM >= 0) G(eqM, k) -= 1.0;

                    // 支路方程（k 行）： i_{n+1} - alpha*(Vp - Vm) = rhs
                    if (eqP >= 0) G(k, eqP) += -alpha;
                    if (eqM >= 0) G(k, eqM) += +alpha;
                    G(k, k) += 1.0;

                    I(k) += rhs;
                }


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

                    const MosCapTrapState& st = mosPrev[m.get()];

                    stampCapBE(eqG, eqS, Cgs, dtTry, st.vgsPrev, G, I);
                    stampCapBE(eqG, eqD, Cgd, dtTry, st.vgdPrev, G, I);
                    stampCapBE(eqS, eqB, CsJ, dtTry, st.vsbPrev, G, I);
                    stampCapBE(eqD, eqB, CdJ, dtTry, st.vdbPrev, G, I);
                }
            } else {
                // ===== TR：电容/电感/MOS 电容 =====
                for (const auto& c : caps) {
                    double Cval = c->getC();
                    if (Cval <= 0.0) continue;

                    const auto& nodes = c->getNodeIds();
                    int n1 = nodes[0];
                    int n2 = nodes[1];
                    int eq1 = ckt.nodes[n1].eqIndex;
                    int eq2 = ckt.nodes[n2].eqIndex;

                    const CapTrapState& st = capState[c.get()];
                    stampCapTR(eq1, eq2, Cval, dtTry, st.vPrev, st.iPrev, G, I);
                }

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

                    const IndTrapState& st = indState[L.get()];

                    double G_eq = dtTry / (2.0 * Lval);  // 等效电导
                    double I_hist = st.iPrev + G_eq * st.vPrev;  // 历史电流源

                    // 在节点方程中直接添加诺顿等效贡献
                    if (eqP >= 0) {
                        G(eqP, eqP) += G_eq;
                        I(eqP) += -I_hist;
                    }
                    if (eqM >= 0) {
                        G(eqM, eqM) += G_eq;
                        I(eqM) += I_hist;
                    }
                    if (eqP >= 0 && eqM >= 0) {
                        G(eqP, eqM) -= G_eq;
                        G(eqM, eqP) -= G_eq;
                    }
                    if (eqP >= 0) G(eqP, k) += 1.0;
                    if (eqM >= 0) G(eqM, k) -= 1.0;

                    G(k, k) += 1.0;
                }


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

                    const MosCapTrapState& st = mosPrev[m.get()];

                    stampCapTR(eqG, eqS, Cgs, dtTry, st.vgsPrev, st.igsPrev, G, I);
                    stampCapTR(eqG, eqD, Cgd, dtTry, st.vgdPrev, st.igdPrev, G, I);
                    stampCapTR(eqS, eqB, CsJ, dtTry, st.vsbPrev, st.isbPrev, G, I);
                    stampCapTR(eqD, eqB, CdJ, dtTry, st.vdbPrev, st.idbPrev, G, I);
                }
            }

            // gmin
            stampGlobalGmin(ckt, G, gmin);

            // 解线性方程
            VectorXd xRaw = Solver::solveLinearSystemLU(G, I);
            if (!xRaw.allFinite()) {
                gmin = std::min(gmin * 10.0, 1e-4);
                continue;
            }

            // ConvController 更新
            auto st = ctrl.update(
                xSolve, xRaw, prevErr, iter,
                alpha, gmin, 1.0,
                tol
            );

            xSolve   = st.xNext;
            alpha    = st.alphaNext;
            gmin     = st.gminNext;
            prevErr  = st.error;

            if (st.converged) return true;
        }

        return false;
    };

    auto passTruncation = [&](const VectorXd& xPred, const VectorXd& xSol) -> bool {
        for (const auto& node : ckt.nodes) {
            if (node.eqIndex < 0) continue;
            int k = node.eqIndex;
            if (k < 0 || k >= xSol.size() || k >= xPred.size()) continue;

            double v  = xSol(k);
            double vp = xPred(k);
            double av = std::abs(v);
            double ap = std::abs(vp);
            double vmax = (av > ap) ? av : ap;

            if (std::abs(v - vp) > bound) return false;
        }
        return true;
    };

    double tNow = 0.0;
    int    step = 0;

    while (tNow + 1e-15 < tstop) {
        bool accepted = false;
        double dtTry  = dt;

        const double dtRemain0 = tstop - tNow;
        const bool allowFinalBelowDtMin = (dtRemain0 < dtMin * (1.0 - 1e-15));

        while (!accepted) {
            if (!allowFinalBelowDtMin && dtTry < dtMin * (1.0 - 1e-15)) {
                throw std::runtime_error("Transient(TR): dt reduced below dtMin, abort.");
            }

            double dtRemain = tstop - tNow;
            if (dtTry > dtRemain) dtTry = dtRemain;

            const bool useBE = (step == 0); // 第一“接受步”用 BE

            VectorXd xBefore = x;

            // 预测：x_pred = x + (dtTry/dtPrev)*(x - xPrevStep)
            VectorXd xPred = xBefore;
            if (step > 0) {
                double ratio = (dtPrev > 0.0) ? (dtTry / dtPrev) : 1.0;
                xPred = xBefore + ratio * (xBefore - xPrevStep);
            }

            x = xPred;

            bool converged = solveOneStep(tNow + dtTry, dtTry, useBE, x);
            if (!converged) {
                if (!allowFinalBelowDtMin && dtTry <= dtMin * (1.0 + 1e-15)) {
                    throw std::runtime_error("Transient(TR): reached dtMin but still not converged, abort.");
                }
                dtTry *= 0.5;
                x = xBefore;
                continue;
            }

            // 第一“接受步”：只要收敛就接受；从第二步开始阈值判断
            if (step > 0) {
                if (!passTruncation(xPred, x)) {
                    if (!allowFinalBelowDtMin && dtTry <= dtMin * (1.0 + 1e-15)) {
                        throw std::runtime_error("Transient(TR): reached dtMin but truncation test still fails, abort.");
                    }
                    dtTry *= 0.5;
                    x = xBefore;
                    continue;
                }
            }

            // ===== 接受该步：更新状态 =====
            if (useBE) {
                // 用 BE 生成 TR 所需 iPrev（关键点：第一步后必须把 iPrev 初始化出来）
                for (auto& c : caps) {
                    CapTrapState& st = capState[c.get()];
                    const auto& nodes = c->getNodeIds();
                    int n1 = nodes[0];
                    int n2 = nodes[1];
                    double v1 = getNodeVoltage(ckt, x, n1);
                    double v2 = getNodeVoltage(ckt, x, n2);
                    double vNow = v1 - v2;

                    double Cval = c->getC();
                    double iNow = (Cval > 0.0) ? (Cval / dtTry) * (vNow - st.vPrev) : 0.0;

                    st.vPrev = vNow;
                    st.iPrev = iNow;
                }

                for (auto& L : inds) {
                    auto it = indState.find(L.get());
                    if (it == indState.end()) continue;

                    IndTrapState& st = it->second;

                    const auto& nodes = L->getNodeIds();
                    int np = nodes[0];
                    int nm = nodes[1];
                    double vP = getNodeVoltage(ckt, x, np);
                    double vM = getNodeVoltage(ckt, x, nm);
                    double vNow = vP - vM;

                    int k = L->getBranchEqIndex();
                    double iNow = (k >= 0 && k < x.size()) ? x(k) : st.iPrev;

                    st.vPrev = vNow;
                    st.iPrev = iNow;
                }

                for (auto& m : mosfets) {
                    MosCapTrapState& st = mosPrev[m.get()];
                    const auto& nodes = m->getNodeIds();
                    int nD = nodes[0];
                    int nG = nodes[1];
                    int nS = nodes[2];
                    int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

                    double vD = getNodeVoltage(ckt, x, nD);
                    double vG = getNodeVoltage(ckt, x, nG);
                    double vS = getNodeVoltage(ckt, x, nS);
                    double vB = getNodeVoltage(ckt, x, nB);

                    double Cj0 = m->getCj0();
                    double Cg0 = m->getCg0();
                    double Cgs = Cg0;
                    double Cgd = Cg0;
                    double CsJ = Cj0;
                    double CdJ = Cj0;

                    double vgsNow = vG - vS;
                    double igsNow = (Cgs > 0.0) ? (Cgs / dtTry) * (vgsNow - st.vgsPrev) : 0.0;
                    st.vgsPrev = vgsNow; st.igsPrev = igsNow;

                    double vgdNow = vG - vD;
                    double igdNow = (Cgd > 0.0) ? (Cgd / dtTry) * (vgdNow - st.vgdPrev) : 0.0;
                    st.vgdPrev = vgdNow; st.igdPrev = igdNow;

                    double vsbNow = vS - vB;
                    double isbNow = (CsJ > 0.0) ? (CsJ / dtTry) * (vsbNow - st.vsbPrev) : 0.0;
                    st.vsbPrev = vsbNow; st.isbPrev = isbNow;

                    double vdbNow = vD - vB;
                    double idbNow = (CdJ > 0.0) ? (CdJ / dtTry) * (vdbNow - st.vdbPrev) : 0.0;
                    st.vdbPrev = vdbNow; st.idbPrev = idbNow;
                }
            } else {
                // 正常 TR 更新
                for (auto& c : caps) {
                    CapTrapState& st = capState[c.get()];
                    const auto& nodes = c->getNodeIds();
                    int n1 = nodes[0];
                    int n2 = nodes[1];
                    double v1 = getNodeVoltage(ckt, x, n1);
                    double v2 = getNodeVoltage(ckt, x, n2);
                    double vNow = v1 - v2;

                    double Cval = c->getC();
                    double iNow = (Cval > 0.0) ? (2.0 * Cval / dtTry) * (vNow - st.vPrev) - st.iPrev : 0.0;

                    st.vPrev = vNow;
                    st.iPrev = iNow;
                }

                for (auto& L : inds) {
                    auto it = indState.find(L.get());
                    if (it == indState.end()) continue;

                    IndTrapState& st = it->second;

                    const auto& nodes = L->getNodeIds();
                    int np = nodes[0];
                    int nm = nodes[1];
                    double vP = getNodeVoltage(ckt, x, np);
                    double vM = getNodeVoltage(ckt, x, nm);
                    double vNow = vP - vM;

                    double Lval = L->getL();
                    double iNow = (Lval > 0.0) ? st.iPrev + (dtTry / (2.0 * Lval)) * (vNow + st.vPrev) : st.iPrev;

                    st.vPrev = vNow;
                    st.iPrev = iNow;
                }

                for (auto& m : mosfets) {
                    MosCapTrapState& st = mosPrev[m.get()];
                    const auto& nodes = m->getNodeIds();
                    int nD = nodes[0];
                    int nG = nodes[1];
                    int nS = nodes[2];
                    int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

                    double vD = getNodeVoltage(ckt, x, nD);
                    double vG = getNodeVoltage(ckt, x, nG);
                    double vS = getNodeVoltage(ckt, x, nS);
                    double vB = getNodeVoltage(ckt, x, nB);

                    double Cj0 = m->getCj0();
                    double Cg0 = m->getCg0();
                    double Cgs = Cg0;
                    double Cgd = Cg0;
                    double CsJ = Cj0;
                    double CdJ = Cj0;

                    double vgsNow = vG - vS;
                    double igsNow = (Cgs > 0.0) ? (2.0 * Cgs / dtTry) * (vgsNow - st.vgsPrev) - st.igsPrev : 0.0;
                    st.vgsPrev = vgsNow; st.igsPrev = igsNow;

                    double vgdNow = vG - vD;
                    double igdNow = (Cgd > 0.0) ? (2.0 * Cgd / dtTry) * (vgdNow - st.vgdPrev) - st.igdPrev : 0.0;
                    st.vgdPrev = vgdNow; st.igdPrev = igdNow;

                    double vsbNow = vS - vB;
                    double isbNow = (CsJ > 0.0) ? (2.0 * CsJ / dtTry) * (vsbNow - st.vsbPrev) - st.isbPrev : 0.0;
                    st.vsbPrev = vsbNow; st.isbPrev = isbNow;

                    double vdbNow = vD - vB;
                    double idbNow = (CdJ > 0.0) ? (2.0 * CdJ / dtTry) * (vdbNow - st.vdbPrev) - st.idbPrev : 0.0;
                    st.vdbPrev = vdbNow; st.idbPrev = idbNow;
                }
            }

            // 接受：推进时间、输出
            tNow += dtTry;
            dumpRow(tNow, x);

            // 维护预测历史
            xPrevStep = xBefore;
            dtPrev    = dtTry;

            // 当前 dt * 2
            dt = dtTry >= 0.5 *dtMax ? dt : dtTry * 2.0;

            ++step;
            accepted = true;
        }
    }
    if (sim.verbose) {
        std::cout << "Transient analysis (Trapezoidal + ConvController) finished. "
                << "Results written to '" << outFile << "'.\n";
    }
}

// ========= 后向欧拉单周期积分（用于周期性稳态分析） =========
VectorXd TransientAnalysis::integrateOnePeriodBE(
    const VectorXd& x0,
    double t0,
    double dt,
    double T,
    const std::function<void(double, const VectorXd&)>& dumpRow)
{
    const int N = ckt.numUnknowns();
    if (N <= 0) {
        throw std::runtime_error("Transient(BE-Period): circuit has no unknowns.");
    }
    if (x0.size() != N) {
        throw std::runtime_error("Transient(BE-Period): x0 size mismatch.");
    }
    if (dt <= 0.0 || T <= 0.0) {
        throw std::runtime_error("Transient(BE-Period): dt and T must be > 0.");
    }

    // 步数检查
    double stepsDouble = T / dt;
    if (!std::isfinite(stepsDouble) || stepsDouble <= 0.0) {
        throw std::runtime_error("Transient(BE-Period): invalid T/dt.");
    }
    int nSteps = static_cast<int>(std::floor(stepsDouble + 1e-12));
    if (nSteps <= 0) {
        throw std::runtime_error("Transient(BE-Period): T/dt too small, no steps.");
    }

    // ===== 1. 收集 C / L / MOS 元件，并用 x0 建立初始状态 =====
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

    // 显式电容：v^0 = V(n1)-V(n2)，从 x0 初始化
    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, x0, n1);
        double v2 = getNodeVoltage(ckt, x0, n2);
        capVprev[c.get()] = v1 - v2;
    }

    // 电感：i_L^0 = x0 里的支路电流
    std::unordered_map<const Inductor*, double> indIprev;
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        double i0 = 0.0;
        if (k >= 0 && k < x0.size()) {
            i0 = x0(k);
        }
        indIprev[L.get()] = i0;
    }

    // MOS 寄生电容状态：从 x0 初始化
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

    // ===== 2. 主时间步循环（完全照 runBackwardEuler，只是初值换成 x0，不写文件） =====
    const int    maxNewtonIters = 50;
    const double tol            = 1e-6;
    const double gmin           = 1e-6;
    const double alpha          = 0.45;

    MatrixXd G(N, N);
    VectorXd I(N);

    VectorXd x = x0;

    // t = t0 的一行（如果需要输出）
    if (dumpRow) {
        dumpRow(t0, x);
    }

    for (int step = 0; step < nSteps; ++step) {
        double tNow = t0 + (step + 1) * dt;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 非 C / 非 L / 非 MOS 元件
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

            // 3) 显式电容（后向欧拉）
            for (const auto& c : caps) {
                double Cval = c->getC();
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c.get()];
                stampCapBE(eq1, eq2, Cval, dt, vPrev, G, I);
            }

            // 4) 电感（后向欧拉 Thévenin）
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
                double G_eq = dt / Lval;
                double I_hist = iPrev;

                if (eqP >= 0) {
                    G(eqP, eqP) += G_eq;
                    I(eqP) += -I_hist;
                }
                if (eqM >= 0) {
                    G(eqM, eqM) += G_eq; 
                    I(eqM) += I_hist;
                }
                if (eqP >= 0 && eqM >= 0) {
                    G(eqP, eqM) -= G_eq;
                    G(eqM, eqP) -= G_eq;
                }
                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;
                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
            }

            // 5) MOS 寄生电容（用 Cj0 粗略近似）
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

                stampCapBE(eqG, eqS, Cgs, dt, stPrev.vgsPrev, G, I);
                stampCapBE(eqG, eqD, Cgd, dt, stPrev.vgdPrev, G, I);
                stampCapBE(eqS, eqB, CsJ, dt, stPrev.vsbPrev, G, I);
                stampCapBE(eqD, eqB, CdJ, dt, stPrev.vdbPrev, G, I);
            }

            // 6) gmin 到地
            stampGlobalGmin(ckt, G, gmin);

            // 7) 解线性方程 + 简单阻尼
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);
            if (!xNew.allFinite()) {
                throw std::runtime_error("Transient(BE-Period): LU produced NaN/Inf.");
            }

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

        // ===== 3. 步进成功：更新电容 / 电感 / MOS 状态 =====
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
