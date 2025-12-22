
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
        xdc = computeDcOperatingPoint();
    } catch (const std::exception& e) {
        std::cerr << "DC operating point failed: " << e.what() << "\n";
        return;
    }

    if (xdc.size() != N) {
        std::cerr << "Transient: DC solution size mismatch.\n";
        return;
    }

    const auto& cache = ckt.elementCache();
    const auto& caps = cache.capacitors;
    const auto& inds = cache.inductors;
    const auto& mosfets = cache.mos;

    // 显式电容：v^{0} = V(n1) - V(n2)
    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (const auto* c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, xdc, n1);
        double v2 = getNodeVoltage(ckt, xdc, n2);
        capVprev[c] = v1 - v2;
    }

    // 电感：i_L^{0} = DC 解中的支路电流
    std::unordered_map<const Inductor*, double> indIprev;
    for (const auto* L : inds) {
        int k = L->getBranchEqIndex();
        double i0 = 0.0;
        if (k >= 0 && k < xdc.size()) {
            i0 = xdc(k);
        }
        indIprev[L] = i0;
    }

    // MOS 寄生电容状态：从 DC 工作点初始化
    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    for (const auto* m : mosfets) {
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
        mosPrev[m] = st;
    }

    // ===== 2. 打开输出文件，写表头 =====
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Cannot open transient output file '" << outFile << "'.\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    struct BranchElemRef {
        const VoltageSource* vs = nullptr;
        const Inductor* ind = nullptr;
    };
    std::vector<BranchElemRef> branchElems;
    branchElems.reserve(cache.voltageSources.size() + cache.inductors.size());
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            branchElems.push_back(BranchElemRef{vs.get(), nullptr});
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            branchElems.push_back(BranchElemRef{nullptr, L.get()});
        }
    }

    ofs << "time";
    // 节点电压
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    // 电压源与电感电流
    for (const auto& b : branchElems) {
        if (b.vs) {
            ofs << ",I(" << b.vs->getName() << ")";
        } else if (b.ind) {
            ofs << ",I(" << b.ind->getName() << ")";
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
        for (const auto& b : branchElems) {
            int k = -1;
            if (b.vs) {
                k = b.vs->getBranchEqIndex();
            } else if (b.ind) {
                k = b.ind->getBranchEqIndex();
            }
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            ofs << "," << Ibr;
        }
        ofs << "\n";
    };

    // ===== 3. 时间步循环：Newton + LU + 后向欧拉 =====
    if (sim.verbose) {
        std::cout << "[TRAN] tstep="  << std::scientific << cfg.tstep
                  << ", tstop="       << cfg.tstop
                  << ", tstart="      << cfg.tstart << "\n";
    }

    int nSteps = static_cast<int>(std::floor(tstop / dt + 1e-12));
    if (sim.verbose) {
        std::cout << "[TRAN] total steps = " << nSteps << "\n";
    }

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
            for (const auto* r : cache.resistors) r->stamp(G, I, ckt, x, ctx);
            for (const auto* vs : cache.voltageSources) vs->stamp(G, I, ckt, x, ctx);
            for (const auto* is : cache.currentSources) is->stamp(G, I, ckt, x, ctx);

            // 2) 线性 MOS 导电部分（Ids、gm、gds 等）
            for (const auto* m : mosfets) {
                m->stamp(G, I, ckt, x, ctx);
            }

            // 3) 显式电容的后向欧拉伴随模型
            for (const auto* c : caps) {
                double Cval = c->getC();
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c]; // V(n1) - V(n2) at previous step
                stampCapBE(eq1, eq2, Cval, dt, vPrev, G, I);
            }

            // 4) 电感的后向欧拉 Thévenin 伴随模型
            for (const auto* L : inds) {
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
                double iPrev = indIprev[L];
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
            for (const auto* m : mosfets) {
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

                const MosCapState& stPrev = mosPrev[m];

                if (Cgs > 0.0) stampCapBE(eqG, eqS, Cgs, dt, stPrev.vgsPrev, G, I);
                if (Cgd > 0.0) stampCapBE(eqG, eqD, Cgd, dt, stPrev.vgdPrev, G, I);
                if (CsJ > 0.0) stampCapBE(eqS, eqB, CsJ, dt, stPrev.vsbPrev, G, I);
                if (CdJ > 0.0) stampCapBE(eqD, eqB, CdJ, dt, stPrev.vdbPrev, G, I);
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
        for (const auto* c : caps) {
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            capVprev[c] = v1 - v2;
        }
        // 电感
        for (const auto* L : inds) {
            int k = L->getBranchEqIndex();
            double iL = 0.0;
            if (k >= 0 && k < x.size()) {
                iL = x(k);
            }
            indIprev[L] = iL;
        }
        // MOS 寄生电容状态
        for (const auto* m : mosfets) {
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
            mosPrev[m] = st;
        }

        dumpRow(tNow, x);
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

    const double dt     = cfg.tstep;
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

    const auto& cache = ckt.elementCache();
    const auto& caps = cache.capacitors;
    const auto& inds = cache.inductors;
    const auto& mosfets = cache.mos;

    // 电容：上一时刻电压/电流（TR 需要记 i^n）
    struct CapTrapState {
        double vPrev = 0.0;  // V(n1)-V(n2) at t^n
        double iPrev = 0.0;  // i 从 n1->n2, at t^n
    };
    std::unordered_map<const CapacitorElement*, CapTrapState> capState;
    for (const auto* c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, xdc, n1);
        double v2 = getNodeVoltage(ckt, xdc, n2);
        CapTrapState st;
        st.vPrev = v1 - v2;
        st.iPrev = 0.0;  // DC 下 dv/dt = 0，i 取 0
        capState[c] = st;
    }

    // 电感：上一时刻电压/电流
    struct IndTrapState {
        double vPrev = 0.0;  // V(np)-V(nm) at t^n
        double iPrev = 0.0;  // I_L at t^n
    };
    std::unordered_map<const Inductor*, IndTrapState> indState;
    for (const auto* L : inds) {
        int k = L->getBranchEqIndex();
        if (k < 0 || k >= N) continue;
        const auto& nodes = L->getNodeIds();
        int np = nodes[0];
        int nm = nodes[1];
        double vP = getNodeVoltage(ckt, xdc, np);
        double vM = getNodeVoltage(ckt, xdc, nm);
        IndTrapState st;
        st.vPrev = vP - vM;
        st.iPrev = xdc(k);   // DC 下的电感电流
        indState[L] = st;
    }

    // MOS 寄生电容：4 个电压/电流
    struct MosCapTrapState {
        double vgsPrev = 0.0, igsPrev = 0.0;
        double vgdPrev = 0.0, igdPrev = 0.0;
        double vsbPrev = 0.0, isbPrev = 0.0;
        double vdbPrev = 0.0, idbPrev = 0.0;
    };
    std::unordered_map<const MosfetBase*, MosCapTrapState> mosPrev;
    for (const auto* m : mosfets) {
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
        // DC 下假定 C 电流为 0
        st.igsPrev = st.igdPrev = st.isbPrev = st.idbPrev = 0.0;
        mosPrev[m] = st;
    }

    // ===== 2. 打开输出文件，写表头 =====
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "Cannot open transient (TR) output file '" << outFile << "'.\n";
        return;
    }
    ofs << std::scientific << std::setprecision(9);

    struct BranchElemRef {
        const VoltageSource* vs = nullptr;
        const Inductor* ind = nullptr;
    };
    std::vector<BranchElemRef> branchElems;
    branchElems.reserve(cache.voltageSources.size() + cache.inductors.size());
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            branchElems.push_back(BranchElemRef{vs.get(), nullptr});
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            branchElems.push_back(BranchElemRef{nullptr, L.get()});
        }
    }

    ofs << "time";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            ofs << ",V(" << node.name << ")";
        }
    }
    for (const auto& b : branchElems) {
        if (b.vs) {
            ofs << ",I(" << b.vs->getName() << ")";
        } else if (b.ind) {
            ofs << ",I(" << b.ind->getName() << ")";
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
        for (const auto& b : branchElems) {
            int k = -1;
            if (b.vs) {
                k = b.vs->getBranchEqIndex();
            } else if (b.ind) {
                k = b.ind->getBranchEqIndex();
            }
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            ofs << "," << Ibr;
        }
        ofs << "\n";
    };

    // ===== 3. 主时间步循环（Trapezoidal + ConvController） =====
    if (sim.verbose) {
        std::cout << "[TRAN-TR] tstep=" << std::scientific << dt
                  << ", tstop="        << tstop
                  << ", tstart="       << tstart << "\n";
    }

    int totalSteps = static_cast<int>(std::floor(tstop / dt + 1e-12));
    if (sim.verbose) {
        std::cout << "[TRAN-TR] total steps = " << totalSteps << "\n";
    }

    VectorXd x         = xdc;   // 当前解
    VectorXd xPrevStep = xdc;   // 上一时间步的解（做线性预测用）

    dumpRow(0.0, xdc);

    ConvController ctrl;           // 复用 DC 里的收敛控制器
    const double tol          = 1e-6;
    const int    maxNewtonIters = 60;
    MatrixXd G(N, N);
    VectorXd I(N);
    VectorXd xRaw(N);

    for (int step = 0; step < totalSteps; ++step) {
        double tNow = (step + 1) * dt;

        // 简单线性预测提高初始点质量：x_pred = 2*x - xPrevStep
        if (step > 0) {
            VectorXd xPred = x + (x - xPrevStep);
            xPrevStep = x;
            x = xPred;
        } else {
            xPrevStep = x;
        }

        double gmin    = ctrl.baseGmin(1.0);  // 这里 rampScale 固定为 1.0
        double alpha   = ctrl.initialAlphaLU();
        double prevErr = std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();

            AnalysisContext ctx;
            ctx.type        = AnalysisType::TRAN;
            ctx.sourceScale = 1.0;
            ctx.time        = tNow;
            ctx.omega       = 0.0;

            // 1) 非 C / 非 L / 非 MOS 元件（含电压源、独立源等）
            for (const auto* r : cache.resistors) r->stamp(G, I, ckt, x, ctx);
            for (const auto* vs : cache.voltageSources) vs->stamp(G, I, ckt, x, ctx);
            for (const auto* is : cache.currentSources) is->stamp(G, I, ckt, x, ctx);

            // 2) MOS 导电部分
            for (const auto* m : mosfets) {
                m->stamp(G, I, ckt, x, ctx);
            }

            // 3) 显式电容 TR 伴随（Norton 等效）
            for (const auto* c : caps) {
                double Cval = c->getC();
                if (Cval <= 0.0) continue;

                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;

                const CapTrapState& st = capState[c];
                stampCapTR(eq1, eq2, Cval, dt, st.vPrev, st.iPrev, G, I);
            }

            // 4) 电感 TR 伴随（Thevenin 等效）
            for (const auto* L : inds) {
                double Lval = L->getL();
                if (Lval <= 0.0) continue;

                const auto& nodes = L->getNodeIds();
                int np = nodes[0];
                int nm = nodes[1];
                int eqP = ckt.nodes[np].eqIndex;
                int eqM = ckt.nodes[nm].eqIndex;
                int k   = L->getBranchEqIndex();
                if (k < 0 || k >= N) continue;

                const IndTrapState& st = indState[L];
                double vPrev = st.vPrev;  // V^n
                double iPrev = st.iPrev;  // i^n

                double R_eq  = 2.0 * Lval / dt;
                // 正确的历史等效：v^{n+1} - R_eq * i^{n+1} = -R_eq * i^n - v^n
                double V_hist = -R_eq * iPrev - vPrev;

                // KCL：节点电流 = I_L
                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;

                // 支路方程：Vp - Vm - R_eq * I_L = V_hist
                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
                G(k, k) += -R_eq;
                I(k)    += V_hist;
            }

            // 5) MOS 寄生电容 TR
            for (const auto* m : mosfets) {
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

                const MosCapTrapState& st = mosPrev[m];

                // Gate-source
                stampCapTR(eqG, eqS, Cgs, dt, st.vgsPrev, st.igsPrev, G, I);
                // Gate-drain
                stampCapTR(eqG, eqD, Cgd, dt, st.vgdPrev, st.igdPrev, G, I);
                // Source-bulk
                stampCapTR(eqS, eqB, CsJ, dt, st.vsbPrev, st.isbPrev, G, I);
                // Drain-bulk
                stampCapTR(eqD, eqB, CdJ, dt, st.vdbPrev, st.idbPrev, G, I);
            }

            // 6) gmin 到地（用 ConvController 控的 gmin）
            stampGlobalGmin(ckt, G, gmin);

            // 7) 解线性方程，使用 LU
            xRaw = Solver::solveLinearSystemLU(G, I);
            if (!xRaw.allFinite()) {
                // 如果有 NaN/Inf，适当增大 gmin 再试一轮
                gmin = std::min(gmin * 10.0, 1e-4);
                continue;
            }

            // 8) 用 ConvController 更新 alpha / gmin 和 x（和 DC 一样）
             auto st = ctrl.update(
                x, xRaw, prevErr, iter,
                alpha, gmin, 1.0,  // rampScale = 1.0
                tol
            );

            x       = st.xNext;
            alpha   = st.alphaNext;
            gmin    = st.gminNext;
            prevErr = st.error;

            if (st.converged) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: transient (TR) Newton did not converge at t="
                          << std::scientific << tNow
                          << " (err=" << st.error
                          << ", alpha=" << alpha
                          << ", gmin=" << gmin << ")\n";
            }
        }

        // ===== 4. 时间步收敛后：更新所有 TR 状态 & 输出 =====
        // 显式电容：i^{n+1} = (2C/dt)*(v^{n+1}-v^n) - i^n
        for (const auto* c : caps) {
            CapTrapState& st = capState[c];
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            double vNow = v1 - v2;
            double Cval = c->getC();
            double iNow = (2.0 * Cval / dt) * (vNow - st.vPrev) - st.iPrev;
            st.vPrev = vNow;
            st.iPrev = iNow;
        }

        // 电感：i^{n+1} = i^n + (dt/(2L))*(v^{n+1}+v^n)
        for (const auto* L : inds) {
            IndTrapState& st = indState[L];
            const auto& nodes = L->getNodeIds();
            int np = nodes[0];
            int nm = nodes[1];
            double vP = getNodeVoltage(ckt, x, np);
            double vM = getNodeVoltage(ckt, x, nm);
            double vNow = vP - vM;
            double Lval = L->getL();
            double iNow = st.iPrev + (dt / (2.0 * Lval)) * (vNow + st.vPrev);
            st.vPrev = vNow;
            st.iPrev = iNow;
        }

        // MOS 寄生电容电压/电流
        for (const auto* m : mosfets) {
            MosCapTrapState& st = mosPrev[m];
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
            double igsNow = (2.0 * Cgs / dt) * (vgsNow - st.vgsPrev) - st.igsPrev;
            st.vgsPrev = vgsNow;
            st.igsPrev = igsNow;

            double vgdNow = vG - vD;
            double igdNow = (2.0 * Cgd / dt) * (vgdNow - st.vgdPrev) - st.igdPrev;
            st.vgdPrev = vgdNow;
            st.igdPrev = igdNow;

            double vsbNow = vS - vB;
            double isbNow = (2.0 * CsJ / dt) * (vsbNow - st.vsbPrev) - st.isbPrev;
            st.vsbPrev = vsbNow;
            st.isbPrev = isbNow;

            double vdbNow = vD - vB;
            double idbNow = (2.0 * CdJ / dt) * (vdbNow - st.vdbPrev) - st.idbPrev;
            st.vdbPrev = vdbNow;
            st.idbPrev = idbNow;
        }

        dumpRow(tNow, x);
    }

    if (sim.verbose) {
        std::cout << "Transient analysis (Trapezoidal + ConvController) finished. "
                  << "Results written to '" << outFile << "'.\n";
    }
}

// ========= 后向欧拉单周期积分（用于周期性稳态分析） =========
// 从给定初值 x0 出发，用后向欧拉在 [t0, t0+T] 上积分，
// 完全复用 runBackwardEuler 的 stamp 逻辑（C/L/MOS + gmin）。
// dumpRow 可选，用于在每个时间点输出波形；返回末态 x(t0+T)。
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

    const auto& cache = ckt.elementCache();
    const auto& caps = cache.capacitors;
    const auto& inds = cache.inductors;
    const auto& mosfets = cache.mos;

    // 显式电容：v^0 = V(n1)-V(n2)，从 x0 初始化
    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (const auto* c : caps) {
        const auto& nodes = c->getNodeIds();
        int n1 = nodes[0];
        int n2 = nodes[1];
        double v1 = getNodeVoltage(ckt, x0, n1);
        double v2 = getNodeVoltage(ckt, x0, n2);
        capVprev[c] = v1 - v2;
    }

    // 电感：i_L^0 = x0 里的支路电流
    std::unordered_map<const Inductor*, double> indIprev;
    for (const auto* L : inds) {
        int k = L->getBranchEqIndex();
        double i0 = 0.0;
        if (k >= 0 && k < x0.size()) {
            i0 = x0(k);
        }
        indIprev[L] = i0;
    }

    // MOS 寄生电容状态：从 x0 初始化
    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    for (const auto* m : mosfets) {
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
        mosPrev[m] = st;
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
            for (const auto* r : cache.resistors) r->stamp(G, I, ckt, x, ctx);
            for (const auto* vs : cache.voltageSources) vs->stamp(G, I, ckt, x, ctx);
            for (const auto* is : cache.currentSources) is->stamp(G, I, ckt, x, ctx);

            // 2) MOS 导电部分
            for (const auto* m : mosfets) {
                m->stamp(G, I, ckt, x, ctx);
            }

            // 3) 显式电容（后向欧拉）
            for (const auto* c : caps) {
                double Cval = c->getC();
                const auto& nodes = c->getNodeIds();
                int n1 = nodes[0];
                int n2 = nodes[1];
                int eq1 = ckt.nodes[n1].eqIndex;
                int eq2 = ckt.nodes[n2].eqIndex;
                double vPrev = capVprev[c];
                stampCapBE(eq1, eq2, Cval, dt, vPrev, G, I);
            }

            // 4) 电感（后向欧拉 Thévenin）
            for (const auto* L : inds) {
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
                double iPrev = indIprev[L];
                double V_hist = -R_eq * iPrev;

                if (eqP >= 0) G(eqP, k) += 1.0;
                if (eqM >= 0) G(eqM, k) -= 1.0;

                if (eqP >= 0) G(k, eqP) += 1.0;
                if (eqM >= 0) G(k, eqM) -= 1.0;
                G(k, k) += -R_eq;
                I(k)    += V_hist;
            }

            // 5) MOS 寄生电容（用 Cj0 粗略近似）
            for (const auto* m : mosfets) {
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

                const MosCapState& stPrev = mosPrev[m];

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
        for (const auto* c : caps) {
            const auto& nodes = c->getNodeIds();
            int n1 = nodes[0];
            int n2 = nodes[1];
            double v1 = getNodeVoltage(ckt, x, n1);
            double v2 = getNodeVoltage(ckt, x, n2);
            capVprev[c] = v1 - v2;
        }
        for (const auto* L : inds) {
            int k = L->getBranchEqIndex();
            double iL = 0.0;
            if (k >= 0 && k < x.size()) {
                iL = x(k);
            }
            indIprev[L] = iL;
        }
        for (const auto* m : mosfets) {
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
            mosPrev[m] = st;
        }

        if (dumpRow) {
            dumpRow(tNow, x);
        }
    }

    return x;
}
