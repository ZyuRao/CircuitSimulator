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

// =======================================================================
// 积分一个周期，并可选地计算 Monodromy 矩阵 (Sensitivity Analysis)
// =======================================================================
VectorXd integrateOnePeriodPssBE(const Circuit& ckt,
                                 const SimulationConfig& sim,
                                 double dt, double T,
                                 const VectorXd& x0,
                                 MatrixXd* outMonodromy, // 如果非空，则计算灵敏度
                                 const std::function<void(double, const VectorXd&)>& dumpRow)
{
    const int N = ckt.numUnknowns();
    double stepsDouble = T / dt;
    int nSteps = (int)std::floor(stepsDouble + 1e-12);
    
    // 1. 初始化历史状态
    std::vector<std::shared_ptr<CapacitorElement>> caps;
    std::vector<std::shared_ptr<Inductor>> inds;
    std::vector<std::shared_ptr<MosfetBase>> mosfets;
    
    for (const auto& e : ckt.elements) {
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) caps.push_back(c);
        else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) inds.push_back(L);
        else if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) mosfets.push_back(m);
    }

    std::unordered_map<const CapacitorElement*, double> capVprev;
    for (auto& c : caps) {
        const auto& nodes = c->getNodeIds();
        capVprev[c.get()] = getNodeVoltage(ckt, x0, nodes[0]) - getNodeVoltage(ckt, x0, nodes[1]);
    }

    std::unordered_map<const Inductor*, double> indIprev;
    for (auto& L : inds) {
        int k = L->getBranchEqIndex();
        indIprev[L.get()] = (k >= 0 && k < x0.size()) ? x0(k) : 0.0;
    }

    std::unordered_map<const MosfetBase*, MosCapState> mosPrev;
    for (auto& m : mosfets) {
        const auto& nodes = m->getNodeIds();
        double vD = getNodeVoltage(ckt, x0, nodes[0]);
        double vG = getNodeVoltage(ckt, x0, nodes[1]);
        double vS = getNodeVoltage(ckt, x0, nodes[2]);
        double vB = getNodeVoltage(ckt, x0, (nodes.size()>3 ? nodes[3] : nodes[2]));
        mosPrev[m.get()] = {vG - vS, vG - vD, vS - vB, vD - vB};
    }

    VectorXd x = x0;
    MatrixXd G(N, N);
    VectorXd I(N);
    MatrixXd B(N, N); // 用于灵敏度分析的 RHS 矩阵 (C/dt)

    // 初始化灵敏度矩阵 S 为单位阵 (S_0 = I)
    MatrixXd S;
    if (outMonodromy) {
        S = MatrixXd::Identity(N, N);
    }

    const int maxNewtonIters = 50;
    const double tol = 1e-6;
    double gmin = 1e-9;
    
    if (dumpRow) dumpRow(0.0, x);

    // === 时间步循环 ===
    for (int step = 0; step < nSteps; ++step) {
        double tNow = (step + 1) * dt;
        
        // --- 1. Newton 求解当前步 x_new ---
        // 为了确保灵敏度矩阵准确，这里不使用激进的 damp，而是标准 Newton
        // 或者保留原有的 damped 逻辑，但收敛后 G 必须是准确的 Jacobian
        
        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            G.setZero();
            I.setZero();
            
            AnalysisContext ctx;
            ctx.type = AnalysisType::TRAN;
            ctx.time = tNow;
            
            // Stamp 静态元件
            for (const auto& e : ckt.elements) {
                if (!std::dynamic_pointer_cast<CapacitorElement>(e) && 
                    !std::dynamic_pointer_cast<Inductor>(e) &&
                    !std::dynamic_pointer_cast<MosfetBase>(e)) {
                    e->stamp(G, I, ckt, x, ctx);
                }
            }
            // Stamp MOS 导电
            for (const auto& m : mosfets) m->stamp(G, I, ckt, x, ctx);
            
            // Stamp 动态元件 (Cap/Ind/MosCap)
            for (const auto& c : caps) {
                double vPrev = capVprev[c.get()];
                const auto& nodes = c->getNodeIds();
                pssStampCapBE(ckt.nodes[nodes[0]].eqIndex, ckt.nodes[nodes[1]].eqIndex, 
                              c->getC(), dt, vPrev, G, I);
            }
            for (const auto& L : inds) {
                double iPrev = indIprev[L.get()];
                const auto& nodes = L->getNodeIds();
                pssStampIndBE(ckt.nodes[nodes[0]].eqIndex, ckt.nodes[nodes[1]].eqIndex,
                              L->getBranchEqIndex(), L->getL(), dt, iPrev, G, I);
            }
            for (const auto& m : mosfets) {
                const auto& st = mosPrev[m.get()];
                const auto& nodes = m->getNodeIds();
                int eqD = ckt.nodes[nodes[0]].eqIndex;
                int eqG = ckt.nodes[nodes[1]].eqIndex;
                int eqS = ckt.nodes[nodes[2]].eqIndex;
                int eqB = ckt.nodes[(nodes.size()>3?nodes[3]:nodes[2])].eqIndex;
                
                double Cj0 = m->getCj0(), Cg0 = m->getCg0();
                pssStampCapBE(eqG, eqS, Cg0, dt, st.vgsPrev, G, I);
                pssStampCapBE(eqG, eqD, Cg0, dt, st.vgdPrev, G, I);
                pssStampCapBE(eqS, eqB, Cj0, dt, st.vsbPrev, G, I);
                pssStampCapBE(eqD, eqB, Cj0, dt, st.vdbPrev, G, I);
            }

            stampGlobalGmin(ckt, G, gmin);

            VectorXd dx = Solver::solveLinearSystemLU(G, I); // I 是残差 (G*x_old - I_src) ? 
            // 注意：Element stamp 约定 I向量存的是源项。Solver 需要解 G*x = I。
            // 实际上这里的 G 已经是 Jacobian 吗？
            // 原有的 element.cpp stamp 是直接填充 G 和 RHS I。
            // 对于非线性器件，G 是线性化 conductance。
            // 所以 solve 出来的是 x_new，而不是 dx。
            // Newton Update: x_new = G^-1 * I
            
            if (!dx.allFinite()) {
                gmin *= 10; continue;
            }

            VectorXd xOld = x;
            x = x + 0.8 * (dx - x); // 简单阻尼
            
            if ((x - xOld).norm() < tol) break;
        }

        // --- 2. 灵敏度传递 (Sensitivity Transfer) ---
        if (outMonodromy) {
            // 方程: d/dt(q(x)) + f(x) = u(t)
            // BE:   (q(x_n) - q(x_{n-1}))/h + f(x_n) = u_n
            // Diff w.r.t x_0:
            //       (1/h * C_n + G_n) * S_n = (1/h * C_n) * S_{n-1}
            // LHS 矩阵 (1/h*C + G) 正是刚刚收敛时的 G (Newton Jacobian)
            // RHS 矩阵需要重新组装只包含 C/h 的部分
            
            B.setZero();
            stampDynamicMatrix(ckt, B, dt);
            
            // 计算 RHS: B * S_{n-1}
            MatrixXd RHS = B * S;
            
            // 解线性方程: G * S_n = RHS
            // 注意：这里 G 必须是最后一次迭代收敛时的 Jacobian。
            // 此时 G 已经包含了 stampDynamicMatrix 的内容(作为 LHS)。
            
            // 使用 LU 分解求解多右端项
            // 自行实现或调用 solver (Solver::solveLinearSystemLU 通常只解 Vector)
            // 这里为了简单，按列求解
            for (int col = 0; col < N; ++col) {
                VectorXd rhsVec = RHS.col(col);
                VectorXd sCol = Solver::solveLinearSystemLU(G, rhsVec);
                S.col(col) = sCol;
            }
        }

        // --- 3. 更新历史状态 ---
        for (const auto& c : caps) {
            const auto& nodes = c->getNodeIds();
            capVprev[c.get()] = getNodeVoltage(ckt, x, nodes[0]) - getNodeVoltage(ckt, x, nodes[1]);
        }
        for (const auto& L : inds) {
            int k = L->getBranchEqIndex();
            indIprev[L.get()] = (k>=0 && k<x.size()) ? x(k) : 0.0;
        }
        for (const auto& m : mosfets) {
            const auto& nodes = m->getNodeIds();
            double vD = getNodeVoltage(ckt, x, nodes[0]);
            double vG = getNodeVoltage(ckt, x, nodes[1]);
            double vS = getNodeVoltage(ckt, x, nodes[2]);
            double vB = getNodeVoltage(ckt, x, (nodes.size()>3?nodes[3]:nodes[2]));
            mosPrev[m.get()] = {vG-vS, vG-vD, vS-vB, vD-vB};
        }
        
        if (dumpRow) dumpRow(tNow, x);
    }

    if (outMonodromy) {
        *outMonodromy = S;
    }
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
    double dt = (cfg.tstep <= 0.0) ? cfg.periodT / 200.0 : cfg.tstep;
    int N = ckt.numUnknowns();
    if (N <= 0) return;

    // 1. DC 初始化
    DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
    VectorXd x0 = dc.run();
    
    std::cout << "PSS-SHOOT: Sensitivity + Broyden Method\n"
              << "  Period=" << cfg.periodT << ", dt=" << dt << "\n";

    VectorXd xInit = x0;
    MatrixXd J_shoot(N, N); // Shooting 雅可比矩阵 M - I
    VectorXd F_curr, F_prev;
    VectorXd dx_prev;
    
    bool useBroyden = false;

    for (int iter = 0; iter < cfg.maxIters; ++iter) {
        VectorXd xT;
        MatrixXd M; // Monodromy Matrix
        
        // 策略: 
        // 第一次迭代 (iter=0) 或 Broyden 失败重启时，计算精确灵敏度
        // 其他时候只积分，不计算灵敏度，用 Broyden 更新 J_shoot
        bool computeSens = (!useBroyden);

        try {
            xT = integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, 
                                         computeSens ? &M : nullptr, nullptr);
        } catch (const std::exception& e) {
            std::cerr << "Integration failed: " << e.what() << "\n";
            return;
        }

        // 残差 F(x) = x(T) - x(0)
        F_curr = xT - xInit;
        double err = F_curr.norm();
        std::cout << "Iter " << iter << " ||F|| = " << std::scientific << err << "\n";

        if (err < cfg.tol) {
            std::cout << "Converged.\n";
            break;
        }

        // --- Jacobian Update ---
        if (computeSens) {
            // 使用精确的 Monodromy 矩阵
            // J_shoot = M - I
            J_shoot = M - MatrixXd::Identity(N, N);
            
            // 下次尝试使用 Broyden 以加速
            useBroyden = true;
        } else {
            // Broyden Update (Rank-1 update)
            // J_k = J_{k-1} + (dF - J_{k-1}*dx) * dx^T / (dx^T * dx)
            // dF = F_curr - F_prev
            // dx = dx_prev (即 xInit_curr - xInit_prev)
            
            VectorXd dF = F_curr - F_prev;
            VectorXd Jdx = J_shoot * dx_prev;
            double dxNorm2 = dx_prev.squaredNorm();
            
            if (dxNorm2 > 1e-20) {
                J_shoot += (dF - Jdx) * dx_prev.transpose() / dxNorm2;
            } else {
                // 如果步长太小，重置为灵敏度分析
                useBroyden = false;
                std::cout << "  (Broyden step too small, resetting to Sensitivity)\n";
                // 回退本次迭代，或者仅仅下一次强制计算灵敏度
                // 这里选择简单策略：继续求解，但标记下一次强制计算
                // 其实最好现在就重算，但为了代码简单，本次用旧 J 继续
            }
        }

        // --- Newton Step ---
        // solve J * dx = -F
        VectorXd dx = Solver::solveLinearSystemLU(J_shoot, -F_curr);
        
        if (!dx.allFinite()) {
            std::cout << "  (Linear solve failed, resetting)\n";
            xInit = x0; // 重置
            useBroyden = false;
            continue;
        }

        // 阻尼更新
        double alpha = cfg.relax;
        xInit += alpha * dx;
        
        // 保存用于下一次 Broyden 的状态
        F_prev = F_curr;
        dx_prev = alpha * dx;
    }

    // 最终输出波形
    std::ofstream ofs(outFile);
    // ... (保留原有的输出代码) ...
    
    // Header
    ofs << "time";
    for (const auto& node : ckt.nodes) if (node.eqIndex >= 0) ofs << ",V(" << node.name << ")";
    // ... Sources/Inductors ...
    ofs << "\n";

    auto dumpRow = [&](double t, const VectorXd& x) {
        ofs << t;
        for (const auto& node : ckt.nodes) 
            if (node.eqIndex >= 0 && node.eqIndex < x.size()) ofs << "," << x(node.eqIndex);
        // ... (其他电流输出逻辑同原文件)
        ofs << "\n";
    };

    dumpRow(0.0, xInit);
    integrateOnePeriodPssBE(ckt, sim, dt, cfg.periodT, xInit, nullptr, dumpRow);
    std::cout << "Results written to " << outFile << "\n";
}
