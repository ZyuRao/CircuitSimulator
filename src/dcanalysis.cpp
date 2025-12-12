#include "analysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "solver.hpp"
#include "sim.hpp"
#include "utils.hpp"

#include <iostream>
#include <memory>
#include <limits>
#include <algorithm>

using Eigen::MatrixXd;
using Eigen::VectorXd;




DcAnalysis::DcAnalysis(const Circuit& ckt_,
                       const SimulationConfig& sim_,
                       DcSolverKind solver_)
    : ckt(ckt_), sim(sim_), solverKind(solver_) {}

AnalysisContext DcAnalysis::makeDcCtx(double sourceScale) {
    AnalysisContext ctx;
    ctx.type        = AnalysisType::DC;
    ctx.sourceScale = sourceScale;
    ctx.time        = 0.0;
    ctx.omega       = 0.0;
    return ctx;
}

Eigen::VectorXd DcAnalysis::run() {
    // 是否包含非线性器件（MOS 等）
    bool nonlinear = hasNonlinearDevices(ckt);

    switch (solverKind) {
    case DcSolverKind::LU:
        return nonlinear ? solveNewtonLU() : solveDirectLU();
    case DcSolverKind::GaussSeidel:
    default:
        return nonlinear ? solveNewtonGS() : solveDirectGS();
    }
}

Eigen::VectorXd DcAnalysis::runLU() {
    bool nonlinear = hasNonlinearDevices(ckt);
    return nonlinear ? solveNewtonLU() : solveDirectLU();
}

Eigen::VectorXd DcAnalysis::runGaussSeidel() {
    bool nonlinear = hasNonlinearDevices(ckt);
    return nonlinear ? solveNewtonGS() : solveDirectGS();
}

// 只处理线性电路，线性方程 G x = I，用 LU 直接解
VectorXd DcAnalysis::solveDirectLU() const {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(N > 0 ? N : 1);

    if (N == 0) {
        std::cerr << "DC solve (LU): no unknowns.\n";
        return x;
    }

    MatrixXd G = MatrixXd::Zero(N, N);
    VectorXd I = VectorXd::Zero(N);

    // 线性电路，不需要迭代；sourceScale = 1.0
    AnalysisContext ctx = makeDcCtx(1.0);

    for (const auto& e : ckt.elements) {
        e->stamp(G, I, ckt, x, ctx);
    }

    // 使用自写 LU 分解
    x = Solver::solveLinearSystemLU(G, I);
    return x;
}

// 只处理线性电路，线性方程 G x = I，用 Gauss-Seidel 迭代
VectorXd DcAnalysis::solveDirectGS() const {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(N > 0 ? N : 1);

    if (N == 0) {
        std::cerr << "DC solve (GS): no unknowns.\n";
        return x;
    }

    MatrixXd G = MatrixXd::Zero(N, N);
    VectorXd I = VectorXd::Zero(N);

    AnalysisContext ctx = makeDcCtx(1.0);

    for (const auto& e : ckt.elements) {
        e->stamp(G, I, ckt, x, ctx);
    }

    // 使用自写 Gauss-Seidel 迭代
    x = Solver::solveLinearSystemGaussSeidel(G, I, 2000, 1e-10);
    return x;
}


// 非线性 DC：外层 Newton，内层用 LU 解线性化后的方程
VectorXd DcAnalysis::solveNewtonLU() const {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(std::max(N, 1));

    if (N == 0) {
        std::cerr << "DC solve (Newton + LU): no unknowns.\n";
        return x;
    }

    const int    rampSteps      = 10;    // 电源从 0 ~ 1 分 10 步 ramp
    const int    maxNewtonIters = 50;    // 每个 ramp 步最长 Newton 迭代次数
    const double tol            = 1e-9;  // Newton 收敛阈值（欧氏范数）

    ConvController ctrl;
    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale   = static_cast<double>(step) / rampSteps;
        double alpha   = ctrl.initialAlphaLU();
        double gmin    = ctrl.baseGmin(scale);
        double prevErr = std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            AnalysisContext ctx = makeDcCtx(scale);

            // 对当前迭代的 x 线性化并 stamp
            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, ctx);
            }

            stampGlobalGmin(ckt, G, gmin);

            // 解 G x_new = I
            VectorXd xRaw = Solver::solveLinearSystemLU(G, I);
            if (!xRaw.allFinite()) {
                gmin = std::min(gmin * 10.0, 1e-4);
                continue;
            }

            auto st = ctrl.update(
                x, xRaw, prevErr, iter,
                alpha, gmin, scale, tol
            );

            x       = st.xNext;
            alpha   = st.alphaNext;
            gmin    = st.gminNext;
            prevErr = st.error;

            if (st.converged) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: Newton (LU) did not converge at ramp step "
                          << step << " (err=" << st.error
                          << ", alpha=" << alpha
                          << ", gmin=" << gmin << ")\n";
            }
        }
    }

    return x;
}


VectorXd DcAnalysis::solveNewtonGS() const {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(std::max(N, 1));
    if (N == 0) return x;

    // --- Newton/continuation 参数（迭代法比 LU 慢，尺度要更现实）---
    const int    maxNewtonIters = 60;     // 每个 scale 的 Newton 迭代上限
    const int    maxLineSearch  = 12;     // 线搜索回退次数
    const double tolF_abs       = 1e-9;   // 残差绝对阈值
    const double tolF_rel       = 1e-7;   // 残差相对阈值
    const double tolX_rel       = 1e-6;   // 步长相对阈值（配合残差一起用）
    const double alphaMin       = 0.02;   // 允许很小，但不是默认一直用它
    const double alphaInit      = 1.0;    // 默认先尝试 full Newton
    const double c1             = 1e-4;   // Armijo 常数

    ConvController ctrl;                 // 只用它的 baseGmin 思路
    x.setZero(N);

    // continuation：从 0 -> 1 的 scale，用自适应步长推进
    double goodScale = 0.0;
    VectorXd xGood = x;

    double dScale = 0.1;                 // 初始步长
    const double minDScale = 1e-3;

    auto buildSystem = [&](const VectorXd& xAt,
                           double scale,
                           double gmin,
                           MatrixXd& G,
                           VectorXd& I)
    {
        G.setZero(N, N);
        I.setZero(N);
        AnalysisContext ctx = makeDcCtx(scale);
        for (const auto& e : ckt.elements) e->stamp(G, I, ckt, xAt, ctx);
        stampGlobalGmin(ckt, G, gmin);
    };

    auto residualNorm = [&](const MatrixXd& G,
                            const VectorXd& I,
                            const VectorXd& xAt) -> double
    {
        VectorXd F = G * xAt - I;
        double fn = F.norm();
        return std::isfinite(fn) ? fn : std::numeric_limits<double>::infinity();
    };

    while (goodScale < 1.0 - 1e-15) {
        double targetScale = std::min(1.0, goodScale + dScale);

        // 这一段是“尝试推进到 targetScale”
        VectorXd xTry = xGood;
        double gmin = ctrl.baseGmin(targetScale);

        bool stepOK = false;

        for (int it = 0; it < maxNewtonIters; ++it) {
            MatrixXd G(N, N);
            VectorXd I(N);

            buildSystem(xTry, targetScale, gmin, G, I);

            double Fnorm = residualNorm(G, I, xTry);
            double Inorm = I.norm();
            double tolF  = tolF_abs + tolF_rel * std::max(1.0, Inorm);

            // 先用残差判断：已经够好就算这个 scale 收敛
            if (Fnorm <= tolF) {
                stepOK = true;
                break;
            }

            // 解线性化方程：G * xRaw = I
            auto ls = Solver::solveLinearSystemGaussSeidelInfo(G, I, xTry,
                                                   /*maxIters=*/8000,
                                                   /*tol=*/1e-10);
            VectorXd xRaw = ls.x;
            if (!xRaw.allFinite()) {
                gmin = std::min(gmin * 10.0, 1e-2);
                continue;
            }

            // 线性残差验收（用原方程 Gx=I）
            VectorXd rlin = G * xRaw - I;
            double linRel = ls.relResidual;

            // 线性解太差才拒绝；阈值别太苛刻（让 line-search 发挥作用）
            if (!std::isfinite(linRel) || linRel > 5e-3) {
                gmin = std::min(gmin * 10.0, 1e-2);
                continue;
            }

            VectorXd dx = xRaw - xTry;
            double dxInf = dx.lpNorm<Eigen::Infinity>();
            double xInf  = xTry.lpNorm<Eigen::Infinity>();

            // 如果步长已经很小，也可以认为收敛（配合残差会更稳）
            if (dxInf <= (1e-12 + tolX_rel * std::max(1.0, xInf)) && Fnorm <= 10 * tolF) {
                stepOK = true;
                break;
            }

            // Armijo 线搜索：优先 full step，失败就缩 alpha
            double alpha = alphaInit;
            VectorXd bestXls = xTry;
            double bestFn = Fnorm;
            bool accepted = false;

            for (int lsit = 0; lsit < maxLineSearch; ++lsit) {
                VectorXd cand = xTry + alpha * dx;

                MatrixXd Gc(N, N);
                VectorXd Ic(N);
                buildSystem(cand, targetScale, gmin, Gc, Ic);

                double Fcand = residualNorm(Gc, Ic, cand);

                // Armijo: ||F(x+αdx)|| <= (1 - c1*α) ||F(x)||
                if (std::isfinite(Fcand) && (Fcand <= (1.0 - c1 * alpha) * Fnorm)) {
                    xTry = cand;
                    accepted = true;
                    break;
                }

                if (Fcand < bestFn) { bestFn = Fcand; bestXls = cand; }

                alpha *= 0.5;
                if (alpha < alphaMin) break;
            }

            if (!accepted) {
                // 线搜索都不行：增大 gmin，并把 xTry 拉回“最不差的候选”
                xTry = bestXls;
                gmin = std::min(gmin * 10.0, 1e-2);
            } else {
                // accepted：gmin 缓慢往 base 拉回（别一直顶在很大）
                double gminBase = ctrl.baseGmin(targetScale);
                gmin = 0.7 * gmin + 0.3 * gminBase;
            }
        }

        if (!stepOK) {
            // 这个 scale 推不动：缩小 dScale 重试
            dScale *= 0.5;
            std::cerr << "WARNING: Newton (GS) failed at scale=" << targetScale
                      << ", reducing dScale to " << dScale << "\n";
            if (dScale < minDScale) {
                std::cerr << "WARNING: dScale too small, stop continuation at scale=" << goodScale << "\n";
                return xGood;
            }
            continue;
        }

        // 成功推进一个 step
        xGood = xTry;
        goodScale = targetScale;

        // 如果走得很顺，稍微放大步长
        dScale = std::min(0.2, dScale * 1.25);
    }

    return xGood;
}




// ====================== 对外接口 ======================

// 显式：DC + LU（内部统一用自写 LU）
// VectorXd dcSolveLU(const Circuit& ckt) {
//     if (hasNonlinearDevices(ckt)) {
//         return dcSolveNewtonLU(ckt);
//     } else {
//         return dcSolveDirectLU(ckt);
//     }
// }

// // 显式：DC + Gauss-Seidel
// VectorXd dcSolveGaussSeidel(const Circuit& ckt) {
//     if (hasNonlinearDevices(ckt)) {
//         return dcSolveNewtonGS(ckt);
//     } else {
//         return dcSolveDirectGS(ckt);
//     }
// }

// // 老接口：现在默认等价于 Gauss-Seidel 版本
// VectorXd dcSolve(const Circuit& ckt) {
//     return dcSolveLU(ckt);
// }

ConvController::ConvController() : alphaMin(0.02), alphaMax(0.5), gminHighBase(1e-6)
            , gminLowBase(3.35e-7), gminAbsMax(1e-4), fastConvRatio(0.7), slowConvRatio(1.05) {}


ConvStatus ConvController::update(
    const Eigen::VectorXd& x, const Eigen::VectorXd& xRaw,
    double prevErr, int iter, double alphaCurrent, 
    double gminCurrent, double rampScale, double tol
) const {
    ConvStatus st;
    double alpha = std::clamp(alphaCurrent, alphaMin, alphaMax);
    Eigen::VectorXd dx = xRaw - x;
    double err = dx.lpNorm<Eigen::Infinity>();
    Eigen::VectorXd xNew = x + alpha * dx;

    double gminBase = baseGmin(rampScale);
    double gminNext = gminBase;

    if (iter == 0 || !std::isfinite(prevErr)) {
        // 第一轮：直接靠近 base
        gminNext = gminBase;
    } else {
        // 从第二轮开始，根据误差变化做一点调节
        if (err > prevErr * slowConvRatio) {
            // 收敛明显变差：减小 alpha、略微增大 gmin
            alpha    = std::max(alpha * 0.7, alphaMin);
            gminNext = std::min(gminCurrent * 2.0, gminAbsMax);
        } else if (err < prevErr * fastConvRatio) {
            // 收敛不错：略微增大 alpha，让牛顿稍微大胆一些，
            // 并把 gmin 轻轻往 base 拉
            alpha    = std::min(alpha * 1.1, alphaMax);
            gminNext = 0.5 * gminCurrent + 0.5 * gminBase;
        } else {
            // 正常收敛：gmin 缓慢往 base 靠近
            gminNext = 0.7 * gminCurrent + 0.3 * gminBase;
        }
        xNew = x + alpha * (xRaw - x);
        err = (xNew - x).norm();
    }

    st.xNext     = std::move(xNew);
    st.alphaNext = alpha;
    st.gminNext  = gminNext;
    st.error     = err;
    st.converged = (err < tol);

    return st;
}
