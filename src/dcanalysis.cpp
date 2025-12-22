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
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct DcAnalysis::DcWorkspace {
    MatrixXd linearG;
    VectorXd linearI;
    MatrixXd G;
    VectorXd I;
    VectorXd F;
    VectorXd xBuf;
};

DcAnalysis::~DcAnalysis() = default;

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

void DcAnalysis::initWorkspace(int N) const {
    workspace = std::make_unique<DcWorkspace>();
    workspace->linearG = MatrixXd::Zero(N, N);
    workspace->linearI = VectorXd::Zero(N);
    workspace->G = MatrixXd::Zero(N, N);
    workspace->I = VectorXd::Zero(N);
    workspace->F = VectorXd::Zero(N);
    workspace->xBuf = VectorXd::Zero(N);
    linearStampScale = std::numeric_limits<double>::quiet_NaN();
}

DcAnalysis::DcWorkspace& DcAnalysis::ws(int N) const {
    if (!workspace || workspace->G.rows() != N || workspace->G.cols() != N) {
        initWorkspace(N);
    }
    return *workspace;
}

void DcAnalysis::buildLinearStamp(int N, double sourceScale) const {
    DcWorkspace& w = ws(N);
    if (!std::isnan(linearStampScale) && std::abs(linearStampScale - sourceScale) < 1e-15) {
        return;
    }
    w.linearG.setZero();
    w.linearI.setZero();
    w.xBuf.setZero();
    AnalysisContext ctx = makeDcCtx(sourceScale);
    const auto& cache = ckt.elementCache();
    for (const auto* r : cache.resistors) r->stamp(w.linearG, w.linearI, ckt, w.xBuf, ctx);
    for (const auto* c : cache.capacitors) c->stamp(w.linearG, w.linearI, ckt, w.xBuf, ctx);
    for (const auto* L : cache.inductors) L->stamp(w.linearG, w.linearI, ckt, w.xBuf, ctx);
    for (const auto* vs : cache.voltageSources) vs->stamp(w.linearG, w.linearI, ckt, w.xBuf, ctx);
    for (const auto* is : cache.currentSources) is->stamp(w.linearG, w.linearI, ckt, w.xBuf, ctx);
    linearStampScale = sourceScale;
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

    DcWorkspace& w = ws(N);
    const auto& cache = ckt.elementCache();

    const int    rampSteps      = 10;    // 电源从 0 ~ 1 分 10 步 ramp
    const int    maxNewtonIters = 50;    // 每个 ramp 步最长 Newton 迭代次数

    auto nonlinearResidual = [&](const VectorXd& xAt, double scale, double gmin) -> double {
        w.G = w.linearG;
        w.I = w.linearI;
        AnalysisContext ctx = makeDcCtx(scale);
        for (const auto* m : cache.mos) m->stamp(w.G, w.I, ckt, xAt, ctx);
        stampGlobalGmin(ckt, w.G, gmin);

        w.F.noalias() = w.G * xAt - w.I;
        double fn = w.F.norm();
        return std::isfinite(fn) ? fn : std::numeric_limits<double>::infinity();
    };

    ConvController ctrl = ConvController::forDc(DcSolverKind::LU);
    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;
        buildLinearStamp(N, scale);

        double alpha   = ctrl.params().alphaInit;
        double gmin    = ctrl.baseGmin(scale);
        double prevStep = std::numeric_limits<double>::infinity();

        bool stepOK = false;
        double lastRes = std::numeric_limits<double>::infinity();
        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            w.G = w.linearG;
            w.I = w.linearI;

            AnalysisContext ctx = makeDcCtx(scale);

            // 对当前迭代的 x 线性化并 stamp
            for (const auto* m : cache.mos) {
                m->stamp(w.G, w.I, ckt, x, ctx);
            }
            stampGlobalGmin(ckt, w.G, gmin);

            // 解 G x_new = I
            w.xBuf = Solver::solveLinearSystemLU(w.G, w.I);
            if (!w.xBuf.allFinite()) {
                gmin = std::min(gmin * 10.0, ctrl.params().gminAbsMax);
                continue;
            }

            auto st = ctrl.updateStep(x, w.xBuf, prevStep, iter, alpha, gmin, scale);
            x        = st.xNext;
            alpha    = st.alphaNext;
            gmin     = st.gminNext;
            prevStep = st.error;

            lastRes = nonlinearResidual(x, scale, gmin);
            double tolF = ctrl.residualTol(/*rhsNorm=*/w.I.norm());
            double tolStep = ctrl.stepTol(scale);

            if (lastRes <= tolF || prevStep <= tolStep) {
                stepOK = true;
                break;
            }
        }
        if (!stepOK && step == rampSteps) {
            std::cerr << "WARNING: Newton (LU) not fully converged at final scale=1"
                      << " (res=" << lastRes
                      << ", alpha=" << alpha
                      << ", gmin=" << gmin << ")\n";
        }
    }
    return x;
}


VectorXd DcAnalysis::solveNewtonGS() const {
    int N = ckt.numUnknowns();
    VectorXd x0 = VectorXd::Zero(std::max(N, 1));
    if (N == 0) return x0;

    DcWorkspace& w = ws(N);
    const auto& cache = ckt.elementCache();

    ConvController ctrl = ConvController::forDc(DcSolverKind::GaussSeidel);
    const auto& p = ctrl.params();

    VectorXd xGood = VectorXd::Zero(N);
    double goodScale = 0.0;

    double dScale = p.dScaleInit;

    auto buildSystem = [&](const VectorXd& xAt, double scale, double gmin,
                           MatrixXd& G, VectorXd& I) {
        G = w.linearG;
        I = w.linearI;
        AnalysisContext ctx = makeDcCtx(scale);
        for (const auto* m : cache.mos) m->stamp(G, I, ckt, xAt, ctx);
        stampGlobalGmin(ckt, G, gmin);
    };

    auto residualNorm = [&](const VectorXd& xAt, double scale, double gmin) -> double {
        buildSystem(xAt, scale, gmin, w.G, w.I);
        w.F.noalias() = w.G * xAt - w.I;
        double fn = w.F.norm();
        return std::isfinite(fn) ? fn : std::numeric_limits<double>::infinity();
    };

    while (goodScale < 1.0 - 1e-15) {
        double targetScale = std::min(1.0, goodScale + dScale);
        VectorXd xTry = xGood;

        // 每个 scale 开始时：gmin 从 base 起步
        double gmin = ctrl.baseGmin(targetScale);
        buildLinearStamp(N, targetScale);

        bool stepOK = false;

        // 外层“Newton”迭代次数（你也可以把它放进 params）
        const int maxNewtonIters = 60;

        for (int it = 0; it < maxNewtonIters; ++it) {
            w.G = w.linearG;
            w.I = w.linearI;
            AnalysisContext ctx = makeDcCtx(targetScale);
            for (const auto* m : cache.mos) m->stamp(w.G, w.I, ckt, xTry, ctx);
            stampGlobalGmin(ckt, w.G, gmin);

            double Fnorm = (w.G * xTry - w.I).norm();
            double tolF  = ctrl.residualTol(w.I.norm());
            double tolStep = ctrl.stepTol(targetScale);

            // 1) 残差已经够小：认为该 scale 收敛
            if (std::isfinite(Fnorm) && Fnorm <= tolF) {
                stepOK = true;
                break;
            }

            // 2) 解线性化方程：G*xRaw = I （内层 GS）
            //    这里用带缩放 + 正规方程兜底的版本，避免零对角的约束行把 GS 卡死
            auto ls = Solver::solveLinearSystemGaussSeidelInfo(
                w.G, w.I, xTry,
                /*maxIters=*/8000,
                /*tol=*/1e-10
            );

            w.xBuf = ls.x;
            double linRel = ls.relResidual;

            // 线性内解太差：统一用 controller 的拒绝策略
            if (!w.xBuf.allFinite() || !std::isfinite(linRel) || linRel > ctrl.linRelTol(targetScale)) {
                double alphaDummy = p.alphaInit; // onLinearReject 需要一个 alpha 变量
                ctrl.onLinearReject(alphaDummy, gmin);
                continue;
            }

            w.F = w.xBuf - xTry;
            double stepInf = w.F.lpNorm<Eigen::Infinity>();
            double xInf = xTry.lpNorm<Eigen::Infinity>();
            double linTol = ctrl.linRelTol(targetScale);

            // 3) 如果步长已经非常小，也可认为收敛（配合残差的宽松倍数）
            double tinyStep = 1e-12 + tolStep * std::max(1.0, xInf);
            if (stepInf <= tinyStep) {
                if (Fnorm <= 15.0 * tolF || linRel <= 5.0 * linTol) {
                    stepOK = true;
                    break;
                }
            }

            // 4) Armijo line-search（全部用 params）
            double alpha = p.alphaInit;
            VectorXd bestCand = xTry;
            double bestFn = Fnorm;
            bool accepted = false;

            for (int lsit = 0; lsit < p.maxLineSearch; ++lsit) {
                VectorXd cand = xTry + alpha * w.F;

                // 注意：这里用 cand 重新 stamp 后的 residual 做判据（非线性意义上）
                double Fcand = residualNorm(cand, targetScale, gmin);

                if (std::isfinite(Fcand) && (Fcand <= (1.0 - p.armijoC1 * alpha) * Fnorm)) {
                    xTry = cand;
                    accepted = true;
                    break;
                }

                if (Fcand < bestFn) { bestFn = Fcand; bestCand = cand; }

                alpha *= 0.5;
                if (alpha < p.lsAlphaMin) break;
            }

            if (!accepted) {
                // 线搜索失败：回退到最不差的点，并加大 gmin
                xTry = bestCand;
                if (bestFn < Fnorm) {
                    double gbase = ctrl.baseGmin(targetScale);
                    gmin = 0.7 * gmin + 0.3 * gbase;
                } else {
                    double alphaDummy = p.alphaInit;
                    ctrl.onLinearReject(alphaDummy, gmin);
                }
            } else {
                // accepted：gmin 往 base 拉回，别一直顶着
                double gbase = ctrl.baseGmin(targetScale);
                gmin = 0.7 * gmin + 0.3 * gbase;
            }
        }

        if (!stepOK) {
            dScale *= 0.5;
            std::cerr << "WARNING: Newton (GS) failed at scale=" << targetScale
                      << ", reducing dScale to " << dScale << "\n";
            if (dScale < p.dScaleMin) {
                std::cerr << "WARNING: dScale too small, stop continuation at scale=" << goodScale << "\n";
                return xGood;
            }
            continue;
        }

        // 成功推进
        xGood = xTry;
        goodScale = targetScale;

        // 顺利就加大步长
        dScale = std::min(p.dScaleMax, dScale * 1.25);
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

ConvController ConvController::forDc(DcSolverKind kind) {
    ConvParams p;

    if (kind == DcSolverKind::LU) {
        // LU：步子要大一点，别跟 GS 共用 alphaMin=0.02
        p.alphaInit = 0.825;
        p.alphaMin  = 0.1;
        p.alphaMax  = 0.825;

        // gmin：允许 early 更大，失败时也允许顶到更高（但最终回到 lowBase）
        p.gminHighBase = 1e-5;
        p.gminLowBase  = 3.325e-7;
        p.gminAbsMax   = 1e-2;

        // DC step/residual 判据
        p.stepTolStart = 1e-6;
        p.stepTolFinal = 1e-9;
        p.fAbsTol      = 1e-12;
        p.fRelTol      = 1e-8;
        return ConvController(NonlinearSolveRole::DC_Newton_LU, p);
    } else {
        // GS：允许更小 alpha（line-search 也会用到）
        p.alphaInit = 0.95;
        p.alphaMin  = 0.02;
        p.alphaMax  = 0.95;

        p.gminHighBase = 1e-5;
        p.gminLowBase  = 3.325e-7;
        p.gminAbsMax   = 1e-2;

        // GS 内线性解验收阈值（你原本写死在 solveNewtonGS）
        p.linRelTolStart = 5e-3;
        p.linRelTolFinal = 1e-6;

        p.maxLineSearch = 12;
        p.armijoC1      = 1e-4;
        p.lsAlphaMin    = 0.02;

        p.dScaleInit = 0.1;
        p.dScaleMin  = 1e-3;
        p.dScaleMax  = 0.2;

        return ConvController(NonlinearSolveRole::DC_Newton_GS, p);
    }
}

ConvController ConvController::forHb() {
    ConvParams p;
    // HB 通常比 DC 更难：alpha 不要太激进
    p.alphaInit = 0.65;
    p.alphaMin  = 0.5;
    p.alphaMax  = 0.9;

    // 减小 gmin 以降低偏置误差，允许适度抬升但最终回到极低值
    p.gminHighBase = 1e-3;
    p.gminLowBase  = 3.325e-7;
    p.gminAbsMax   = 1e-2;

    p.stepTolStart = 1e-6;
    p.stepTolFinal = 1e-9;
    // HB 收敛标准：保持严格，但略高于 DC 以兼顾稳定性
    p.fAbsTol      = 1e-11;
    p.fRelTol      = 5e-7;

    return ConvController(NonlinearSolveRole::HB_Newton_LU, p);
}

// 关键：update 里不要 inf-norm / 2-norm 混着用
ConvStatus ConvController::updateStep(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& xRaw,
    double prevStep,
    int iter,
    double alphaCurrent,
    double gminCurrent,
    double rampScale
) const {
    ConvStatus st;

    double alpha = std::clamp(alphaCurrent, p_.alphaMin, p_.alphaMax);
    Eigen::VectorXd dx = xRaw - x;
    Eigen::VectorXd xNew = x + alpha * dx;

    // 用一致的 metric：建议 inf-norm
    double stepNow = (xNew - x).lpNorm<Eigen::Infinity>();

    double gminBase = baseGmin(rampScale);
    double gminNext = gminBase;

    if (iter == 0 || !std::isfinite(prevStep)) {
        gminNext = gminBase;
    } else {
        if (stepNow > prevStep * p_.slowConvRatio) {
            alpha    = std::max(alpha * 0.7, p_.alphaMin);
            gminNext = std::min(gminCurrent * 2.0, p_.gminAbsMax);
        } else if (stepNow < prevStep * p_.fastConvRatio) {
            alpha    = std::min(alpha * 1.1, p_.alphaMax);
            gminNext = 0.5 * gminCurrent + 0.5 * gminBase;
        } else {
            gminNext = 0.7 * gminCurrent + 0.3 * gminBase;
        }

        xNew = x + alpha * dx;
        stepNow = (xNew - x).lpNorm<Eigen::Infinity>();
    }

    st.xNext     = std::move(xNew);
    st.alphaNext = alpha;
    st.gminNext  = gminNext;
    st.error     = stepNow;
    // “step 收敛”只是一条腿，LU 警告修复要靠 residual（下面 solveNewtonLU 会加）
    st.converged = false;
    return st;
}
