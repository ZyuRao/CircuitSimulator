#include "dcanalysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "solver.hpp"
#include "sim.hpp"

#include <iostream>
#include <memory>
#include <limits>
#include <algorithm>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// 小工具：构造一个 DC 用的 AnalysisContext
static AnalysisContext makeDcCtx(double sourceScale) {
    AnalysisContext ctx;
    ctx.type        = AnalysisType::DC;
    ctx.sourceScale = sourceScale;
    ctx.time        = 0.0;
    ctx.omega       = 0.0;
    return ctx;
}

// 判定电路中是否存在非线性器件（目前只看 MOS）
static bool hasNonlinearDevices(const Circuit& ckt) {
    for (const auto& e : ckt.elements) {
        if (std::dynamic_pointer_cast<MosfetBase>(e)) {
            return true;
        }
    }
    return false;
}

//全局gmin-to-ground，后续可能有用
static void stampGlobalGmin(const Circuit& ckt, MatrixXd& G, double gmin) {
    for(const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if(eq >= 0 && eq < G.rows()) {
            G(eq, eq) += gmin;
        }
    }
}

// 只处理线性电路，线性方程 G x = I，用 LU 直接解
static VectorXd dcSolveDirectLU(const Circuit& ckt) {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(N);

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
static VectorXd dcSolveDirectGS(const Circuit& ckt) {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(N);

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
static VectorXd dcSolveNewtonLU(const Circuit& ckt) {
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
        double scale = static_cast<double>(step) / rampSteps;

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
                gmin = std::min(gmin * 10.0, 1e-2);
                continue;
            }

            auto st = ctrl.update(
                x, xRaw, prevErr, iter, alpha,
                gmin, scale, tol
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

// 非线性 DC：外层 Newton，内层用 Gauss-Seidel 解线性化后的方程
static VectorXd dcSolveNewtonGS(const Circuit& ckt) {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(std::max(N, 1));

    if (N == 0) {
        std::cerr << "DC solve (Newton + GS): no unknowns.\n";
        return x;
    }

    const int    rampSteps      = 10;
    const int    maxNewtonIters = 60;
    const double tol            = 1e-9;

    ConvController ctrl;

    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;
        double alpha   = ctrl.initialAlphaGS();
        double gmin    = ctrl.baseGmin(scale);
        double prevErr = std::numeric_limits<double>::infinity();
        int maxIterThisStep = maxNewtonIters;
        if (step == rampSteps) {
            maxIterThisStep = maxNewtonIters * 2;  // 最后一步给多一点机会
        }
        for (int iter = 0; iter < maxIterThisStep; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            AnalysisContext ctx = makeDcCtx(scale);

            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, ctx);
            }

            stampGlobalGmin(ckt, G, gmin);

            // 这里用上一轮的 x 当作 Gauss-Seidel 的初值，利用 warm start
            VectorXd xRaw =
                Solver::solveLinearSystemGaussSeidel(G, I, x, 2000, 1e-10);

            if (!xRaw.allFinite()) {
                gmin = std::min(gmin * 10.0, 1e-2);
                std::cerr << "WARNING: GS produced non-finite x, increasing gmin to "
                          << gmin << " at ramp step " << step
                          << ", iter " << iter << "\n";
                continue;
            }
            auto st = ctrl.update(
                        x, xRaw, prevErr, iter,alpha,
                        gmin, scale, tol );

            x       = st.xNext;
            alpha   = st.alphaNext;
            gmin    = st.gminNext;
            prevErr = st.error;

            if (st.converged) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: Newton (GS) did not converge at ramp step "
                          << step << " (err=" << st.error
                          << ", alpha=" << alpha
                          << ", gmin=" << gmin << ")\n";
            }
        }
    }

    return x;
}

// ====================== 对外接口 ======================

// 显式：DC + LU（内部统一用自写 LU）
VectorXd dcSolveLU(const Circuit& ckt) {
    if (hasNonlinearDevices(ckt)) {
        return dcSolveNewtonLU(ckt);
    } else {
        return dcSolveDirectLU(ckt);
    }
}

// 显式：DC + Gauss-Seidel
VectorXd dcSolveGaussSeidel(const Circuit& ckt) {
    if (hasNonlinearDevices(ckt)) {
        return dcSolveNewtonGS(ckt);
    } else {
        return dcSolveDirectGS(ckt);
    }
}

// 老接口：现在默认等价于 Gauss-Seidel 版本
VectorXd dcSolve(const Circuit& ckt) {
    return dcSolveLU(ckt);
}

ConvController::ConvController() : alphaMin(0.1), alphaMax(0.5), gminHighBase(1e-6)
            , gminLowBase(3.35e-7), gminAbsMax(1e-4), fastConvRatio(0.7), slowConvRatio(1.05) {}


ConvStatus ConvController::update(
    const Eigen::VectorXd& x, const Eigen::VectorXd& xRaw,
    double prevErr, int iter, double alphaCurrent, 
    double gminCurrent, double rampScale, double tol
) const {
    ConvStatus st;
    double alpha = std::clamp(0.35, alphaMin, alphaMax);
    Eigen::VectorXd xNew = x + alpha * (xRaw - x);
    double err = (xNew - x).norm();
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
    }

    st.xNext     = std::move(xNew);
    st.alphaNext = alpha;
    st.gminNext  = gminNext;
    st.error     = err;
    st.converged = (err < tol);

    return st;
}
