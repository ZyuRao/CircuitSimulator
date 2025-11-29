#include "dcanalysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "solver.hpp"
#include "sim.hpp"

#include <iostream>
#include <memory>

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

    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            AnalysisContext ctx = makeDcCtx(scale);

            // 对当前迭代的 x 线性化并 stamp
            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, ctx);
            }

            // 解 G x_new = I
            VectorXd xNew = Solver::solveLinearSystemLU(G, I);

            double err = (xNew - x).norm();
            x = xNew;

            if (err < tol) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: Newton (LU) did not converge at ramp step "
                          << step << "\n";
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
    const int    maxNewtonIters = 50;
    const double tol            = 1e-9;

    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            AnalysisContext ctx = makeDcCtx(scale);

            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, ctx);
            }

            // 这里用上一轮的 x 当作 Gauss-Seidel 的初值，利用 warm start
            VectorXd xNew =
                Solver::solveLinearSystemGaussSeidel(G, I, x, 2000, 1e-10);

            double err = (xNew - x).norm();
            x = xNew;

            if (err < tol) {
                break;
            }
            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: Newton (GS) did not converge at ramp step "
                          << step << "\n";
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
    return dcSolveGaussSeidel(ckt);
}
