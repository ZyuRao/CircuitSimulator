#include "dcanalysis.hpp"
#include "circuit.hpp"
#include "element.hpp"
#include "solver.hpp"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

    for (const auto& e : ckt.elements) {
        e->stamp(G, I, ckt, x, 1.0);
    }

    x = Solver::solveLinearSystemLU(G, I);  // 手写 LU :contentReference[oaicite:1]{index=1}
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

    for (const auto& e : ckt.elements) {
        e->stamp(G, I, ckt, x, 1.0);
    }

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

    const int    rampSteps      = 10;
    const int    maxNewtonIters = 50;
    const double tol            = 1e-9;

    x.setZero(N);

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, scale);
            }

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

            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, scale);
            }

            // 这里用上一轮的 x 当作 Gauss-Seidel 的初值，利用 warm start
            VectorXd xNew = Solver::solveLinearSystemGaussSeidel(G, I, x, 2000, 1e-10);

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

// 默认 DC：统一用 LU（方便在项目里“随手就用”）
VectorXd dcSolveLU(const Circuit& ckt) {
    if (hasNonlinearDevices(ckt)) {
        return dcSolveNewtonLU(ckt);
    } else {
        return dcSolveDirectLU(ckt);
    }
}

// 明确指定：DC + Gauss-Seidel
VectorXd dcSolveGaussSeidel(const Circuit& ckt) {
    if (hasNonlinearDevices(ckt)) {
        return dcSolveNewtonGS(ckt);
    } else {
        return dcSolveDirectGS(ckt);
    }
}

// 老接口：默认等价于 dcSolveLU
VectorXd dcSolve(const Circuit& ckt) {
    return dcSolveLU(ckt);
}
