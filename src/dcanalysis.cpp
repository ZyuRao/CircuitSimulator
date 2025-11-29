#include "circuit.hpp"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// =============== 牛顿迭代 DC（带斜坡源） ===============

VectorXd dcSolveNewton(const Circuit& ckt) {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(std::max(N, 1)); // N=0 时给个长度1的 dummy

    if (N == 0) {
        std::cerr << "DC solve: no unknowns.\n";
        return x;
    }

    const int rampSteps = 10;     // 斜坡分成 10 步
    const int maxNewtonIters = 50;
    const double tol = 1e-9;

    x.setZero();

    for (int step = 1; step <= rampSteps; ++step) {
        double scale = static_cast<double>(step) / rampSteps;

        for (int iter = 0; iter < maxNewtonIters; ++iter) {
            MatrixXd G = MatrixXd::Zero(N, N);
            VectorXd I = VectorXd::Zero(N);

            // 根据当前猜测解 x，构建线性化后的 G & I
            AnalysisContext ctx;
            ctx.type        = AnalysisType::DC;
            ctx.sourceScale = scale;

            for (const auto& e : ckt.elements) {
                e->stamp(G, I, ckt, x, ctx);
            }

            // 解 G * x_new = I
            Eigen::ColPivHouseholderQR<MatrixXd> solver(G);
            if (solver.info() != Eigen::Success) {
                std::cerr << "WARNING: matrix decomposition failed at ramp step "
                          << step << ", iter " << iter << "\n";
                break;
            }
            VectorXd xNew = solver.solve(I);
            if (solver.info() != Eigen::Success) {
                std::cerr << "WARNING: linear solve failed at ramp step "
                          << step << ", iter " << iter << "\n";
                break;
            }

            double err = (xNew - x).norm();
            x = xNew;
            if (err < tol) {
                // 收敛
                // std::cout << "Converged: step " << step << ", iter " << iter << "\n";
                break;
            }

            if (iter == maxNewtonIters - 1) {
                std::cerr << "WARNING: Newton did not converge at ramp step "
                          << step << "after" << maxNewtonIters << "iters\n";
            }
        }
    }

    return x;
}

VectorXd dcSolveDirect(const Circuit& ckt) {
    int N = ckt.numUnknowns();
    VectorXd x = VectorXd::Zero(N);

    if (N == 0) {
        std::cerr << "DC solve: no unknowns.\n";
        return x;
    }

    MatrixXd G = MatrixXd::Zero(N, N);
    VectorXd I = VectorXd::Zero(N);

    // 构建线性方程组 G * x = I
    AnalysisContext ctx;
    ctx.type        = AnalysisType::OP;
    ctx.sourceScale = 1.0;

    for (const auto& e : ckt.elements) {
        e->stamp(G, I, ckt, x, ctx);
    }

    // 直接解线性方程组
    Eigen::ColPivHouseholderQR<MatrixXd> solver(G);
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERROR: matrix decomposition failed in direct DC solve.\n";
        return x;
    }
    x = solver.solve(I);
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERROR: linear solve failed in direct DC solve.\n";
        return x;
    }

    return x;
}

VectorXd dcSolve(const Circuit& ckt){
    for (const auto& e : ckt.elements) {
        if (std::dynamic_pointer_cast<MosfetBase>(e)) {
            // 存在非线性器件，使用牛顿迭代
            return dcSolveNewton(ckt);
        }
    }
    return dcSolveDirect(ckt);
}