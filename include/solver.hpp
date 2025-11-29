#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>

// 专门放线性方程 G x = I 的求解算法
// - 直接法：自写 LU 分解（部分主元）
// - 迭代法：Gauss-Seidel

namespace Solver {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// 线性求解器类型（给 DC 分析选择用）
enum class LinearSolver {
    DirectLU,
    GaussSeidel
};

// ================ 自写 LU 分解（Doolittle + 部分主元） =================
//
// A 被分解为 P * A = L * U
// 这里：LU 矩阵同时存 L、U：
//   - 下三角（含对角线下方）：L 的非对角元素（对角线隐含为 1）
//   - 上三角（含对角线）：U
// perm 记录行置换：b_perm[i] = b[perm[i]]
inline bool luDecompose(const MatrixXd& A, MatrixXd& LU, std::vector<int>& perm) {
    const double eps = 1e-15;

    int n = static_cast<int>(A.rows());
    if (n == 0) return false;
    if (A.cols() != n) {
        std::cerr << "LU: matrix is not square.\n";
        return false;
    }

    LU = A;
    perm.resize(n);
    for (int i = 0; i < n; ++i) {
        perm[i] = i;
    }

    for (int k = 0; k < n; ++k) {
        // 选主元：找当前列 k 中绝对值最大的元素
        int pivot = k;
        double maxAbs = std::fabs(LU(k, k));
        for (int i = k + 1; i < n; ++i) {
            double val = std::fabs(LU(i, k));
            if (val > maxAbs) {
                maxAbs = val;
                pivot = i;
            }
        }

        if (maxAbs < eps) {
            std::cerr << "LU: zero (or tiny) pivot at column " << k << ".\n";
            return false;
        }

        // 行交换
        if (pivot != k) {
            LU.row(k).swap(LU.row(pivot));
            std::swap(perm[k], perm[pivot]);
        }

        // 消元，填 L 的本列（对角线下方）和更新 U 的剩余元素
        for (int i = k + 1; i < n; ++i) {
            double factor = LU(i, k) / LU(k, k);
            LU(i, k) = factor;  // 存在 L(i,k)
            for (int j = k + 1; j < n; ++j) {
                LU(i, j) -= factor * LU(k, j);
            }
        }
    }

    return true;
}

// 利用 LU 分解求解 A x = b
inline VectorXd solveLinearSystemLU(const MatrixXd& A, const VectorXd& b) {
    int n = static_cast<int>(A.rows());
    VectorXd x = VectorXd::Zero(n);
    if (n == 0) return x;
    if (A.cols() != n || b.size() != n) {
        std::cerr << "LU solve: dimension mismatch.\n";
        return x;
    }

    MatrixXd LU;
    std::vector<int> perm;
    if (!luDecompose(A, LU, perm)) {
        std::cerr << "LU solve: decomposition failed.\n";
        return x;
    }

    // b_perm = P * b
    VectorXd b_perm(n);
    for (int i = 0; i < n; ++i) {
        b_perm(i) = b(perm[i]);
    }

    // 前代：L y = b_perm （L 对角线为 1）
    VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        double sum = b_perm(i);
        for (int j = 0; j < i; ++j) {
            sum -= LU(i, j) * y(j);
        }
        y(i) = sum;  // 因为 L(i,i) = 1
    }

    // 回代：U x = y
    for (int i = n - 1; i >= 0; --i) {
        double sum = y(i);
        for (int j = i + 1; j < n; ++j) {
            sum -= LU(i, j) * x(j);
        }
        double diag = LU(i, i);
        if (std::fabs(diag) < 1e-15) {
            std::cerr << "LU solve: zero diagonal at row " << i << ".\n";
            x(i) = 0.0;
        } else {
            x(i) = sum / diag;
        }
    }

    return x;
}

// ================ Gauss-Seidel 迭代法 =================
//
// 解 A x = b
// x0 为初值（可以用上一轮牛顿迭代的 x，当作 warm start）
// 遇到对角元接近 0 时，增加一个极小量进行正则化，避免直接失败。
// 返回最后一次迭代结果（不保证一定收敛）
inline VectorXd solveLinearSystemGaussSeidel(const MatrixXd& A,
                                             const VectorXd& b,
                                             const VectorXd& x0,
                                             int maxIters = 1000,
                                             double tol = 1e-10) {
    int n = static_cast<int>(A.rows());
    VectorXd x = x0;
    if (n == 0) return x;

    if (A.cols() != n || b.size() != n) {
        std::cerr << "Gauss-Seidel: dimension mismatch.\n";
        return VectorXd::Zero(std::max(n, 0));
    }

    if (x.size() != n) {
        // 如果给的初值维度不对，就从 0 向量开始
        x = VectorXd::Zero(n);
    }

    VectorXd xOld = x;
    const double diagEps = 1e-12; // 对角正则化的极小量

    for (int iter = 0; iter < maxIters; ++iter) {
        xOld = x;

        for (int i = 0; i < n; ++i) {
            double diag = A(i, i);

            // 对角元过小，用一个极小量替代，避免除以 0
            if (std::fabs(diag) < diagEps) {
                // 保持原来的符号（如果有），否则默认为正
                double sign = (diag >= 0.0 ? 1.0 : -1.0);
                diag = sign * diagEps;
            }

            double sum = b(i);

            // j < i 用最新的 x(j)，j > i 用上一轮的 xOld(j)
            for (int j = 0; j < i; ++j) {
                sum -= A(i, j) * x(j);
            }
            for (int j = i + 1; j < n; ++j) {
                sum -= A(i, j) * xOld(j);
            }

            x(i) = sum / diag;
        }

        double err = (x - xOld).norm();
        if (err < tol) {
            // 收敛
            break;
        }
    }

    return x;
}

// 不给初值时，从 0 开始
inline VectorXd solveLinearSystemGaussSeidel(const MatrixXd& A,
                                             const VectorXd& b,
                                             int maxIters = 1000,
                                             double tol = 1e-10) {
    VectorXd x0 = VectorXd::Zero(b.size());
    return solveLinearSystemGaussSeidel(A, b, x0, maxIters, tol);
}

} // namespace Solver
