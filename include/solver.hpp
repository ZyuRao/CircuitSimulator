#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>

// 专门放线性方程 G x = I 的求解算法
// - 直接法：自写 LU 分解（部分主元）
// - 迭代法：Gauss-Seidel

namespace Solver {

    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using MatrixXcd = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXcd = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;


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

    //通用LU（部分主元）
    template <typename MatrixType>
    inline bool luDecompose(const MatrixType& A,
                    MatrixType& LU,
                    std::vector<int>& perm) {
        using std::abs;
        const int n = static_cast<int>(A.rows());
        if (A.cols() != n) {
            std::cerr << "LU: matrix is not square.\n";
            return false;
        }

        LU = A;
        perm.resize(n);
        for (int i = 0; i < n; ++i) perm[i] = i;

        for (int k = 0; k < n; ++k) {
            // 选主元
            int pivot = k;
            auto maxVal = abs(LU(k, k));
            for (int i = k + 1; i < n; ++i) {
                auto v = abs(LU(i, k));
                if (v > maxVal) {
                    maxVal = v;
                    pivot = i;
                }
            }

            if (maxVal == typename MatrixType::Scalar(0)) {
                std::cerr << "LU: singular matrix at column " << k << "\n";
                return false;
            }

            if (pivot != k) {
                LU.row(k).swap(LU.row(pivot));
                std::swap(perm[k], perm[pivot]);
            }

            // 消元
            for (int i = k + 1; i < n; ++i) {
                LU(i, k) /= LU(k, k);
                for (int j = k + 1; j < n; ++j) {
                    LU(i, j) -= LU(i, k) * LU(k, j);
                }
            }
        }
        return true;
    }

    // 通用LU求解
    template<typename MatrixType, typename VectorType>
    VectorType luSolve(
        const MatrixType& LU,
        const std::vector<int>& perm,
        const VectorType& b
    ) {
        const int n = static_cast<int>(LU.rows());
        VectorType x = VectorType::Zero(n);
        VectorType y = VectorType::Zero(n);

        // 应用置换：b_p[i] = b[perm[i]]
        for (int i = 0; i < n; ++i) {
            y(i) = b(perm[i]);
        }

        // 前代：Ly = b_p （L 对角为 1）
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                y(i) -= LU(i, j) * y(j);
            }
        }

        // 回代：Ux = y
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) {
                y(i) -= LU(i, j) * x(j);
            }
            auto diag = LU(i, i);
            if (diag == typename MatrixType::Scalar(0)) {
                std::cerr << "LU solve: zero diagonal at row " << i << "\n";
                x(i) = typename MatrixType::Scalar(0);
            } else {
                x(i) = y(i) / diag;
            }
        }
        return x;
    }

    //实数版
    inline VectorXd solveLinearSystemLU(const MatrixXd& A, const VectorXd& b) {
        MatrixXd LU;
        std::vector<int> perm;
        if (!luDecompose(A, LU, perm)) {
            std::cerr << "solveLinearSystemLU (real): LU failed, "
                        "fallback to zero solution.\n";
            return VectorXd::Zero(b.size());
        }
        return luSolve<MatrixXd, VectorXd>(LU, perm, b);
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


    inline VectorXcd solveLinearSystemLU(const MatrixXcd& A, const VectorXcd& b) {
        MatrixXcd LU;
        std::vector<int> perm;
        if (!luDecompose(A, LU, perm)) {
            std::cerr << "solveLinearSystemLU (complex): LU failed, "
                        "fallback to zero solution.\n";
            return VectorXcd::Zero(b.size());
        }
        return luSolve<MatrixXcd, VectorXcd>(LU, perm, b);
    }

} // namespace Solver
