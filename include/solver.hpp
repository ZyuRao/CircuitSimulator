#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <random>

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
    struct IterSolveInfo {
        VectorXd x;
        bool     converged = false;
        double   relResidual = std::numeric_limits<double>::infinity(); // 用原方程 A*x-b 评估
        int      iters = 0;
    };

    inline void equilibrationScalesMaxAbs(const MatrixXd& A, VectorXd& rScale, VectorXd& cScale) {
        const int n = static_cast<int>(A.rows());
        rScale = VectorXd::Ones(n);
        cScale = VectorXd::Ones(n);

        for (int i = 0; i < n; ++i) {
            double mx = A.row(i).cwiseAbs().maxCoeff();
            if (!std::isfinite(mx) || mx < 1e-18) mx = 1.0;
            rScale(i) = 1.0 / mx;
        }
        for (int j = 0; j < n; ++j) {
            double mx = A.col(j).cwiseAbs().maxCoeff();
            if (!std::isfinite(mx) || mx < 1e-18) mx = 1.0;
            cScale(j) = 1.0 / mx;
        }
    }

    // SPD 上的 GS（可选 under-relaxation）
    inline IterSolveInfo gaussSeidelSPDInfo(const MatrixXd& M,
                                        const VectorXd& rhs,
                                        const VectorXd& x0,
                                        int maxIters,
                                        double tol,
                                        double omega = 1.0)
    {
        IterSolveInfo info;
        const int n = static_cast<int>(M.rows());
        info.x = (x0.size() == n) ? x0 : VectorXd::Zero(n);
        if (n == 0) { info.converged = true; info.relResidual = 0.0; return info; }

        VectorXd x = info.x, xOld = x;

        auto relResSPD = [&](const VectorXd& xx) -> double {
            VectorXd r = M * xx - rhs;
            double rn = r.norm();
            double bn = rhs.norm();
            if (!std::isfinite(rn)) return std::numeric_limits<double>::infinity();
            return (bn > 0.0) ? (rn / bn) : rn;
        };

        double best = std::numeric_limits<double>::infinity();
        VectorXd bestX = x;

        const double diagFloor = 1e-18;

        for (int it = 0; it < maxIters; ++it) {
            xOld = x;

            for (int i = 0; i < n; ++i) {
                double diag = M(i, i);
                if (!std::isfinite(diag) || std::fabs(diag) < diagFloor) {
                    diag = (diag >= 0.0 ? diagFloor : -diagFloor);
                }
                double sum = rhs(i);
                for (int j = 0; j < i; ++j)     sum -= M(i, j) * x(j);
                for (int j = i + 1; j < n; ++j) sum -= M(i, j) * xOld(j);

                double xGS = sum / diag;
                double xi  = (1.0 - omega) * xOld(i) + omega * xGS;

                if (!std::isfinite(xi)) {
                    // 爆了：回退到 best，并收紧 omega
                    x = bestX;
                    omega = std::max(0.05, omega * 0.5);
                    break;
                }
                x(i) = xi;
            }

            double rr = relResSPD(x);
            if (rr < best) { best = rr; bestX = x; }

            if (rr <= tol) {
                info.x = x;
                info.converged = true;
                info.relResidual = rr;
                info.iters = it + 1;
                return info;
            }

            // 残差明显变差：收紧 omega，避免发散
            double rrOld = relResSPD(xOld);
            if (rr > rrOld * 1.5 && it > 3) {
                omega = std::max(0.05, omega * 0.7);
                x = bestX;
            }
        }

        info.x = bestX;
        info.converged = false;
        info.relResidual = best;
        info.iters = maxIters;
        return info;
    }

    inline Solver::IterSolveInfo solveLinearSystemGaussSeidelInfo(const MatrixXd& A,
                                                              const VectorXd& b,
                                                              const VectorXd& x0,
                                                              int maxIters = 8000,
                                                              double tol = 1e-10)
    {
        IterSolveInfo info;
        const int n = static_cast<int>(A.rows());
        if (n == 0) { info.x = VectorXd(); info.converged = true; info.relResidual = 0.0; info.iters = 0; return info; }
        if (A.cols() != n || b.size() != n) {
            std::cerr << "Iter solver (GS): dimension mismatch.\n";
            info.x = VectorXd::Zero(std::max(n, 1));
            info.converged = false;
            info.relResidual = std::numeric_limits<double>::infinity();
            info.iters = 0;
            return info;
        }

        auto relResidualOriginal = [&](const VectorXd& xCand) -> double {
            VectorXd r = A * xCand - b;
            double rn = r.norm();
            double bn = b.norm();
            if (!std::isfinite(rn)) return std::numeric_limits<double>::infinity();
            return (bn > 0.0) ? (rn / bn) : rn;
        };

        // ================== 0) 行列缩放（非常关键） ==================
        VectorXd rS, cS;
        equilibrationScalesMaxAbs(A, rS, cS);

        MatrixXd As = A;
        VectorXd bs = b;
        for (int i = 0; i < n; ++i) { As.row(i) *= rS(i); bs(i) *= rS(i); }
        for (int j = 0; j < n; ++j) { As.col(j) *= cS(j); }

        // x = Dc * y  (Dc = diag(cS))
        VectorXd y = VectorXd::Zero(n);
        if (x0.size() == n) {
            for (int j = 0; j < n; ++j) {
                double s = cS(j);
                if (!std::isfinite(s) || std::fabs(s) < 1e-18) s = 1.0;
                y(j) = x0(j) / s;
            }
        }

        auto recoverX = [&](const VectorXd& yy) -> VectorXd {
            VectorXd x = VectorXd::Zero(n);
            for (int j = 0; j < n; ++j) {
                double s = cS(j);
                if (!std::isfinite(s) || std::fabs(s) < 1e-18) s = 1.0;
                x(j) = yy(j) * s;
            }
            return x;
        };

        // ================== 1) 在缩放后的系统上跑 SOR/GS ==================
        const double diagFloor = 1e-18;
        double omega = 1.0;                 // 先从 1.0 开始，稳定优先
        VectorXd yOld = y;

        VectorXd xBest = recoverX(y);
        double bestRR  = relResidualOriginal(xBest);

        // 早停/判定“迭代没进展”的阈值
        const int   stallWindow = 50;
        const double stallRatio = 0.9995;   // 50 轮内残差下降不到 0.05% 就算停滞

        int stallCnt = 0;
        double rrLast = bestRR;

        for (int it = 0; it < maxIters; ++it) {
            yOld = y;

            for (int i = 0; i < n; ++i) {
                double diag = As(i, i);
                if (!std::isfinite(diag) || std::fabs(diag) < diagFloor) {
                    diag = (diag >= 0.0 ? diagFloor : -diagFloor);
                }

                double sum = bs(i);
                for (int j = 0; j < i; ++j)     sum -= As(i, j) * y(j);
                for (int j = i + 1; j < n; ++j) sum -= As(i, j) * yOld(j);

                double yGS = sum / diag;
                double yi  = (1.0 - omega) * yOld(i) + omega * yGS;

                if (!std::isfinite(yi)) {
                    y = yOld;
                    omega = std::max(0.05, omega * 0.5);
                    break;
                }
                y(i) = yi;
            }

            VectorXd xCand = recoverX(y);
            double rr = relResidualOriginal(xCand);

            if (rr < bestRR) { bestRR = rr; xBest = xCand; }
            else omega = std::max(0.05, omega * 0.9);

            if (std::isfinite(bestRR) && bestRR <= tol) {
                info.x = xBest;
                info.converged = true;
                info.relResidual = bestRR;
                info.iters = it + 1;
                return info;
            }

            // 停滞检测：如果一直几乎不降，就别在原始系统上耗命了
            if (it > 5) {
                if (rr >= rrLast * stallRatio) stallCnt++;
                else stallCnt = 0;
                rrLast = std::min(rrLast, rr);

                if (stallCnt >= stallWindow) {
                    break;
                }
            }
        }

        // ================== 2) 兜底：正规方程 (A^T A) x = (A^T b) ==================
        // 在 MNA 含理想源/电感约束时，这一步经常是“救命稻草”
        MatrixXd AtA = A.transpose() * A;
        VectorXd Atb = A.transpose() * b;

        // 极小正则，防止 AtA 病态到对角接近 0
        double diagMax = AtA.diagonal().cwiseAbs().maxCoeff();
        if (!std::isfinite(diagMax) || diagMax < 1e-18) diagMax = 1.0;
        AtA.diagonal().array() += 1e-14 * diagMax;

        // 用你已有的 SPD-GS 解正规方程
        auto spd = gaussSeidelSPDInfo(AtA, Atb,
                                    /*x0=*/xBest,
                                    /*maxIters=*/6000,
                                    /*tol=*/std::max(tol, 1e-10),
                                    /*omega=*/1.0);

        VectorXd x2 = spd.x;
        double rr2 = relResidualOriginal(x2);

        // 选更好的那个
        if (rr2 < bestRR) { xBest = x2; bestRR = rr2; }

        info.x = xBest;
        info.converged = (std::isfinite(bestRR) && bestRR <= tol);
        info.relResidual = bestRR;
        info.iters = maxIters;
        return info;
    }



    // 兼容旧接口：只返回 x
    inline VectorXd solveLinearSystemGaussSeidel(const MatrixXd& A,
                                                const VectorXd& b,
                                                const VectorXd& x0,
                                                int maxIters = 5000,
                                                double tol = 1e-10)
    {
        return solveLinearSystemGaussSeidelInfo(A, b, x0, maxIters, tol).x;
    }

    inline VectorXd solveLinearSystemGaussSeidel(const MatrixXd& A,
                                                const VectorXd& b,
                                                int maxIters = 5000,
                                                double tol = 1e-10)
    {
        VectorXd x0 = VectorXd::Zero(A.cols());
        return solveLinearSystemGaussSeidel(A, b, x0, maxIters, tol);
    }

} // namespace Solver
