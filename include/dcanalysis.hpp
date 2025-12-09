#pragma once

#include <Eigen/Dense>

class Circuit;

// DC 工作点（默认使用 Gauss-Seidel 作为线性求解器）
Eigen::VectorXd dcSolve(const Circuit& ckt);

// DC 工作点：内部统一用手写 LU
Eigen::VectorXd dcSolveLU(const Circuit& ckt);

// DC 工作点：内部统一用 Gauss-Seidel
Eigen::VectorXd dcSolveGaussSeidel(const Circuit& ckt);

struct ConvStatus {
    Eigen::VectorXd xNext;
    double alphaNext;
    double gminNext;
    double error;
    bool converged;
};

class ConvController {
private:
    double alphaMin, alphaMax;
    double gminHighBase, gminLowBase, gminAbsMax;

    double fastConvRatio;
    double slowConvRatio;
public:
    ConvController();

    ConvStatus update(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& xRaw,
        double prevErr,
        int    iter,
        double alphaCurrent,
        double gminCurrent,
        double rampScale,
        double tol
    ) const;
    // 给定 ramp 进度，计算当前步的基础 gmin
    double baseGmin(double rampScale) const {
        rampScale = std::clamp(rampScale, 0.0, 1.0);
        return gminHighBase * (1.0 - rampScale) + gminLowBase * rampScale;
    }

    // LU 版的初始 alpha
    double initialAlphaLU() const { return 0.5; }

    // GS 版可以稍微激进一点
    double initialAlphaGS() const { return 0.7; }
};
