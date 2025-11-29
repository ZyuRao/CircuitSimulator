#pragma once

#include <Eigen/Dense>

class Circuit;

// DC 工作点（默认使用 Gauss-Seidel 作为线性求解器）
Eigen::VectorXd dcSolve(const Circuit& ckt);

// DC 工作点：内部统一用手写 LU
Eigen::VectorXd dcSolveLU(const Circuit& ckt);

// DC 工作点：内部统一用 Gauss-Seidel
Eigen::VectorXd dcSolveGaussSeidel(const Circuit& ckt);
