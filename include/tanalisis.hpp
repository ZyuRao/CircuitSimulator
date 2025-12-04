#pragma once

#include <string>
#include <Eigen/Dense>
#include "circuit.hpp"
#include "sim.hpp"

// 统一 DC 工作点接口（封装 dcanalysis 的 dcSolve）
Eigen::VectorXd computeDcOperatingPoint(const Circuit& ckt);

// 后向欧拉瞬态仿真：
// - 使用 DC 工作点作为 t=0 初值
// - 使用 .TRAN 里的 tstep / tstop / tstart
// - 输出 CSV 文件：time, 所有节点电压, 所有电压源/电感支路电流
void runTransientAnalysisBackwardEuler(const Circuit& ckt,
                                       const SimulationConfig& sim,
                                       const std::string& outFile);
