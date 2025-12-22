#pragma once

#include <string>
#include "circuit.hpp"
#include "sim.hpp"

// 高层调度：解析网表 -> 运行各个分析 -> 生成 CSV / 触发画图
class Runner {
public:
    explicit Runner(const std::string& outDir = "out", bool verbose = false);

    // 返回 0 表示成功，非 0 表示错误码
    int run(const std::string& netlistPath);

private:
    std::string outDir;
    bool verbose = false;
};
