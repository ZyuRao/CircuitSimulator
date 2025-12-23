#include "parser.hpp"
#include "utils.hpp"
#include "element.hpp"
#include "runner.hpp"
#include "analysis.hpp"
#include <future>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include <cstdlib>

namespace {

using ColumnMap = std::unordered_map<std::string, int>;

std::string analysisTag(AnalysisType a) {
    switch (a) {
        case AnalysisType::OP:   return "op";
        case AnalysisType::DC:   return "dc";
        case AnalysisType::AC:   return "ac";
        case AnalysisType::TRAN: return "tran";
        case AnalysisType::HB:   return "hb";
        default: return "unknown";
    }
}

std::string probeHeaderName(const ProbeSpec& p) {
    if (p.kind == ProbeKind::NodeVoltage) {
        if (p.node2.empty() || isGroundName(p.node2)) {
            return "V(" + p.node1 + ")";
        }
        return "V(" + p.node1 + "," + p.node2 + ")";
    }
    if (p.kind == ProbeKind::DiffVoltage) {
        return "V(" + p.node1 + "," + p.node2 + ")";
    }
    if (p.kind == ProbeKind::BranchCurrent) {
        return "I(" + p.eleName + ")";
    }
    return p.expr;
}

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : line) {
        if (ch == ',') {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

double parseNum(const std::string& s) {
    try {
        return std::stod(s);
    } catch (...) {
        return 0.0;
    }
}

std::vector<std::string> readHeader(const std::string& path) {
    std::ifstream ifs(path);
    std::string line;
    if (ifs && std::getline(ifs, line)) {
        return splitCsv(line);
    }
    return {};
}

struct ColumnRequest {
    enum class Kind { Direct, Diff } kind = Kind::Direct;
    std::string name;
    int idxA = -1;
    int idxB = -1; // diff: idxA - idxB，idxB=-1 表示地
};

std::vector<ProbeSpec> collectCsvProbes(const SimulationConfig& sim, AnalysisType a) {
    std::vector<ProbeSpec> out;
    std::unordered_set<std::string> seen;

    auto add = [&](const ProbeSpec& p) {
        std::string key = probeHeaderName(p);
        if (seen.insert(key).second) out.push_back(p);
    };

    for (const auto& pc : sim.printCommands) {
        if (!matchAnalysis(pc, a)) continue;
        for (const auto& p : pc.probes) add(p);
    }
    for (const auto& p : sim.plotProbes) add(p);
    return out;
}

std::vector<ProbeSpec> collectPlotProbes(const SimulationConfig& sim, AnalysisType a) {
    std::vector<ProbeSpec> out;
    for (const auto& pc : sim.probeCommands) {
        if (!(pc.analysis == AnalysisType::NONE || pc.analysis == a)) continue;
        out.insert(out.end(), pc.probes.begin(), pc.probes.end());
    }
    return out;
}

std::vector<ProbeSpec> mergeProbes(std::vector<ProbeSpec> lhs,
                                   const std::vector<ProbeSpec>& rhs) {
    std::unordered_set<std::string> seen;
    for (const auto& p : lhs) {
        seen.insert(probeHeaderName(p));
    }
    for (const auto& p : rhs) {
        std::string key = probeHeaderName(p);
        if (seen.insert(key).second) lhs.push_back(p);
    }
    return lhs;
}

ColumnMap buildColumnMap(const std::vector<std::string>& header) {
    ColumnMap m;
    for (size_t i = 0; i < header.size(); ++i) {
        m[header[i]] = static_cast<int>(i);
    }
    return m;
}

bool filterCsv(const std::string& rawPath,
               const std::string& dstPath,
               const std::vector<ProbeSpec>& requested) {
    std::ifstream ifs(rawPath);
    if (!ifs) {
        std::cerr << "Cannot open raw CSV: " << rawPath << "\n";
        return false;
    }
    std::string headerLine;
    if (!std::getline(ifs, headerLine)) {
        std::cerr << "Empty CSV: " << rawPath << "\n";
        return false;
    }
    auto header = splitCsv(headerLine);
    ColumnMap colMap = buildColumnMap(header);

    int timeIdx = -1;
    auto itTime = colMap.find("time");
    if (itTime != colMap.end()) timeIdx = itTime->second;
    auto itT = colMap.find("t");
    if (timeIdx < 0 && itT != colMap.end()) timeIdx = itT->second;

    // 构建列请求
    std::vector<ColumnRequest> cols;
    std::unordered_set<std::string> seenNames;
    bool includeAll = requested.empty();
    if (includeAll) {
        for (size_t i = 0; i < header.size(); ++i) {
            if ((int)i == timeIdx) continue;
            ColumnRequest req;
            req.name = header[i];
            req.kind = ColumnRequest::Kind::Direct;
            req.idxA = static_cast<int>(i);
            cols.push_back(req);
        }
    } else {
        for (const auto& p : requested) {
            ColumnRequest req;
            req.name = probeHeaderName(p);

            if (!seenNames.insert(req.name).second) continue;

            if (p.kind == ProbeKind::NodeVoltage || p.kind == ProbeKind::DiffVoltage) {
                std::string colA = "V(" + p.node1 + ")";
                auto itA = colMap.find(colA);
                if (itA == colMap.end()) {
                    std::cerr << "[CSV] Missing column " << colA << " for probe " << req.name << "\n";
                    continue;
                }
                req.idxA = itA->second;
                if (!p.node2.empty() && !isGroundName(p.node2)) {
                    std::string colB = "V(" + p.node2 + ")";
                    auto itB = colMap.find(colB);
                    if (itB == colMap.end()) {
                        std::cerr << "[CSV] Missing column " << colB << " for probe " << req.name << "\n";
                        continue;
                    }
                    req.kind = ColumnRequest::Kind::Diff;
                    req.idxB = itB->second;
                } else {
                    req.kind = ColumnRequest::Kind::Direct;
                    req.idxB = -1;
                }
            } else if (p.kind == ProbeKind::BranchCurrent) {
                std::string col = "I(" + p.eleName + ")";
                auto it = colMap.find(col);
                if (it == colMap.end()) {
                    std::cerr << "[CSV] Missing column " << col << " for probe " << req.name << "\n";
                    continue;
                }
                req.kind = ColumnRequest::Kind::Direct;
                req.idxA = it->second;
            } else {
                continue;
            }
            cols.push_back(req);
        }
    }

    if (cols.empty()) {
        std::cerr << "[CSV] No valid probes for " << dstPath << ".\n";
        return false;
    }

    std::filesystem::path outPath(dstPath);
    std::filesystem::create_directories(outPath.parent_path());

    std::string tmpPath = dstPath;
    bool samePath = (rawPath == dstPath);
    if (samePath) {
        tmpPath = dstPath + ".tmp";
    }

    std::ofstream ofs(tmpPath);
    if (!ofs) {
        std::cerr << "Cannot open output CSV: " << dstPath << "\n";
        return false;
    }
    ofs << std::scientific << std::setprecision(9);

    // header
    ofs << "time";
    for (const auto& c : cols) {
        ofs << "," << c.name;
    }
    ofs << "\n";

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        auto fields = splitCsv(line);
        double tval = 0.0;
        if (timeIdx >= 0 && timeIdx < (int)fields.size()) {
            tval = parseNum(fields[timeIdx]);
        }
        ofs << tval;
        for (const auto& c : cols) {
            double va = 0.0;
            if (c.idxA >=0 && c.idxA < (int)fields.size()) va = parseNum(fields[c.idxA]);
            double vb = 0.0;
            if (c.kind == ColumnRequest::Kind::Diff) {
                if (c.idxB >=0 && c.idxB < (int)fields.size()) vb = parseNum(fields[c.idxB]);
            }
            double v = (c.kind == ColumnRequest::Kind::Direct) ? va : (va - vb);
            ofs << "," << v;
        }
        ofs << "\n";
    }

    if (samePath) {
        std::filesystem::rename(tmpPath, dstPath);
    }
    return true;
}

bool runPythonPlot(const std::string& csvPath,
                   const std::vector<std::string>& columns,
                   const std::string& outPng) {
    if (columns.empty()) return true;
    
    // 构建Python命令参数列表（避免命令行字符串拼接）
    std::vector<std::string> args;
    args.push_back("plot_tran.py");
    args.push_back(csvPath);
    
    for (const auto& c : columns) {
        args.push_back(c);
    }
    
    if (!outPng.empty()) {
        args.push_back("--out");
        args.push_back(outPng);
    }
    
    // 使用更可靠的方式执行Python
    std::string command;
#ifdef _WIN32
    command = "python";
#else
    command = "python3";
#endif
    
    // 构建完整的命令行
    std::ostringstream fullCmd;
    fullCmd << command;
    for (const auto& arg : args) {
        fullCmd << " \"" << arg << "\"";
    }
    
    // 在Windows上不重定向输出，以便看到错误信息
#ifndef _WIN32
    fullCmd << " > /dev/null 2>&1";
#endif
    
    std::cout << "Executing: " << fullCmd.str() << std::endl;
    
    int rc = std::system(fullCmd.str().c_str());
    if (rc != 0) {
        std::cerr << "Plot script failed with code: " << rc << std::endl;
        
        // 检查文件是否存在
        if (!std::filesystem::exists(csvPath)) {
            std::cerr << "CSV file does not exist: " << csvPath << std::endl;
        } else {
            std::cerr << "CSV file exists, size: " << std::filesystem::file_size(csvPath) << " bytes" << std::endl;
        }
        return false;
    }
    return true;
}

// 写出“全量” DC/OP CSV（包含所有节点电压与 Vsrc/L 支路电流）
bool writeFullStaticCsv(const Circuit& ckt,
                        const Eigen::VectorXd& x,
                        const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream ofs(path);
    if (!ofs) return false;
    ofs << std::scientific << std::setprecision(9);

    ofs << "t";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) ofs << ",V(" << node.name << ")";
    }
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            ofs << ",I(" << vs->getName() << ")";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            ofs << ",I(" << L->getName() << ")";
        }
    }
    ofs << "\n";

    ofs << 0.0;
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0) {
            double v = (node.eqIndex < x.size()) ? x(node.eqIndex) : 0.0;
            ofs << "," << v;
        }
    }
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            int k = vs->getBranchEqIndex();
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            ofs << "," << Ibr;
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            int k = L->getBranchEqIndex();
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            ofs << "," << Ibr;
        }
    }
    ofs << "\n";
    return true;
}

} // namespace

Runner::Runner(const std::string& outDir_, bool verbose_)
    : outDir(outDir_), verbose(verbose_) {}

int Runner::run(const std::string& netlistPath) {
    Circuit ckt;
    SimulationConfig sim;

    if (verbose) {
        std::cout << "Reading netlist: " << netlistPath << "\n";
    }
    if (!parseNetlist(netlistPath, ckt, sim)) {
        std::cerr << "Parse netlist failed.\n";
        return 1;
    }
    ckt.assignEquationIndices();
    sim.verbose = verbose;

    std::filesystem::create_directories(outDir);
    std::filesystem::path caseStem = std::filesystem::path(netlistPath).stem();
    auto csvPath = [&](const std::string& tag) {
        return (std::filesystem::path(outDir) / (caseStem.string() + "_" + tag + ".csv")).string();
    };

    auto rawPath = [&](const std::string& tag) {
        return (std::filesystem::path(outDir) / (caseStem.string() + "_" + tag + "_raw.csv")).string();
    };

    Eigen::VectorXd opSolution;
    bool opReady = false;
    std::mutex opMutex;
    auto ensureOp = [&]() -> const Eigen::VectorXd& {
        std::lock_guard<std::mutex> lock(opMutex);
        if (!opReady) {
            DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
            opSolution = dc.run();
            opReady = true;
        }
        return opSolution;
    };

    // OP
    if (sim.doOp) {
        if (sim.verbose) {
            std::cout << "[OP] Running operating point\n";
        }
        const Eigen::VectorXd& xop = ensureOp();
        std::string raw = rawPath("op");
        std::string finalCsv = csvPath("op");
        if (!writeFullStaticCsv(ckt, xop, raw)) {
            std::cerr << "[OP] Failed to write raw CSV\n";
            return 2;
        }
        auto probes = mergeProbes(collectCsvProbes(sim, AnalysisType::OP),
                                  collectPlotProbes(sim, AnalysisType::OP));
        bool ok = filterCsv(raw, finalCsv, probes);
        if (!ok) return 3;
    }

    // DC sweep 未实现
    if (!sim.dcSweeps.empty()) {
        std::cerr << "[DC] .DC sweep is not implemented; skipped.\n";
    }

    // 并行运行可相互独立的分析（目前 TRAN 与 HB）
    std::vector<std::future<int>> tasks;

    auto runTranTask = [&]() -> int {
        if (!sim.tran.enabled) return 0;
        if (sim.verbose) {
            std::cout << "[TRAN] Running transient (TR)\n";
        }
        std::string raw = rawPath("tran");
        std::string finalCsv = csvPath("tran");
        try {
            TransientAnalysis tran(ckt, sim, raw);
            tran.runTrapezoidal();
        } catch (const std::exception& e) {
            std::cerr << "[TRAN] Exception: " << e.what() << "\n";
            return 4;
        }
        auto plotProbes = collectPlotProbes(sim, AnalysisType::TRAN);
        auto probes = mergeProbes(collectCsvProbes(sim, AnalysisType::TRAN), plotProbes);
        bool ok = filterCsv(raw, finalCsv, probes);
        if (!ok) return 5;

        std::vector<std::string> plotCols;
        auto header = readHeader(finalCsv);
        std::unordered_set<std::string> avail(header.begin(), header.end());
        for (const auto& p : plotProbes) {
            std::string name = probeHeaderName(p);
            if (avail.count(name)) {
                plotCols.push_back(name);
            } else {
                std::cerr << "[TRAN] Probe column not in CSV: " << name << "\n";
            }
        }
        if (!plotCols.empty() && sim.enablePlot) {
            std::string pngPath = (std::filesystem::path(outDir) / (caseStem.string() + "_tran_probe.png")).string();
            runPythonPlot(finalCsv, plotCols, pngPath);
        }
        return 0;
    };

    auto runHbTask = [&]() -> int {
        if (!sim.hb.enabled) return 0;
        if (sim.verbose) {
            std::cout << "[HB] Running harmonic balance\n";
        }
        std::string raw = rawPath("hb");
        std::string finalCsv = csvPath("hb");
        try {
            const Eigen::VectorXd& xdc = ensureOp();
            HbAnalysis hb(ckt, sim, xdc);
            Eigen::VectorXd xhb;
            bool ok = hb.run(xhb, raw);
            if (!ok) {
                std::cerr << "[HB] Convergence failed.\n";
                return 6;
            }
        } catch (const std::exception& e) {
            std::cerr << "[HB] Exception: " << e.what() << "\n";
            return 7;
        }
        auto plotProbes = collectPlotProbes(sim, AnalysisType::HB);
        auto probes = mergeProbes(collectCsvProbes(sim, AnalysisType::HB), plotProbes);
        bool ok = filterCsv(raw, finalCsv, probes);
        if (!ok) return 8;

        std::vector<std::string> plotCols;
        auto header = readHeader(finalCsv);
        std::unordered_set<std::string> avail(header.begin(), header.end());
        for (const auto& p : plotProbes) {
            std::string name = probeHeaderName(p);
            if (avail.count(name)) {
                plotCols.push_back(name);
            } else {
                std::cerr << "[HB] Probe column not in CSV: " << name << "\n";
            }
        }
        if (!plotCols.empty() && sim.enablePlot) {
            std::string pngPath = (std::filesystem::path(outDir) / (caseStem.string() + "_hb_probe.png")).string();
            runPythonPlot(finalCsv, plotCols, pngPath);
        }
        return 0;
    };

    if (sim.tran.enabled) {
        tasks.push_back(std::async(std::launch::async, runTranTask));
    }
    if (sim.hb.enabled) {
        tasks.push_back(std::async(std::launch::async, runHbTask));
    }

    for (auto& f : tasks) {
        int rc = f.get();
        if (rc != 0) return rc;
    }

    if (sim.ac.enabled) {
        std::cerr << "[AC] AC analysis not implemented; skipped.\n";
    }

    if (sim.verbose) {
        std::cout << "All requested analyses finished. Outputs in " << outDir << "\n";
    }
    return 0;
}
