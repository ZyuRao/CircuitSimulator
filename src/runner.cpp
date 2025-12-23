#include "parser.hpp"
#include "utils.hpp"
#include "element.hpp"
#include "runner.hpp"
#include "analysis.hpp"

#include <future>
#include <chrono>
#include <algorithm>
#include <cctype>
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
    std::ostringstream cmd;
    cmd << "python3 plot_tran.py '" << csvPath << "'";
    for (const auto& c : columns) {
        cmd << " '" << c << "'";
    }
    if (!outPng.empty()) {
        cmd << " --out='" << outPng << "'";
    }
    cmd << " > /dev/null";
    int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        std::cerr << "Plot script failed for " << csvPath << " (code " << rc << ")\n";
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

Runner::Runner(const std::string& outDir_) : outDir(outDir_) {}


namespace {

struct RunnerOptions {
    bool manualRunList = false;      // 若 true：忽略网表 enable/doOp，只按 runOp/runTran/runHb 运行
    bool runOp   = false;
    bool runTran = false;
    bool runHb   = false;

    DcSolverKind dcSolver = DcSolverKind::GaussSeidel;

    enum class TranMethod { TR, BE };
    TranMethod tranMethod = TranMethod::TR;

    bool parallel = true;
    bool timing   = true;

    bool printHelp = false;
};

static std::string toLowerCopy(std::string s) {
    for (char& ch : s) ch = (char)std::tolower((unsigned char)ch);
    return s;
}

static bool startsWith(const std::string& s, const std::string& pfx) {
    return s.size() >= pfx.size() && s.compare(0, pfx.size(), pfx) == 0;
}

static std::vector<std::string> splitCommaList(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == ',') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

static void applyRunList(RunnerOptions& opt, const std::string& listRaw) {
    std::string list = toLowerCopy(listRaw);
    auto items = splitCommaList(list);
    if (items.empty()) return;

    opt.manualRunList = true;
    opt.runOp = opt.runTran = opt.runHb = false;

    for (const auto& it : items) {
        if (it == "all") {
            opt.runOp = opt.runTran = opt.runHb = true;
            continue;
        }
        if (it == "auto") {
            opt.manualRunList = false; // 退回网表驱动
            continue;
        }
        if (it == "op" || it == "dc") opt.runOp = true;
        else if (it == "tran" || it == "transient") opt.runTran = true;
        else if (it == "hb") opt.runHb = true;
    }
}

static RunnerOptions parseRunnerArgs(int argc, char** argv) {
    RunnerOptions opt;
    if (argc <= 2) return opt;

    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        std::string al = toLowerCopy(a);

        if (al == "--help" || al == "-h") { opt.printHelp = true; continue; }
        if (al == "--no-parallel") { opt.parallel = false; continue; }
        if (al == "--parallel") { opt.parallel = true; continue; }
        if (al == "--no-timing") { opt.timing = false; continue; }
        if (al == "--timing") { opt.timing = true; continue; }

        // 位置参数：op/tran/hb/all
        if (al == "op" || al == "dc") { applyRunList(opt, "op"); continue; }
        if (al == "tran" || al == "transient") { applyRunList(opt, "tran"); continue; }
        if (al == "hb") { applyRunList(opt, "hb"); continue; }
        if (al == "all") { applyRunList(opt, "all"); continue; }

        if (startsWith(al, "--run=")) {
            applyRunList(opt, a.substr(std::string("--run=").size()));
            continue;
        }

        if (startsWith(al, "--dc-solver=")) {
            std::string v = toLowerCopy(a.substr(std::string("--dc-solver=").size()));
            if (v == "lu") opt.dcSolver = DcSolverKind::LU;
            else opt.dcSolver = DcSolverKind::GaussSeidel;
            continue;
        }

        if (startsWith(al, "--tran-method=")) {
            std::string v = toLowerCopy(a.substr(std::string("--tran-method=").size()));
            if (v == "be" || v == "backwardeuler") opt.tranMethod = RunnerOptions::TranMethod::BE;
            else opt.tranMethod = RunnerOptions::TranMethod::TR;
            continue;
        }

        // 未知参数：不失败，仅警告
        std::cerr << "[Runner] WARNING: unknown option '" << a << "' ignored.\n";
    }

    return opt;
}

static void printRunnerHelp(const char* argv0) {
    std::cout
        << "Usage:\n"
        << "  " << argv0 << " <netlist.sp> [op|tran|hb|all] [options]\n\n"
        << "Options:\n"
        << "  --run=op,tran,hb|all     Manually choose analyses to run (overrides netlist enable flags)\n"
        << "  --dc-solver=gs|lu        Choose DC/OP solver (default: gs)\n"
        << "  --tran-method=tr|be      Choose transient integrator (default: tr)\n"
        << "  --no-parallel            Run analyses sequentially\n"
        << "  --no-timing              Disable timing prints\n"
        << "  -h, --help               Show this help\n\n"
        << "Notes:\n"
        << "  - To enable CLI options, main.cpp should call: runner.run(argc, argv)\n"
        << "  - If main.cpp still calls runner.run(netlistPath), you can use env vars:\n"
        << "      SIM_RUN=op,tran,hb|all, SIM_DCSOLVER=gs|lu, SIM_TRAN=tr|be, SIM_TIMING=0/1, SIM_PARALLEL=0/1\n";
}

struct Timer {
    using Clock = std::chrono::steady_clock;
    Clock::time_point t0;
    explicit Timer() : t0(Clock::now()) {}
    double ms() const {
        auto t1 = Clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        return (double)us / 1000.0;
    }
};

static const char* dcSolverName(DcSolverKind k) {
    switch (k) {
    case DcSolverKind::LU: return "LU";
    case DcSolverKind::GaussSeidel:
    default: return "GaussSeidel";
    }
}

static const char* tranMethodName(RunnerOptions::TranMethod m) {
    return (m == RunnerOptions::TranMethod::BE) ? "BE" : "TR";
}

// 核心实现：把“旧 Runner::run()”逻辑搬到这里，只增加：手动选择 + 计时 + 方法选择
static int runImpl(const std::string& outDir, const std::string& netlistPath, RunnerOptions opt) {
    Timer totalTimer;

    Circuit ckt;
    SimulationConfig sim;

    Timer parseTimer;
    std::cout << "Reading netlist: " << netlistPath << "\n";
    if (!parseNetlist(netlistPath, ckt, sim)) {
        std::cerr << "Parse netlist failed.\n";
        return 1;
    }
    ckt.assignEquationIndices();
    if (opt.timing) {
        std::cout << "[TIME] Netlist parse + index: " << std::fixed << std::setprecision(3) << parseTimer.ms() << " ms\n";
    }

    // 运行计划：默认按网表；若 opt.manualRunList=true 则按用户选择
    bool wantOp   = opt.manualRunList ? opt.runOp   : sim.doOp;
    bool wantTran = opt.manualRunList ? opt.runTran : sim.tran.enabled;
    bool wantHb   = opt.manualRunList ? opt.runHb   : sim.hb.enabled;

    // “手动选择模式”的提示
    if (opt.manualRunList) {
        std::cout << "[Runner] Manual run list enabled: "
                  << (wantOp ? "OP " : "")
                  << (wantTran ? "TRAN " : "")
                  << (wantHb ? "HB " : "")
                  << "\n";
    }

    std::cout << "[Runner] DC solver = " << dcSolverName(opt.dcSolver)
              << ", TRAN method = " << tranMethodName(opt.tranMethod)
              << ", parallel = " << (opt.parallel ? "on" : "off")
              << ", timing = " << (opt.timing ? "on" : "off") << "\n";

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
    double opSolveMs = 0.0;

    auto ensureOp = [&]() -> const Eigen::VectorXd& {
        if (!opReady) {
            Timer t;
            DcAnalysis dc(ckt, sim, opt.dcSolver);
            opSolution = dc.run();
            opReady = true;
            opSolveMs = t.ms();
            if (opt.timing) {
                std::cout << "[TIME] DC/OP (" << dcSolverName(opt.dcSolver) << "): "
                          << std::fixed << std::setprecision(3) << opSolveMs << " ms\n";
            }
        }
        return opSolution;
    };

    // OP
    if (wantOp) {
        Timer tTotal;
        std::cout << "[OP] Running operating point\n";
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

        if (opt.timing) {
            std::cout << "[TIME] OP total: " << std::fixed << std::setprecision(3) << tTotal.ms() << " ms\n";
        }
    }

    // DC sweep 未实现
    if (!sim.dcSweeps.empty()) {
        std::cerr << "[DC] .DC sweep is not implemented; skipped.\n";
    }

    // 并行运行可相互独立的分析（目前 TRAN 与 HB）
    std::vector<std::future<int>> tasks;

    auto runTranTask = [&]() -> int {
        if (!wantTran || !sim.tran.enabled) {
            if (wantTran && !sim.tran.enabled) {
                std::cerr << "[TRAN] Requested but .TRAN missing; skipped.\n";
            }
            return 0;
        }

        Timer tAll;
        std::cout << "[TRAN] Running transient (" << tranMethodName(opt.tranMethod) << ")\n";
        std::string raw = rawPath("tran");
        std::string finalCsv = csvPath("tran");

        double solveMs = 0.0;
        try {
            Timer tSolve;
            TransientAnalysis tran(ckt, sim, raw);
            if (opt.tranMethod == RunnerOptions::TranMethod::BE) tran.runBackwardEuler();
            else tran.runTrapezoidal();
            solveMs = tSolve.ms();
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
        if (!plotCols.empty()) {
            std::string pngPath = (std::filesystem::path(outDir) / (caseStem.string() + "_tran_probe.png")).string();
            runPythonPlot(finalCsv, plotCols, pngPath);
        }

        if (opt.timing) {
            std::cout << "[TIME] TRAN solve: " << std::fixed << std::setprecision(3) << solveMs
                      << " ms, total: " << tAll.ms() << " ms\n";
        }
        return 0;
    };

    auto runHbTask = [&]() -> int {
        if (!wantHb || !sim.hb.enabled) {
            if (wantHb && !sim.hb.enabled) {
                std::cerr << "[HB] Requested but .HB missing; skipped.\n";
            }
            return 0;
        }

        Timer tAll;
        std::cout << "[HB] Running harmonic balance\n";
        std::string raw = rawPath("hb");
        std::string finalCsv = csvPath("hb");

        double solveMs = 0.0;
        try {
            const Eigen::VectorXd& xdc = ensureOp();
            Timer tSolve;
            HbAnalysis hb(ckt, sim, xdc);
            Eigen::VectorXd xhb;
            bool okRun = hb.run(xhb, raw);
            solveMs = tSolve.ms();
            if (!okRun) {
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
        if (!plotCols.empty()) {
            std::string pngPath = (std::filesystem::path(outDir) / (caseStem.string() + "_hb_probe.png")).string();
            runPythonPlot(finalCsv, plotCols, pngPath);
        }

        if (opt.timing) {
            std::cout << "[TIME] HB solve: " << std::fixed << std::setprecision(3) << solveMs
                      << " ms, total: " << tAll.ms() << " ms\n";
        }
        return 0;
    };

    if (opt.parallel) {
        if (wantTran && sim.tran.enabled) tasks.push_back(std::async(std::launch::async, runTranTask));
        if (wantHb && sim.hb.enabled)     tasks.push_back(std::async(std::launch::async, runHbTask));
    } else {
        int rc = runTranTask();
        if (rc != 0) return rc;
        rc = runHbTask();
        if (rc != 0) return rc;
    }

    for (auto& f : tasks) {
        int rc = f.get();
        if (rc != 0) return rc;
    }

    if (sim.ac.enabled) {
        std::cerr << "[AC] AC analysis not implemented; skipped.\n";
    }

    if (opt.timing) {
        std::cout << "[TIME] Total: " << std::fixed << std::setprecision(3) << totalTimer.ms() << " ms\n";
    }

    std::cout << "All requested analyses finished. Outputs in " << outDir << "\n";
    return 0;
}

} // namespace

int Runner::run(int argc, char** argv) {
    if (argc < 2 || argv == nullptr) {
        std::cerr << "Usage: sim <netlist.sp> [op|tran|hb|all] [options]\n";
        return 1;
    }

    RunnerOptions opt = parseRunnerArgs(argc, argv);
    if (opt.printHelp) {
        printRunnerHelp(argv[0]);
        return 0;
    }
    return runImpl(this->outDir, std::string(argv[1]), opt);
}

int Runner::run(const std::string& netlistPath) {
    // 兼容旧入口：保留“网表驱动”行为；同时支持用环境变量临时手动选择/切换方法
    RunnerOptions opt;

    if (const char* v = std::getenv("SIM_RUN")) {
        applyRunList(opt, v);
    }
    if (const char* v = std::getenv("SIM_DCSOLVER")) {
        std::string s = toLowerCopy(v);
        opt.dcSolver = (s == "lu") ? DcSolverKind::LU : DcSolverKind::GaussSeidel;
    }
    if (const char* v = std::getenv("SIM_TRAN")) {
        std::string s = toLowerCopy(v);
        opt.tranMethod = (s == "be") ? RunnerOptions::TranMethod::BE : RunnerOptions::TranMethod::TR;
    }
    if (const char* v = std::getenv("SIM_TIMING")) {
        std::string s = toLowerCopy(v);
        opt.timing = !(s == "0" || s == "false" || s == "off");
    }
    if (const char* v = std::getenv("SIM_PARALLEL")) {
        std::string s = toLowerCopy(v);
        opt.parallel = !(s == "0" || s == "false" || s == "off");
    }

    return runImpl(this->outDir, netlistPath, opt);
}


