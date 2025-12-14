// 编译、运行命令
// g++ -std=c++17 -O2 src\*.cpp -I include -I E:\eigen-5.0.0 -o CircuitSimulator.exe
// .\CircuitSimulator.exe tests\RLC_s3.sp dc dc_RLC_s3.csv --dc-solver=gs
// .\CircuitSimulator.exe tests\RLC_s3.sp dc dc_RLC_s3.csv --dc-solver=lu
// .\CircuitSimulator.exe tests\RLC_s3.sp tran tran_RLC_s3.csv --tran-method=tr
// .\CircuitSimulator.exe tests\RLC_s3.sp tran tran_RLC_s3.csv --tran-method=be
// .\CircuitSimulator.exe tests\RLC_s3.sp hb hb_RLC_s3.csv
// .\CircuitSimulator.exe tests\RLC_s3.sp shoot shoot_RLC_s3.csv --period=1e-6 --dt=5e-9




#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <cstdio>

#include "parser.hpp"
#include "circuit.hpp"
#include "analysis.hpp"
#include "sim.hpp"
#include "element.hpp"


static inline bool startsWith(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && std::equal(p.begin(), p.end(), s.begin());
}

static std::vector<std::string> splitCsvSimple(const std::string& line) {
    // 简化版 CSV：不支持引号/转义（你现在的输出不含引号，足够）
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

static void printDcLikeSolution(const Circuit& ckt, const Eigen::VectorXd& x) {
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n==== Node voltages ====\n";
    for (const auto& node : ckt.nodes) {
        if (node.eqIndex >= 0 && node.eqIndex < x.size()) {
            std::cout << "V(" << node.name << ") = " << x(node.eqIndex) << " V"
                      << "   [eqIndex=" << node.eqIndex << "]\n";
        } else {
            std::cout << "V(" << node.name << ") = 0.000000 V   [GND]\n";
        }
    }

    std::cout << "\n==== Branch currents (V sources / L) ====\n";
    for (const auto& e : ckt.elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            int k = vs->getBranchEqIndex();
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            std::cout << "I(" << vs->getName() << ") = " << Ibr << " A"
                      << "   [eqIndex=" << k << "]\n";
        } else if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            int k = L->getBranchEqIndex();
            double Ibr = (k >= 0 && k < x.size()) ? x(k) : 0.0;
            std::cout << "I(" << L->getName() << ") = " << Ibr << " A"
                      << "   [eqIndex=" << k << "]\n";
        }
    }
}

static bool writeDcCsv(const Circuit& ckt, const Eigen::VectorXd& x, const std::string& outCsv) {
    std::ofstream ofs(outCsv);
    if (!ofs) return false;
    ofs << std::scientific << std::setprecision(9);

    // header
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

    // single row, t=0
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

static bool readLastCsvRow(const std::string& csv,
                           std::vector<std::string>& headerOut,
                           std::vector<double>& valuesOut) {
    std::ifstream ifs(csv);
    if (!ifs) return false;

    std::string header;
    if (!std::getline(ifs, header)) return false;
    headerOut = splitCsvSimple(header);

    std::string line, last;
    while (std::getline(ifs, line)) {
        // 跳过空行
        bool allSpace = true;
        for (char c : line) if (!std::isspace(static_cast<unsigned char>(c))) { allSpace = false; break; }
        if (!allSpace) last = line;
    }
    if (last.empty()) return false;

    auto cols = splitCsvSimple(last);
    valuesOut.clear();
    valuesOut.reserve(cols.size());
    for (const auto& s : cols) {
        try {
            valuesOut.push_back(std::stod(s));
        } catch (...) {
            valuesOut.push_back(0.0);
        }
    }
    return true;
}

static void printLastTimepointVoltages(const std::string& csv) {
    std::vector<std::string> header;
    std::vector<double> vals;
    if (!readLastCsvRow(csv, header, vals)) {
        std::cerr << "[WARN] Cannot read last row from " << csv << "\n";
        return;
    }
    if (header.empty() || vals.empty()) return;

    double t = (vals.size() > 0) ? vals[0] : 0.0;
    std::cout << "\n==== Last timepoint (from CSV) ====\n";
    std::cout << "t = " << std::scientific << std::setprecision(9) << t << "\n";

    // 只打印 V(...) 列
    for (size_t i = 1; i < header.size() && i < vals.size(); ++i) {
        const std::string& name = header[i];
        if (startsWith(name, "V(")) {
            std::cout << name << " = " << vals[i] << " V\n";
        }
    }
}

// ----------------- CLI parsing -----------------
struct Cmd {
    std::string netlist;
    std::string mode;     // dc / tran / hb / shoot
    std::string outCsv;

    std::string dcSolver = "gs";   // gs / lu
    std::string tranMethod = "tr"; // tr / be

    // Optional overrides
    bool noCsv = false;

    // HB overrides
    double hbF0 = -1.0;
    int    hbNHarm = -1;

    // Shooting config
    double shootPeriodT = -1.0;
    double shootDt = -1.0;
    int    shootMaxIters = -1;
    double shootTol = -1.0;
    double shootRelax = -1.0;
};

static void printHelp() {
    std::cout <<
R"(Usage:
  CircuitSimulator <netlist.sp> [mode] [out.csv] [options]

Modes:
  dc       DC operating point (supports --dc-solver=gs|lu)
  tran     transient (uses .TRAN in netlist; supports --tran-method=tr|be)
  hb       harmonic balance steady-state (requires .hb; supports --f0, --nharm override)
  shoot    shooting steady-state (needs period; supports --period, --dt, --max-iters, --tol, --relax)

Options:
  --mode=dc|tran|hb|shoot
  --out=FILE.csv
  --no-csv

  --dc-solver=gs|lu
  --tran-method=tr|be

  --f0=Hz          (override .hb f0)
  --nharm=K        (override .hb nHarm)

  --period=T       (shooting period)
  --dt=dt          (shooting timestep; default T/200)
  --max-iters=N    (shooting outer iters)
  --tol=val
  --relax=val

Examples:
  CircuitSimulator tests/test.sp dc dc.csv --dc-solver=gs
  CircuitSimulator tests/test.sp tran tran.csv --tran-method=tr
  CircuitSimulator tests/test.sp hb hb.csv
  CircuitSimulator tests/test.sp shoot shoot.csv --period=1e-7 --dt=5e-10
)";
}

static bool parseCmd(int argc, char** argv, Cmd& cmd) {
    if (argc < 2) {
        printHelp();
        return false;
    }
    cmd.netlist = argv[1];

    int i = 2;
    // allow positional: [mode] [out.csv]
    if (i < argc && argv[i][0] != '-') {
        cmd.mode = toLower(argv[i]);
        ++i;
    }
    if (i < argc && argv[i][0] != '-') {
        cmd.outCsv = argv[i];
        ++i;
    }

    auto readValue = [&](const std::string& token, const std::string& key, std::string& valOut) -> bool {
        if (startsWith(token, key + "=")) {
            valOut = token.substr(key.size() + 1);
            return true;
        }
        return false;
    };

    for (; i < argc; ++i) {
        std::string tok = argv[i];
        if (tok == "--help" || tok == "-h") {
            printHelp();
            return false;
        }
        if (tok == "--no-csv") {
            cmd.noCsv = true;
            continue;
        }

        std::string v;
        if (readValue(tok, "--mode", v)) { cmd.mode = toLower(v); continue; }
        if (readValue(tok, "--out", v))  { cmd.outCsv = v; continue; }

        if (readValue(tok, "--dc-solver", v)) { cmd.dcSolver = toLower(v); continue; }
        if (readValue(tok, "--tran-method", v)) { cmd.tranMethod = toLower(v); continue; }

        if (readValue(tok, "--f0", v)) {
            cmd.hbF0 = std::stod(v);
            continue;
        }
        if (readValue(tok, "--nharm", v)) {
            cmd.hbNHarm = std::stoi(v);
            continue;
        }

        if (readValue(tok, "--period", v)) { cmd.shootPeriodT = std::stod(v); continue; }
        if (readValue(tok, "--dt", v))     { cmd.shootDt = std::stod(v); continue; }
        if (readValue(tok, "--max-iters", v)) { cmd.shootMaxIters = std::stoi(v); continue; }
        if (readValue(tok, "--tol", v))    { cmd.shootTol = std::stod(v); continue; }
        if (readValue(tok, "--relax", v))  { cmd.shootRelax = std::stod(v); continue; }

        // ignore unknowns but warn
        if (startsWith(tok, "--")) {
            std::cerr << "[WARN] Unknown option: " << tok << "\n";
        }
    }

    if (cmd.mode.empty()) cmd.mode = "dc";
    if (cmd.outCsv.empty() && !cmd.noCsv) {
        if (cmd.mode == "dc") cmd.outCsv = "dc.csv";
        else if (cmd.mode == "tran") cmd.outCsv = "tran.csv";
        else if (cmd.mode == "hb") cmd.outCsv = "hb.csv";
        else if (cmd.mode == "shoot") cmd.outCsv = "shoot.csv";
        else cmd.outCsv = "out.csv";
    }
    return true;
}

// ----------------- main -----------------
int main(int argc, char** argv) {
    Cmd cmd;
    if (!parseCmd(argc, argv, cmd)) return 1;

    Circuit ckt;
    SimulationConfig sim;

    std::cout << "Reading netlist: " << cmd.netlist << "\n";
    if (!parseNetlist(cmd.netlist, ckt, sim)) {
        std::cerr << "ERROR: parseNetlist() failed.\n";
        return 2;
    }

    // assign equation indices
    ckt.assignEquationIndices();

    std::cout << "\n==== Circuit summary ====\n";
    std::cout << "Nodes        : " << ckt.nodes.size() << "\n";
    std::cout << "Elements     : " << ckt.elements.size() << "\n";
    std::cout << "Unknowns     : " << ckt.numUnknowns()
              << "  (nodeEq=" << ckt.numNodeEquations()
              << ", branchEq=" << ckt.numVoltageBranches() << ")\n";

    // ------------- mode dispatch -------------
    const std::string mode = toLower(cmd.mode);

    if (mode == "dc") {
        DcSolverKind kind = DcSolverKind::GaussSeidel;
        if (cmd.dcSolver == "lu") kind = DcSolverKind::LU;

        std::cout << "\n[DC] solver = " << (kind == DcSolverKind::LU ? "LU" : "GaussSeidel") << "\n";

        Eigen::VectorXd xdc;
        try {
            DcAnalysis dc(ckt, sim, kind);
            xdc = dc.run();
        } catch (const std::exception& e) {
            std::cerr << "DC solve failed: " << e.what() << "\n";
            return 3;
        }

        printDcLikeSolution(ckt, xdc);

        if (!cmd.noCsv) {
            if (!writeDcCsv(ckt, xdc, cmd.outCsv)) {
                std::cerr << "[DC] Cannot write csv: " << cmd.outCsv << "\n";
                return 4;
            }
            std::cout << "[DC] CSV written: " << cmd.outCsv << "\n";
        }
        return 0;
    }

    if (mode == "tran") {
        // allow method select
        std::cout << "\n[TRAN] method = " << cmd.tranMethod << "\n";

        if (cmd.noCsv) {
            std::cerr << "[TRAN] --no-csv given, but current TransientAnalysis writes waveform to CSV.\n";
            std::cerr << "       Please provide --out=xxx.csv or omit --no-csv.\n";
            return 5;
        }

        try {
            TransientAnalysis tran(ckt, sim, cmd.outCsv);
            if (cmd.tranMethod == "be") tran.runBackwardEuler();
            else tran.runTrapezoidal();
        } catch (const std::exception& e) {
            std::cerr << "TRAN exception: " << e.what() << "\n";
            return 6;
        }

        std::cout << "[TRAN] CSV written: " << cmd.outCsv << "\n";

        // print last timepoint node voltages for accuracy check
        printLastTimepointVoltages(cmd.outCsv);
        return 0;
    }

    if (mode == "hb") {
        // override .hb if needed
        if (cmd.hbF0 > 0.0) {
            sim.hb.f0 = cmd.hbF0;
            sim.hb.enabled = true;
        }
        if (cmd.hbNHarm >= 0) {
            sim.hb.nHarm = cmd.hbNHarm;
            sim.hb.enabled = true;
        }

        if (!sim.hb.enabled) {
            std::cerr << "[HB] No .hb card found (and no override). Abort.\n";
            return 7;
        }
        if (cmd.noCsv) {
            std::cerr << "[HB] --no-csv given, but HB time waveform is written to CSV by design.\n";
            std::cerr << "     Please provide --out=xxx.csv or omit --no-csv.\n";
            return 8;
        }

        // DC op for initial
        Eigen::VectorXd xdc;
        try {
            DcAnalysis dc(ckt, sim, DcSolverKind::GaussSeidel);
            xdc = dc.run();
        } catch (const std::exception& e) {
            std::cerr << "[HB] DC operating point failed: " << e.what() << "\n";
            return 9;
        }

        std::cout << "\n[HB] f0=" << std::scientific << sim.hb.f0
                  << " Hz, nHarm=" << sim.hb.nHarm
                  << ", out=" << cmd.outCsv << "\n";

        Eigen::VectorXd xhb;
        try {
            HbAnalysis hb(ckt, sim, xdc);
            bool ok = hb.run(xhb, cmd.outCsv);
            if (!ok) {
                std::cerr << "[HB] Failed to converge.\n";
                return 10;
            }
        } catch (const std::exception& e) {
            std::cerr << "[HB] exception: " << e.what() << "\n";
            return 11;
        }


        std::cout << "[HB] CSV written: " << cmd.outCsv << "\n";
        return 0;
    }

    if (mode == "shoot") {
        if (cmd.noCsv) {
            std::cerr << "[SHOOT] --no-csv given, but shooting outputs one-period waveform to CSV.\n";
            std::cerr << "        Please provide --out=xxx.csv or omit --no-csv.\n";
            return 12;
        }

        ShootingPssConfig cfg;
        // infer period from --period, else from .hb if present
        if (cmd.shootPeriodT > 0.0) {
            cfg.periodT = cmd.shootPeriodT;
        } else if (sim.hb.enabled && sim.hb.f0 > 0.0) {
            cfg.periodT = 1.0 / sim.hb.f0;
            std::cout << "[SHOOT] period not given, inferred from .hb: T=" << std::scientific << cfg.periodT << "\n";
        } else {
            std::cerr << "[SHOOT] Need --period=T, or provide a valid .hb f0 to infer T=1/f0.\n";
            return 13;
        }

        cfg.tstep = (cmd.shootDt > 0.0) ? cmd.shootDt : 0.0;
        if (cmd.shootMaxIters > 0) cfg.maxIters = cmd.shootMaxIters;
        if (cmd.shootTol > 0.0) cfg.tol = cmd.shootTol;
        if (cmd.shootRelax > 0.0) cfg.relax = cmd.shootRelax;

        std::cout << "\n[SHOOT] out=" << cmd.outCsv << "\n";
        try {
            runPssShootingAnalysis(ckt, sim, cfg, cmd.outCsv);
        } catch (const std::exception& e) {
            std::cerr << "[SHOOT] exception: " << e.what() << "\n";
            return 14;
        }


        std::cout << "[SHOOT] CSV written: " << cmd.outCsv << "\n";
        return 0;
    }

    std::cerr << "Unknown mode: " << cmd.mode << "\n";
    printHelp();
    return 20;
}
