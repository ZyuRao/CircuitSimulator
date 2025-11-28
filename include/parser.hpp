#pragma once

#include <string>
#include <vector>
#include "circuit.hpp"

enum class AnalysisType {
    NONE,
    OP,
    DC,
    AC,
    TRAN
};

enum class AcSweepType {
    LIN,
    DEC,
    OCT
};

enum class ProbeKind {
    NodeVoltage,
    DiffVoltage,
    BranchCurrent
};

struct DCSweepConfig {
    std::string sourceName;
    double start = 0.0;
    double stop = 0.0;
    double step = 0.0;
};

struct TranConfig {
    bool enabled = false;
    double tstep = 0.0;
    double tstop = 0.0;
    double tstart = 0.0;
};

struct AcConfig {
    bool enabled = false;
    AcSweepType sweepType = AcSweepType::DEC;
    int nPoints = 0;
    double fstart = 0.0;
    double fstop = 0.0;
};

struct ProbeSpec {
    ProbeKind kind = ProbeKind::NodeVoltage;
    std::string expr;

    std::string node1;
    std::string node2;

    std::string eleName;
};

struct PrintCommand {
    AnalysisType analysis = AnalysisType::NONE;
    std::vector<ProbeSpec> probes;
};

class SimulationConfig {
public:
    bool doOp = false;
    std::vector<DCSweepConfig> dcSweeps;
    TranConfig tran;
    AcConfig ac;
    std::vector<PrintCommand> printCommands;


    bool hasAnyAnalysis() const {
        return doOp || !dcSweeps.empty() || tran.enabled || ac.enabled;
    }

    void ensureDefaultOp() {
        doOp = !hasAnyAnalysis();
    }
};

class NetlistParser {
public:
    NetlistParser(Circuit& circuit, SimulationConfig& simConfig);

    bool parseFile(const std::string& filename);

    bool parseStream(std::istream& in, const std::string& originName = "<stream>");

private:
    struct Statement {
        int lineNo = 0;
        std::string raw;
        std::vector<std::string> tokens;
    };

    Circuit& ckt;
    SimulationConfig& sim;
    std::string sourceName;
    std::vector<Statement> stmts;

    void lex(std::istream& in);

    void parseStatements();

     // ---- 工具函数 ----
    static std::string rtrimLocal(const std::string& s);
    static std::string stripInlineComment(const std::string& s);
    static bool isFullLineComment(const std::string& s);

    ProbeSpec parseProbeToken(const std::string& token) const;
    AnalysisType parseAnalysisTypeFromToken(const std::string& tok) const;
    AcSweepType  parseAcSweepTypeFromToken(const std::string& tok) const;

    void parseDeviceStmt(const Statement& st);
    void parseDotCard (const Statement& st);
    void parseTitle(const Statement& st, bool& titleConsumed);

    // 器件
    void parseResistor     (const Statement& st);
    void parseCapacitor    (const Statement& st);
    void parseInductor     (const Statement& st);
    void parseVoltageSource(const Statement& st);
    void parseCurrentSource(const Statement& st);
    void parseMosfet       (const Statement& st);

    // 控制卡
    void parseOpCard   (const Statement& st);
    void parseDcCard   (const Statement& st);
    void parseTranCard (const Statement& st);
    void parseAcCard   (const Statement& st);
    void parsePrintCard(const Statement& st);
    void parseModelCard(const Statement& st);

};

inline bool parseNetlist(
    const std::string& filename, Circuit& ckt,
    SimulationConfig& sim
) {
    NetlistParser parser(ckt, sim);
    bool ok = parser.parseFile(filename);
    sim.ensureDefaultOp();
    return ok;
}



