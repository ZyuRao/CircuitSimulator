#include "parser.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <stdexcept>

NetlistParser::NetlistParser(Circuit& circuit, SimulationConfig& simConfig)
    : ckt(circuit), sim(simConfig) {}
    
std::string NetlistParser::rtrimLocal(const std::string& s)
{
    auto pos = s.find_last_not_of(" \t\r\n");
    if(pos == std::string::npos) return "";

    return s.substr(0, pos + 1);
}

std::string NetlistParser::stripInlineComment(const std::string& s)
{
    auto pos = s.find('$');
    if(pos == std::string::npos) return s;

    return s.substr(0,pos);
}

bool NetlistParser::isFullLineComment(const std::string& s) {
    std::size_t i = 0;
    while(i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;

    if(i >= s.size()) return false;
    char c = s[i];
    return (c == '*' || c == ';');
}


// ---- 词法分析：合并续行 + 按空白切 token ----
bool NetlistParser::parseFile(const std::string& filename) {
    std::ifstream fin(filename);
    if(!fin) {
        std::cerr << "无法打开网表文件" << filename << "\n";
        return false;
    }

    sourceName = filename;
    return parseStream(fin, filename);
}

bool NetlistParser::parseStream(std::istream& in, const std::string& originName)
{
    sourceName = originName;
    lex(in);
    parseStatements();
    return true;
}

void NetlistParser::lex(std::istream& in) {
    stmts.clear();

    std::string physical;
    std::string logical;
    int logicalStartLine = 0;
    int lineNo = 0;

    auto flushlogical =[&]() {
        if(logical.empty()) return;

        std::string s = stripInlineComment(logical);
        s = rtrimLocal(ltrim(s));
        if(s.empty()) {
            logical.clear();
            return;
        }

        Statement st;
        st.lineNo = logicalStartLine;
        st.raw = s;

        std::istringstream iss(s);
        std::string tok;
        while(iss >> tok) {
            st.tokens.push_back(tok);
        }
        if(!st.tokens.empty()) {
            stmts.push_back(std::move(st));
        }
        logical.clear();
    };

    while(std::getline(in, physical)) {
        ++lineNo;

        if(!physical.empty() && physical.back() == '\r') {
            physical.pop_back();
        }
        
        std::string s = stripInlineComment(physical);
        s = rtrimLocal(ltrim(s));
        if(s.empty()) continue;
        
        if(isFullLineComment(s)) continue;

        std::size_t i = 0;
        while(i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;

        bool isCont = (i < s.size() && s[i] == '+');

        if(isCont) {
            std::string rest = s.substr(i + 1);
            rest = ltrim(rest);

            if(!logical.empty()) {
                logical += " ";
                logical += rest;
            } else {
                // 没有前一行却遇到 '+'
                logicalStartLine = lineNo;
                logical = rest;
            } 
        } else {
            // 新逻辑行开始
            if (!logical.empty()) {
                flushlogical();
            }
            logicalStartLine = lineNo;
            logical = s;
        }
    }

    if(!logical.empty()) {
        flushlogical();
    }
}

void NetlistParser::parseStatements() {
    bool titleConsumed = false;

    for(const auto& st : stmts) {
        if(st.tokens.empty()) continue;

        const std::string& head = st.tokens[0];
        if(head.empty()) continue;

        if(head[0] == '.') {
            parseDotCard(st);
            continue;
        }

        parseTitle(st, titleConsumed);
        if(!titleConsumed && st.tokens.empty()) {
            continue;
        }

        if(!st.tokens.empty()) {
            parseDeviceStmt(st);
        }
    }

    sim.ensureDefaultOp();
}

void NetlistParser::parseTitle(const Statement& st, bool& titleConsumed) {
    if(!titleConsumed || st.tokens.empty()) return;

    const std::string& head = st.tokens[0];
    char c0 = static_cast<char>(
        std::toupper(static_cast<unsigned char>(head[0]))
    );

    if(c0 != 'R' && c0 != 'C' && c0 != 'L' &&
        c0 != 'V' && c0 != 'I' && c0 != 'M' &&
        head[0] != '.') {
            titleConsumed = true;
    }

}

void NetlistParser::parseDeviceStmt(const Statement& st) {
    if(st.tokens.empty()) return;

    char c0 = static_cast<char>(
        std::toupper(static_cast<unsigned char>(st.tokens[0][0]))
    );

    switch (c0) {
        case 'R': parseResistor(st);       break;
        case 'C': parseCapacitor(st);      break;
        case 'L': parseInductor(st);       break;
        case 'V': parseVoltageSource(st);  break;
        case 'I': parseCurrentSource(st);  break;
        case 'M': parseMosfet(st);         break;
        default:
            std::cerr << "Line " << st.lineNo
                      << ": unsupported element or syntax: " << st.raw << "\n";
            break;
    }
}

void NetlistParser::parseResistor(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo << ": invalid resistor: " << st.raw << "\n";
        return;
    }
    double val = 0.0;
    try {
        val = parseSpiceNumber(t[3]);
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo << ": cannot parse R value: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    ckt.addResistor(t[0], t[1], t[2], val);
}

void NetlistParser::parseCapacitor(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo << ": invalid capacitor: " << st.raw << "\n";
        return;
    }
    double val = 0.0;
    try {
        val = parseSpiceNumber(t[3]);
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo << ": cannot parse C value: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }
    ckt.addCapacitor(t[0], t[1], t[2], val);
}

void NetlistParser::parseInductor(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo << ": invalid inductor: " << st.raw << "\n";
        return;
    }
    double val = 0.0;
    try {
        val = parseSpiceNumber(t[3]);
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo << ": cannot parse L value: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }
    ckt.addInductor(t[0], t[1], t[2], val);
}
//支持：Vname np nm value / Vname np nm DC value
// 支持：
//  Vname np nm value                （等价于 DC value）
//  Vname np nm DC value
//  Vname np nm SIN voff vamp freq [phase]
void NetlistParser::parseVoltageSource(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo << ": invalid voltage source: "
                  << st.raw << "\n";
        return;
    }
<<<<<<< HEAD

    const std::string& name = t[0];
    const std::string& np   = t[1];
    const std::string& nm   = t[2];

    // 关键字在 t[3]
    std::string kw = toLower(t[3]);

    auto parseNum = [&](const std::string& s, double& out) -> bool {
        try {
            out = parseSpiceNumber(s);
            return true;
        } catch (...) {
            return false;
        }
    };

    if (kw == "dc") {
        if (t.size() < 5) {
            std::cerr << "Line " << st.lineNo << ": DC needs a value: "
                      << st.raw << "\n";
            return;
        }
        double dc = 0.0;
        if (!parseNum(t[4], dc)) {
            std::cerr << "Line " << st.lineNo
                      << ": cannot parse DC value in voltage source: "
                      << st.raw << "\n";
            return;
        }
        ckt.addVoltageSource(name, np, nm, dc);
        return;
    }

    if (kw == "sin") {
        if (t.size() < 7) {
            std::cerr << "Line " << st.lineNo
                      << ": SIN needs at least voff, vamp, freq: "
                      << st.raw << "\n";
            return;
        }
        double voff  = 0.0;
        double vamp  = 0.0;
        double freq  = 0.0;
        double phase = 0.0;

        if (!parseNum(t[4], voff) ||
            !parseNum(t[5], vamp) ||
            !parseNum(t[6], freq)) {
            std::cerr << "Line " << st.lineNo
                      << ": cannot parse SIN parameters: " << st.raw << "\n";
            return;
        }
        if (t.size() >= 8) {
            parseNum(t[7], phase); // 失败就保留 0
=======
    SourceSpec spec;
    try {
        if(t.size() >= 5 && toLower(t[3]) == "dc") {
            spec.dcValue = parseSpiceNumber(t[4]);
        } else {
            spec.dcValue = parseSpiceNumber(t[3]);
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
        }

        ckt.addVoltageSourceSin(name, np, nm, voff, vamp, freq, phase);
        return;
    }

    // 兼容简写：Vname np nm value
    double dc = 0.0;
    if (!parseNum(t[3], dc)) {
        std::cerr << "Line " << st.lineNo
                  << ": cannot parse V value: " << t[3]
                  << " in '" << st.raw << "'\n";
        return;
    }
<<<<<<< HEAD
    ckt.addVoltageSource(name, np, nm, dc);
=======

    ckt.addVoltageSource(t[0], t[1], t[2], spec);
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
}



void NetlistParser::parseCurrentSource(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo << ": invalid current source: " << st.raw << "\n";
        return;
    }

    SourceSpec spec;
    try {
        if (t.size() >= 5 && toLower(t[3]) == "dc") {
            spec.dcValue = parseSpiceNumber(t[4]);
        } else {
            spec.dcValue = parseSpiceNumber(t[3]);
        }
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo << ": cannot parse I value: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    ckt.addCurrentSource(t[0], t[1], t[2], spec);
}

void NetlistParser::parseMosfet(const Statement& st) {
    const auto& t = st.tokens;

    // 原支持格式：  M name nd ng ns model W L  （7 个 token）
    // 助教给的格式：M name nd ng ns type W L modelId  （8 个 token）
    if (t.size() != 7 && t.size() != 8) {
        std::cerr << "Line " << st.lineNo
                  << ": invalid MOSFET: " << st.raw << "\n";
        return;
    }

    const std::string& name = t[0];
    const std::string& nd   = t[1];
    const std::string& ng   = t[2];
    const std::string& ns   = t[3];

    std::string modelId;
    if (t.size() == 7) {
        // 教材原始语法：M name nd ng ns model W L
        modelId = t[4];
    } else {
        // buffer.sp 这种：M name nd ng ns p/n W L modelId
        // 忽略 t[4] 的 p/n，把最后一个 token 当模型 ID（1 / 2）
        modelId = t.back();
    }

    double W = 0.0;
    double L = 0.0;
    try {
        W = parseSpiceNumber(t[5]);
        L = parseSpiceNumber(t[6]);
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo
                  << ": cannot parse MOS W/L: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    ckt.addMosfet(name, nd, ng, ns, modelId, W, L);
}


AcSweepType NetlistParser::parseAcSweepTypeFromToken(const std::string& tok) const {
    std::string s = toLower(tok);
    if (s == "lin") return AcSweepType::LIN;
    if (s == "oct") return AcSweepType::OCT;
    return AcSweepType::DEC;  // 默认 DEC
}

void NetlistParser::parseDotCard(const Statement& st) {
    if (st.tokens.empty()) return;
    std::string head = toLower(st.tokens[0]);

    if (head == ".op") {
        parseOpCard(st);
    } else if (head == ".dc") {
        parseDcCard(st);
    } else if (head == ".tran") {
        parseTranCard(st);
    } else if (head == ".ac") {
        parseAcCard(st);
    } else if (head == ".print") {
        parsePrintCard(st);
    } else if (head == ".model") {
        parseModelCard(st);
    } else {
        std::cerr << "Line " << st.lineNo
                  << ": unsupported control card: " << st.raw << "\n";
    }
}


void NetlistParser::parseOpCard(const Statement& st) {
    (void)st;
    sim.doOp = true;
}

void NetlistParser::parseDcCard(const Statement& st) {
    const auto& t = st.tokens;
    if(t.size() < 5) {
        std::cerr << "Line " << st.lineNo << ": invalid .DC syntax: " << st.raw << "\n";
        return;
    }

    DCSweepConfig dc;
    dc.sourceName = t[1];
    try {
        dc.start = parseSpiceNumber(t[2]);
        dc.stop = parseSpiceNumber(t[3]);
        dc.step = parseSpiceNumber(t[4]);
    } catch(const std::exception& e) {
        std::cerr << "Line " << st.lineNo << ": cannot parse .DC numbers: "
                << e.what() << " in '" << st.raw << "'\n";
        return;
    }
    sim.dcSweeps.push_back(dc);
}

void NetlistParser::parseTranCard(const Statement& st) {
    const auto& t = st.tokens;
    if(t.size() < 3) {
        std::cerr << "Line " << st.lineNo << ": invalid .TRAN syntax: "
                  << st.raw << "\n";
        return;
    }

    TranConfig cfg;

    try {
        cfg.tstep = parseSpiceNumber(t[1]);
        cfg.tstop = parseSpiceNumber(t[2]);
        if(t.size() >= 4) {
            cfg.tstart = parseSpiceNumber(t[3]);
        } else {
            cfg.tstart = 0.0;
        }
    } catch(const std::exception& e) {
        std::cerr << "Line " << st.lineNo
                  << ": cannot parse .TRAN numbers: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    cfg.enabled = true;
    sim.tran = cfg;
}

void NetlistParser::parseAcCard(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 5) {
        std::cerr << "Line " << st.lineNo << ": invalid .AC syntax: " << st.raw << "\n";
        return;
    }

    AcConfig cfg;
    cfg.sweepType = parseAcSweepTypeFromToken(t[1]);

    try {
        cfg.nPoints = std::stoi(t[2]);
        cfg.fstart  = parseSpiceNumber(t[3]);
        cfg.fstop   = parseSpiceNumber(t[4]);
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo
                  << ": cannot parse .AC arguments: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    cfg.enabled = true;
    sim.ac = cfg;
}

//.PRINT解析

ProbeSpec NetlistParser::parseProbeToken(const std::string& token) const {
    ProbeSpec p;
    p.expr = token;
    
    if(token.empty()) return p;

    char c0 = static_cast<char>(
        std::toupper(static_cast<unsigned char>(token[0]))
    );

    auto findParen = [&](const std::string& s) -> std::pair<int,int> {
        int l = -1, r = -1;
        for(int i = 0; i <(int)s.size(); ++i) {
            if(s[i] == '(' && l == -1) l = i;
            if(s[i] == ')') r = i;
        }
        return {l, r};
    };

    auto trimSpaces = [&](std::string s) {
        s = rtrimLocal(ltrim(s));
        return s;
    };

    if(c0 == 'V') {
        p.kind = ProbeKind::NodeVoltage;
        auto lr = findParen(token);
        if(lr.first >= 0 && lr.second > lr.first + 1) {
            std::string inside = token.substr(lr.first + 1, lr.second - lr.first - 1);
            auto commaPos = inside.find(',');
            if(commaPos == std::string::npos) {
                p.node1 = trimSpaces(inside);
                p.node2.clear();
                p.kind = ProbeKind::NodeVoltage;
            } else {
                p.node1 = trimSpaces(inside.substr(0, commaPos));
                p.node2 = trimSpaces(inside.substr(commaPos + 1));
                p.kind = ProbeKind::DiffVoltage;
            }
        }
    }

    else if(c0 == 'I') {
        p.kind = ProbeKind::BranchCurrent;
        auto lr = findParen(token);
        if(lr.first >= 0 && lr.second > lr.first + 1) {
            std::string inside = token.substr(lr.first + 1, lr.second - lr.first - 1);
            p.eleName = trimSpaces(inside);
        }
    }

    return p;
}

void NetlistParser::parsePrintCard(const Statement& st) {
    const auto& t = st.tokens;
    if(t.size() < 3) {
        std::cerr << "Line " << st.lineNo << ": invalid .PRINT: " << st.raw << "\n";
        return;
    }

    PrintCommand pc;

    pc.analysis = parseAnalysisTypeFromToken(t[1]);
    if(pc.analysis == AnalysisType::NONE) {
        std::cerr << "Line " << st.lineNo << ": unknown analysis type in .PRINT: "
                  << t[1] << " in '" << st.raw << "'\n";
        return;
    }

    for(std::size_t i = 2; i < t.size(); ++i) {
        ProbeSpec p = parseProbeToken(t[i]);
        pc.probes.push_back(std::move(p));
    }

    sim.printCommands.push_back(std::move(pc));
}

void NetlistParser::parseModelCard(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo
                  << ": invalid .MODEL: " << st.raw << "\n";
        return;
    }

    MosModel m;
    m.name = t[1];

    for (std::size_t i = 2; i + 1 < t.size(); i += 2) {
        std::string key = toLower(t[i]);
        double val = 0.0;
        try {
            val = parseSpiceNumber(t[i + 1]);
        } catch (const std::exception& e) {
            std::cerr << "Line " << st.lineNo
                      << ": cannot parse .MODEL param value '"
                      << t[i + 1] << "' in '" << st.raw
                      << "': " << e.what() << "\n";
            continue;
        }

        if      (key == "vt")      m.VT      = val;
        else if (key == "mu")      m.MU      = val;
        else if (key == "cox")     m.COX     = val;
        else if (key == "lambda")  m.LAMBDA  = val;
        else if (key == "cjo"
              || key == "cj0")     m.CJO     = val;
    }

    m.isP = (m.VT < 0.0);

    ckt.addMosModel(m);
}
