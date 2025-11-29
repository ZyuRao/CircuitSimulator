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
//新增支持：Vname np nm SIN v0 va freq [td [phi]]
void NetlistParser::parseVoltageSource(const Statement& st) {
    const auto& t = st.tokens;
    if (t.size() < 4) {
        std::cerr << "Line " << st.lineNo
                  << ": invalid voltage source: " << st.raw << "\n";
        return;
    }

    SourceSpec spec;
    int idx = 3; // 当前要解析的 token 下标

    // 解析 DC 部分（支持三种写法）:
    //   Vname np nm value
    //   Vname np nm DC value
    //   Vname np nm SIN ...
    try {
        if (t.size() >= 5 && toLower(t[3]) == "dc") {
            // V ... DC value [后面可再跟 SIN ...]
            spec.dcValue = parseSpiceNumber(t[4]);
            idx = 5;
        } else {
            std::string low3 = toLower(t[3]);
            if (low3 == "sin") {
                // V ... SIN v0 va freq [td [phi]]
                spec.dcValue = 0.0;
                idx = 3; // 让后面统一从 "SIN" 开始解析
            } else {
                // V ... <number> [可能后面再跟 SIN ...]
                spec.dcValue = parseSpiceNumber(t[3]);
                idx = 4;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Line " << st.lineNo
                  << ": cannot parse V DC value: " << e.what()
                  << " in '" << st.raw << "'\n";
        return;
    }

    // 小工具：解析 "SIN v0 va freq [td [phi]]"
    auto parseSIN = [&](int sinIdx) {
        if (toLower(t[sinIdx]) != "sin") return;

        int need = sinIdx + 3; // 至少 SIN + 3 个参数
        if (t.size() < (std::size_t)need + 1) {
            std::cerr << "Line " << st.lineNo
                      << ": SIN needs at least 3 parameters (v0 va freq): "
                      << st.raw << "\n";
            return;
        }

        SinSpec sin;
        try {
            sin.v0   = parseSpiceNumber(t[sinIdx + 1]);
            sin.va   = parseSpiceNumber(t[sinIdx + 2]);
            sin.freq = parseSpiceNumber(t[sinIdx + 3]);

            if ((int)t.size() > sinIdx + 4) {
                sin.td = parseSpiceNumber(t[sinIdx + 4]);
            }
            if ((int)t.size() > sinIdx + 5) {
                sin.phi = parseSpiceNumber(t[sinIdx + 5]);
            }
        } catch (const std::exception& e) {
            std::cerr << "Line " << st.lineNo
                      << ": cannot parse SIN parameters: " << e.what()
                      << " in '" << st.raw << "'\n";
            return;
        }

        spec.tran.type = WaveformType::SIN;
        spec.tran.sine = sin;
    };

    // 继续解析瞬态波形（只考虑 SIN）
    if (idx < (int)t.size()) {
        if (toLower(t[idx]) == "sin") {
            parseSIN(idx);
        }
    }

    ckt.addVoltageSource(t[0], t[1], t[2], spec);
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

AnalysisType NetlistParser::parseAnalysisTypeFromToken(
    const std::string& tok
) const {
    std::string t = toLower(tok);

    if (t == "op")   return AnalysisType::OP;
    if (t == "dc")   return AnalysisType::DC;
    if (t == "ac")   return AnalysisType::AC;
    if (t == "tran") return AnalysisType::TRAN;
    return AnalysisType::NONE;
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
    m.name = t[1]; // 对应 .MODEL 1 / .MODEL 2 里的 "1" / "2"

    for (size_t i = 2; i + 1 < t.size(); i += 2) {
        std::string key = toLower(t[i]);
        double val = 0.0;
        try {
            val = parseSpiceNumber(t[i + 1]);
        } catch (const std::exception& e) {
            std::cerr << "Line " << st.lineNo
                      << ": cannot parse .MODEL param " << t[i]
                      << " = " << t[i+1] << " : " << e.what() << "\n";
            return;
        }

        if      (key == "vt")     m.VT      = val;
        else if (key == "mu")     m.MU      = val;
        else if (key == "cox")    m.COX     = val;
        else if (key == "lambda") m.LAMBDA  = val;
        else if (key == "cj0" || key == "cjo") m.CJO = val;
    }

    // VT 的符号决定 NMOS / PMOS，内部存成正数
    if (m.VT < 0.0) {
        m.isP = true;
        m.VT  = -m.VT;
    } else {
        m.isP = false;
    }

    ckt.addMosModel(m);
}