#pragma once
#include <string>
#include <vector>

struct Node {
    int id;
    std::string name;
    int eqIndex;                       // 对应 MNA 方程号（-1 为 GND）
    std::vector<int> attachedElements; // elements 中的下标

    Node(int i, const std::string& n)
        : id(i), name(n), eqIndex(-1) {}
};

std::string ltrim(const std::string& s) {
    std::size_t pos = s.find_first_not_of(" \t\r\n");
    if (pos == std::string::npos) return "";
    return s.substr(pos);
}

std::string toLower(const std::string& s) {
    std::string t = s;
    for (char &c : t) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return t;
}

// 支持科学计数法 + SPICE 后缀：10k, 1u, 3e12, 3.3meg 等
double parseSpiceNumber(const std::string& token) {
    std::string s = toLower(token);

    // 先尝试直接用 stod 解析（可以处理 3e12、-1.2e-3 等）
    try {
        std::size_t pos = 0;
        double base = std::stod(s, &pos);
        if (pos == s.size()) {
            // 完全解析，没有后缀
            return base;
        }
        // 有后缀，剩余部分当作 SPICE 后缀
        std::string suf = s.substr(pos);
        double factor = 1.0;
        if      (suf == "f")   factor = 1e-15;
        else if (suf == "p")   factor = 1e-12;
        else if (suf == "n")   factor = 1e-9;
        else if (suf == "u")   factor = 1e-6;
        else if (suf == "m")   factor = 1e-3;
        else if (suf == "k")   factor = 1e3;
        else if (suf == "meg") factor = 1e6;
        else if (suf == "g")   factor = 1e9;
        else if (suf == "t")   factor = 1e12;
        else                   factor = 1.0;
        return base * factor;
    } catch (...) {
        // stod 不认，比如 "10k" 这种，从头开始找后缀
        std::size_t pos = std::string::npos;
        for (std::size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if (std::isalpha(static_cast<unsigned char>(c))) {
                pos = i;
                break;
            }
        }
        if (pos == std::string::npos) {
            // 实在不行，返回 0
            return 0.0;
        }
        double base = std::stod(s.substr(0, pos));
        std::string suf = s.substr(pos);
        double factor = 1.0;
        if      (suf == "f")   factor = 1e-15;
        else if (suf == "p")   factor = 1e-12;
        else if (suf == "n")   factor = 1e-9;
        else if (suf == "u")   factor = 1e-6;
        else if (suf == "m")   factor = 1e-3;
        else if (suf == "k")   factor = 1e3;
        else if (suf == "meg") factor = 1e6;
        else if (suf == "g")   factor = 1e9;
        else if (suf == "t")   factor = 1e12;
        else                   factor = 1.0;
        return base * factor;
    }
}

bool isGroundName(const std::string& n) {
    std::string lower = toLower(n);
    return (lower == "0" || lower == "gnd");
}

