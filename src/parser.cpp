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


