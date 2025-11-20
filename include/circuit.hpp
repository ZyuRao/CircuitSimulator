#pragma once

#include <vector>
#include <unordered_map>
#include <memory>
#include "element.hpp"
#include "utils.hpp"

class Circuit {
public:
    std::vector<Node> nodes;
    std::vector<std::shared_ptr<Element>> elements;
    std::unordered_map<std::string, int> nodeNameToId;

    // 获取或创建节点
    int getOrCreateNode(const std::string& name);

    int numNodeEquations() const;
    int numVoltageBranches() const;
    int numUnknowns() const;
    void assignEquationIndices();

    // 添加元件的工厂方法
    void addResistor(const std::string& name, const std::string& n1, const std::string& n2, double value);
    void addCapacitor(const std::string& name, const std::string& n1, const std::string& n2, double value);
    void addInductor(const std::string& name, const std::string& n1, const std::string& n2, double value);
    void addCurrentSource(const std::string& name, const std::string& np, 
                         const std::string& nm, double value);
    void addVoltageSource(const std::string& name, const std::string& np,
                         const std::string& nm, double value);
    void addMosfet(const std::string& name, const std::string& nd, 
                    const std::string& ng, const std::string& ns, const std::string& nb, const std::string& model);

    // 打印电路连接信息
    void printConnectivity() const;
};