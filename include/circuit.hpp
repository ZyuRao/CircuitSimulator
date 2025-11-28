#pragma once

#include <vector>
#include <unordered_map>
#include <memory>
#include <string>
#include <unordered_map>
#include "element.hpp"
#include "utils.hpp"

struct Node {
    int id;
    std::string name;
    int eqIndex;                       // 对应 MNA 方程号（-1 为 GND）
    std::vector<int> attachedElements; // elements 中的下标

    Node(int i, const std::string& n)
        : id(i), name(n), eqIndex(-1) {}
};

//后续可以修改
struct MosModel {
    std::string name; //ModelId
    bool isP = false;

    double VT     = 0.7;   
    double MU     = 1e-3;
    double COX    = 1e-3;
    double LAMBDA = 0.0;
    double CJO    = 0.0;
};

class Circuit {
public:
    std::vector<Node> nodes;
    std::vector<std::shared_ptr<Element>> elements;
    std::unordered_map<std::string, int> nodeNameToId;

    std::unordered_map<std::string, MosModel> mosModels;
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
                    const std::string& ng, const std::string& ns, const std::string& modelId,
                    double W, double L);

    void addMosModel(const MosModel& m);

    const MosModel* findMosModel(const std::string& id) const;

    // 打印电路连接信息
    void printConnectivity() const;
};