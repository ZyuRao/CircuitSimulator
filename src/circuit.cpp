#include "circuit.hpp"
#include <iostream>
#include <memory>

int Circuit::getOrCreateNode(const std::string& name) {
    auto it = nodeNameToId.find(name);
    if (it != nodeNameToId.end()) {
        return it->second;
    }
    int id = static_cast<int>(nodes.size());
    nodes.emplace_back(id, name);
    nodeNameToId[name] = id;
    return id;
}

int Circuit::numNodeEquations() const {
    int count = 0;
    for (const auto& node : nodes) {
        if (!isGroundName(node.name)) {
            ++count;
        }
    }
    return count;
}

// 只有电压源和电感引入支路电流未知量
int Circuit::numVoltageBranches() const {
    int count = 0;
    for (const auto& e : elements) {
        if (std::dynamic_pointer_cast<VoltageSource>(e) ||
            std::dynamic_pointer_cast<Inductor>(e)) {
            ++count;
        }
    }
    return count;
}

int Circuit::numUnknowns() const {
    return numNodeEquations() + numVoltageBranches();
}

void Circuit::assignEquationIndices() {
    int eq = 0;
    // 节点电压未知量
    for (auto& node : nodes) {
        if (isGroundName(node.name)) {
            node.eqIndex = -1;
        } else {
            node.eqIndex = eq++;
        }
    }

    // 电压源 / 电感 电流未知量
    for (auto& e : elements) {
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            vs->setBranchEqIndex(eq++);
        } else if (auto ind = std::dynamic_pointer_cast<Inductor>(e)) {
            ind->setBranchEqIndex(eq++);
        }
    }
}

void Circuit::addResistor(const std::string& name,
                          const std::string& n1,
                          const std::string& n2,
                          double value) {
    int id1 = getOrCreateNode(n1);
    int id2 = getOrCreateNode(n2);
    auto e = std::make_shared<Resistor>(name, id1, id2, value);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[id1].attachedElements.push_back(idx);
    nodes[id2].attachedElements.push_back(idx);
}

void Circuit::addCapacitor(const std::string& name,
                           const std::string& n1,
                           const std::string& n2,
                           double value) {
    int id1 = getOrCreateNode(n1);
    int id2 = getOrCreateNode(n2);
    auto e = std::make_shared<CapacitorElement>(name, id1, id2, value);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[id1].attachedElements.push_back(idx);
    nodes[id2].attachedElements.push_back(idx);
}

void Circuit::addInductor(const std::string& name,
                          const std::string& n1,
                          const std::string& n2,
                          double value) {
    int id1 = getOrCreateNode(n1);
    int id2 = getOrCreateNode(n2);
    auto e = std::make_shared<Inductor>(name, id1, id2, value);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[id1].attachedElements.push_back(idx);
    nodes[id2].attachedElements.push_back(idx);
}

void Circuit::addCurrentSource(const std::string& name,
                               const std::string& np,
                               const std::string& nm,
                               const SourceSpec& spec) {
    int idp = getOrCreateNode(np);
    int idm = getOrCreateNode(nm);
    auto e = std::make_shared<CurrentSource>(name, idp, idm, spec);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[idp].attachedElements.push_back(idx);
    nodes[idm].attachedElements.push_back(idx);
}

void Circuit::addVoltageSource(const std::string& name,
                               const std::string& np,
                               const std::string& nm,
                               const SourceSpec& spec) {
    int idp = getOrCreateNode(np);
    int idm = getOrCreateNode(nm);
    auto e = std::make_shared<VoltageSource>(name, idp, idm, spec);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[idp].attachedElements.push_back(idx);
    nodes[idm].attachedElements.push_back(idx);
}

void Circuit::addVoltageSourceSin(const std::string& name,
                                  const std::string& np,
                                  const std::string& nm,
                                  double voff, double vamp,
                                  double freq, double phaseDeg) {
    int idp = getOrCreateNode(np);
    int idm = getOrCreateNode(nm);
    auto e  = std::make_shared<VoltageSource>(name, idp, idm,
                                              voff, vamp, freq, phaseDeg);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[idp].attachedElements.push_back(idx);
    nodes[idm].attachedElements.push_back(idx);
}

void Circuit::addMosfet(
    const std::string& name, const std::string& nd, 
    const std::string& ng, const std::string& ns,
    const std::string& modelId, double W, double L
) {
    const MosModel* m = findMosModel(modelId);
    if(!m) {
        std::cerr << "Unknown MOS model: " << modelId << "\n";
        return;
    }

    int idd = getOrCreateNode(nd);
    int idg = getOrCreateNode(ng);
    int ids = getOrCreateNode(ns);
    int idb = getOrCreateNode("0");

    double K = m->MU * m->COX * (W/L);
    double Vth_mag = std::abs(m->VT);
    double lambda = m->LAMBDA;
    double CJo = m->CJO;

    std::shared_ptr<Element> e;
    if (m->isP) {
        e = std::make_shared<PMosElement>(
            name, idd, idg, ids, idb, Vth_mag, K, 
            lambda, CJo
        );
    } else {
        e = std::make_shared<NMosElement>(
            name, idd, idg, ids, idb, Vth_mag, K, 
            lambda, CJo
        );
    }

    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[idd].attachedElements.push_back(idx);
    nodes[idg].attachedElements.push_back(idx);
    nodes[ids].attachedElements.push_back(idx);
    nodes[idb].attachedElements.push_back(idx);
}

<<<<<<< HEAD
void Circuit::addMosModel(const MosModel& m)
{
    std::string key = toLower(m.name);
    mosModels[key] = m;
}

const MosModel* Circuit::findMosModel(const std::string& id) const
{
    std::string key = toLower(id);
    auto it = mosModels.find(key);
    if (it == mosModels.end()) {
        return nullptr;
    }
    return &it->second;
=======
void Circuit::addMosModel(const MosModel& m) {
    mosModels[m.name] = m;
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
}

void Circuit::printConnectivity() const {
    std::cout << "========== 节点与连接关系 ==========\n";
    for (const auto& node : nodes) {
        std::cout << "Node " << node.name
                  << " (id=" << node.id
                  << ", eqIndex=" << node.eqIndex
                  << "): ";
        for (int ei : node.attachedElements) {
            std::cout << elements[ei]->getName() << " ";
        }
        std::cout << "\n";
    }
}

const MosModel* Circuit::findMosModel(const std::string& id) const {
    auto it = mosModels.find(id);
    if (it == mosModels.end()) return nullptr;
    return &it->second;
}