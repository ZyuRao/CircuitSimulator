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

void Circuit::invalidateCache() {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    cache_.reset();
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
    invalidateCache();
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
    invalidateCache();
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
    invalidateCache();
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
    invalidateCache();
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
    invalidateCache();
    int idp = getOrCreateNode(np);
    int idm = getOrCreateNode(nm);
    auto e = std::make_shared<VoltageSource>(name, idp, idm, spec);
    int idx = static_cast<int>(elements.size());
    elements.push_back(e);
    nodes[idp].attachedElements.push_back(idx);
    nodes[idm].attachedElements.push_back(idx);
}

void Circuit::addMosfet(const std::string& name,
                        const std::string& drain,
                        const std::string& gate,
                        const std::string& source,
                        const std::string& modelId,
                        double W,
                        double L)
{
    invalidateCache();
    const MosModel* m = findMosModel(modelId);
    if (!m) {
        throw std::runtime_error("Unknown MOS model: " + modelId);
    }

    int idd = getOrCreateNode(drain);
    int idg = getOrCreateNode(gate);
    int ids = getOrCreateNode(source);

    // 你的网表是 3 节点 MOS：默认 bulk = source
    int idb = ids;

    double Vth_mag = std::abs(m->VT);

    double K = 0.0;
    if (L > 0.0) {
        K = m->MU * m->COX * (W / L);
    }

    double lambda = m->LAMBDA;

    // 按你给的公式：结电容直接用 CJ0，不乘 W*L
    double Cj0 = m->CJO;

    // 按你给的公式：Cg0 = Cox * W * L
    double Cg0 = 0.0;
    if (m->COX > 0.0 && W > 0.0 && L > 0.0) {
        Cg0 = 0.5 * m->COX * W * L;
    }

    std::shared_ptr<Element> e;
    if (m->isP) {
        e = std::make_shared<PMosElement>(name, idd, idg, ids, idb, Vth_mag, K, lambda, Cg0, Cj0);
    } else {
        e = std::make_shared<NMosElement>(name, idd, idg, ids, idb, Vth_mag, K, lambda, Cg0, Cj0);
    }

    int idx = (int)elements.size();
    elements.push_back(e);

    // 附着关系：bulk=source 时不额外重复挂载
    nodes[idd].attachedElements.push_back(idx);
    nodes[idg].attachedElements.push_back(idx);
    nodes[ids].attachedElements.push_back(idx);
    if (idb != ids) nodes[idb].attachedElements.push_back(idx);
}


void Circuit::addMosModel(const MosModel& m) {
    mosModels[m.name] = m;
}

const Circuit::ElementCache& Circuit::elementCache() const {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    if (cache_) return *cache_;

    cache_.emplace();
    for (const auto& e : elements) {
        if (auto m = std::dynamic_pointer_cast<MosfetBase>(e)) {
            cache_->mos.push_back(m.get());
            continue;
        }
        if (auto r = std::dynamic_pointer_cast<Resistor>(e)) {
            cache_->resistors.push_back(r.get());
            continue;
        }
        if (auto c = std::dynamic_pointer_cast<CapacitorElement>(e)) {
            cache_->capacitors.push_back(c.get());
            continue;
        }
        if (auto L = std::dynamic_pointer_cast<Inductor>(e)) {
            cache_->inductors.push_back(L.get());
            continue;
        }
        if (auto vs = std::dynamic_pointer_cast<VoltageSource>(e)) {
            cache_->voltageSources.push_back(vs.get());
            continue;
        }
        if (auto is = std::dynamic_pointer_cast<CurrentSource>(e)) {
            cache_->currentSources.push_back(is.get());
            continue;
        }
    }
    return *cache_;
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
