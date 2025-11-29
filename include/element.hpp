#pragma once

#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include "sim.hpp"

class Circuit;

class Element {
protected:
    std::string name;
    std::vector<int> nodeIds;  // 存储节点的索引

public:
    Element(const std::string& n, const std::vector<int>& nodes)
        : name(n), nodeIds(nodes) {}

    virtual ~Element() {}

    const std::string& getName() const { return name; }
    const std::vector<int>& getNodeIds() const { return nodeIds; }

    // 在 MNA 方程中插入元件的影响
    virtual void stamp(
        Eigen::MatrixXd& G, Eigen::VectorXd& I, const Circuit& ckt,
        const Eigen::VectorXd& x, const AnalysisContext& ctx
    ) const = 0;

    virtual void stampAC(Eigen::MatrixXcd& /*Y*/, Eigen::VectorXcd& /*J*/,
                         const Circuit& /*ckt*/,
                         double /*omega*/) const
    {
        // 默认啥也不做，只有真正 AC 需要的元件去重载
    }
};

class Resistor : public Element {
private:
    double R;
public:
    Resistor(const std::string& n, int n1, int n2, double r)
        : Element(n, {n1, n2}), R(r) {}

    void stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
               const Circuit& ckt,
               const Eigen::VectorXd& x,
               const AnalysisContext& ctx) const override;
};

// 电流源 I（从 nodeIds[0] -> nodeIds[1]，值为 value）
class CurrentSource : public Element {
    SourceSpec spec;
public:
    CurrentSource(const std::string& n, int np, int nm, const SourceSpec& s)
        : Element(n, {np, nm}), spec(s) {}

    const SourceSpec& setSpec() const { return spec; }

    void stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
               const Circuit& ckt,
               const Eigen::VectorXd& x,
               const AnalysisContext& ctx) const override;
    void stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                const Circuit& ckt, double omega) const override;
};

class VoltageSource : public Element {
<<<<<<< HEAD
public:
    enum class Waveform { DC, SIN };

private:
    Waveform waveform_;
    double   dcValue_;      // 对 DC：直流值；对 SIN：VOFF
    double   sinAmp_;       // VAMP
    double   sinFreq_;      // 频率 Hz
    double   sinPhaseDeg_;  // 相位（度）
    int      branchEqIndex_;

public:
    // 直流源：V = dc
    VoltageSource(const std::string& n, int np, int nm, double dc)
        : Element(n, {np, nm}),
          waveform_(Waveform::DC),
          dcValue_(dc),
          sinAmp_(0.0),
          sinFreq_(0.0),
          sinPhaseDeg_(0.0),
          branchEqIndex_( -1 ) {}

    // 正弦源：V(t) = VOFF + VAMP * sin(2π f t + phase)
    VoltageSource(const std::string& n, int np, int nm,
                  double voff, double vamp,
                  double freq, double phaseDeg)
        : Element(n, {np, nm}),
          waveform_(Waveform::SIN),
          dcValue_(voff),     // DC 工作点使用 VOFF
          sinAmp_(vamp),
          sinFreq_(freq),
          sinPhaseDeg_(phaseDeg),
          branchEqIndex_( -1 ) {}

    void setBranchEqIndex(int idx) { branchEqIndex_ = idx; }
    int  getBranchEqIndex() const { return branchEqIndex_; }

    Waveform getWaveform()   const { return waveform_; }
    double   getDcValue()    const { return dcValue_; }      // DC / OP 用这个
    double   getSinAmp()     const { return sinAmp_; }
    double   getSinFreq()    const { return sinFreq_; }
    double   getSinPhaseDeg() const { return sinPhaseDeg_; }
=======
    SourceSpec spec;
    int branchEqIndex;
public:
    VoltageSource(const std::string& n, int np, int nm, const SourceSpec& s)
        : Element(n, {np, nm}), spec(s), branchEqIndex(-1) {}

    void setBranchEqIndex(int idx) { branchEqIndex = idx; }
    int  getBranchEqIndex() const { return branchEqIndex; }
    const SourceSpec& getSpec() const { return spec; }
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27

    void stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
               const Circuit& ckt,
               const Eigen::VectorXd& x,
               const AnalysisContext& ctx) const override;

    void stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                const Circuit& ckt, double omega) const override;
};


// 电容（DC 中视为开路）
class CapacitorElement : public Element {
    double C;
public:
    CapacitorElement(const std::string& n, int n1, int n2, double c)
        : Element(n, {n1, n2}), C(c) {}

    void stamp(Eigen::MatrixXd& /*G*/, Eigen::VectorXd& /*I*/,
               const Circuit& /*ckt*/,
               const Eigen::VectorXd& /*x*/,
               const AnalysisContext& ctx) const override {
        // DC 中视为开路，不 stamp
    }
};

// 电感（DC 中视为 0V 电压源）
class Inductor : public Element {
    double L;
    int branchEqIndex;
public:
    Inductor(const std::string& n, int n1, int n2, double l)
        : Element(n, {n1, n2}), L(l), branchEqIndex(-1) {}

    void setBranchEqIndex(int idx) { branchEqIndex = idx; }
    int  getBranchEqIndex() const { return branchEqIndex; }

    void stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
               const Circuit& ckt,
               const Eigen::VectorXd& /*x*/,
               const AnalysisContext& ctx) const override;
};

// =============== MOSFET：NMOS / PMOS ===============
//
// 使用非常简化的 Level-1 (Shichman-Hodges) 模型：
//   Ids (D->S) = 0, Vgs <= Vth
//   Triode:    Ids = K * [ (Vgs-Vth)*Vds - 0.5*Vds^2 ]
//              (Vgs > Vth, Vds < Vgs - Vth)
//   Saturation: Ids = 0.5 * K * (Vgs - Vth)^2  (Vds >= Vgs-Vth)
//
// 忽略体效应、λ 等，为了演示 Newton 线性化与非线性支持。
//
// 对 PMOS，使用极性参数 p = -1：
//   Vgs_eff = p*(Vg - Vs), Vds_eff = p*(Vd - Vs)
//   用上式算出 Ids_eff，再令实际 Ids = p * Ids_eff
//   导数对 Vd/Vg/Vs 的形式与 NMOS 相同。
class MosfetBase : public Element {
protected:
    bool isP;       // false: NMOS, true: PMOS
    double Vth;     // 阈值（正数），NMOS: Vgs > Vth 打开；PMOS: Vsg > Vth 打开
    double K;
    double lambda;
    double Cj0;       // K (A/V^2)，越大越“强”

public:
    MosfetBase(const std::string& n,
               int nd, int ng, int ns, int nb,
               bool isPchannel,
               double Vth_ , double K_, double lambda_, double Cj0_   )
        : Element(n, {nd, ng, ns, nb}),
          isP(isPchannel),
          Vth(Vth_),
          K(K_) , lambda(lambda_), Cj0(Cj0_) {}

    void stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
               const Circuit& ckt,
               const Eigen::VectorXd& x,
<<<<<<< HEAD
               double sourceScale) const override;
=======
               const AnalysisContext& ctx) const override;
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
};

class NMosElement : public MosfetBase {
public:
    NMosElement(const std::string& n, int nd, int ng, int ns, int nb,
        double Vth_ , double K_, double lambda_, double Cj0_)
        : MosfetBase(n, nd, ng, ns, nb, false, Vth_, K_, lambda_, Cj0_) {}
};

class PMosElement : public MosfetBase {
public:
    PMosElement(const std::string& n, int nd, int ng, int ns, int nb,
        double Vth_ , double K_, double lambda_, double Cj0_)
        : MosfetBase(n, nd, ng, ns, nb, true, Vth_, K_, lambda_, Cj0_) {}
};

