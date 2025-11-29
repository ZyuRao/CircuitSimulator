#define _USW_MATH_DEFINES
#include "element.hpp"
#include "circuit.hpp"
#include <iostream>
#include <complex>

// =============== 各元件 stamp 实现 ===============

void Resistor::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                     const Circuit& ckt,
                     const Eigen::VectorXd& /*x*/,
                     const AnalysisContext& ) const {
    (void)I; // 未使用

    int n1 = nodeIds[0];
    int n2 = nodeIds[1];
    int eq1 = ckt.nodes[n1].eqIndex;
    int eq2 = ckt.nodes[n2].eqIndex;

    if (R == 0.0) {
        std::cerr << "Warning: resistor " << name << " has zero resistance.\n";
        return;
    }
    double g = 1.0 / R;

    if (eq1 >= 0) G(eq1, eq1) += g;
    if (eq2 >= 0) G(eq2, eq2) += g;
    if (eq1 >= 0 && eq2 >= 0) {
        G(eq1, eq2) -= g;
        G(eq2, eq1) -= g;
    }
}

void CurrentSource::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                          const Circuit& ckt,
                          const Eigen::VectorXd& /*x*/,
                          const AnalysisContext& ctx) const {
    (void)G;
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;

    double Ival = 0.0; // 斜坡源：逐步放大

    switch(ctx.type) {
        case AnalysisType::OP:
        case AnalysisType::DC:
            Ival = spec.evalDC(ctx.sourceScale);
            break;
        case AnalysisType::TRAN:
            Ival = spec.evalTran(ctx.time);
            break;
        case AnalysisType::AC:
            return;
        case AnalysisType::NONE:
            return;
    }

    // 定义为从 p -> m 的电流源，所以“离开 p 的电流”为 +Ival
    // 在 Gv = I 形式下，我们把独立电流源的贡献放入 I 向量
    // 约定：I(row) 存的是“-离开节点的独立电流”，
    // 这样最后解出来的方程是 sum(branch currents leaving) + I(source) = 0
    if (eqP >= 0) I(eqP) -= Ival;  // 离开 p 的电流 +Ival => I -= Ival
    if (eqM >= 0) I(eqM) += Ival;  // 离开 m 的电流 -Ival => I += Ival
}

void CurrentSource::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                            const Circuit& ckt, double) const {
    using cd = std::complex<double>;
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;

    double phaseRad = spec.acPhaseDeg * M_PI / 180.0;
    cd Iac = std::polar(spec.acMag, phaseRad);

    if(eqP >= 0) J(eqP) -= Iac;
    if(eqM >= 0) J(eqM) += Iac;
}

void VoltageSource::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                          const Circuit& ckt,
                          const Eigen::VectorXd& /*x*/,
<<<<<<< HEAD
                          double sourceScale) const {
    int np  = nodeIds[0];
    int nm  = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex_;
=======
                          const AnalysisContext& ctx) const {
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex;

>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27

    if (k < 0 || k >= G.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for " << name << "\n";
        return;
    }

<<<<<<< HEAD
    // DC/OP 分析使用等效直流值（DC 源直接是 dcValue，SIN 源用 VOFF）
    double Vdc = getDcValue();
    double V   = sourceScale * Vdc;

    // 节点方程中的电流未知量 I_v
=======
    double Vval = 0.0;

    switch(ctx.type) {
        case AnalysisType::OP:
        case AnalysisType::DC:
            Vval = spec.evalDC(ctx.sourceScale);
            break;
        case AnalysisType::TRAN:
            Vval = spec.evalTran(ctx.time);
            break;
        case AnalysisType::AC:
            return;
        case AnalysisType::NONE:
            return;
    }
    // 节点方程中的电压源电流 I_v
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;

    // 电压源方程：V(p) - V(m) = V
    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;

    I(k) += Vval;
}

void VoltageSource::stampAC(
    Eigen::MatrixXcd& Y, Eigen::VectorXcd& J, 
    const Circuit& ckt, double 
) const {
    using cd = std::complex<double>;
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex;

    if (k < 0 || k >= Y.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for " << name << "\n";
        return;
    }

    double phaseRad = spec.acPhaseDeg * M_PI / 180.0;
    cd Vac = std::polar(spec.acMag, phaseRad);

    if (eqP >= 0) Y(eqP, k) += cd(1.0, 0.0);
    if (eqM >= 0) Y(eqM, k) -= cd(1.0, 0.0);

    if (eqP >= 0) Y(k, eqP) += cd(1.0, 0.0);
    if (eqM >= 0) Y(k, eqM) -= cd(1.0, 0.0);

    J(k) += Vac;
}



// 电感在 DC 中视为 0V 电压源（短路 + 支路电流未知量）
void Inductor::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                     const Circuit& ckt,
                     const Eigen::VectorXd& /*x*/,
                     const AnalysisContext& ) const {
    (void)I;
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex;

    if (k < 0 || k >= G.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for inductor " << name << "\n";
        return;
    }

    // 和 0V 电压源完全一样，只是 V = 0（不往 I(k) 里加东西）
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;

    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;
}

// MOSFET 牛顿线性化 stamp（修正符号问题）
void MosfetBase::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                       const Circuit& ckt,
                       const Eigen::VectorXd& x,
<<<<<<< HEAD
                       double sourceScale) const {
    (void)sourceScale; // MOS 导纳不参与源斜坡，只依赖当前解 x

=======
                       const AnalysisContext&) const {
>>>>>>> 2edfa30d876afc48a5c4ddd7f6e3757c7a097b27
    // 节点：D G S B
    int nD = nodeIds[0];
    int nG = nodeIds[1];
    int nS = nodeIds[2];
    int nB = nodeIds.size() > 3 ? nodeIds[3] : nodeIds[2];

    int eqD = ckt.nodes[nD].eqIndex;
    int eqG = ckt.nodes[nG].eqIndex;
    int eqS = ckt.nodes[nS].eqIndex;
    int eqB = ckt.nodes[nB].eqIndex;

    auto getV = [&](int eq) -> double {
        if (eq >= 0 && eq < x.size()) return x(eq);
        return 0.0;
    };

    double Vd = getV(eqD);
    double Vg = getV(eqG);
    double Vs = getV(eqS);
    (void)eqB; // 目前忽略体效应
    // double Vb = getV(eqB);

    double p = isP ? -1.0 : 1.0;

    // 有效电压：将 PMOS 映射到等效 NMOS
    double Vgs_eff = p * (Vg - Vs);
    double Vds_eff = p * (Vd - Vs);

    double Ids_eff = 0.0;
    double dId_dVds_eff = 0.0;
    double dId_dVgs_eff = 0.0;

    // 简单 Level-1 模型 + channel-length modulation (λ)
    if (Vgs_eff > Vth) {
        double Vov = Vgs_eff - Vth; // overdrive

        // 先计算不含 λ 的 Ids0 及其导数
        double Id0        = 0.0;
        double dId0_dVds  = 0.0;
        double dId0_dVgs  = 0.0;

        if (Vds_eff < Vov) {
            // Triode 区
            Id0       = K * (Vov * Vds_eff - 0.5 * Vds_eff * Vds_eff);
            dId0_dVds = K * (Vov - Vds_eff);
            dId0_dVgs = K * Vds_eff;
        } else {
            // Saturation 区
            Id0       = 0.5 * K * Vov * Vov;
            dId0_dVds = 0.0;
            dId0_dVgs = K * Vov;
        }

        // 引入 λ： Ids = Id0 * (1 + λ Vds)
        double onePlusLambdaVds = 1.0 + lambda * Vds_eff;

        Ids_eff      = Id0 * onePlusLambdaVds;
        dId_dVds_eff = dId0_dVds * onePlusLambdaVds + lambda * Id0;
        dId_dVgs_eff = dId0_dVgs * onePlusLambdaVds;
    } else {
        // 截止区：电流近似为 0
        Ids_eff      = 0.0;
        dId_dVds_eff = 0.0;
        dId_dVgs_eff = 0.0;
    }


    // 映射回实际器件：Ids 为从 D -> S 的电流
    double Ids = p * Ids_eff;

    // 对物理节点电压的偏导（链式法则）
    double gd = dId_dVds_eff;                         // ∂Ids/∂Vd
    double gg = dId_dVgs_eff;                         // ∂Ids/∂Vg
    double gs = -(dId_dVds_eff + dId_dVgs_eff);       // ∂Ids/∂Vs

    // 把 Ids 线性化为：Ids ≈ gd*Vd + gg*Vg + gs*Vs + cst
    double cst = Ids - gd * Vd - gg * Vg - gs * Vs;

    // 关键：这里的“支路电流”是从 D->S，且我们在 KCL 中使用“离开节点的电流”为正。
    //
    // 在 D 节点：离开节点的电流 = +Ids
    // 在 S 节点：离开节点的电流 = -Ids
    //
    // 所以线性化后：
    //   D 节点的离开电流 ≈ gd*Vd + gg*Vg + gs*Vs + cst
    //   S 节点的离开电流 ≈ -gd*Vd - gg*Vg - gs*Vs - cst
    //
    // 而我们的 MNA 形式是  G * v = I，等价于  (所有支路电流之和) + I(source) = 0，
    // 其中 I(source) 是独立源的“负”号。对非线性支路的常数项 cst，我们
    // 等价看成一个独立电流源，于是要对 I(row) 加上“负的 cst”。

    // D 节点：G(D,D)+=gd, G(D,G)+=gg, G(D,S)+=gs, I(D) -= cst
    if (eqD >= 0) {
        if (eqD >= 0) G(eqD, eqD) += gd;
        if (eqG >= 0) G(eqD, eqG) += gg;
        if (eqS >= 0) G(eqD, eqS) += gs;
        I(eqD) -= cst;
    }

    // S 节点：离开电流是 -Ids = (-gd)*Vd + (-gg)*Vg + (-gs)*Vs - cst
    //        => G(S,D)+=-gd, G(S,G)+=-gg, G(S,S)+=-gs, I(S) += cst
    if (eqS >= 0) {
        if (eqD >= 0) G(eqS, eqD) += -gd;
        if (eqG >= 0) G(eqS, eqG) += -gg;
        if (eqS >= 0) G(eqS, eqS) += -gs;
        I(eqS) += cst;
    }

    // Gate / Bulk 节点理想 DC 不导通，这里不对其 KCL stamp（即 Ig=Ib=0）
}