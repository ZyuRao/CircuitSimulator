#define _USE_MATH_DEFINES
#include "element.hpp"
#include "circuit.hpp"
#include <iostream>
#include <complex>

// =============== 各元件 stamp 实现 ===============


static std::complex<double> evalSourcePhasor(const SourceSpec& spec,
                                             const AnalysisContext& ctx)
{
    using cd = std::complex<double>;
    cd j(0.0, 1.0);

    // AC 小信号：沿用你们原逻辑
    if (ctx.type == AnalysisType::AC) {
        double phaseRad = spec.acPhaseDeg * M_PI / 180.0;
        return std::polar(spec.acMag, phaseRad);
    }

    // HB：把“单音 SIN”映射到对应谐波系数
    if (ctx.type == AnalysisType::HB) {
        // DC 分量（h=0）：dcValue + sine.v0，并乘 sourceScale
        if (ctx.hbHarm == 0) {
            return cd(spec.evalDC(ctx.sourceScale), 0.0);
        }

        // 只支持 SIN（PULSE/PWL 先不支持，HB 里返回 0）
        if (spec.tran.type != WaveformType::SIN) return cd(0.0, 0.0);

        const auto& s = spec.tran.sine;
        if (s.freq <= 0.0 || ctx.hbF0 <= 0.0) return cd(0.0, 0.0);

        // 单音 HB：只在 s.freq ≈ h*f0 时给该谐波注入
        double target = ctx.hbHarm * ctx.hbF0;
        double tol = 1e-6 * std::max(1.0, s.freq);
        if (std::abs(s.freq - target) > tol) return cd(0.0, 0.0);

        // v(t)=v0+va*sin(ω(t-td)+phi) => 正频率系数 Vk = (-j)*(va/2)*exp(j*(phi-ωtd))
        double omega = 2.0 * M_PI * s.freq;
        double phase = s.phi - omega * s.td;

        return (-j) * (0.5 * s.va * ctx.sourceScale) * std::exp(j * phase);
    }

    return cd(0.0, 0.0);
}


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

void CurrentSource::stampAC(Eigen::MatrixXcd&, Eigen::VectorXcd& J,
                            const Circuit& ckt,
                            const AnalysisContext& ctx) const
{
    using cd = std::complex<double>;
    int np  = nodeIds[0];
    int nm  = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;

    cd I = evalSourcePhasor(spec, ctx);

    if (eqP >= 0) J(eqP) -= I;
    if (eqM >= 0) J(eqM) += I;
}


void VoltageSource::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                          const Circuit& ckt,
                          const Eigen::VectorXd& /*x*/,
                          const AnalysisContext& ctx) const {
    int np = nodeIds[0];
    int nm = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex;


    if (k < 0 || k >= G.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for " << name << "\n";
        return;
    }

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
    if (eqP >= 0) G(eqP, k) += 1.0;
    if (eqM >= 0) G(eqM, k) -= 1.0;

    // 电压源方程：V(p) - V(m) = V
    if (eqP >= 0) G(k, eqP) += 1.0;
    if (eqM >= 0) G(k, eqM) -= 1.0;

    I(k) += Vval;
}

void VoltageSource::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                            const Circuit& ckt,
                            const AnalysisContext& ctx) const
{
    using cd = std::complex<double>;
    int np  = nodeIds[0];
    int nm  = nodeIds[1];
    int eqP = ckt.nodes[np].eqIndex;
    int eqM = ckt.nodes[nm].eqIndex;
    int k   = branchEqIndex;

    if (k < 0 || k >= Y.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for " << name << "\n";
        return;
    }

    // MNA stamping（和 AC/HB 都一样）
    if (eqP >= 0) Y(eqP, k) += cd(1.0, 0.0);
    if (eqM >= 0) Y(eqM, k) -= cd(1.0, 0.0);
    if (eqP >= 0) Y(k, eqP) += cd(1.0, 0.0);
    if (eqM >= 0) Y(k, eqM) -= cd(1.0, 0.0);

    cd V = evalSourcePhasor(spec, ctx);
    J(k) += V;
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

// MOSFET 牛顿线性化 stamp（修正符号问题 + 动态源漏互换）
void MosfetBase::stamp(Eigen::MatrixXd& G, Eigen::VectorXd& I,
                       const Circuit& ckt,
                       const Eigen::VectorXd& x,
                       const AnalysisContext&) const {
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

    // ======== 根据端电压动态决定源/漏 ========
    {
        if (!isP) {
            // NMOS：source 取电势更低的一端
            if (Vd < Vs) {
                // 互换 D/S 电压
                double tmpV = Vd;
                Vd = Vs;
                Vs = tmpV;
                // 互换 D/S 对应的方程索引
                int tmpEq = eqD;
                eqD = eqS;
                eqS = tmpEq;
            }
        } else {
            // PMOS：source 取电势更高的一端
            if (Vd > Vs) {
                double tmpV = Vd;
                Vd = Vs;
                Vs = tmpV;
                int tmpEq = eqD;
                eqD = eqS;
                eqS = tmpEq;
            }
        }
    }

    double p = isP ? -1.0 : 1.0;

    // 有效电压：将 PMOS 映射到等效 NMOS
    double Vgs_eff = p * (Vg - Vs);
    double Vds_eff = p * (Vd - Vs);

    double Ids_eff = 0.0;
    double dId_dVds_eff = 0.0;
    double dId_dVgs_eff = 0.0;

    // 先算“不带 λ”的 Ids0 / gds0 / gm0
    double Ids0  = 0.0;
    double gds0  = 0.0;
    double gm0   = 0.0;

    bool on = false;
    if (Vgs_eff > Vth && Vds_eff >= 0) {
        on = true;
        double Vov = Vgs_eff - Vth; // overdrive

        // 先计算不含 λ 的 Ids0 及其导数
        double Id0        = 0.0;
        double dId0_dVds  = 0.0;
        double dId0_dVgs  = 0.0;

        if (Vds_eff < Vov) {
            // Triode 区
            Ids0 = K * (Vov * Vds_eff - 0.5 * Vds_eff * Vds_eff);
            gds0 = K * (Vov - Vds_eff);  // ∂Ids0/∂Vds
            gm0  = K * Vds_eff;          // ∂Ids0/∂Vgs
        } else {
            // Saturation 区
            Ids0 = 0.5 * K * Vov * Vov;
            gds0 = 0.0;
            gm0  = K * Vov;
        }
    }

    {
        const double gmin = 1e-12;
        if(!on) {
            Ids0 = 0.0;
            gm0  = 0.0;
            gds0 = gmin;
        }
    }

    // 统一加上沟道长度调制：Ids = Ids0 * (1 + λVds)
    double factor = 1.0 + lambda * Vds_eff;
    if(factor < 0.0) factor = 0.0;
    Ids_eff       = Ids0 * factor;

    // ∂Ids/∂Vds = gds0 * (1 + λVds) + Ids0 * λ
    dId_dVds_eff  = gds0 * factor + Ids0 * lambda;

    // ∂Ids/∂Vgs = gm0 * (1 + λVds)
    dId_dVgs_eff  = gm0 * factor;

    // 映射回实际器件：Ids 为从 D -> S 的电流
    double Ids = p * Ids_eff;

    // 对物理节点电压的偏导（链式法则）
    double gd = dId_dVds_eff;                         // ∂Ids/∂Vd
    double gg = dId_dVgs_eff;                         // ∂Ids/∂Vg
    double gs = -(dId_dVds_eff + dId_dVgs_eff);       // ∂Ids/∂Vs

    // 把 Ids 线性化为：Ids ≈ gd*Vd + gg*Vg + gs*Vs + cst
    double cst = Ids - gd * Vd - gg * Vg - gs * Vs;

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

void Resistor::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                       const Circuit& ckt, const AnalysisContext& ctx) const {
    (void)J;
    int n1 = nodeIds[0];
    int n2 = nodeIds[1];
    int eq1 = ckt.nodes[n1].eqIndex;
    int eq2 = ckt.nodes[n2].eqIndex;

    if (R <= 0.0) {
        std::cerr << "Warning: resistor " << name << " has non-positive R in AC.\n";
        return;
    }
    std::complex<double> g(1.0 / R, 0.0);

    if (eq1 >= 0) Y(eq1, eq1) += g;
    if (eq2 >= 0) Y(eq2, eq2) += g;
    if (eq1 >= 0 && eq2 >= 0) {
        Y(eq1, eq2) -= g;
        Y(eq2, eq1) -= g;
    }
} 

void CapacitorElement::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                               const Circuit& ckt,
                               const AnalysisContext& ctx) const
{
    (void)J;
    int eq1 = ckt.nodes[nodeIds[0]].eqIndex;
    int eq2 = ckt.nodes[nodeIds[1]].eqIndex;

    if (C <= 0.0) return;
    if (ctx.omega == 0.0) return;

    std::complex<double> y(0.0, ctx.omega * C); // jωC
    if (eq1 >= 0) Y(eq1, eq1) += y;
    if (eq2 >= 0) Y(eq2, eq2) += y;
    if (eq1 >= 0 && eq2 >= 0) {
        Y(eq1, eq2) -= y;
        Y(eq2, eq1) -= y;
    }
}


void Inductor::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                       const Circuit& ckt,
                       const AnalysisContext& ctx) const
{
    (void)J;
    int eqP = ckt.nodes[nodeIds[0]].eqIndex;
    int eqM = ckt.nodes[nodeIds[1]].eqIndex;
    int k   = branchEqIndex;

    if (k < 0 || k >= Y.rows()) {
        std::cerr << "Internal error: invalid branchEqIndex for inductor "
                  << name << " in AC/HB.\n";
        return;
    }

    std::complex<double> z(0.0, ctx.omega * L); // jωL

    if (eqP >= 0) Y(eqP, k) += std::complex<double>(1.0, 0.0);
    if (eqM >= 0) Y(eqM, k) -= std::complex<double>(1.0, 0.0);

    if (eqP >= 0) Y(k, eqP) += std::complex<double>(1.0, 0.0);
    if (eqM >= 0) Y(k, eqM) -= std::complex<double>(1.0, 0.0);

    // V(p)-V(m) - jωL*I = 0
    Y(k, k) += -z;
}

void MosfetBase::stampAC(Eigen::MatrixXcd& Y, Eigen::VectorXcd& J,
                         const Circuit& ckt,
                         const AnalysisContext& ctx) const
{
    (void)J;
    if (ctx.omega == 0.0) return; // DC 下电容开路

    using cd = std::complex<double>;

    auto stampCap = [&](int eq1, int eq2, double C) {
        if (C <= 0.0) return;
        cd y(0.0, ctx.omega * C);
        if (eq1 >= 0) Y(eq1, eq1) += y;
        if (eq2 >= 0) Y(eq2, eq2) += y;
        if (eq1 >= 0 && eq2 >= 0) {
            Y(eq1, eq2) -= y;
            Y(eq2, eq1) -= y;
        }
    };

    int nD = nodeIds[0], nG = nodeIds[1], nS = nodeIds[2];
    int nB = (nodeIds.size() > 3) ? nodeIds[3] : nodeIds[2];

    int eqD = ckt.nodes[nD].eqIndex;
    int eqG = ckt.nodes[nG].eqIndex;
    int eqS = ckt.nodes[nS].eqIndex;
    int eqB = ckt.nodes[nB].eqIndex;

    // 你 HB 里现在用的模型：Cj0 常数拆分
    double Cgs = 1.0 * Cg0, Cgd = 1.0 * Cg0;
    double Csb = 1.0 * Cj0, Cdb = 1.0 * Cj0;

    stampCap(eqG, eqS, Cgs);
    stampCap(eqG, eqD, Cgd);
    stampCap(eqS, -1, Csb);
    stampCap(eqD, -1, Cdb);
}

void MosfetBase::evalIdsGmGds(
    const Circuit& ckt, 
    const Eigen::VectorXd& nodeVoltages,
    double& Ids, double& gm, double& gds
) const {
    // 网表给的 D/G/S/B 结点编号
    int nD = nodeIds[0];
    int nG = nodeIds[1];
    int nS = nodeIds[2];
    int nB = nodeIds.size() > 3 ? nodeIds[3] : nodeIds[2];

    int eqD = ckt.nodes[nD].eqIndex;
    int eqG = ckt.nodes[nG].eqIndex;
    int eqS = ckt.nodes[nS].eqIndex;
    int eqB = ckt.nodes[nB].eqIndex;
    (void)eqB;

    auto getV = [&](int eq) -> double {
        if (eq < 0 || eq >= nodeVoltages.size()) return 0.0;
        return nodeVoltages(eq);
    };

    // 原始结点电压（保持一份）
    double Vd0 = getV(eqD);
    double Vg0 = getV(eqG);
    double Vs0 = getV(eqS);
    double Vb0 = getV(eqB);
    (void)Vb0; // 目前没用到体效应

    // p = +1: NMOS, p = -1: PMOS（统一成“等效 NMOS”思路）
    double p = isP ? -1.0 : 1.0;

    // ========= 1. 先根据电压决定“有效的 D/S”方向 =========
    double Vds_eff_A = p * (Vd0 - Vs0);
    bool swapped = (Vds_eff_A < 0.0);

    double Vd_eff = Vd0;
    double Vs_eff = Vs0;
    if (swapped) {
        double tmp = Vd_eff;
        Vd_eff = Vs_eff;
        Vs_eff = tmp;
    }

    // ========= 2. 在“有效 D/S”坐标系下计算等效 NMOS 电流和导数 =========
    double Vgs_eff = p * (Vg0 - Vs_eff);
    double Vds_eff = p * (Vd_eff - Vs_eff);  // 经过上面的选择，正常情况下 >= 0

    // 先算不带 λ 的 Ids0 / gm0 / gds0
    double Ids0 = 0.0;
    double gm0  = 0.0;
    double gds0 = 0.0;
    bool on     = false;

    if (Vgs_eff > Vth && Vds_eff >= 0.0) {
        on = true;
        double Vov = Vgs_eff - Vth;  // overdrive
        if (Vov <= 0.0) {
            on = false;
        } else {
            if (Vds_eff < Vov) {
                // Triode 区
                Ids0 = K * (Vov * Vds_eff - 0.5 * Vds_eff * Vds_eff);
                gds0 = K * (Vov - Vds_eff);  // ∂Ids0/∂Vds_eff
                gm0  = K * Vds_eff;          // ∂Ids0/∂Vgs_eff
            } else {
                // Saturation 区
                Ids0 = 0.5 * K * Vov * Vov;
                gds0 = 0.0;
                gm0  = K * Vov;
            }
        }
    }

    // 关断时给一个很小的 gmin，保证 HB 线性化矩阵不会完全断路
    const double gmin = 1e-12;
    if (!on) {
        Ids0 = 0.0;
        gm0  = 0.0;
        gds0 = gmin;
    }

    // 加上沟道长度调制：Ids_eff = Ids0 * (1 + λ Vds_eff)
    double factor = 1.0 + lambda * Vds_eff;
    if (factor < 0.0) factor = 0.0;

    double Ids_eff       = Ids0 * factor;
    double dId_dVds_eff  = gds0 * factor + Ids0 * lambda;  // ∂Ids_eff/∂Vds_eff
    double dId_dVgs_eff  = gm0  * factor;                  // ∂Ids_eff/∂Vgs_eff

    // ========= 3. 从“等效 NMOS + 有效 D/S”映射回“原网表 D/S” =========
    if (!swapped) {
        Ids = p * Ids_eff;
        gm  = dId_dVgs_eff;
        gds = dId_dVds_eff;
    } else {
        Ids = -p * Ids_eff;
        gm  = -dId_dVgs_eff;
        gds =  dId_dVgs_eff + dId_dVds_eff;
    }
}

void MosfetBase::stampNonlinearConductance(
    Eigen::MatrixXd& Gnl,
    const Circuit& ckt,
    const Eigen::VectorXd& nodeVoltages
) const {
    double Ids, gm, gds;
    evalIdsGmGds(ckt, nodeVoltages, Ids, gm, gds);

    int nD = nodeIds[0];
    int nG = nodeIds[1];
    int nS = nodeIds[2];
    int nB = nodeIds.size() > 3 ? nodeIds[3] : nodeIds[2];

    int eqD = ckt.nodes[nD].eqIndex;
    int eqG = ckt.nodes[nG].eqIndex;
    int eqS = ckt.nodes[nS].eqIndex;
    int eqB = ckt.nodes[nB].eqIndex;

    auto add = [&](int r, int c, double val) {
        if (r >= 0 && r < Gnl.rows() && c >= 0 && c < Gnl.cols()) {
            Gnl(r, c) += val;
        }
    };

     // 对 vD 的导数
    add(eqD, eqD,  gds);          // ∂I_D/∂vD
    add(eqS, eqD, -gds);          // ∂I_S/∂vD

    // 对 vG 的导数
    add(eqD, eqG,  gm);           // ∂I_D/∂vG
    add(eqS, eqG, -gm);           // ∂I_S/∂vG

    // 对 vS 的导数
    double dId_dvS = -gm - gds;
    add(eqD, eqS,  dId_dvS);      // ∂I_D/∂vS
    add(eqS, eqS, -dId_dvS);      // ∂I_S/∂vS

    (void)eqB;
}