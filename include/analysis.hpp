#pragma once

#include<Eigen/Dense>
#include "circuit.hpp"
#include "sim.hpp"
#include "solver.hpp"
#include <string>
#include <algorithm>
#include <limits>


enum class DcSolverKind {
    LU,
    GaussSeidel
};

enum class NonlinearSolveRole {
    DC_Newton_LU,
    DC_Newton_GS,
    HB_Newton_LU
};

struct ConvParams {
    // damping / step control
    double alphaInit = 1.0;
    double alphaMin  = 0.1;
    double alphaMax  = 1.0;

    // gmin schedule
    double gminHighBase = 1e-4;   // scale=0 附近
    double gminLowBase  = 1e-12;  // scale=1 附近
    double gminAbsMax   = 1e-2;   // 允许在失败时顶上去（但最终会拉回 base）

    // heuristics
    double fastConvRatio = 0.7;
    double slowConvRatio = 1.05;

    // DC: convergence check
    double stepTolStart = 1e-6;   // scale=0 的 step 容忍
    double stepTolFinal = 1e-9;   // scale=1 的 step 容忍
    double fAbsTol      = 1e-12;  // residual: abs
    double fRelTol      = 1e-8;   // residual: rel * max(1, ||rhs||)

    // GS inner linear-solve acceptance (你原来写在 solveNewtonGS 里的 acceptTol)
    double linRelTolStart = 5e-3; // scale=0
    double linRelTolFinal = 1e-6; // scale=1

    // continuation / line-search (你原来写在 solveNewtonGS 里的参数)
    int    maxLineSearch = 12;
    double armijoC1      = 1e-4;
    double lsAlphaMin    = 0.02;

    double dScaleInit = 0.1;
    double dScaleMin  = 1e-3;
    double dScaleMax  = 0.2;
};
struct ConvStatus {
    Eigen::VectorXd xNext;
    double alphaNext = 1.0;
    double gminNext = 0.0;
    double error = std::numeric_limits<double>::infinity();
    bool converged = false;
};

class ConvController {
public:
    static ConvController forDc(DcSolverKind kind);
    static ConvController forHb();

    ConvController() = default; // 保留以免你其它地方一堆报错
    ConvController(NonlinearSolveRole role, const ConvParams& p)
        : role_(role), p_(p) {}

    const ConvParams& params() const { return p_; }

    // 关键：baseGmin 你原公式写错了，这里修正为线性插值
    double baseGmin(double rampScale) const {
        double s = std::clamp(rampScale, 0.0, 1.0);
        double g = p_.gminHighBase * (1.0 - s) + p_.gminLowBase * s;
        return std::max(g, 1e-9);   // DC 建议别低于 1e-9，除非你做了“浮动节点检测”
    }

    double stepTol(double rampScale) const {
        double s = std::clamp(rampScale, 0.0, 1.0);
        return p_.stepTolStart * (1.0 - s) + p_.stepTolFinal * s;
    }

    double linRelTol(double rampScale) const {
        double s = std::clamp(rampScale, 0.0, 1.0);
        return p_.linRelTolStart * (1.0 - s) + p_.linRelTolFinal * s;
    }

    double residualTol(double rhsNorm) const {
        return p_.fAbsTol + p_.fRelTol * std::max(1.0, rhsNorm);
    }

    // 线性内解拒绝时的统一策略（你原来散落在 solveNewtonGS 里）
    void onLinearReject(double& alpha, double& gmin) const {
        gmin  = std::min(gmin * 10.0, p_.gminAbsMax);
        alpha = std::max(alpha * 0.5, p_.alphaMin);
    }

    double initialAlphaLU() const { return std::clamp(p_.alphaInit, p_.alphaMin, p_.alphaMax); }
    double initialAlphaGS() const { return std::clamp(p_.alphaInit, p_.alphaMin, p_.alphaMax); }

    ConvStatus update(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& xRaw,
        double prevErr,
        int iter,
        double alphaCurrent,
        double gminCurrent,
        double rampScale,
        double tol
    ) const {
        // 这里 prevErr 在你旧代码里其实是“上一次步长/误差”
        auto st = updateStep(x, xRaw, prevErr, iter, alphaCurrent, gminCurrent, rampScale);
        st.converged = std::isfinite(st.error) && (st.error < tol);
        return st;
    }

    // LU/GS 都能用的“基于 step 的 update”
    ConvStatus updateStep(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& xRaw,
        double prevStep,
        int iter,
        double alphaCurrent,
        double gminCurrent,
        double rampScale
    ) const;

private:
    NonlinearSolveRole role_ = NonlinearSolveRole::DC_Newton_LU;
    ConvParams p_;
};

class DcAnalysis {
private:
    const Circuit& ckt;
    const SimulationConfig& sim;
    DcSolverKind solverKind;

    // 小工具：构造 DC 的 AnalysisContext
    static AnalysisContext makeDcCtx(double sourceScale);  

     // 只处理纯线性电路
    Eigen::VectorXd solveDirectLU() const;
    Eigen::VectorXd solveDirectGS() const;

    Eigen::VectorXd solveNewtonLU() const;
    Eigen::VectorXd solveNewtonGS() const;

public:
    explicit DcAnalysis(const Circuit& ckt_,
                        const SimulationConfig& SimulationConfig,
                        DcSolverKind solver_ = DcSolverKind::GaussSeidel);

    // 对外唯一入口：自动判断是否有非线性器件，调用合适的路径
    Eigen::VectorXd run();

    // 如果强制指定：
    Eigen::VectorXd runLU();          
    Eigen::VectorXd runGaussSeidel(); 
};

class TransientAnalysis {
public:
    TransientAnalysis(const Circuit& ckt,
                      const SimulationConfig& sim,
                      const std::string& outFile = "tran_out.csv");

    void runBackwardEuler();

    void runTrapezoidal();

    Eigen::VectorXd integrateOnePeriodBE(
        const Eigen::VectorXd& x0,
        double t0,
        double dt,
        double T,
        const std::function<void(double, const Eigen::VectorXd&)>& dumpRow = nullptr
    );
private:
    const Circuit& ckt;
    const SimulationConfig& sim;
    std::string outFile;
        
    friend class DcAnalysis;
     // 从 DC 求工作点：内部直接用 DcAnalysis
    Eigen::VectorXd computeDcOperatingPoint() const;

    // MOS 寄生电容状态
    struct MosCapState {
        double vgsPrev = 0.0;
        double vgdPrev = 0.0;
        double vsbPrev = 0.0;
        double vdbPrev = 0.0;
    };

    // 后向欧拉的“电容 + 历史电流源”等效 stamp
    static void stampCapBE(int eq1, int eq2,
                           double C, double dt,
                           double vPrev,
                           Eigen::MatrixXd& G,
                           Eigen::VectorXd& I);
    // 梯形法的“电容 + 历史电流源”等效 stamp（Norton 等效）
    static void stampCapTR(int eq1, int eq2,
                           double C, double dt,
                           double vPrev, double iPrev,
                           Eigen::MatrixXd& G,
                           Eigen::VectorXd& I);
};


class HbAnalysis {
private:
    const Circuit& ckt;
    const SimulationConfig& sim;
    Eigen::VectorXd xdc;   // DC 工作点

    int    N;             // MNA unknown 数
    int    K;             // 谐波阶数（sim.hb.nHarm）
    int    nTimeSamples;  // 一个周期内的时域采样点数（2 的幂）
    double f0;            // 基波频率
    double omega0;        // 2πf0
    double T;             // 1/f0
    std::vector<std::vector<double>> dvRealTable; // 预计算 cos/sin 表
    std::vector<std::vector<double>> dvImagTable;

    // int numUnknowns()  const { return N; }
    // int numHarmonics() const { return K; }

    int numRealVars() const { return 2 * (K + 1) * N; }

    // ---- 帮助函数：实数向量 <-> 复数谐波系数 ----
    using CMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using CVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    struct HbVarIndex {
        int h;
        int node;
        bool isReal;
    };

    HbVarIndex decodeVarIndex(int p) const;
    // --- 实数 -> 谐波系数 ---
    void unpackRealToHarmonics(const Eigen::VectorXd& x,
                               std::vector<CVector>& Vk) const;
    // ---- 时域采样 / 频域变换 ----
    void harmonicsToTimeDomain(const std::vector<CVector>& Vk,
                               std::vector<Eigen::VectorXd>& x_t) const;
    // --- 时域非线性电流 -> 谐波电流谱（DFT）---
    void timeDomainToHarmonics(
        const std::vector<Eigen::VectorXd>& Inl_t,
        std::vector<CVector>& Inl_k) const;

        // 在某个时刻，根据 v(t) 评估所有 MOS 的非线性电流，组成 Inl
    void evalNonlinearCurrentsAtTime(const Eigen::VectorXd& v_t,
        Eigen::VectorXd& Inl_t,
        Eigen::MatrixXd& Gnl_t) const;

    void evalMosCapChargeAtTime(
        const Eigen::VectorXd& v_t,
        Eigen::VectorXd& Qcap_t
    ) const;

      // --- 线性部分：频域 Y_k / J_k ---
    void buildLinearYJ(int h, double omega_k, double gmin, double sourceScale,
                       CMatrix& Yk, CVector& Jk) const;

    // --- 残差 F(X) ---
    void computeResidualAndTimeDomainJacobian(
        const Eigen::VectorXd& x, double gmin, double sourceScale,
        Eigen::VectorXd& F,
        std::vector<Eigen::MatrixXd>& Gnl_t_vec) const;

    // 用解析方式构建 Jacobian：
    //   J = ∂F/∂x = (线性 Y_k 部分) + (非线性 gm/gds + Fourier + FFT 部分)
    void buildJacobianAnalytic(
        const Eigen::VectorXd& x, double gmin, double sourceScale,
        const std::vector<Eigen::MatrixXd>& Gnl_t_vec,
        Eigen::MatrixXd& J) const;

    // Newton 求解 F(x)=0，使用 ConvController 做阻尼 + gmin stepping
    bool newtonSolve(Eigen::VectorXd& x, double rampScale, double& bestResidual, bool allowRelax = false) const;

    void writeHbTimeCsv(const Eigen::VectorXd& x, const std::string& outFile) const;
public:
    HbAnalysis(const Circuit& ckt,
               const SimulationConfig& sim,
               const Eigen::VectorXd& dcOp);
    bool run(Eigen::VectorXd& xOut, const std::string& outFile) const;
    bool run(Eigen::VectorXd& xOut) const;
};

class AcAnalysis {
private:
    const Circuit& ckt;
    const SimulationConfig& sim;
    Eigen::VectorXd xdc;

    std::vector<double> buildFreqs() const;

public:
    AcAnalysis(const Circuit& ckt, const SimulationConfig& sim,
               const Eigen::VectorXd& dcOp);

    void run() const;
};

// 判定电路中是否存在非线性器件（目前只看 MOS）
inline static bool hasNonlinearDevices(const Circuit& ckt) {
    for (const auto& e : ckt.elements) {
        if (std::dynamic_pointer_cast<MosfetBase>(e)) {
            return true;
        }
    }
    return false;
}

// ========= 小工具：取节点电压 =========
inline static double getNodeVoltage(const Circuit& ckt,
                             const Eigen::VectorXd& x,
                             int nodeId)
{
    int eq = ckt.nodes[nodeId].eqIndex;
    if (eq >= 0 && eq < x.size()) return x(eq);
    return 0.0;
}

// ========= 小工具：全局 gmin-to-ground =========
template<typename MatrixType>
inline static void stampGlobalGmin(const Circuit& ckt,
                            MatrixType& G,
                            double gmin)
{
    int N = G.rows();
    for (const auto& node : ckt.nodes) {
        int eq = node.eqIndex;
        if (eq >= 0 && eq < N) {
            G(eq, eq) += gmin;
        }
    }
}

// Shooting PSS 配置结构体
struct ShootingPssConfig {
    double periodT = 1e-6;   // 周期 T
    double tstep   = 1e-9;   // 时间步长 Δt
    int    maxIters = 50;    // Shooting 外层最大迭代次数
    double tol      = 1e-6;  // 收敛阈值 ||x(T) - x(0)||
    double relax    = 0.5;   // 松弛因子 x0_{k+1} = x0_k + relax * (xT - x0_k)
};

// Shooting 方法主分析入口
void runPssShootingAnalysis(const Circuit& ckt,
                            const SimulationConfig& sim,
                            const ShootingPssConfig& cfg,
                            const std::string& outFile);

// 积分一个周期 T（后向欧拉），并可选输出每步结果
Eigen::VectorXd integrateOnePeriodPssBE(
    const Circuit& ckt,
    const SimulationConfig& sim,
    double dt, double T,
    const Eigen::VectorXd& x0,
    const std::function<void(double, const Eigen::VectorXd&)>& dumpRow
);
