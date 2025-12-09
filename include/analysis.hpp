#pragma once

#include<Eigen/Dense>
#include "circuit.hpp"
#include "sim.hpp"
#include "solver.hpp"
#include <string>


enum class DcSolverKind {
    LU,
    GaussSeidel
};

struct ConvStatus {
    Eigen::VectorXd xNext;
    double alphaNext;
    double gminNext;
    double error;
    bool converged;
};

class ConvController {
private:
    double alphaMin, alphaMax;
    double gminHighBase, gminLowBase, gminAbsMax;

    double fastConvRatio;
    double slowConvRatio;
public:
    ConvController();

    ConvStatus update(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& xRaw,
        double prevErr,
        int    iter,
        double alphaCurrent,
        double gminCurrent,
        double rampScale,
        double tolStep
    ) const;
    // 给定 ramp 进度，计算当前步的基础 gmin
    double baseGmin(double rampScale) const {
        rampScale = std::clamp(rampScale, 0.0, 1.0);
        return (gminHighBase - gminLowBase) * (1.0 - rampScale) + gminLowBase * rampScale;
    }

    // LU 版的初始 alpha
    double initialAlphaLU() const { return 0.35; }

    // GS 版可以稍微激进一点
    double initialAlphaGS() const { return 0.45; }
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
    // --- 实数 <-> 谐波系数 ---
    void unpackRealToHarmonics(const Eigen::VectorXd& x,
                               std::vector<CVector>& Vk) const;
    // ---- 时域采样 / 频域变换 ----
    void harmonicsToTimeDomain(const std::vector<CVector>& Vk,
                               std::vector<Eigen::VectorXd>& x_t) const;
    // --- 时域非线性电流 -> 谐波电流谱（DFT）---
    void timeDomainCurrentsToHarmonics(
        const std::vector<Eigen::VectorXd>& Inl_t,
        std::vector<CVector>& Inl_k) const;

        // 在某个时刻，根据 v(t) 评估所有 MOS 的非线性电流，组成 Inl
    void evalNonlinearCurrentsAtTime(const Eigen::VectorXd& v_t,
        Eigen::VectorXd& Inl_t,
        Eigen::MatrixXd& Gnl_t) const;

      // --- 线性部分：频域 Y_k / J_k ---
    void buildLinearYJ(double omega_k, double gmin,
                       CMatrix& Yk,
                       CVector& Jk) const;

    // --- 残差 F(X) ---
    void computeResidualAndTimeDomainJacobian(
        const Eigen::VectorXd& x, double gmin,
        Eigen::VectorXd& F,
        std::vector<Eigen::MatrixXd>& Gnl_t_vec) const;

    // 用解析方式构建 Jacobian：
    //   J = ∂F/∂x = (线性 Y_k 部分) + (非线性 gm/gds + Fourier + FFT 部分)
    void buildJacobianAnalytic(
        const Eigen::VectorXd& x, double gmin,
        const std::vector<Eigen::MatrixXd>& Gnl_t_vec,
        Eigen::MatrixXd& J) const;

    // Newton 求解 F(x)=0，使用 ConvController 做阻尼 + gmin stepping
    bool newtonSolve(Eigen::VectorXd& x) const;

public:
    HbAnalysis(const Circuit& ckt,
               const SimulationConfig& sim,
               const Eigen::VectorXd& dcOp);
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
