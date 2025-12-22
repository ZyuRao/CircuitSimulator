#define _USE_MATH_DEFINES
#include "analysis.hpp"
#include "runtime.hpp"
#include <unsupported/Eigen/FFT>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <thread>

struct HbJacobianScratch {
    Eigen::VectorXd deltaCol;
    std::vector<std::complex<double>> fftFreq;
    std::vector<std::complex<double>> fftTime;
    Eigen::FFT<double> fft;
};

struct HbAnalysis::HbWorkspace {
    std::vector<HbAnalysis::CVector> Vk;
    std::vector<HbAnalysis::CVector> linearJk;
    std::vector<HbAnalysis::CMatrix> linearYk;
    std::vector<HbAnalysis::CVector> Ilin_k;
    std::vector<HbAnalysis::CVector> Inl_k;
    std::vector<Eigen::VectorXd> v_t;
    std::vector<Eigen::VectorXd> Inl_t;
    std::vector<Eigen::MatrixXd> Gnl_t;
    Eigen::VectorXd F;
    Eigen::VectorXd Ftry;
    Eigen::VectorXd rhs;
    Eigen::VectorXd dx;
    Eigen::MatrixXd J;
    std::vector<int> qH;
    std::vector<int> qNode;
    std::vector<unsigned char> qIsReal;
    std::vector<unsigned char> qActive;
    std::vector<int> qActiveList;
    std::vector<int> qNodeActive;
    std::vector<const double*> qDvActive;
    std::vector<double*> qColPtr;
    std::vector<const double*> qDvPtr;
    std::vector<int> rowBaseByH;
    std::vector<int> rowBaseImByH;
    Eigen::VectorXd deltaCol;
    std::vector<HbJacobianScratch> jacScratch;
    std::vector<std::complex<double>> fftFreq;
    std::vector<std::complex<double>> fftTime;
    Eigen::FFT<double> fft;
    std::vector<int> luPerm;
};

HbAnalysis::~HbAnalysis() = default;

HbAnalysis::HbAnalysis(const Circuit&          ckt_,
                       const SimulationConfig& sim_,
                       const Eigen::VectorXd&  dcOp_)
    : ckt(ckt_), sim(sim_), xdc(dcOp_) {

    N = ckt.numUnknowns();
    if (!sim.hb.enabled) {
        throw std::runtime_error("HbAnalysis: HB not enabled in SimulationConfig.");
    }

    f0 = sim.hb.f0;    
    K  = sim.hb.nHarm;
    omega0 = 2.0 * M_PI * f0;
    T      = (f0 > 0.0) ? (1.0 / f0) : 0.0;

    // 采样点数：至少满足 Nyquist（2K+1），再乘一个适度安全系数，向上取 2 的幂
    // 8x 过于保守会拖慢大阶谐波；4x 仍能提供足够抗混叠裕度
    int minSamples = std::max(64, 4 * (2 * K + 1));
    nTimeSamples = 1;
    while (nTimeSamples < minSamples) {
        nTimeSamples <<= 1;
    }

    if (N <= 0) {
        std::cerr << "[HB] Warning: circuit has no unknowns.\n";
    }
    if (K < 0) {
        std::cerr << "[HB] Warning: numHarmonics < 0, clamp to 0.\n";
        K = 0;
    }
    if (nTimeSamples <= 0) {
        std::cerr << "[HB] Warning: numTimeSamples <= 0, set to 128.\n";
        nTimeSamples = 128;
    }

    // 预计算时域基函数表，避免每次 Jacobian 重复生成 sin/cos
    dvRealTable.assign(K + 1, std::vector<double>(nTimeSamples, 0.0));
    dvImagTable.assign(K + 1, std::vector<double>(nTimeSamples, 0.0));
    for (int h = 0; h <= K; ++h) {
        for (int n = 0; n < nTimeSamples; ++n) {
            double t   = (static_cast<double>(n) / nTimeSamples) * T;
            double ang = h * omega0 * t;
            if (h == 0) {
                dvRealTable[h][n] = 1.0;
                dvImagTable[h][n] = 0.0;
            } else {
                dvRealTable[h][n] = 2.0 * std::cos(ang);
                dvImagTable[h][n] = -2.0 * std::sin(ang);
            }
        }
    }
    initWorkspace();
}

void HbAnalysis::initWorkspace() {
    workspace = std::make_unique<HbWorkspace>();
    auto& w = *workspace;
    int H = K + 1;
    int nVars = numRealVars();

    w.Vk.assign(H, CVector::Zero(N));
    w.linearJk.assign(H, CVector::Zero(N));
    w.linearYk.assign(H, CMatrix::Zero(N, N));
    w.Ilin_k.assign(H, CVector::Zero(N));
    w.Inl_k.assign(H, CVector::Zero(N));
    w.v_t.assign(nTimeSamples, Eigen::VectorXd::Zero(N));
    w.Inl_t.assign(nTimeSamples, Eigen::VectorXd::Zero(N));
    w.Gnl_t.assign(nTimeSamples, Eigen::MatrixXd::Zero(N, N));
    w.F = Eigen::VectorXd::Zero(nVars);
    w.Ftry = Eigen::VectorXd::Zero(nVars);
    w.rhs = Eigen::VectorXd::Zero(nVars);
    w.dx = Eigen::VectorXd::Zero(nVars);
    w.J = Eigen::MatrixXd::Zero(nVars, nVars);
    w.qH.assign(nVars, 0);
    w.qNode.assign(nVars, 0);
    w.qIsReal.assign(nVars, 0);
    w.qActive.assign(nVars, 1);
    w.qActiveList.clear();
    w.qNodeActive.clear();
    w.qDvActive.clear();
    w.qColPtr.clear();
    w.qDvPtr.assign(nVars, nullptr);
    w.rowBaseByH.assign(H, 0);
    w.rowBaseImByH.assign(H, 0);
    w.deltaCol = Eigen::VectorXd::Zero(nVars);
    w.jacScratch.clear();
    w.fftFreq.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    w.fftTime.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    w.luPerm.clear();
    linearCacheScale = std::numeric_limits<double>::quiet_NaN();

    for (int h = 0; h < H; ++h) {
        w.rowBaseByH[h] = h * 2 * N;
        w.rowBaseImByH[h] = w.rowBaseByH[h] + N;
    }
    if (N > 0) {
        const int block = 2 * N;
        for (int q = 0; q < nVars; ++q) {
            int h = q / block;
            int offset = q - h * block;
            bool isReal = (offset < N);
            int node = isReal ? offset : (offset - N);
            w.qH[q] = h;
            w.qNode[q] = node;
            w.qIsReal[q] = static_cast<unsigned char>(isReal ? 1 : 0);
            w.qActive[q] = static_cast<unsigned char>((h == 0 && !isReal) ? 0 : 1);
            if (w.qActive[q]) {
                if (isReal) {
                    w.qDvPtr[q] = dvRealTable[h].data();
                } else {
                    w.qDvPtr[q] = dvImagTable[h].data();
                }
                w.qActiveList.push_back(q);
                w.qNodeActive.push_back(node);
                w.qDvActive.push_back(w.qDvPtr[q]);
            }
        }
    }
}

HbAnalysis::HbWorkspace& HbAnalysis::ws() const {
    if (!workspace) {
        throw std::runtime_error("HB workspace not initialized");
    }
    return *workspace;
}

HbAnalysis::HbVarIndex HbAnalysis::decodeVarIndex(int p) const {
    HbVarIndex idx;
    int H = K + 1;
    int block = 2 * N;

    idx.h = p / block;
    int offset = p % block;

    if (offset < N) {
        idx.node   = offset;
        idx.isReal = true;
    } else {
        idx.node   = offset - N;
        idx.isReal = false;
    }
    return idx;
}

void HbAnalysis::unpackRealToHarmonics(
    const Eigen::VectorXd& x,
    std::vector<CVector>& Vk
) const {
    int H = K + 1;
    if (static_cast<int>(Vk.size()) != H) {
        Vk.assign(H, CVector::Zero(N));
    } else {
        for (auto& v : Vk) {
            if (v.size() != N) {
                v = CVector::Zero(N);
            } else {
                v.setZero();
            }
        }
    }

    for(int h = 0; h < H; ++h) {
        int base = h * 2 * N;
        for(int i = 0; i < N; i++) {
            double vr = x(base + i);
            double vi = x(base + N + i);
            Vk[h](i) = std::complex<double>(vr, vi);
        }
    }
}


//谐波 时域采样

void HbAnalysis::harmonicsToTimeDomain(
    const std::vector<CVector>& Vk,
    std::vector<Eigen::VectorXd>& v_t
) const {
    HbWorkspace& w = ws();
    int H = K + 1;

    if (static_cast<int>(v_t.size()) != nTimeSamples) {
        v_t.assign(nTimeSamples, Eigen::VectorXd::Zero(N));
    }
    for (auto& vt : v_t) {
        if (vt.size() != N) {
            vt = Eigen::VectorXd::Zero(N);
        } else {
            vt.setZero();
        }
    }

    if (static_cast<int>(w.fftFreq.size()) != nTimeSamples) {
        w.fftFreq.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }
    if (static_cast<int>(w.fftTime.size()) != nTimeSamples) {
        w.fftTime.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }

    const double Nfft = static_cast<double>(nTimeSamples);

    for (int i = 0; i < N; ++i) {
        // 构造「双边谱」：0, 1..K, 其余为 0，负频率用共轭补齐
        std::fill(w.fftFreq.begin(), w.fftFreq.end(), std::complex<double>(0.0, 0.0));
        w.fftFreq[0] = Vk[0](i) * Nfft; // DC
        for (int h = 1; h < H && h < nTimeSamples; ++h) {
            w.fftFreq[h] = Vk[h](i) * Nfft;
            // 负频率位置 N - h，用共轭
            int negIdx = nTimeSamples - h;
            if (negIdx >= 0 && negIdx < nTimeSamples) {
                w.fftFreq[negIdx] = std::conj(Vk[h](i)) * Nfft;
            }
        }
        w.fft.inv(w.fftTime, w.fftFreq);
        
        for (int n = 0; n < nTimeSamples; ++n) {
            v_t[n](i) = w.fftTime[n].real();
        }
    }
}

//时域非线性电流or电荷->FFT频谱域
void HbAnalysis::timeDomainToHarmonics(
    const std::vector<Eigen::VectorXd>& Inl_t, 
    std::vector<CVector>& Inl_k
) const {
    HbWorkspace& w = ws();
    int H = K + 1;
    
    if (static_cast<int>(Inl_k.size()) != H) {
        Inl_k.assign(H, CVector::Zero(N));
    }
    for (auto& v : Inl_k) {
        if (v.size() != N) {
            v = CVector::Zero(N);
        } else {
            v.setZero();
        }
    }

    if (static_cast<int>(w.fftTime.size()) != nTimeSamples) {
        w.fftTime.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }
    if (static_cast<int>(w.fftFreq.size()) != nTimeSamples) {
        w.fftFreq.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }

    for(int i = 0; i < N; ++i) {
        for(int n = 0; n < nTimeSamples; n++) {
            double val = (i < Inl_t[n].size()) ? Inl_t[n](i) : 0.0;
            w.fftTime[n] = std::complex<double>(val, 0.0);
        }

        w.fft.fwd(w.fftFreq, w.fftTime);

        for (int h = 0; h < H && h < static_cast<int>(w.fftFreq.size()); ++h) {
            // 统一用 1/N 缩放，使逆/正变换一致
            Inl_k[h](i) = w.fftFreq[h] / static_cast<double>(nTimeSamples);
        }
    }
}

void HbAnalysis::evalMosCapChargeAtTime(
    const Eigen::VectorXd& v_t,
    Eigen::VectorXd& Qcap_t
) const {
    Qcap_t.setZero(N);

    auto addCapCharge = [&](int eq1, int eq2, double C) {
        if(C <= 0.0) return;
        double v1 = (eq1 >= 0 && eq1 < N) ? v_t(eq1) : 0.0;
        double v2 = (eq2 >= 0 && eq2 < N) ? v_t(eq2) : 0.0;
        double dv = v1 - v2;
        if (eq1 >= 0 && eq1 < N) Qcap_t(eq1) += C * dv;
        if (eq2 >= 0 && eq2 < N) Qcap_t(eq2) -= C * dv; 
    };

    for (const auto* m : ckt.elementCache().mos) {
        const auto& nodes = m->getNodeIds();
        int nD = nodes[0];
        int nG = nodes[1];
        int nS = nodes[2];
        int nB = (nodes.size() > 3) ? nodes[3] : nodes[2];

        int eqD = ckt.nodes[nD].eqIndex;
        int eqG = ckt.nodes[nG].eqIndex;
        int eqS = ckt.nodes[nS].eqIndex;
        int eqB = ckt.nodes[nB].eqIndex;
        //MOS的寄生电容模型
        double Cj0 = m->getCj0();
        double Cg0 = m->getCg0();
        double Cgs = Cg0, Cgd = Cg0,
               CsJ = Cj0, CdJ = Cj0;

        addCapCharge(eqG, eqS, Cgs);//Gate-source
        addCapCharge(eqG, eqD, Cgd);//gate-drain
        addCapCharge(eqS, eqB, CsJ);//source-bulk
        addCapCharge(eqD, eqB, CdJ);//drain-bulk

    }
}

void HbAnalysis::prepareLinearCache(double sourceScale) const {
    HbWorkspace& w = ws();
    int H = K + 1;
    if (static_cast<int>(w.linearYk.size()) != H) {
        w.linearYk.assign(H, CMatrix::Zero(N, N));
    }
    if (static_cast<int>(w.linearJk.size()) != H) {
        w.linearJk.assign(H, CVector::Zero(N));
    }
    if (!std::isnan(linearCacheScale) && std::abs(linearCacheScale - sourceScale) < 1e-15) {
        return;
    }
    for (int h = 0; h < H; ++h) {
        double omega_k = h * omega0;
        buildLinearYJ(h, omega_k, /*gmin=*/0.0, sourceScale, w.linearYk[h], w.linearJk[h]);
    }
    linearCacheScale = sourceScale;
}

void HbAnalysis::evalNonlinearCurrentsAtTime(
    const Eigen::VectorXd& v_t,
    Eigen::VectorXd& Inl_t, Eigen::MatrixXd* Gnl_t
) const {
    Inl_t.setZero(N);
    if (Gnl_t) {
        Gnl_t->setZero(N, N);
    }
    auto addI = [&](int eq, double val) {
        if (eq >= 0 && eq < N) Inl_t(eq) += val;
    };

    for (const auto* m : ckt.elementCache().mos) {
        double Ids = 0.0, gm = 0.0, gds = 0.0;
        m->evalIdsGmGds(ckt, v_t, Ids, gm, gds);

        const auto& nodes = m->getNodeIds();
        int nD = nodes[0];
        int nS = nodes[2];

        int eqD = ckt.nodes[nD].eqIndex;
        int eqS = ckt.nodes[nS].eqIndex;

        addI(eqD,  Ids);
        addI(eqS, -Ids);

        if (Gnl_t) {
            m->stampNonlinearConductance(*Gnl_t, ckt, v_t);
        }
    }
}

void HbAnalysis::buildLinearYJ(
    int h, double omega_k, double gmin, double sourceScale,
    CMatrix& Yk, CVector& Jk
) const {
    Yk = CMatrix::Zero(N, N);
    Jk = CVector::Zero(N);

    AnalysisContext ctx;
    ctx.type        = AnalysisType::HB;
    ctx.omega       = omega_k;
    ctx.sourceScale = sourceScale;
    ctx.hbHarm      = h;
    ctx.hbF0        = f0;

    for (const auto& e : ckt.elements) {
        // 现在 MOS 也可以进来：它的 stampAC 只 stamp 寄生电容
        e->stampAC(Yk, Jk, ckt, ctx);
    }

    stampGlobalGmin(ckt, Yk, gmin);
}


void HbAnalysis::computeResidualAndTimeDomainJacobian(
    const Eigen::VectorXd& x,
    double gmin, double sourceScale,
    Eigen::VectorXd& F,
    bool needGnl
) const {
    int H = K + 1;
    if (x.size() != numRealVars()) {
        throw std::runtime_error("HbAnalysis::computeResidual: x size mismatch.");
    }

    HbWorkspace& w = ws();

    unpackRealToHarmonics(x, w.Vk);
    harmonicsToTimeDomain(w.Vk, w.v_t);

    if (static_cast<int>(w.Inl_t.size()) != nTimeSamples) {
        w.Inl_t.assign(nTimeSamples, Eigen::VectorXd::Zero(N));
    }
    if (needGnl && static_cast<int>(w.Gnl_t.size()) != nTimeSamples) {
        w.Gnl_t.assign(nTimeSamples, Eigen::MatrixXd::Zero(N, N));
    }
    for (int n = 0; n < nTimeSamples; ++n) {
        evalNonlinearCurrentsAtTime(
            w.v_t[n],
            w.Inl_t[n],
            needGnl ? &w.Gnl_t[n] : nullptr
        );
    }

    timeDomainToHarmonics(w.Inl_t, w.Inl_k);

    // 线性部分
    for (int h = 0; h < H; ++h) {
        auto& lin = w.Ilin_k[h];
        if (lin.size() != N) lin = CVector::Zero(N);
        lin.noalias() = w.linearYk[h] * w.Vk[h];
        if (gmin != 0.0) {
            for (int i = 0; i < N; ++i) {
                lin(i) += gmin * w.Vk[h](i);
            }
        }
    }

    F.setZero(numRealVars());
    for (int h = 0; h < H; ++h) {
        int base = h * 2 * N;
        for (int i = 0; i < N; ++i) {
            std::complex<double> val = w.Ilin_k[h](i) - w.linearJk[h](i);
            F(base + i) = val.real();
            F(base + i + N) = val.imag();
        }
    }

    for (int h = 0; h < H; ++h) {
        int base = h * 2 * N;
        for (int i = 0; i < N; ++i) {
            std::complex<double> val = w.Inl_k[h](i);
            F(base + i)     += val.real();
            F(base + i + N) += val.imag();
        }
    }
}

void HbAnalysis::buildJacobianAnalytic(
    const Eigen::VectorXd& x, double gmin, double sourceScale,
    const std::vector<Eigen::MatrixXd>& Gnl_t_vec,
    Eigen::MatrixXd& J
) const {
    const int n = numRealVars();
    const int H = K + 1;
    (void)sourceScale;
    J.setZero(n, n);

    HbWorkspace& w = ws();
    const int hpLimit = std::min(H, static_cast<int>(w.fftFreq.size()));

    // ==== 1) 线性部分：逐谐波块填充，避免 O(n^2) 遍历 ====
    for (int h = 0; h < H; ++h) {
        int rowBase = h * 2 * N;
        const auto& Yk = w.linearYk[h];
        for (int ip = 0; ip < N; ++ip) {
            for (int iq = 0; iq < N; ++iq) {
                std::complex<double> Y = Yk(ip, iq);
                double yr = Y.real();
                double yi = Y.imag();

                int colRe = rowBase + iq;
                int colIm = rowBase + N + iq;

                J(rowBase + ip,     colRe) += yr;
                J(rowBase + ip,     colIm) += -yi;
                J(rowBase + N + ip, colRe) += yi;
                J(rowBase + N + ip, colIm) += yr;
            }
        }
    }
    if (gmin != 0.0) {
        for (int h = 0; h < H; ++h) {
            int rowBase = h * 2 * N;
            for (int ip = 0; ip < N; ++ip) {
                int rowRe = rowBase + ip;
                int rowIm = rowBase + N + ip;
                int colRe = rowBase + ip;
                int colIm = rowBase + N + ip;
                J(rowRe, colRe) += gmin;
                J(rowIm, colIm) += gmin;
            }
        }
    }

    // ==== 2) 非线性部分：用 FFT 同时累积所有 hp，复杂度 O(H*N^2*Nfft log Nfft) ====
    if (static_cast<int>(w.fftTime.size()) != nTimeSamples) {
        w.fftTime.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }
    if (static_cast<int>(w.fftFreq.size()) != nTimeSamples) {
        w.fftFreq.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
    }

    const double fftScale = 1.0 / static_cast<double>(nTimeSamples);

    const int activeCount = static_cast<int>(w.qActiveList.size());
    if (static_cast<int>(w.qColPtr.size()) != activeCount) {
        w.qColPtr.resize(activeCount, nullptr);
    }
    double* jData = J.data();
    for (int qi = 0; qi < activeCount; ++qi) {
        int q = w.qActiveList[qi];
        w.qColPtr[qi] = jData + q * J.rows();
    }
    const int nRows = J.rows();
    auto processRange = [&](int qStart, int qEnd,
                            Eigen::VectorXd& deltaCol,
                            std::vector<std::complex<double>>& fftTime,
                            std::vector<std::complex<double>>& fftFreq,
                            Eigen::FFT<double>& fft) {
        if (deltaCol.size() != nRows) {
            deltaCol = Eigen::VectorXd::Zero(nRows);
        }
        if (static_cast<int>(fftTime.size()) != nTimeSamples) {
            fftTime.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
        }
        if (static_cast<int>(fftFreq.size()) != nTimeSamples) {
            fftFreq.assign(nTimeSamples, std::complex<double>(0.0, 0.0));
        }

        for (int qi = qStart; qi < qEnd; ++qi) {
            int iq = w.qNodeActive[qi];
            const double* dv_dx = w.qDvActive[qi];
            double* Jcol = w.qColPtr[qi];
            if (!dv_dx || !Jcol) continue;
            deltaCol.setZero();
            double* deltaData = deltaCol.data();

            for (int ip = 0; ip < N; ++ip) {
                // wTime = Gnl(ip,iq,t) * dv_dx(t)
                for (int n_t = 0; n_t < nTimeSamples; ++n_t) {
                    const auto& Gnl = Gnl_t_vec[n_t];
                    const double* gcol = Gnl.data() + static_cast<std::ptrdiff_t>(iq) * N;
                    double g = gcol[ip];
                    fftTime[n_t] = std::complex<double>(g * dv_dx[n_t], 0.0);
                }

                fft.fwd(fftFreq, fftTime);

                for (int hp = 0; hp < hpLimit; ++hp) {
                    std::complex<double> dInl = fftFreq[hp] * fftScale;

                    int rowBase = w.rowBaseByH[hp];
                    int rowBaseIm = w.rowBaseImByH[hp];
                    int rowRe     = rowBase + ip;
                    int rowIm     = rowBaseIm + ip;

                    deltaData[rowRe] += dInl.real();
                    deltaData[rowIm] += dInl.imag();
                }
            }
            for (int row = 0; row < nRows; ++row) {
                Jcol[row] += deltaData[row];
            }
        }
    };

    const auto& opts = runtimeOptions();
    int numThreads = 1;
    if (opts.parallel && activeCount > 1) {
        numThreads = (opts.threads > 0)
            ? opts.threads
            : static_cast<int>(std::thread::hardware_concurrency());
        if (numThreads <= 0) numThreads = 1;
        if (numThreads > activeCount) numThreads = activeCount;
    }

    if (numThreads <= 1) {
        processRange(0, activeCount, w.deltaCol, w.fftTime, w.fftFreq, w.fft);
    } else {
        if (static_cast<int>(w.jacScratch.size()) != numThreads) {
            w.jacScratch.resize(numThreads);
        }
        const int chunk = (activeCount + numThreads - 1) / numThreads;
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for (int t = 0; t < numThreads; ++t) {
            int start = t * chunk;
            int end = std::min(activeCount, start + chunk);
            if (start >= end) break;
            threads.emplace_back([&, start, end, t]() {
                auto& scratch = w.jacScratch[t];
                processRange(start, end,
                             scratch.deltaCol,
                             scratch.fftTime,
                             scratch.fftFreq,
                             scratch.fft);
            });
        }
        for (auto& th : threads) {
            th.join();
        }
    }
}


bool HbAnalysis::newtonSolve(Eigen::VectorXd& x, double rampScale, double& bestResidual, bool allowRelax) const {
    const int n        = numRealVars();
    const int maxIters = 40;
    const int logEvery = 5;

    ConvController ctrl = ConvController::forHb();
    HbWorkspace& w = ws();

    double alpha   = ctrl.params().alphaInit;
    double gmin    = ctrl.baseGmin(rampScale);
    double prevStep= std::numeric_limits<double>::infinity();
    bool loggedSolveSize = false;

    Eigen::VectorXd& F    = w.F;
    Eigen::VectorXd& Ftry = w.Ftry;
    Eigen::MatrixXd& J    = w.J;
    Eigen::VectorXd& dx   = w.dx;
    Eigen::VectorXd bestX = x;
    double bestNorm = std::numeric_limits<double>::infinity();
    bestResidual = bestNorm;

    prepareLinearCache(rampScale);

    for (int it = 0; it < maxIters; ++it) {
        // ===== 1) 残差 F(x) + 时域 Gnl(t) =====
        computeResidualAndTimeDomainJacobian(x, gmin, rampScale, F, /*needGnl=*/true);

        const double normF = F.norm();
        if (sim.verbose && (it == 0 || ((it + 1) % logEvery) == 0)) {
            std::cout << "[HB] iter " << (it + 1) << "/" << maxIters
                      << " scale=" << rampScale
                      << " ||F||=" << normF
                      << " alpha=" << alpha
                      << " gmin=" << gmin << "\n";
        }
        if (normF < bestNorm) {
            bestNorm = normF;
            bestX    = x;
            bestResidual = bestNorm;
        }

        // residualTol 的 rhsScale：用 ||x||∞ 做一个稳定的尺度（比用 ||F|| 自己更靠谱）
        const double rhsScale = std::max(1.0, x.lpNorm<Eigen::Infinity>());
        const double tolF     = ctrl.residualTol(rhsScale);
        const double tolStep  = ctrl.stepTol(rampScale);

        if (normF < tolF) {
            // std::cout << "[HB] Converged by residual.\n";
            return true;
        }

        // ===== 2) 构造 Jacobian =====
        buildJacobianAnalytic(x, gmin, rampScale, w.Gnl_t, J);

        // ===== 3) 解线性方程 J dx = -F =====
        if (sim.verbose && !loggedSolveSize) {
            std::cout << "[HB] solve_linear n=" << J.rows() << "\n";
            loggedSolveSize = true;
        }
        w.rhs = -F;
        bool luOk = Solver::luFactorizeInPlace(J, w.luPerm);
        if (!luOk) {
            std::cerr << "[HB] LU failed, fallback to zero solution.\n";
        }
        if (!luOk) {
            dx.setZero();
        } else {
            Solver::luSolveInPlace(J, w.luPerm, w.rhs, dx);
        }
        if (!dx.allFinite()) {
            // 线性解拒绝：提高 gmin，减小 alpha，重来
            ctrl.onLinearReject(alpha, gmin);
            continue;
        }

        // ===== 4) controller 给一个“基于 step 的推荐 alpha/gmin” =====
        Eigen::VectorXd xRaw = x + dx; // raw Newton step
        ConvStatus st = ctrl.updateStep(x, xRaw, prevStep, it, alpha, gmin, rampScale);

        // 我们用 st 的 alpha/gmin 作为 first-try
        double alphaTry = st.alphaNext;
        double gminTry  = st.gminNext;

        // ===== 5) 残差验收（轻量 backtracking，最多 3 次）=====
        bool accepted = false;
        Eigen::VectorXd xNext;

        const int maxBT = 3; // 控制 FFT 代价，不要开到 12
        for (int bt = 0; bt < maxBT; ++bt) {
            xNext = x + alphaTry * dx;

            computeResidualAndTimeDomainJacobian(xNext, gminTry, rampScale, Ftry, /*needGnl=*/false);
            const double normFtry = Ftry.norm();

            // 只要 residual 下降就接受（HB 这里比 Armijo 更实用、更稳）
            if (normFtry < normF) {
                accepted = true;
                if (normFtry < bestNorm) {
                    bestNorm = normFtry;
                    bestX    = xNext;
                }
                break;
            }

            // 不下降就减小步长
            alphaTry *= 0.5;
            if (alphaTry < ctrl.params().lsAlphaMin) break;
        }

        if (!accepted) {
            // residual 没变好：按统一策略惩罚 alpha/gmin，然后重来
            ctrl.onLinearReject(alpha, gmin);
            continue;
        }

        // ===== 6) 接受步长 =====
        prevStep = (xNext - x).lpNorm<Eigen::Infinity>();
        x = std::move(xNext);
        alpha = std::clamp(alphaTry, ctrl.params().alphaMin, ctrl.params().alphaMax);
        // 接受后把 gmin 拉回基线，避免过大 gmin 造成直流偏差
        double gbase = ctrl.baseGmin(rampScale);
        gmin = 0.5 * gminTry + 0.5 * gbase;

        // ===== 7) step + residual 双判据 =====
        if (prevStep < tolStep) {
            // Ftry 是接受点的残差（上一段已经算过）
            if (Ftry.norm() < tolF) {
                return true;
            }
        }

        if (it == maxIters - 1) {
            std::cerr << "[HB] Newton did not converge: "
                      << "||F||=" << normF
                      << ", step=" << prevStep
                      << ", alpha=" << alpha
                      << ", gmin="  << gmin
                      << ", scale=" << rampScale << "\n";
        }
    }

    // 兜底：如果严格收敛失败，但残差已经足够小，则放宽接受，避免“无输出”
    (void)allowRelax; // 暂不使用松弛验收

    return false;
}

void HbAnalysis::writeHbTimeCsv(
    const Eigen::VectorXd& x, 
    const std::string& outFile
) const {
    std::vector<CVector> Vk;
    std::vector<Eigen::VectorXd> v_t;
    unpackRealToHarmonics(x, Vk);
    harmonicsToTimeDomain(Vk, v_t);
    
    std::ofstream ofs(outFile);
    if (!ofs) {
        std::cerr << "[HB] Cannot open '" << outFile << "' for write.\n";
        return;
    }

    ofs << std::scientific << std::setprecision(9);

    // eqIndex -> node 名字的映射
    std::vector<int> eq2node(N, -1);
    for (int nid = 0; nid < (int)ckt.nodes.size(); ++nid) {
        int eq = ckt.nodes[nid].eqIndex;
        if (eq >= 0 && eq < N) eq2node[eq] = nid;
    }
    // 表头
    ofs << "time";
    for (int eq = 0; eq < N; ++eq) {
        int nid = eq2node[eq];
        if (nid < 0) continue;
        ofs << ",V(" << ckt.nodes[nid].name << ")";
    }
    ofs << "\n";
    // 每个时间采样点一行
    for (int n = 0; n < nTimeSamples; ++n) {
        double t = (static_cast<double>(n) / nTimeSamples) * T;
        ofs << t;
        const auto& vt = v_t[n];
        for (int eq = 0; eq < N; ++eq) {
            int nid = eq2node[eq];
            if (nid < 0) continue;
            double v = (eq >= 0 && eq < vt.size()) ? vt(eq) : 0.0;
            ofs << "," << v;
        }
        ofs << "\n";
    }
    if (sim.verbose) {
        std::cout << "[HB] Time-domain waveform written to '" << outFile << "'.\n";
    }
}

bool HbAnalysis::run(Eigen::VectorXd& xOut) const {
    if (N <= 0 || K < 0 || nTimeSamples <= 0 || f0 <= 0.0) {
        std::cerr << "[HB] Invalid configuration.\n";
        return false;
    }
    if (sim.verbose) {
        std::cout << "[HB] Config: f0=" << f0
                  << " nHarm=" << K
                  << " nTimeSamples=" << nTimeSamples
                  << " nVars=" << numRealVars() << "\n";
    }

    Eigen::VectorXd x(numRealVars());
    x.setZero();
     // 初值：DC 工作点作为 k=0 的实部
    int H = K + 1;
    if (xdc.size() == N) {
        for (int i = 0; i < N; ++i) {
            x(i) = xdc(i); // h=0, 实部
        }
    } else {
        std::cerr << "[HB] Warning: DC operating point size mismatch.\n";
    }

    // 用一个周期的瞬态积分构造初始波形，再做 FFT 得到谐波初值
    Eigen::VectorXd xInitFull;
    bool haveInit = false;
    try {
        const double dt = T / static_cast<double>(nTimeSamples);
        std::vector<Eigen::VectorXd> samples;
        samples.reserve(nTimeSamples);

        TransientAnalysis tran(ckt, sim, "hb_init_dummy.csv");
        if (sim.verbose) {
            std::cout << "[HB] Init: integrating one period (BE), dt=" << dt
                      << " steps=" << nTimeSamples << "\n";
        }
        auto dump = [&](double /*t*/, const Eigen::VectorXd& xt) {
            if (static_cast<int>(samples.size()) < nTimeSamples) {
                samples.push_back(xt);
            }
        };
        tran.integrateOnePeriodBE(xdc, 0.0, dt, T, dump);

        if (static_cast<int>(samples.size()) >= nTimeSamples) {
            samples.resize(nTimeSamples);
            std::vector<CVector> VkGuess;
            timeDomainToHarmonics(samples, VkGuess);

            xInitFull.resize(numRealVars());
            for (int h = 0; h < H; ++h) {
                int base = h * 2 * N;
                for (int i = 0; i < N; ++i) {
                    std::complex<double> v = VkGuess[h](i);
                    xInitFull(base + i)      = v.real();
                    xInitFull(base + N + i)  = v.imag();
                }
            }
            haveInit = true;
        }
    } catch (const std::exception& e) {
        std::cerr << "[HB] Transient init failed: " << e.what() << "\n";
    }

    // 粗阶预解：先用更低的谐波阶求解，再映射到目标阶作为初值
    if (K > 8) {
        int coarseK = std::max(3, std::min(8, K / 2));
        SimulationConfig simCoarse = sim;
        simCoarse.hb.nHarm = coarseK;
        try {
            HbAnalysis coarse(ckt, simCoarse, xdc);
            Eigen::VectorXd xCoarse;
            if (coarse.run(xCoarse)) {
                // 覆盖初始猜测
                x.setZero();
                for (int i = 0; i < N && i < xCoarse.size(); ++i) {
                    x(i) = xCoarse(i);
                }
                for (int h = 1; h <= coarseK; ++h) {
                    int baseSrc = h * 2 * N;
                    int baseDst = h * 2 * N;
                    for (int i = 0; i < N; ++i) {
                        x(baseDst + i)      = xCoarse(baseSrc + i);
                        x(baseDst + N + i)  = xCoarse(baseSrc + N + i);
                    }
                }
                haveInit = true;
            }
        } catch (const std::exception& e) {
            std::cerr << "[HB] Coarse HB init failed: " << e.what() << "\n";
        }
    }

    // 如果拿到了瞬态初值，按最小 ramp 比例缩放作为第一步的起点
    std::vector<double> rampTargets = {0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0};
    if (haveInit) {
        x = xInitFull * rampTargets.front();
    }

    // 连续延拓：如果直接跳到目标 scale 失败，就把步长减半继续逼近
    double prevScale = 0.0;
    Eigen::VectorXd xBestScale = x;
    double bestScaleReached = 0.0;
    double bestResidualSeen = std::numeric_limits<double>::infinity();
    for (double target : rampTargets) {
        if (sim.verbose) {
            std::cout << "[HB] Continuation target scale=" << target << "\n";
        }
        double current = prevScale;
        double step    = target - current;
        int    guard   = 0;

        while (current < target - 1e-6) {
            double next = std::min(target, current + step);
            Eigen::VectorXd xTry = x;
            if (prevScale > 0.0) {
                double factor = next / prevScale;
                xTry *= factor;
            }

            double res = std::numeric_limits<double>::infinity();
            bool ok = newtonSolve(xTry, next, res, /*allowRelax=*/false);
            if (ok) {
                x        = std::move(xTry);
                current   = next;
                prevScale = current;
                if (sim.verbose) {
                    std::cout << "[HB] scale=" << current
                              << " converged (residual=" << res << ")\n";
                }
                if (res < bestResidualSeen) {
                    bestResidualSeen = res;
                    xBestScale = x;
                    bestScaleReached = current;
                }
                // 成功后适当放大步长，加速逼近目标
                step = std::min(step * 1.5, target - current);
                guard = 0;
            } else {
                step *= 0.5;
                guard++;
                if (step < 1e-5 || guard > 20) {
                    std::cerr << "[HB] Continuation failed at scale="
                              << next << " (target=" << target << ")\n";
                    break;
                }
            }
        }
    }

    // 尝试强行推到 scale=1：以当前最佳解为起点，严格验收
    if (prevScale < 1.0 - 1e-6) {
        if (bestScaleReached > 0.0) {
            x = xBestScale;
        }
        double res1 = std::numeric_limits<double>::infinity();
        bool ok1 = newtonSolve(x, 1.0, res1, /*allowRelax=*/false);
        if (!ok1) {
            std::cerr << "[HB] Final push to scale=1 failed, best residual="
                      << res1 << "\n";
            if (bestScaleReached > 0.0) {
                x = xBestScale;
                std::cerr << "[HB] Using best reached scale="
                          << bestScaleReached
                          << " with residual=" << bestResidualSeen << "\n";
            } else {
                return false;
            }
        }
    }

    xOut = x;
    return true;
} 

bool HbAnalysis::run(Eigen::VectorXd& xOut, const std::string& outFile) const {
    bool ok = run(xOut);
    if(!ok) return false;
    writeHbTimeCsv(xOut, outFile);
    return true;
}
