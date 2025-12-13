#define _USE_MATH_DEFINES
#include "analysis.hpp"
#include <unsupported/Eigen/FFT>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <fstream>

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

    // 采样点数：至少满足 Nyquist（2K+1），再乘一个安全系数，向上取 2 的幂
    int minSamples = std::max(64, 8 * (2 * K + 1));
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
    Vk.assign(H, CVector::Zero(N));

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
    int H = K + 1;

    v_t.assign(nTimeSamples, Eigen::VectorXd::Zero(N));

    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> freq(nTimeSamples);
    std::vector<std::complex<double>> time;

    const double Nfft = static_cast<double>(nTimeSamples);

    for (int i = 0; i < N; ++i) {
        // 构造「双边谱」：0, 1..K, 其余为 0，负频率用共轭补齐
        std::fill(freq.begin(), freq.end(), std::complex<double>(0.0, 0.0));
        freq[0] = Vk[0](i) * Nfft; // DC
        for (int h = 1; h < H && h < nTimeSamples; ++h) {
            freq[h] = Vk[h](i) * Nfft;
            // 负频率位置 N - h，用共轭
            int negIdx = nTimeSamples - h;
            if (negIdx >= 0 && negIdx < nTimeSamples) {
                freq[negIdx] = std::conj(Vk[h](i)) * Nfft;
            }
        }
        fft.inv(time, freq);
        
        for (int n = 0; n < nTimeSamples; ++n) {
            v_t[n](i) = time[n].real();
        }
    }
}

//时域非线性电流or电荷->FFT频谱域
void HbAnalysis::timeDomainToHarmonics(
    const std::vector<Eigen::VectorXd>& Inl_t, 
    std::vector<CVector>& Inl_k
) const {
    int H = K + 1;
    
    Inl_k.assign(H, CVector::Zero(N));
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> in(nTimeSamples);
    std::vector<std::complex<double>> out;

    for(int i = 0; i < N; ++i) {
        for(int n = 0; n < nTimeSamples; n++) {
            double val = (i < Inl_t[n].size()) ? Inl_t[n](i) : 0.0;
            in[n] = std::complex<double>(val, 0.0);
        }

        fft.fwd(out, in);

        for (int h = 0; h < H && h < static_cast<int>(out.size()); ++h) {
            // 统一用 1/N 缩放，使逆/正变换一致
            Inl_k[h](i) = out[h] / static_cast<double>(nTimeSamples);
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

    for(const auto& e : ckt.elements) {
        auto m = std::dynamic_pointer_cast<MosfetBase>(e);
        if(!m) continue;

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
        double Cgs = 0.5 * Cj0, Cgd = 0.5 * Cj0,
               CsJ = Cj0, CdJ = Cj0;

        addCapCharge(eqG, eqS, Cgs);//Gate-source
        addCapCharge(eqG, eqD, Cgd);//gate-drain
        addCapCharge(eqS, eqB, CsJ);//source-bulk
        addCapCharge(eqD, eqB, CdJ);//drain-bulk

    }
}

void HbAnalysis::evalNonlinearCurrentsAtTime(
    const Eigen::VectorXd& v_t,
    Eigen::VectorXd& Inl_t, Eigen::MatrixXd& Gnl_t
) const {
    Inl_t.setZero(N);
    Gnl_t.setZero(N, N);
    auto addI = [&](int eq, double val) {
        if (eq >= 0 && eq < N) Inl_t(eq) += val;
    };

    for(const auto& e : ckt.elements) {
        auto m = std::dynamic_pointer_cast<MosfetBase>(e);
        if(!m) continue;

        double Ids = 0.0, gm = 0.0, gds = 0.0;
        m->evalIdsGmGds(ckt, v_t, Ids, gm, gds);

        const auto& nodes = m->getNodeIds();
        int nD = nodes[0];
        int nS = nodes[2];

        int eqD = ckt.nodes[nD].eqIndex;
        int eqS = ckt.nodes[nS].eqIndex;

        addI(eqD,  Ids);
        addI(eqS, -Ids);

        m->stampNonlinearConductance(Gnl_t, ckt, v_t);
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
    std::vector<Eigen::MatrixXd>& Gnl_t_vec
) const {
    int H = K + 1;
    if (x.size() != numRealVars()) {
        throw std::runtime_error("HbAnalysis::computeResidual: x size mismatch.");
    }

    //1.x->Vk
    std::vector<CVector> Vk;
    unpackRealToHarmonics(x, Vk);

    //2.Linear Part
    std::vector<CMatrix> Yk(H);
    std::vector<CVector> Jk(H), Ilin_k(H);
    for(int h = 0; h < H; h++) {
        double omega_k = h * omega0;
        buildLinearYJ(h, omega_k, gmin, sourceScale, Yk[h], Jk[h]);
        Ilin_k[h] = Yk[h] * Vk[h];
    }

    //3.Vk->时域
    std::vector<Eigen::VectorXd> v_t;
    harmonicsToTimeDomain(Vk, v_t);

    //4.时域：非线性电流 i_nl(t)，以及 d i / d v;MOS电容电荷Qcap(t_n)
    Gnl_t_vec.clear();
    Gnl_t_vec.reserve(nTimeSamples);

    std::vector<Eigen::VectorXd> Inl_t(nTimeSamples, Eigen::VectorXd::Zero(N));
    for (int n = 0; n < nTimeSamples; ++n) {
        Eigen::MatrixXd Gnl_t(N, N);
        evalNonlinearCurrentsAtTime(v_t[n], Inl_t[n], Gnl_t);
        Gnl_t_vec.push_back(Gnl_t);
    }

    // std::vector<Eigen::VectorXd> Qcap_t(nTimeSamples, Eigen::VectorXd::Zero(N));
    // for(int n = 0; n < nTimeSamples; ++n) {
    //     evalMosCapChargeAtTime(v_t[n], Qcap_t[n]);
    // }

    //5.时域 -> 频域：非线性电流 I_nl(V); Qcap(t_n) -> Qk -> IQ_k
    std::vector<CVector> Inl_k, Qk;
    timeDomainToHarmonics(Inl_t, Inl_k);
    // timeDomainToHarmonics(Qcap_t, Qk);


    //6.非线性电荷 Q(V) 与 ΩQ(V)
    // std::vector<CVector> Iq_k(H, CVector::Zero(N));
    // for(int h = 0; h < H; ++h) {
    //     double omega_k = h * omega0;
    //     std::complex<double> jw(0.0, omega_k);
    //     Iq_k[h] = jw * Qk[h];
    // }

    
    //7.频域KCL
    F.resize(numRealVars());
    F.setZero();

    for(int h = 0;h < H; ++h) {
        CVector rk = Ilin_k[h] + Inl_k[h] - Jk[h];
        int base = h * 2 * N;
        for(int i = 0; i < N; i++) {
            std::complex<double> val = rk[i];
            F(base + i) = val.real();
            F(base + i + N) = val.imag();
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
    J.setZero(n, n);

    //线性部分
    std::vector<CMatrix> Yk(H);
    std::vector<CVector> Jk(H);

    for (int h = 0; h < H; ++h) {
        double omega_k = h * omega0;
        Yk[h] = CMatrix::Zero(N, N);
        Jk[h] = CVector::Zero(N);

        buildLinearYJ(h, omega_k, gmin, sourceScale, Yk[h], Jk[h]);
    }

    for (int p = 0; p < n; ++p) {
        HbVarIndex vp = decodeVarIndex(p);
        int hp = vp.h;
        int ip = vp.node;
        bool isReP = vp.isReal;

        if (hp >= H || ip >= N) continue;

        for (int q = 0; q < n; ++q) {
            HbVarIndex vq = decodeVarIndex(q);
            int hq    = vq.h;
            int iq    = vq.node;
            bool isReQ = vq.isReal;

            if (hq != hp || iq >= N) continue;

            std::complex<double> Y = Yk[hp](ip, iq);

            double dRe_dRe =  Y.real();
            double dRe_dIm = -Y.imag();
            double dIm_dRe =  Y.imag();
            double dIm_dIm =  Y.real();

            double contrib = 0.0;
            if (isReP && isReQ) contrib = dRe_dRe;
            else if (isReP && !isReQ) contrib = dRe_dIm;
            else if (!isReP && isReQ) contrib = dIm_dRe;
            else                      contrib = dIm_dIm;

            J(p, q) += contrib;
        }
    }

    //非线性
    std::vector<std::vector<double>> cosTable(H, std::vector<double>(nTimeSamples));
    std::vector<std::vector<double>> sinTable(H, std::vector<double>(nTimeSamples));

    for (int h = 0; h < H; ++h) {
        for (int n_t = 0; n_t < nTimeSamples; ++n_t) {
            double t = (static_cast<double>(n_t) / nTimeSamples) * T;
            double ang = h * omega0 * t;
            cosTable[h][n_t] = std::cos(ang);
            sinTable[h][n_t] = std::sin(ang);
        }
    }
    double fftScale = 1.0 /static_cast<double>(nTimeSamples);
    for (int p = 0; p < n; ++p) {
        HbVarIndex vp = decodeVarIndex(p);
        int hp    = vp.h;
        int ip    = vp.node;
        bool isReP = vp.isReal;

        if (hp >= H || ip >= N) continue;

        for (int q = 0; q < n; ++q) {
            HbVarIndex vq = decodeVarIndex(q);
            int hq    = vq.h;
            int iq    = vq.node;
            bool isReQ = vq.isReal;

            if (iq >= N) continue;

            double sum = 0.0;

            for (int n_t = 0; n_t < nTimeSamples; ++n_t) {
                const Eigen::MatrixXd& Gnl_n = Gnl_t_vec[n_t];

                double dInl_dv = Gnl_n(ip, iq); // ∂I_nl(ip,t_n)/∂v(iq,t_n)

                double dv_dx = 0.0;
                if (hq == 0) {
                    // DC 分量：只依赖 Re(V0)
                    dv_dx = isReQ ? 1.0 : 0.0;
                } else {
                    double c = cosTable[hq][n_t];
                    double s = sinTable[hq][n_t];
                    if (isReQ) dv_dx = 2.0 * c;
                    else       dv_dx = -2.0 * s;
                }

                double c_hp = cosTable[hp][n_t];
                double s_hp = sinTable[hp][n_t];

                double dRe_dInl = fftScale * c_hp;
                double dIm_dInl = -fftScale * s_hp;

                double dF_dInl = isReP ? dRe_dInl : dIm_dInl;

                sum += dF_dInl * dInl_dv * dv_dx;
            }

            J(p, q) += sum;
        }
    }
}


bool HbAnalysis::newtonSolve(Eigen::VectorXd& x, double rampScale) const {
    const int n        = numRealVars();
    const int maxIters = 25;

    ConvController ctrl = ConvController::forHb();

    double alpha   = ctrl.params().alphaInit;
    double gmin    = ctrl.baseGmin(rampScale);
    double prevStep= std::numeric_limits<double>::infinity();

    Eigen::VectorXd F(n), Ftry(n);

    for (int it = 0; it < maxIters; ++it) {
        // ===== 1) 残差 F(x) + 时域 Gnl(t) =====
        std::vector<Eigen::MatrixXd> Gnl_t_vec;
        computeResidualAndTimeDomainJacobian(x, gmin, rampScale, F, Gnl_t_vec);

        const double normF = F.norm();

        // residualTol 的 rhsScale：用 ||x||∞ 做一个稳定的尺度（比用 ||F|| 自己更靠谱）
        const double rhsScale = std::max(1.0, x.lpNorm<Eigen::Infinity>());
        const double tolF     = ctrl.residualTol(rhsScale);
        const double tolStep  = ctrl.stepTol(rampScale);

        if (normF < tolF) {
            // std::cout << "[HB] Converged by residual.\n";
            return true;
        }

        // ===== 2) 构造 Jacobian =====
        Eigen::MatrixXd J(n, n);
        buildJacobianAnalytic(x, gmin, rampScale, Gnl_t_vec, J);

        // ===== 3) 解线性方程 J dx = -F =====
        Eigen::VectorXd dx = Solver::solveLinearSystemLU(J, -F);
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

            std::vector<Eigen::MatrixXd> dummy;
            computeResidualAndTimeDomainJacobian(xNext, gminTry, rampScale, Ftry, dummy);
            const double normFtry = Ftry.norm();

            // 只要 residual 下降就接受（HB 这里比 Armijo 更实用、更稳）
            if (normFtry < normF) {
                accepted = true;
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
        gmin  = gminTry;

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
                      << ", gmin="  << gmin << "\n";
        }
    }

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
    std::cout << "[HB] Time-domain waveform written to '" << outFile << "'.\n";
}

bool HbAnalysis::run(Eigen::VectorXd& xOut) const {
    if (N <= 0 || K < 0 || nTimeSamples <= 0 || f0 <= 0.0) {
        std::cerr << "[HB] Invalid configuration.\n";
        return false;
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
    std::vector<double> ramps = {0.05, 0.1, 0.2, 0.4, 0.7, 1.0};
    for (double s : ramps) {
        bool ok = newtonSolve(x, s);
        if (!ok) return false;
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