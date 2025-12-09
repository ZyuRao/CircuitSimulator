#define _USE_MATH_DEFINES
#include "analysis.hpp"
#include <unsupported/Eigen/FFT>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

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

    for (int n = 0; n < nTimeSamples; ++n) {
        double t = (static_cast<double>(n) / nTimeSamples) * T;

        for (int i = 0; i < N; ++i) {
            double v = 0.0;

            // DC 分量
            if (i < Vk[0].size()) {
                v = Vk[0](i).real();
            }

            // 高次谐波
            for (int h = 1; h < H; ++h) {
                double ang = h * omega0 * t;
                double c = std::cos(ang);
                double s = std::sin(ang);
                const std::complex<double>& Vh = Vk[h](i);
                v += 2.0 * ( Vh.real() * c - Vh.imag() * s );
            }

            v_t[n](i) = v;
        }
    }
}

//时域非线性电流->FFT频谱域
void HbAnalysis::timeDomainCurrentsToHarmonics(
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
    double omega_k, double gmin,
    CMatrix& Yk,
    CVector& Jk
) const {
    Yk = CMatrix::Zero(N, N);
    Jk = CVector::Zero(N);
    for (const auto& e : ckt.elements) {
        if (std::dynamic_pointer_cast<MosfetBase>(e)) {
            continue;
        }
        e->stampAC(Yk, Jk, ckt, omega_k);
    }

    stampGlobalGmin(ckt, Yk, gmin);
}


void HbAnalysis::computeResidualAndTimeDomainJacobian(
    const Eigen::VectorXd& x,
    double gmin,
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
        buildLinearYJ(omega_k, gmin, Yk[h], Jk[h]);
        Ilin_k[h] = Yk[h] * Vk[h];
    }

    //3.Vk->时域
    std::vector<Eigen::VectorXd> v_t;
    harmonicsToTimeDomain(Vk, v_t);

    //4.nonlinear cuurent Inl(t_n)
    std::vector<Eigen::VectorXd> Inl_t(nTimeSamples, Eigen::VectorXd::Zero(N));
    Gnl_t_vec.assign(nTimeSamples, Eigen::VectorXd::Zero(N, N));
    for(int n = 0; n < nTimeSamples; ++n) {
        evalNonlinearCurrentsAtTime(v_t[n], Inl_t[n], Gnl_t_vec[n]);
    }

    // 5.Inl(t_n) -> Inl_k via FFT
    std::vector<CVector> Inl_k;
    timeDomainCurrentsToHarmonics(Inl_t, Inl_k);

    //6.频域KCL
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
    const Eigen::VectorXd& x, double gmin,
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

        buildLinearYJ(omega_k, gmin, Yk[h], Jk[h]);
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


bool HbAnalysis::newtonSolve(Eigen::VectorXd& x) const {
    const int n = numRealVars();
    const int maxIters = 20;
    const double tolF = 1e-6;
    const double tolStep = 1e-6;
    
    Eigen::VectorXd F(n);
    ConvController ctrl;
    double alpha = ctrl.initialAlphaLU();
    double gmin = ctrl.baseGmin(1.0);
    double prevErr = std::numeric_limits<double>::infinity();

    for(int it = 0; it < maxIters; ++it) {
        std::vector<Eigen::MatrixXd> Gnl_t_vec;
        //1.当前残差F(x)
        computeResidualAndTimeDomainJacobian(x, gmin, F, Gnl_t_vec);
        double normF = F.norm();
        // std::cout << "[HB] Newton iter " << it
        //           << ", ||F||=" << normF
        //           << ", gmin="  << gmin
        //           << ", alpha=" << alpha << "\n";
        if (normF < tolF) {
            std::cout << "[HB] Converged by residual.\n";
            return true;
        }

        //2.构造解析Jacobian
        Eigen::MatrixXd J(n, n);
        buildJacobianAnalytic(x, gmin, Gnl_t_vec, J);

        //3.解J dx = -F
        Eigen::VectorXd dx = Solver::solveLinearSystemLU(J, -F);
        // if (!dx.allFinite()) {
        //     std::cerr << "[HB] Newton: LU produced NaN/Inf.\n";
        //     return false;
        // }
        auto st = ctrl.update(
            x, x + dx, prevErr, it, alpha, gmin, 1.0, tolStep
        );

        x       = st.xNext;
        alpha   = st.alphaNext;
        gmin    = st.gminNext;
        prevErr = st.error;

        if (st.converged) {
            // 步长收敛后，再检查残差
            computeResidualAndTimeDomainJacobian(x, gmin, F, Gnl_t_vec);
            if (F.norm() < tolF) {
                // std::cout << "[HB] Converged by step & residual.\n";
                return true;
            }
        }

        if (it == maxIters - 1) {
            std::cerr << "[HB] Newton did not converge: "
                      << "||F||=" << normF
                      << ", stepErr=" << st.error
                      << ", alpha=" << alpha
                      << ", gmin="  << gmin << "\n";
        }
    }

    return false;    
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
    bool ok = newtonSolve(x);
    xOut = x;
    return ok;
} 