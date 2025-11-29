#pragma once

#define _USE_MATH_DEFINES
#include <string>
#include <vector>
#include <cmath>
#include "utils.hpp"
#define M_PI 3.14159265358979323846

enum class AnalysisType {
    NONE,
    OP,
    DC,
    AC,
    TRAN
};

enum class AcSweepType {
    LIN,
    DEC,
    OCT
};

enum class WaveformType {
    NONE,
    PULSE,
    SIN,
    PWL
};

enum class ProbeKind {
    NodeVoltage,
    DiffVoltage,
    BranchCurrent
};

struct AnalysisContext {
    AnalysisType type = AnalysisType::OP;

    double sourceScale = 1.0;
    double time = 0.0;
    double omega = 0.0;
};

struct PulseSpec {
    double v1 = 0.0;//初始电平
    double v2 = 0.0;//脉冲电平
    double td = 0.0;//延时
    double tr = 0.0;//上升时间
    double tf = 0.0;//下降时间
    double ton = 0.0;
    double per = 0.0; //周期  =0表示单次脉冲
};

struct SinSpec {
    double v0 = 0.0;//DC offset
    double va = 0.0;//幅度
    double freq = 0.0;//Hz
    double td = 0.0;//延时
    double phi = 0.0;//初相位
};

struct PwlSpec {
    std::vector<double> t;//时间点
    std::vector<double> v;//对应电压
};

struct TranWaveform {
    WaveformType type = WaveformType::NONE;
    PulseSpec pulse;
    SinSpec sine;
    PwlSpec pwl;

    double eval(double t) const {
        switch(type) {
            case WaveformType::NONE:
                return 0.0;
            
            case WaveformType::PULSE: {
                if(pulse.per <= 0.0) {
                    //单次脉冲
                    double tau = t - pulse.td;
                    if(tau <= 0.0) return pulse.v1;

                    if(tau < pulse.tr) {
                        double k = clamp01(tau / pulse.tr);
                        return pulse.v1 + k * (pulse.v2 - pulse.v1);
                    } else if (tau < pulse.tr + pulse.ton) {
                        return pulse.v2;
                    } else {
                        double tfall = tau - (pulse.tr + pulse.ton);
                        double k = clamp01(tfall / pulse.tf);
                        return pulse.v2 + k * (pulse.v1 - pulse.v2);
                    }
                } else {
                    //周期性脉冲
                     if (t < pulse.td) return pulse.v1;
                        double tau = std::fmod(t - pulse.td, pulse.per);
                        if (tau < 0.0) tau += pulse.per;

                        if (tau < pulse.tr) {
                            double k = clamp01(tau / pulse.tr);
                            return pulse.v1 + (pulse.v2 - pulse.v1) * k;
                        } else if (tau < pulse.tr + pulse.ton) {
                            return pulse.v2;
                        } else if (tau < pulse.tr + pulse.ton + pulse.tf) {
                            double tfall = tau - (pulse.tr + pulse.ton);
                            double k = clamp01(tfall / pulse.tf);
                            return pulse.v2 + (pulse.v1 - pulse.v2) * k;
                        } else {
                            return pulse.v1;
                        }
                }
            }
                
            case WaveformType::SIN: {
                if(t < sine.td) return sine.v0;
                double tau = t - sine.td;
                double w = 2.0 * M_PI * sine.freq;
                return sine.v0 + sine.va * std::sin(w * tau + sine.phi);
            }

            case WaveformType::PWL: {
                const auto& tt = pwl.t;
                const auto& vv = pwl.v;
                if(tt.empty()) return 0.0;
                if(t <= tt.front()) return vv.front();
                if(t >= tt.back()) return vv.back();
                //线性插值
                for(std::size_t i = 0; i + 1 < tt.size(); ++i) {
                    if( t > tt[i] && t <= tt[i + 1]) {
                        double k = (t - tt[i]) / (tt[i+1] - tt[i]);
                        return vv[i] + (vv[i+1] - vv[i]) * k;
                    }
                }
                return vv.back();
            }
            
            default:
                return 0.0;
        }
    }
};

struct SourceSpec {
    double dcValue = 0.0;
    double acMag = 0.0;
    double acPhaseDeg = 0.0;
    TranWaveform tran;

    double evalDC(double scale) const {
        return dcValue * scale;
    }

    double evalTran(double t) const {
        return dcValue + tran.eval(t);
    }
};

struct DCSweepConfig {
    std::string sourceName;
    double start = 0.0;
    double stop = 0.0;
    double step = 0.0;
};

struct TranConfig {
    bool enabled = false;
    double tstep = 0.0;
    double tstop = 0.0;
    double tstart = 0.0;
};

struct AcConfig {
    bool enabled = false;
    AcSweepType sweepType = AcSweepType::DEC;
    int nPoints = 0;
    double fstart = 0.0;
    double fstop = 0.0;
};

struct ProbeSpec {
    ProbeKind kind = ProbeKind::NodeVoltage;
    std::string expr;

    std::string node1;
    std::string node2;

    std::string eleName;
};

struct PrintCommand {
    AnalysisType analysis = AnalysisType::NONE;
    std::vector<ProbeSpec> probes;
};

class SimulationConfig {
public:
    bool doOp = false;
    std::vector<DCSweepConfig> dcSweeps;
    TranConfig tran;
    AcConfig ac;
    std::vector<PrintCommand> printCommands;


    bool hasAnyAnalysis() const {
        return doOp || !dcSweeps.empty() || tran.enabled || ac.enabled;
    }

    void ensureDefaultOp() {
        doOp = !hasAnyAnalysis();
    }
};