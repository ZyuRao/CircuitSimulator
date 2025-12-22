#include "runtime.hpp"
#include <cstdlib>
#include <string>
#include <cctype>
#include <thread>

namespace {

int parseIntEnv(const char* name, int defaultValue, bool& wasSet) {
    const char* v = std::getenv(name);
    if (!v) {
        wasSet = false;
        return defaultValue;
    }
    wasSet = true;
    char* endp = nullptr;
    long out = std::strtol(v, &endp, 10);
    if (endp == v) return defaultValue;
    if (out < 0) return defaultValue;
    return static_cast<int>(out);
}

bool parseBoolEnv(const char* name, bool defaultValue, bool& wasSet) {
    const char* v = std::getenv(name);
    if (!v) {
        wasSet = false;
        return defaultValue;
    }
    wasSet = true;
    char* endp = nullptr;
    long out = std::strtol(v, &endp, 10);
    if (endp != v) {
        return out != 0;
    }
    std::string s(v);
    for (auto& ch : s) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    if (s == "false" || s == "off" || s == "no") return false;
    if (s == "true" || s == "on" || s == "yes") return true;
    return defaultValue;
}

void initRuntimeOptions(RuntimeOptions& opts) {
    bool setParallel = false;
    opts.parallel = parseBoolEnv("CSIM_PARALLEL", opts.parallel, setParallel);
    bool setThreads = false;
    opts.threads = parseIntEnv("CSIM_THREADS", opts.threads, setThreads);
    bool setBlock = false;
    opts.luBlockSize = parseIntEnv("CSIM_LU_BLOCK", opts.luBlockSize, setBlock);
    bool setThreshold = false;
    opts.luBlockThreshold = parseIntEnv("CSIM_LU_THRESHOLD", opts.luBlockThreshold, setThreshold);
    bool setGemm = false;
    opts.luGemmBlockSize = parseIntEnv("CSIM_LU_GEMM_BLOCK", opts.luGemmBlockSize, setGemm);
    bool setGemmExp = false;
    opts.luGemmExperiment = parseBoolEnv("CSIM_LU_GEMM_EXPERIMENT",
                                         opts.luGemmExperiment,
                                         setGemmExp);
    bool setTrailPar = false;
    opts.luTrailParallel = parseBoolEnv("CSIM_LU_TRAIL_PAR",
                                        opts.luTrailParallel,
                                        setTrailPar);

    if (opts.luBlockSize <= 0) opts.luBlockSize = 16;
    if (opts.luBlockThreshold <= 0) opts.luBlockThreshold = opts.luBlockSize * 4;
    if (!setGemm) opts.luGemmBlockSize = 0;
    if (opts.luGemmBlockSize < 0) opts.luGemmBlockSize = 0;
    if (opts.threads < 0) opts.threads = 0;
    if (opts.threads == 0) {
        unsigned int hc = std::thread::hardware_concurrency();
        int hw = hc > 0 ? static_cast<int>(hc) : 16;
        int cap = 16;
        opts.threads = (hw < cap) ? hw : cap;
    }
}

} // namespace

RuntimeOptions& runtimeOptions() {
    static RuntimeOptions opts;
    static bool inited = false;
    if (!inited) {
        initRuntimeOptions(opts);
        inited = true;
    }
    return opts;
}
