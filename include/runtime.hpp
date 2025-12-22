#pragma once

struct RuntimeOptions {
    bool parallel = true;   // enable internal parallel loops
    int  threads  = 0;      // 0 = auto (set in init to min(16, hw))
    bool bench    = false;  // run serial+parallel and compare
    double cmpAbsTol = 1e-9;
    double cmpRelTol = 1e-6;
    int  luBlockSize = 16;
    int  luBlockThreshold = 64;
    int  luGemmBlockSize = 0;
    bool luGemmExperiment = false;
    bool luTrailParallel = false;
};

RuntimeOptions& runtimeOptions();
