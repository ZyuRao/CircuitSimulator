#pragma once

#include <chrono>
#include <mutex>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>

struct TimingEntry {
    double totalMs = 0.0;
    int count = 0;
};

bool profilingEnabled();

class TimingRegistry {
public:
    void add(const std::string& label, double ms);
    std::vector<std::pair<std::string, TimingEntry>> snapshot() const;
    void print(std::ostream& os) const;
    bool writeCsv(const std::string& path) const;

private:
#ifdef ENABLE_PROF
    mutable std::mutex mu_;
    std::unordered_map<std::string, TimingEntry> entries_;
    std::vector<std::string> order_;
#endif
};

class Timer {
public:
    Timer();
    void reset();
    double elapsedMs() const;

private:
    std::chrono::steady_clock::time_point start_;
};

class ScopedTimer {
public:
    ScopedTimer(TimingRegistry& registry, const std::string& label);
    ScopedTimer(TimingRegistry* registry, const std::string& label);
    ~ScopedTimer();

private:
    TimingRegistry* registry_;
    std::string label_;
    std::chrono::steady_clock::time_point start_;
};
