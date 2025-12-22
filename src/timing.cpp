#include "timing.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace {

using Clock = std::chrono::steady_clock;

double toMs(const Clock::duration& d) {
    return std::chrono::duration<double, std::milli>(d).count();
}

} // namespace

bool profilingEnabled() {
#ifdef ENABLE_PROF
    return true;
#else
    return false;
#endif
}

#ifdef ENABLE_PROF

void TimingRegistry::add(const std::string& label, double ms) {
    std::lock_guard<std::mutex> lock(mu_);
    auto& entry = entries_[label];
    if (entry.count == 0) {
        order_.push_back(label);
    }
    entry.totalMs += ms;
    entry.count += 1;
}

std::vector<std::pair<std::string, TimingEntry>> TimingRegistry::snapshot() const {
    std::lock_guard<std::mutex> lock(mu_);
    std::vector<std::pair<std::string, TimingEntry>> out;
    out.reserve(order_.size());
    for (const auto& label : order_) {
        auto it = entries_.find(label);
        if (it != entries_.end()) {
            out.emplace_back(label, it->second);
        }
    }
    return out;
}

void TimingRegistry::print(std::ostream& os) const {
    auto data = snapshot();
    if (data.empty()) return;

    auto flags = os.flags();
    auto prec = os.precision();

    os.setf(std::ios::fixed);
    os.precision(3);
    os << "Timing summary (ms):\n";
    double totalMs = 0.0;
    for (const auto& item : data) {
        if (item.first == "total_wall") {
            totalMs = item.second.totalMs;
            break;
        }
    }
    if (totalMs <= 0.0) {
        for (const auto& item : data) {
            totalMs += item.second.totalMs;
        }
    }
    os << "  total: " << totalMs << "\n";
    for (const auto& item : data) {
        const auto& entry = item.second;
        double avg = (entry.count > 0) ? (entry.totalMs / entry.count) : 0.0;
        double ratio = (totalMs > 0.0) ? (entry.totalMs / totalMs * 100.0) : 0.0;
        os << "  " << item.first << ": " << entry.totalMs
           << " (n=" << entry.count;
        if (entry.count > 1) {
            os << ", avg=" << avg;
        }
        os << ", " << ratio << "%)\n";
    }

    os.flags(flags);
    os.precision(prec);
}

bool TimingRegistry::writeCsv(const std::string& path) const {
    auto data = snapshot();
    if (data.empty()) return false;

    std::filesystem::path p(path);
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path());
    }

    std::ofstream ofs(path);
    if (!ofs) return false;

    double totalMs = 0.0;
    for (const auto& item : data) {
        if (item.first == "total_wall") {
            totalMs = item.second.totalMs;
            break;
        }
    }
    if (totalMs <= 0.0) {
        for (const auto& item : data) totalMs += item.second.totalMs;
    }

    ofs << "stage,total_ms,count,avg_ms,ratio\n";
    ofs << std::fixed << std::setprecision(6);
    for (const auto& item : data) {
        const auto& entry = item.second;
        double avg = (entry.count > 0) ? (entry.totalMs / entry.count) : 0.0;
        double ratio = (totalMs > 0.0) ? (entry.totalMs / totalMs * 100.0) : 0.0;
        ofs << item.first << "," << entry.totalMs << "," << entry.count << "," << avg << "," << ratio << "\n";
    }
    return true;
}

#else

void TimingRegistry::add(const std::string&, double) {}

std::vector<std::pair<std::string, TimingEntry>> TimingRegistry::snapshot() const {
    return {};
}

void TimingRegistry::print(std::ostream&) const {}

bool TimingRegistry::writeCsv(const std::string&) const { return false; }

#endif

Timer::Timer() {
    reset();
}

void Timer::reset() {
    start_ = Clock::now();
}

double Timer::elapsedMs() const {
    return toMs(Clock::now() - start_);
}

ScopedTimer::ScopedTimer(TimingRegistry& registry, const std::string& label)
    : ScopedTimer(&registry, label) {}

ScopedTimer::ScopedTimer(TimingRegistry* registry, const std::string& label)
    : registry_(registry), label_(label), start_(Clock::now()) {}

ScopedTimer::~ScopedTimer() {
    if (!registry_) return;
#ifdef ENABLE_PROF
    registry_->add(label_, toMs(Clock::now() - start_));
#endif
}
