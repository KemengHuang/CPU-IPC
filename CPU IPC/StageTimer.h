#pragma once

#include <chrono>

class StageTimer final {
public:
    void start()
    {
        startTime_ = Clock::now();
    }

    void stop()
    {
        endTime_ = Clock::now();
    }

    double milliseconds() const
    {
        return std::chrono::duration<double, std::milli>(
            endTime_ - startTime_).count();
    }

private:
    using Clock = std::chrono::steady_clock;
    Clock::time_point startTime_{};
    Clock::time_point endTime_{};
};
