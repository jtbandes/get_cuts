#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <chrono>
#include <cstdio>
#include <string>

// Helper class to display a progress bar on stderr.
class Progress {
    static const int PROGRESS_WIDTH = 60;
    using Clock = std::chrono::steady_clock;

    std::string _name;
    size_t _totalBytes = 0;
    size_t _bytesRead = 0;
    size_t _bytesReadAtLastReport = 0;
    Clock::time_point _startTime;
    Clock::time_point _lastReportTime;

    double secondsSince(Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    }

    void report() {
        double percentRead = _bytesRead / double(_totalBytes);
        int filledWidth = percentRead * PROGRESS_WIDTH;

        double rate = double(_bytesRead - _bytesReadAtLastReport) / 1024 / 1024 / secondsSince(_lastReportTime);

        // "\r" returns to beginning of line, "esc [ K" clears the line
        // https://en.wikipedia.org/wiki/ANSI_escape_code#CSI_sequences
        std::fprintf(stderr, "\r\x1b[K%s [%-*s] %2.1lf%% (%2.1lf MB/s)",
            _name.c_str(), PROGRESS_WIDTH, std::string(filledWidth, '=').c_str(), percentRead * 100, rate);
        std::fflush(stderr);

        _bytesReadAtLastReport = _bytesRead;
        _lastReportTime = Clock::now();
    }

public:
    Progress(std::string name, size_t totalBytes) : _name(std::move(name)), _totalBytes(totalBytes) {
        _lastReportTime = _startTime = Clock::now();
    }

    void addBytesRead(size_t bytesRead) {
        _bytesRead += bytesRead;
        if (_bytesRead > _bytesReadAtLastReport + 10'000'000) {
            report();
        }
    }

    void finish() {
        double totalElapsed = secondsSince(_startTime);
        std::fprintf(stderr, "\r\x1b[K%s [%s] Done in %2.1lfs (%2.1lf MB/s avg)\n",
            _name.c_str(), std::string(PROGRESS_WIDTH, '=').c_str(), totalElapsed, double(_bytesRead) / 1024 / 1024 / totalElapsed);
    }
};
