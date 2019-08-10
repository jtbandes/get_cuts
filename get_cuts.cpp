// compile:
//   clang++ -std=c++17 -O3 get_cuts.cpp -o get_cuts

// run:
//   ./get_cuts newerHIJ.txt

#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <string>
#include <string_view>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iomanip>
// #include <charconv>
#include <sstream>

enum class Var : size_t {
    VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_N_SD,
// VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_C11, VAR_C10, VAR_ANG1, VAR_ANG05, VAR_N_SD, VAR_C11_SD, VAR_C10_SD, VAR_ANG1_SD,
    NUM_VARS
};

const size_t WEIGHT_INSERT_POINT = 1;
const size_t Z_INSERT_POINT = 8;
const size_t FLAG_INSERT_POINT = 13;



using Jet = std::vector<double>;

struct CutClause {
    Var var;
    double min;
    double max;

    bool matches(const Jet& jet) const {
        return min <= jet[size_t(var)] && jet[size_t(var)] <= max;
    }
};
using Cut = std::vector<CutClause>;


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

        // "\r" returns to beginning of line, "ESC [ K" clears the line
        std::printf("\r\x1b[K%s [%-*s] %2.1lf%% (%2.1lf MB/s)",
            _name.c_str(), PROGRESS_WIDTH, std::string(filledWidth, '=').c_str(), percentRead * 100, rate);
        std::fflush(stdout);

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
        std::printf("\r\x1b[K%s [%s] Done in %2.1lfs (%2.1lf MB/s avg)\n",
            _name.c_str(), std::string(PROGRESS_WIDTH, '=').c_str(), totalElapsed, double(_bytesRead) / 1024 / 1024 / totalElapsed);
    }
};

class LineReader {
    static const size_t MAX_LINE_LENGTH = 1024;

    char* _p = nullptr;  // current position in line
    char* _end = nullptr;  // end of line

    char _buf[MAX_LINE_LENGTH];
    std::ifstream _file;
    Progress _progress; // must be after _file since we use _file during initialization

    void checkEnd() {
        if (_p == _end) {
            throw std::out_of_range("Read past end of line");
        }
    }
    
public:
    LineReader(const char* filename)
        : _file(filename, std::ios::in | std::ios::ate /* open at end to get file size */)
        , _progress(filename, _file.tellg())
    {
        _file.seekg(0);
    }
    bool advance() {
        _file.getline(_buf, MAX_LINE_LENGTH - 1);
        size_t len = std::strlen(_buf);

        // Update progress with strlen instead of _file.tellg() which is slow because it actually seeks the file :(
        _progress.addBytesRead(len + 1);

        _p = _buf;
        _end = _p + len;
        if (_file.eof()) {
            _progress.finish();
        }
        return bool(_file);
    }
    void skip(const char* str) {
        checkEnd();
        if (std::strcmp(_p, str) != 0) {
            throw std::runtime_error(std::string("Expected ") + str);
        }
        _p += std::strlen(str);
    }
    void skip(char c) {
        checkEnd();
        if (*_p != c) {
            throw std::runtime_error(std::string("Expected ") + c);
        }
        ++_p;
    }
    char peek() {
        checkEnd();
        return *_p;
    }
    template<size_t N> void skipDouble() {
        if constexpr (N > 0) {
            readDouble();
            skipDouble<N-1>();
        }
    }
    double readDouble() {
        char* end = _end;
        double val = std::strtod(_p, &end);
        if (end == _p) {
            throw std::runtime_error("Unable to read double");
        }
        _p = end;
        return val;
    }
    bool atEnd() const {
        return _p == _end;
    }
    bool atEOF() const {
        return _file.eof();
    }
};

// In order from slowest to fastest:
// istream >> val;
// sscanf();
// val = strtod(...);


// getline(file, std::string& line);
// file.getline(char* line, 1024);

double getCutJets(const char* filename, size_t takeNum, const std::vector<Cut>& cuts, size_t skipNum, bool strict) {
    LineReader lineReader{filename};

    double totalWeight = 0;

    std::vector<std::vector<Jet>> jetsList(cuts.size());

    lineReader.advance(); // skip header line

    lineReader.advance();
    while (true) {
        lineReader.skip("New Event");
        lineReader.advance();

        double weight = lineReader.readDouble();
        lineReader.skip(',');
        double crossSection = lineReader.readDouble();
        assert(lineReader.atEnd());

        if (!lineReader.advance()) {
            std::cout<<"Done, "<<jetsList[0].size()<<std::endl;
            return crossSection / totalWeight;
        }

        totalWeight += weight;

        int isGluon1 = 2;
        int isGluon2 = 2;
        if (lineReader.peek() == 'H') {
            lineReader.skip('H');
            lineReader.skipDouble<6>();
            isGluon1 = lineReader.readDouble();
            isGluon2 = lineReader.readDouble();
        
            assert(isGluon1 == 0 || isGluon1 == 1 || isGluon1 == 2 /* ??? */);
            assert(isGluon2 == 0 || isGluon2 == 1 || isGluon2 == 2 /* ??? */);
            if (!lineReader.advance()) throw std::runtime_error("Ended after H");
        }

        double zData[5];
        std::fill_n(zData, 5, INFINITY);
        if (lineReader.peek() == 'M') {
            double muData1[4];
            double muData2[4];
            // assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData1[0], &muData1[1], &muData1[2], &muData1[3], &readChars) == 4 && readChars == line.size());
            lineReader.skip('M');
            muData1[0] = lineReader.readDouble();
            muData1[1] = lineReader.readDouble();
            muData1[2] = lineReader.readDouble();
            muData1[3] = lineReader.readDouble();
            if (!lineReader.advance()) throw std::runtime_error("Ended after M1");
            // assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData2[0], &muData2[1], &muData2[2], &muData2[3], &readChars) == 4 && readChars == line.size());
            lineReader.skip('M');
            muData2[0] = lineReader.readDouble();
            muData2[1] = lineReader.readDouble();
            muData2[2] = lineReader.readDouble();
            muData2[3] = lineReader.readDouble();
            if (!lineReader.advance()) throw std::runtime_error("Ended after M2");

            std::transform(std::begin(muData1), std::end(muData1), std::begin(muData2), std::begin(zData), std::plus{});
            zData[4] = std::log((zData[3] + zData[2]) / (zData[3] - zData[2])) / 2.0;
        }

        if (lineReader.peek() == 'N') {  // new event
            continue;
        }

        size_t jetsSeen = 0;
        std::vector<size_t> jetsTaken(cuts.size(), 0);
        do {
            if (lineReader.peek() == 'N') {  // new event
                break;
            }
            jetsSeen++;
            if (jetsSeen <= skipNum) {
                continue;
            }
            if (std::all_of(jetsTaken.begin(), jetsTaken.end(), [&](auto taken) { return taken >= takeNum; })) {
                continue;
            }
            if (strict && jetsSeen > skipNum + takeNum) {
                continue;
            }

            Jet jet;
            jet.reserve(size_t(Var::NUM_VARS));
            while (true) {
                jet.push_back(lineReader.readDouble());
                if (lineReader.atEnd()) {
                    break;
                } else {
                    lineReader.skip(',');
                }
            }
            assert(jet.size() == size_t(Var::NUM_VARS) - 8);

            jet.insert(jet.begin() + WEIGHT_INSERT_POINT, weight);
            jet.insert(jet.begin() + Z_INSERT_POINT, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + FLAG_INSERT_POINT, isGluon1);
            jet.insert(jet.begin() + FLAG_INSERT_POINT+1, isGluon2);

            assert(jet.size() == size_t(Var::NUM_VARS));

            for (size_t i = 0; i < cuts.size(); i++) {
                if (jetsTaken[i] >= takeNum) {
                    continue;
                }

                bool matches = std::all_of(cuts[i].begin(), cuts[i].end(), [&](auto& clause) {
                    return clause.matches(jet);
                });
                if (matches) {
                    jetsTaken[i]++;
                    jetsList[i].push_back(std::move(jet));
                }
            }
        } while (lineReader.advance());
        
        if (lineReader.atEOF()) {
            std::cout<<"Done, "<<jetsList[0].size()<<std::endl;
            return crossSection / totalWeight;
        }
    }
    return -1;
}

int main(int argc, char** argv) {
    double csOnW = getCutJets(argv[1], 2, {
        Cut{
            {Var::VAR_PT, 150, 175},
            {Var::VAR_RAP, -2, 2},
            {Var::VAR_CONST, 1, 100},
            {Var::VAR_M, 0, 60},
        },
    }, 0, true);
    std::cout<<"csOnW " <<csOnW<<std::endl;
    return 0;
}