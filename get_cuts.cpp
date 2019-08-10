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

// char* readDouble(char* first, char* last, double* out) {
//     auto [ptr, ec] = std::from_chars(first, last, static_cast<double&>(*out));
//     if (ec != {}) {
//         throw std::runtime_error("Invalid format: " + first);
//     }
//     return ptr;
// }

struct StrStreambuf : public std::streambuf {
    void reset(std::string& str) {
        setg(str.data(), str.data(), str.data() + str.length());
    }
};

template<typename T>
std::istream& ignore(std::istream& is) {
    T t;
    return is >> t;
}

template<size_t N, typename T>
std::istream& ignore(std::istream& is) {
    if constexpr (N == 0) {
        return is;
    } else {
        return is >> ignore<T> >> ignore<N-1, T>;
    }
}

class LineReader {
    char* _p;
    char* _end;
public:
    void reset(std::string& str) {
        _p = str.data();
        _end = _p + str.length();
    }
    void skip(char c) {
        if (_p == _end || *_p != c) {
            throw std::invalid_argument(std::string("Unable to read ") + c);
        }
        ++_p;
    }
    double readDouble() {
        char* end = _end;
        double val = std::strtod(_p, &end);
        if (end == _p) {
            throw std::invalid_argument("Unable to read double");
        }
        _p = end;
        return val;
    }
    bool atEnd() const {
        return _p == _end;
    }
};

// In order from slowest to fastest:
// istream >> val;
// sscanf();
// val = strtod(...);


// getline(file, std::string& line);
// file.getline(char* line, 1024);

static size_t MAX_LINE_LENGTH = 1024;
static int PROGRESS_WIDTH = 60;

double getCutJets(const char* filename, size_t takeNum, const std::vector<Cut>& cuts, size_t skipNum, bool strict) {
    std::ifstream file{filename};

    // Get total file size
    file.seekg(0, std::ifstream::end);
    std::streampos fileSize = file.tellg();
    file.seekg(0);

    auto startTime = std::chrono::steady_clock::now();

    double totalWeight = 0;

    size_t lineNum = 0;
    std::string line;

    StrStreambuf lineStreambuf;
    std::istream lineStream{&lineStreambuf};

    LineReader lineReader;
    
    size_t bytesReadAtLastReport = 0;
    auto timeAtLastReport = startTime;
    auto readNextLine = [&] {
        lineNum++;

        // std::getline(file, line);
        line.resize(MAX_LINE_LENGTH);
        file.getline(line.data(), line.size() - 1);
        line.resize(std::strlen(line.data()));
        // line.clear();
        // while (true) {
        //     char c = file.rdbuf()->sbumpc();
        //     if (c == '\n') {
        //         break;
        //     }
        //     line.push_back(c);
        // }
        lineReader.reset(line);
        // lineStreambuf.reset(line);
        // lineStream.sync();
        // lineStream.clear();

        if (lineNum % 200000 == 0) {
            auto bytesRead = file.tellg();
            double percentRead = bytesRead / double(fileSize);
            // double percentRead = std::round(bytesRead / double(fileSize) * 100'00) / 100;
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - timeAtLastReport;
            double rate = double(size_t(bytesRead) - bytesReadAtLastReport) / 1024 / 1024 / elapsed.count();
            // std::printf("Read %zu lines (%.1f%%) in %.1fs (%.1fs MB/s)\n", lineNum, percentRead, totalElapsed.count(), );

            int filledWidth = percentRead * PROGRESS_WIDTH;
            // "\r" returns to beginning of line, "ESC [ K" clears the line
            std::printf("\r\x1b[K%s [%-*s] %2.1lf%% (%2.1lf MB/s)",
                filename, PROGRESS_WIDTH, std::string(filledWidth, '=').c_str(), percentRead * 100, rate);
            std::fflush(stdout);

            bytesReadAtLastReport = bytesRead;
            timeAtLastReport = now;
        } else if (!file) {
            double totalElapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
            std::printf("\r\x1b[K%s [%s] Done in %2.1lfs (%2.1lf MB/s avg)\n",
                filename, std::string(PROGRESS_WIDTH, '=').c_str(), totalElapsed, double(fileSize) / 1024 / 1024 / totalElapsed);
        }

        return bool(file);
    };

    std::vector<std::vector<Jet>> jetsList(cuts.size());

    readNextLine(); // skip header line

    readNextLine();
    while (true) {
        assert(line == "New Event");
        readNextLine();

        // lineStreambuf.reset(line);
        // lineStream >> weight;

        // char* p = line.begin();
        // p = readDouble(p, line.end(), &weight);

        int readChars;

        double weight;
        double crossSection;
        // assert(std::sscanf(line.data(), "%lf, %lf%n", &weight, &crossSection, &readChars) == 2 && readChars == line.size());
        weight = lineReader.readDouble();
        lineReader.skip(',');
        crossSection = lineReader.readDouble();

        // lineStream >> weight;
        // assert(lineStream.peek() == ',');
        // lineStream.get();
        // lineStream >> crossSection;
        // assert(lineStream.eof());

        if (!readNextLine()) {
            std::cout<<"Done, "<<jetsList[0].size()<<std::endl;
            return crossSection / totalWeight;
        }

        totalWeight += weight;

        int isGluon1 = 2;
        int isGluon2 = 2;
        if (line[0] == 'H') {
            // assert(std::sscanf(line.data(), "H %*d %*d %*d %*d %*d %*d %d %d", &isGluon1, &isGluon2) == 2);
            lineReader.skip('H');
            lineReader.readDouble();
            lineReader.readDouble();
            lineReader.readDouble();
            lineReader.readDouble();
            lineReader.readDouble();
            lineReader.readDouble();
            isGluon1 = lineReader.readDouble();
            isGluon2 = lineReader.readDouble();
            // lineStream.get();
            // lineStream >> ignore<6, int> >> isGluon1 >> isGluon2;
            // std::cout<<"read gluons " << isGluon1 << ' ' << isGluon2 << std::endl;
            // assert(lineStream.eof());
        
            assert(isGluon1 == 0 || isGluon1 == 1 || isGluon1 == 2 /* ??? */);
            assert(isGluon2 == 0 || isGluon2 == 1 || isGluon2 == 2 /* ??? */);
            if (!readNextLine()) throw std::runtime_error("Ended after H");
        }

// if line[0] == 'M':
//     z_data = [float(x) + float(y) for x, y in zip(line.split()[MU_DATA_START:MU_DATA_END], f.readline().split()[MU_DATA_START:MU_DATA_END])]
//     line_num += 1
//     z_data.append(np.log((z_data[3] + z_data[2]) / (z_data[3] - z_data[2])) / 2)
// else:
//     z_data = [np.inf] * 5
//     f.seek(old_place)
//     line_num -= 1

        double zData[5];
        std::fill_n(zData, 5, std::numeric_limits<double>::infinity());
        if (line[0] == 'M') {
            double muData1[4];
            double muData2[4];
            // assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData1[0], &muData1[1], &muData1[2], &muData1[3], &readChars) == 4 && readChars == line.size());
            lineReader.skip('M');
            muData1[0] = lineReader.readDouble();
            muData1[1] = lineReader.readDouble();
            muData1[2] = lineReader.readDouble();
            muData1[3] = lineReader.readDouble();
            if (!readNextLine()) throw std::runtime_error("Ended after M1");
            // assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData2[0], &muData2[1], &muData2[2], &muData2[3], &readChars) == 4 && readChars == line.size());
            lineReader.skip('M');
            muData2[0] = lineReader.readDouble();
            muData2[1] = lineReader.readDouble();
            muData2[2] = lineReader.readDouble();
            muData2[3] = lineReader.readDouble();
            if (!readNextLine()) throw std::runtime_error("Ended after M2");

            std::transform(std::begin(muData1), std::end(muData1), std::begin(muData2), std::begin(zData), std::plus{});
            zData[4] = std::log((zData[3] + zData[2]) / (zData[3] - zData[2])) / 2.0;
        }

        if (line == "New Event") {
            continue;
        }

        size_t jetsSeen = 0;
        std::vector<size_t> jetsTaken(cuts.size(), 0);
        do {
            if (line == "New Event") {
                break;
            }
            jetsSeen++;
            if (jetsSeen <= skipNum) {
                continue;
            }
            if (std::all_of(jetsTaken.begin(), jetsTaken.end(), [&](auto taken) { return taken >= takeNum; })) {
                // std::cout<<"skip (already taken) "<<jetsTaken[0]<<std::endl;
                continue;
            }
            if (strict && jetsSeen > skipNum + takeNum) {
                // std::cout<<"skip (already seen all)"<<std::endl;
                continue;
            }

            Jet jet;
            jet.reserve(size_t(Var::NUM_VARS));
            char* ptr = line.data();
            double val;
            // std::cout<<"trying "<<line<<std::endl;
            while (true) {
                jet.push_back(lineReader.readDouble());
                if (lineReader.atEnd()) {
                    break;
                } else {
                    lineReader.skip(',');
                }
            }
            // while (std::sscanf(ptr, " %lf%n", &val, &readChars) == 1) {
            //     // std::cout<<"read "<<readChars<<" chars " << val<<std::endl;
            //     ptr += readChars;
            //     if (ptr != line.data() + line.size()) {
            //         ptr += 1; // skip comma
            //     }
            //     jet.push_back(val);
            // }
            // std::cout<<"got jet "<<jet.size()<< "   " << line<<std::endl;
            assert(jet.size() == size_t(Var::NUM_VARS) - 8);

            jet.insert(jet.begin() + WEIGHT_INSERT_POINT, weight);
            jet.insert(jet.begin() + Z_INSERT_POINT, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + FLAG_INSERT_POINT, isGluon1);
            jet.insert(jet.begin() + FLAG_INSERT_POINT+1, isGluon2);

            assert(jet.size() == size_t(Var::NUM_VARS));
            // for (auto f : jet) {
            //     std::cout << f << ", ";
            // }
            // std::cout << std::endl;

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
                    // std::cout<<"TOOK"<<std::endl;
                }
            }
        } while (readNextLine());
        if (jetsSeen>3) {
            // std::cout<<"REad event "<<jetsTaken[0]<<std::endl;
            // throw std::runtime_error("done");
        }
        if (file.eof()) {
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