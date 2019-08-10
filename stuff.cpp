// compile:
//   clang++ -stdlib=libc++ -std=c++17 -O3 get_cuts.cpp -o get_cuts

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
// #include <charconv>
#include <sstream>

enum class Var : size_t {
    VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_N_SD,
// VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_C11, VAR_C10, VAR_ANG1, VAR_ANG05, VAR_N_SD, VAR_C11_SD, VAR_C10_SD, VAR_ANG1_SD
    NUM_VARS
};

const size_t WEIGHT_INSERT_POINT = 1;
const size_t Z_INSERT_POINT = 8;
const size_t FLAG_INSERT_POINT = 13;


struct CutClause {
    Var var;
    double min;
    double max;

    bool matches(const std::vector<double>& vals) const {
        return min <= vals[size_t(var)] && vals[size_t(var)] <= max;
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

void getCutJets(const char* filename, size_t takeNum, const std::vector<Cut>& cuts, size_t skipNum, bool strict) {
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
    
    auto readNextLine = [&] {
        lineNum++;
        return !!std::getline(file, line);
    };

    using Jet = std::vector<double>;

    std::vector<std::vector<Jet>> jetsTaken{cuts.size()};

    readNextLine(); // skip header line

    readNextLine();
    while (true) {
        if (lineNum % 100000 == 0) {
            double amountRead = std::round(file.tellg() / double(fileSize) * 100'00) / 100;
            std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTime;
            std::printf("Read %zu lines (%.1f%%) in %.1fs\n", lineNum, amountRead, elapsed.count());
        }

        assert(line == "New Event");
        readNextLine();

        // lineStreambuf.reset(line);
        // std::istream lineStream{&lineStreambuf};
        // lineStream >> weight;

        // char* p = line.begin();
        // p = readDouble(p, line.end(), &weight);

        int readChars;

        double weight;
        double crossSection;
        assert(std::sscanf(line.data(), "%lf, %lf%n", &weight, &crossSection, &readChars) == 2 && readChars == line.size());
        if (!readNextLine()) throw std::runtime_error("Ended after header");

        totalWeight += weight;

        int isGluon1 = 2;
        int isGluon2 = 2;
        if (line[0] == 'H') {
            assert(std::sscanf(line.data(), "H %*d %*d %*d %*d %*d %*d %d %d", &isGluon1, &isGluon2) == 2);
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
            assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData1[0], &muData1[1], &muData1[2], &muData1[3], &readChars) == 4 && readChars == line.size());
            if (!readNextLine()) throw std::runtime_error("Ended after M1");
            assert(std::sscanf(line.data(), "M %lf %lf %lf %lf %*d%n", &muData2[0], &muData2[1], &muData2[2], &muData2[3], &readChars) == 4 && readChars == line.size());
            if (!readNextLine()) throw std::runtime_error("Ended after M2");

            std::transform(std::begin(muData1), std::end(muData1), std::begin(muData2), std::begin(zData), std::plus{});
            zData[4] = std::log((zData[3] + zData[2]) / (zData[3] - zData[2])) / 2.0;
        }

        if (line == "New Event") {
            continue;
        }
        while (readNextLine()) {
            if (line == "New Event") {
                break;
            }
            Jet jet;
            char* ptr = line.data();
            double val;
            // std::cout<<"trying "<<line<<std::endl;
            while (std::sscanf(ptr, " %lf%n", &val, &readChars) == 1) {
                // std::cout<<"read "<<readChars<<" chars " << val<<std::endl;
                ptr += readChars;
                if (ptr != line.data() + line.size()) {
                    ptr += 1; // skip comma
                }
                jet.push_back(val);
            }
            // std::cout<<"got jet "<<jet.size()<< "   " << line<<std::endl;
            assert(jet.size() == size_t(Var::NUM_VARS) - 8);

            jet.insert(jet.begin() + WEIGHT_INSERT_POINT, weight);
            jet.insert(jet.begin() + Z_INSERT_POINT, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + FLAG_INSERT_POINT, isGluon1);
            jet.insert(jet.begin() + FLAG_INSERT_POINT+1, isGluon2);

            assert(jet.size() == size_t(Var::NUM_VARS));

            for (size_t i = 0; i < cuts.size(); i++) {
                if (jetsTaken[i].size() >= takeNum) {
                    continue;
                }

                bool matches = std::all_of(cuts[i].begin(), cuts[i].end(), [&](auto& clause) {
                    return clause.matches(jet);
                });
                if (matches) {
                    jetsTaken[i].push_back(std::move(jet));
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    getCutJets(argv[1], 2, {
        Cut{
            {Var::VAR_PT, 150, 175},
            {Var::VAR_RAP, -2, 2},
            {Var::VAR_CONST, 1, 100},
            {Var::VAR_M, 0, 60},
        },
    }, 0, true);
    return 0;
}