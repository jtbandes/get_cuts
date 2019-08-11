#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "LineReader.h"

struct NewerFormat {
    enum class Vars : size_t {
        VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_N_SD,
        NUM_VARS
    };

    static const size_t WEIGHT_INSERT_POINT = 1;
    static const size_t Z_INSERT_POINT = 8;
    static const size_t FLAG_INSERT_POINT = 13;
};

struct NewFormat {
    enum class Vars : size_t {
        VAR_NUM, VAR_WEIGHT, VAR_PT, VAR_PSEUDORAP, VAR_PHI, VAR_M, VAR_CONST, VAR_RAP, Z_PX, Z_PY, Z_PZ, Z_E, Z_RAP, GLUON_FLAG_1, GLUON_FLAG_2, VAR_C11, VAR_C10, VAR_ANG1, VAR_ANG05, VAR_N_SD, VAR_C11_SD, VAR_C10_SD, VAR_ANG1_SD,
        NUM_VARS
    };
    static const size_t WEIGHT_INSERT_POINT = 1;
    static const size_t Z_INSERT_POINT = 8;
    static const size_t FLAG_INSERT_POINT = 13;
};

using Jet = std::vector<double>;

template<class Format>
struct CutClause {
    typename Format::Vars var;
    double min;
    double max;

    bool matches(const Jet& jet) const {
        return min <= jet[size_t(var)] && jet[size_t(var)] <= max;
    }
};

template<class Format>
using Cut = std::vector<CutClause<Format>>;

struct CutJetsResult {
    double csOnW = 0;
    std::vector<std::vector<Jet>> jetsList;
};

template<class Format>
CutJetsResult getCutJets(const char* filename, size_t takeNum, const std::vector<Cut<Format>>& cuts, size_t skipNum, bool strict) {
    LineReader reader{filename};

    double totalWeight = 0;
    double crossSection = NAN;  // keep this outside the loop so we can return the last value

    CutJetsResult result;
    result.jetsList.resize(cuts.size());

    reader.nextLine(); // skip header line

    reader.nextLine();
    while (!reader.atEOF()) {
        reader.skip("New Event");

        // Read weight and cross section
        reader.nextLine();
        double weight = reader.readDouble();
        reader.skip(',');
        crossSection = reader.readDouble();
        assert(reader.usedWholeLine());

        totalWeight += weight;

        if (!reader.nextLine()) break;

        // Data which get inserted into each jet in the event
        int isGluon1 = 2;
        int isGluon2 = 2;
        double zData[5];
        std::fill_n(std::begin(zData), 5, INFINITY);

        // Read gluon flag line if present
        if (reader.peek() == 'H') {
            reader.skip('H');
            reader.skipDouble<6>();
            isGluon1 = reader.readDouble();
            isGluon2 = reader.readDouble();
            assert(isGluon1 == 0 || isGluon1 == 1 || isGluon1 == 2 /* ??? */);
            assert(isGluon2 == 0 || isGluon2 == 1 || isGluon2 == 2 /* ??? */);

            if (!reader.nextLine()) break;
        }

        // Read muon data if present
        if (reader.peek() == 'M') {
            double muData1[4];
            double muData2[4];
            
            reader.skip('M');
            std::generate_n(std::begin(muData1), 4, [&] { return reader.readDouble(); });
            if (!reader.nextLine()) throw std::runtime_error("Ended after first M line");

            reader.skip('M');
            std::generate_n(std::begin(muData2), 4, [&] { return reader.readDouble(); });

            std::transform(std::begin(muData1), std::end(muData1), std::begin(muData2), std::begin(zData), std::plus{});
            zData[4] = std::log((zData[3] + zData[2]) / (zData[3] - zData[2])) / 2.0;

            if (!reader.nextLine()) break;
        }

        // Read all jets until the next new event
        size_t jetsSeen = 0;
        std::vector<size_t> jetsTaken(cuts.size(), 0);
        do {
            if (reader.peek() == 'N') {  // new event
                break;
            }
            jetsSeen++;
            if (jetsSeen <= skipNum) {
                // skip jets until skipNum is satisfied
                continue;
            }
            if (std::all_of(jetsTaken.begin(), jetsTaken.end(), [&](auto taken) { return taken >= takeNum; })) {
                // skip all remaining jets once takeNum has been satisfied across all cuts
                continue;
            }
            if (strict && jetsSeen > skipNum + takeNum) {
                // in strict mode, skip all remaining jets if takeNum jets have been considered
                continue;
            }

            Jet jet;
            jet.reserve(size_t(Format::Vars::NUM_VARS));
            reader.readCommaSeparatedDoubles(&jet);

            jet.insert(jet.begin() + Format::WEIGHT_INSERT_POINT, weight);
            jet.insert(jet.begin() + Format::Z_INSERT_POINT, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + Format::FLAG_INSERT_POINT, isGluon1);
            jet.insert(jet.begin() + Format::FLAG_INSERT_POINT+1, isGluon2);

            if (jet.size() != size_t(Format::Vars::NUM_VARS)) {
                throw std::length_error(
                    std::string("Expected jet to have ") + std::to_string(size_t(Format::Vars::NUM_VARS)) +
                    " values, but encountered " + std::to_string(jet.size()));
            }

            for (size_t i = 0; i < cuts.size(); i++) {
                if (jetsTaken[i] >= takeNum) {
                    continue;
                }

                bool matches = std::all_of(cuts[i].begin(), cuts[i].end(), [&](const auto& clause) {
                    return clause.matches(jet);
                });
                if (matches) {
                    jetsTaken[i]++;
                    result.jetsList[i].push_back(std::move(jet));
                }
            }
        } while (reader.nextLine());
    }

    result.csOnW = crossSection / totalWeight;
    return result;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: get_cuts input.txt" << std::endl;
        return 1;
    }

    CutJetsResult result = getCutJets<NewFormat>(argv[1], 2, {
        Cut<NewFormat>{
            {NewFormat::Vars::VAR_PT, 150, 175},
            {NewFormat::Vars::VAR_RAP, -2, 2},
            {NewFormat::Vars::VAR_CONST, 1, 100},
            {NewFormat::Vars::VAR_M, 0, 60},
        },
    }, 0, true);
    
    // Naive output as CSV
    std::printf("# crossSection / totalWeight = %lf\n", result.csOnW);
    for (size_t i = 0; i < result.jetsList.size(); i++) {
        std::printf("# Cut %zu: %zu jets taken\n", i, result.jetsList[i].size());
        for (const auto& jet : result.jetsList[i]) {
            for (double val : jet) {
                std::printf("%lf,", val);
            }
            std::printf("\n");
        }
        std::printf("\n\n");
    }

    return 0;
}
