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

size_t indexOf(const std::vector<std::string>& v, const std::string& x) {
    if (auto found = std::find(v.begin(), v.end(), x); found != v.end()) {
        return found - v.begin();
    }
    throw std::runtime_error(x + " not found");
}

struct Format {
    const std::vector<std::string> vars;
    const size_t weightInsertPoint;
    const size_t zInsertPoint;
    const size_t flagInsertPoint;

    Format(const std::vector<std::string>& vars)
        : vars(vars)
        , weightInsertPoint(indexOf(vars, "VAR_WEIGHT"))
        , zInsertPoint(indexOf(vars, "Z_PX"))
        , flagInsertPoint(indexOf(vars, "GLUON_FLAG_1"))
    {}

    size_t numVars() const {
        return vars.size();
    }

    size_t var(const std::string& name) const {
        return indexOf(vars, name);
    }
};

Format NewerFormat({
    "VAR_NUM", "VAR_WEIGHT", "VAR_PT", "VAR_PSEUDORAP", "VAR_PHI", "VAR_M", "VAR_CONST", "VAR_RAP", "Z_PX", "Z_PY", "Z_PZ", "Z_E", "Z_RAP", "GLUON_FLAG_1", "GLUON_FLAG_2", "VAR_N_SD",
});

Format NewFormat({
    "VAR_NUM", "VAR_WEIGHT", "VAR_PT", "VAR_PSEUDORAP", "VAR_PHI", "VAR_M", "VAR_CONST", "VAR_RAP", "Z_PX", "Z_PY", "Z_PZ", "Z_E", "Z_RAP", "GLUON_FLAG_1", "GLUON_FLAG_2", "VAR_C11", "VAR_C10", "VAR_ANG1", "VAR_ANG05", "VAR_N_SD", "VAR_C11_SD", "VAR_C10_SD", "VAR_ANG1_SD",
});

using Jet = std::vector<double>;

struct CutClause {
    size_t varIndex;
    double min;
    double max;

    bool matches(const Jet& jet) const {
        if (varIndex >= jet.size()) {
            throw std::out_of_range("Variable " + std::to_string(varIndex) + " out of range");
        }
        return min <= jet[varIndex] && jet[varIndex] <= max;
    }
};

struct HistogramSpec {
    size_t varIndex;
    std::vector<double> binEndpoints;
};

struct Histogram {
    HistogramSpec spec;
    std::vector<double> binSums;

    Histogram(const HistogramSpec& spec) {
        if (spec.binEndpoints.size() < 2) {
            throw std::invalid_argument("HistogramSpec must have at least 2 bin endpoints");
        }
        binSums.resize(spec.binEndpoints.size() - 1, 0);
    }

    void add(double weight, const Jet& jet) {
        auto iter = std::lower_bound(spec.binEndpoints.begin(), spec.binEndpoints.end(), jet[spec.varIndex]);
        size_t binIdx = iter - spec.binEndpoints.begin();
        if (binIdx == 0 || binIdx + 1 >= spec.binEndpoints.size()) {
            // value falls outside all bins
            return;
        }
        binSums[binIdx] += weight;
    }
};

struct Cut {
    std::vector<CutClause> clauses;
    std::vector<HistogramSpec> histogramSpecs;
};

struct CutJetsResult {
    double csOnW = 0;
    std::vector<std::vector<Histogram>> histograms;
    // std::vector<std::vector<Jet>> jetsList;
};

CutJetsResult getCutJets(const Format& format, const char* filename, size_t takeNum, const std::vector<Cut>& cuts, size_t skipNum, bool strict) {
    LineReader reader{filename};

    double totalWeight = 0;
    double crossSection = NAN;  // keep this outside the loop so we can return the last value

    CutJetsResult result;
    // result.jetsList.resize(cuts.size());
    // Initialize histograms based on the specs for each cut
    for (size_t cutIdx = 0; cutIdx < cuts.size(); cutIdx++) {
        for (const auto& spec : cuts[cutIdx].histogramSpecs) {
            result.histograms[cutIdx].emplace_back(spec);
        }
    }

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
            jet.reserve(format.numVars());
            reader.readCommaSeparatedDoubles(&jet);

            jet.insert(jet.begin() + format.weightInsertPoint, weight);
            jet.insert(jet.begin() + format.zInsertPoint, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + format.flagInsertPoint, isGluon1);
            jet.insert(jet.begin() + format.flagInsertPoint + 1, isGluon2);

            if (jet.size() != format.numVars()) {
                throw std::length_error(
                    std::string("Expected jet to have ") + std::to_string(format.numVars()) +
                    " values, but encountered " + std::to_string(jet.size()));
            }

            for (size_t i = 0; i < cuts.size(); i++) {
                if (jetsTaken[i] >= takeNum) {
                    continue;
                }

                bool matches = std::all_of(cuts[i].clauses.begin(), cuts[i].clauses.end(), [&](const auto& clause) {
                    return clause.matches(jet);
                });
                if (matches) {
                    jetsTaken[i]++;
                    // result.jetsList[i].push_back(jet);
                }
            }
        } while (reader.nextLine());
    }

    result.csOnW = crossSection / totalWeight;
    return result;
}

int main(int argc, char** argv) {
    std::vector<std::string> args(argv+1, argv+argc);
    if (args.size() < 5) {
        std::cerr << "Usage: get_cuts [--new|--newer] input.txt VAR_1 min1 max1 VAR_2 min2 max2 ..." << std::endl;
        return 1;
    }

    // iterator to help consume arguments one by one
    auto currentArg = args.begin();

    const Format* format;
    {
        const auto& formatArg = *currentArg++;
        if (formatArg == "--new") {
            format = &NewFormat;
        } else if (formatArg == "--newer") {
            format = &NewerFormat;
        } else {
            throw std::runtime_error("Expected --new or --newer");
        }
    }

    const auto& filename = *currentArg++;
    if ((args.end() - currentArg) % 3 != 0) {
        throw std::runtime_error(
            "Wrong number of arguments after filename; expected multiple of 3 but have "
            + std::to_string(args.end() - currentArg));
    }

    // only support a single cut for now; could change this,
    // but would need a way to delineate the cuts in the list of arguments
    Cut cut;
    while (currentArg != args.end()) {
        size_t varIndex = format->var(*currentArg++);
        double min = std::atof((*currentArg++).c_str());
        double max = std::atof((*currentArg++).c_str());
        cut.clauses.push_back({varIndex, min, max});
    }

    CutJetsResult result = getCutJets(*format, filename.c_str(), 2, {cut}, 0, true);
    
    // Naive output as CSV
    // std::printf("# crossSection / totalWeight = %lf\n", result.csOnW);
    // for (size_t i = 0; i < result.jetsList.size(); i++) {
    //     std::printf("# Cut %zu: %zu jets taken\n", i, result.jetsList[i].size());
    //     for (const auto& jet : result.jetsList[i]) {
    //         for (double val : jet) {
    //             std::printf("%lf,", val);
    //         }
    //         std::printf("\n");
    //     }
    //     std::printf("\n\n");
    // }

    return 0;
}
