#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "LineReader.h"
#include "get_cuts.h"


CutJetsResult getCutJets(const Format& format, const char* filename, const GetCutJetsSpec& spec) {
    CutJetsResult result;
    LineReader reader{filename};

    bool useEventProbability = !std::isnan(spec.eventProbabilityMultiplier);
    std::uniform_real_distribution randDouble(0.0, 1.0);
    std::seed_seq seed({spec.randomSeed});
    std::mt19937_64 randEngine(seed);

    double crossSection = NAN;  // keep this outside the loop so we can return the last value

    // result.jetsList.resize(cuts.size());
    // Initialize output histograms based on the specs for each cut
    for (const auto& cut : spec.cuts) {
        result.cutResults.push_back(CutResult{
            .intHistograms = cut.intHistograms,
            .binHistograms = cut.binHistograms,
        });
    }

    reader.nextLine(); // skip header line

    reader.nextLine();
    while (!reader.atEOF()) {
        reader.skip("New Event");
        reader.nextLine();

        double weight = reader.readDouble();
        reader.skip(',');

        bool keepEvent = !useEventProbability || randDouble(randEngine) < weight * spec.eventProbabilityMultiplier;

        if (keepEvent) {
            ++result.numEvents;
            result.totalWeight += weight;
            crossSection = reader.readDouble();
        } else {
            reader.readDouble();
        }

        assert(reader.usedWholeLine());

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
        std::vector<size_t> jetsTaken(spec.cuts.size(), 0);
        do {
            if (reader.peek() == 'N') {  // new event
                break;
            }
            if (!keepEvent) {
                continue;
            }
            jetsSeen++;
            if (jetsSeen <= spec.skipNum) {
                // skip jets until skipNum is satisfied
                continue;
            }
            if (std::all_of(jetsTaken.begin(), jetsTaken.end(), [&](auto taken) { return taken >= spec.takeNum; })) {
                // skip all remaining jets once takeNum has been satisfied across all cuts
                continue;
            }
            if (spec.strict && jetsSeen > spec.skipNum + spec.takeNum) {
                // in strict mode, skip all remaining jets if takeNum jets have been considered
                continue;
            }

            Jet jet;
            jet.reserve(format.numVars());
            reader.readCommaSeparatedDoubles(&jet);

            jet.insert(jet.begin() + format.weightInsertPoint, useEventProbability ? 1.0 : weight);
            jet.insert(jet.begin() + format.zInsertPoint, std::begin(zData), std::end(zData));
            jet.insert(jet.begin() + format.flagInsertPoint, isGluon1);
            jet.insert(jet.begin() + format.flagInsertPoint + 1, isGluon2);

            if (jet.size() != format.numVars()) {
                throw std::length_error(
                    std::string("Expected jet to have ") + std::to_string(format.numVars()) +
                    " values, but encountered " + std::to_string(jet.size()));
            }

            for (size_t i = 0; i < spec.cuts.size(); i++) {
                if (jetsTaken[i] >= spec.takeNum) {
                    continue;
                }

                if (spec.cuts[i].matches(jet)) {
                    jetsTaken[i]++;
                    result.cutResults[i].add(useEventProbability ? 1.0 : weight, jet);
                }
            }
        } while (reader.nextLine());
    }

    result.csOnW = crossSection / result.totalWeight;
    result.finish();
    return result;
}
