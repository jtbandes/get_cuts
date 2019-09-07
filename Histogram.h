#pragma once

#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <cstdint>
#include <string>
#include <unordered_map>

#include "Jet.h"

struct IntHistogram {
    size_t varIndex;
    std::unordered_map<intmax_t, double> binSums;

    IntHistogram(size_t varIndex) : varIndex(varIndex) {}

    void add(double weight, const Jet& jet) {
        double val = jet[varIndex];
        if (std::fmod(val, 1.0) != 0) {
            throw std::runtime_error("Used integer binning, but encountered non-integer " + std::to_string(val));
        }
        intmax_t key = intmax_t(val);
        binSums[key] += weight;
    }
};

struct BinHistogram {
    size_t varIndex;
    std::vector<double> binEndpoints;
    std::vector<double> binSums;

    BinHistogram(size_t varIndex, double min, double max, size_t nBins) : varIndex(varIndex) {
        if (nBins == 0) {
            throw std::invalid_argument("Histogram must have at least 1 bin");
        }
        for (size_t i = 0; i <= nBins; i++) {
            binEndpoints.push_back(min + (max - min) * i / nBins);
        }
        binSums.resize(nBins, 0);
    }

    void add(double weight, const Jet& jet) {
        double val = jet[varIndex];
        auto iter = std::upper_bound(binEndpoints.begin(), binEndpoints.end(), val);
        size_t binIdx = iter - binEndpoints.begin();
        if (binIdx == 0) {
            // value is less than all bins
            return;
        } else if (binIdx <= binSums.size()) {
            binSums[binIdx - 1] += weight;
        } else if (val == binEndpoints.back()) {
            // value falls in last bin (inclusive)
            binSums.back() += weight;
        }
    }
};
