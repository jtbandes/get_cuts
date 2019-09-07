#pragma once

#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <cmath>
#include <cstdint>
#include <string>
#include <map>

#include "Jet.h"

struct IntHistogram {
    size_t varIndex;
    double totalWeight;
    double totalErr;
    std::map<intmax_t, double> binSums;
    std::map<intmax_t, double> binErrs;

    IntHistogram(size_t varIndex) : varIndex(varIndex) {}

    void add(double weight, const Jet& jet) {
        double val = jet[varIndex];
        if (std::fmod(val, 1.0) != 0) {
            throw std::runtime_error("Used integer binning, but encountered non-integer " + std::to_string(val));
        }
        intmax_t key = intmax_t(val);
        binSums[key] += weight;
        binErrs[key] += weight * weight;
        totalWeight += weight;
        totalErr += weight * weight;
    }

    void finish() {
        for (auto& [k, v] : binSums) {
            v /= totalWeight;
        }
        for (auto& [k, v] : binErrs) {
            v = std::sqrt(v) / totalWeight;
        }
    }
};

struct BinHistogram {
    const size_t varIndex;
    const double binWidth;
    double totalWeight = 0;
    double totalErr = 0;
    std::vector<double> binEndpoints;
    std::vector<double> binSums;
    std::vector<double> binErrs;

    BinHistogram(size_t varIndex, double min, double max, size_t nBins)
        : varIndex(varIndex)
        , binWidth((max - min) / nBins)
    {
        if (nBins == 0) {
            throw std::invalid_argument("Histogram must have at least 1 bin");
        }
        for (size_t i = 0; i <= nBins; i++) {
            binEndpoints.push_back(min + (max - min) * i / nBins);
        }
        binSums.resize(nBins, 0);
        binErrs.resize(nBins, 0);
    }

    void add(double weight, const Jet& jet) {
        double val = jet[varIndex];
        auto iter = std::upper_bound(binEndpoints.begin(), binEndpoints.end(), val);
        size_t binIdx = iter - binEndpoints.begin();
        if (binIdx == 0) {
            // value is less than all bins
        } else if (binIdx <= binSums.size()) {
            binSums[binIdx - 1] += weight;
            binErrs[binIdx - 1] += weight * weight;
            totalWeight += weight;
            totalErr += weight * weight;
        } else if (val == binEndpoints.back()) {
            // value falls in last bin (inclusive)
            binSums.back() += weight;
            binErrs.back() += weight * weight;
            totalWeight += weight;
            totalErr += weight * weight;
        } else {
            // value is greater than all bins
        }
    }

    void finish() {
        for (auto& v : binSums) {
            v = v / binWidth / totalWeight;
        }
        for (auto& v : binErrs) {
            v = std::sqrt(v) / binWidth / totalWeight;
        }
    }
};
