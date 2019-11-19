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
    const std::string varName;
    const size_t varIndex;
    double totalWeight;
    double totalErr;
    std::map<intmax_t, double> binSums;
    std::map<intmax_t, double> binErrs;

    IntHistogram(const std::string& varName, size_t varIndex)
        : varName(varName)
        , varIndex(varIndex) {}

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
    const std::string varName;
    const size_t varIndex;
    double totalWeight = 0;
    double totalErr = 0;
    std::vector<double> binEndpoints;
    std::vector<double> binSums;
    std::vector<double> binErrs;

    BinHistogram(const std::string& varName, size_t varIndex, std::vector<double>&& binEndpoints)
        : varName(varName)
        , varIndex(varIndex)
        , binEndpoints(binEndpoints)
    {
        if (binEndpoints.size() < 2) {
            throw std::invalid_argument("Histogram must have at least 1 bin");
        }
        binSums.resize(binEndpoints.size() - 1, 0);
        binErrs.resize(binEndpoints.size() - 1, 0);
    }

    BinHistogram(const std::string& varName, size_t varIndex, double min, double max, size_t nBins)
        : varName(varName)
        , varIndex(varIndex)
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
        for (size_t i = 0; i < binSums.size(); i++) {
            double binWidth = binEndpoints[i + 1] - binEndpoints[i];
            binSums[i] = binSums[i] / binWidth / totalWeight;
            binErrs[i] = std::sqrt(binErrs[i]) / binWidth / totalWeight;
        }
    }
};
