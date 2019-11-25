#pragma once

#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <cstdint>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

#include "Histogram.h"
#include "Jet.h"

inline size_t indexOf(const std::vector<std::string>& v, const std::string& x) {
    if (auto found = std::find(v.begin(), v.end(), x); found != v.end()) {
        return found - v.begin();
    }
    throw std::runtime_error("unrecognized variable " + x);
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

    bool operator==(const CutClause& other) const {
        return varIndex == other.varIndex && min == other.min && max == other.max;
    }
};

struct Cut {
    std::vector<CutClause> clauses;
    std::vector<IntHistogram> intHistograms;
    std::vector<BinHistogram> binHistograms;

    bool matches(const Jet& jet) const {
        return std::all_of(clauses.begin(), clauses.end(), [&](const auto& clause) {
            return clause.matches(jet);
        });
    }
};

struct CutResult {
    size_t totalJetsTaken = 0;
    std::vector<IntHistogram> intHistograms;
    std::vector<BinHistogram> binHistograms;

    void add(double weight, const Jet& jet) {
        ++totalJetsTaken;
        for (auto& hist : intHistograms) {
            hist.add(weight, jet);
        }
        for (auto& hist : binHistograms) {
            hist.add(weight, jet);
        }
    }

    void finish() {
        for (auto& hist : intHistograms) {
            hist.finish();
        }
        for (auto& hist : binHistograms) {
            hist.finish();
        }
    }
};

struct CutJetsResult {
    double csOnW = 0;
    double totalWeight = 0;
    size_t numEvents = 0;
    std::vector<CutResult> cutResults;

    void finish() {
        for (auto& cutResult : cutResults) {
            cutResult.finish();
        }
    }
};

struct GetCutJetsSpec {
    size_t takeNum;
    size_t skipNum;
    bool strict;
    double eventProbabilityMultiplier = NAN;
    long long randomSeed;
    std::vector<Cut> cuts;

    // Initialize by reading from a specification file (or stdin)
    GetCutJetsSpec(const Format& format, std::string&& str) : GetCutJetsSpec(format, std::istringstream(str)) { }
    GetCutJetsSpec(const Format& format, std::istream&& stream) : GetCutJetsSpec(format, stream) { }
    GetCutJetsSpec(const Format& format, std::istream& stream) {
        auto nextWord = [&](const std::string& description) {
            if (stream.eof()) {
                throw std::runtime_error("Expected " + description + " in spec");
            }
            std::string str;
            stream >> str;
            if (!stream) {
                if (str.empty()) {
                    throw std::runtime_error("Expected " + description + " in spec");
                } else if (!stream.eof()) {
                    throw std::runtime_error("Error reading spec; expected " + description);
                }
            }
            return str;
        };
        auto consumeWord = [&](const std::string& expected) {
            std::string actual = nextWord("'" + expected + "'");
            if (actual != expected) {
                throw std::runtime_error("Expected '" + expected + "' but found " + actual);
            }
        };
        auto spaceSeparatedDoublesToEndOfLine = [&]() {
            std::vector<double> result;
            do {
                double value;
                stream >> value;
                if (!stream) {
                    throw std::runtime_error("Error reading spec; expected doubles");
                }
                result.push_back(value);
            } while (stream.peek() != '\n');
            return result;
        };

        consumeWord("takeNum:");
        takeNum = std::atoi(nextWord("integer").c_str());

        consumeWord("skipNum:");
        skipNum = std::atoi(nextWord("integer").c_str());

        consumeWord("strict:");
        {
            std::string strictStr = nextWord("boolean");
            if (strictStr == "true") {
                strict = true;
            } else if (strictStr == "false") {
                strict = false;
            } else {
                throw std::runtime_error("Expected strict: true or strict: false; found " + strictStr);
            }
        }

        consumeWord("eventProbabilityMultiplier:");
        eventProbabilityMultiplier = std::atof(nextWord("double").c_str());

        consumeWord("randomSeed:");
        randomSeed = std::atoll(nextWord("integer").c_str());

        Cut cut;

        auto finishCut = [&] {
            if (!cut.clauses.empty() || !cut.intHistograms.empty() || !cut.binHistograms.empty()) {
                if (cut.clauses.empty()) {
                    throw std::runtime_error("Cut didn't have any clauses");
                } else if (cut.intHistograms.empty() && cut.binHistograms.empty()) {
                    throw std::runtime_error("Cut didn't have any histograms");
                }
                cuts.push_back(cut);
                cut = {};
            }
        };

        // This looks horrible, but actually it is. Check whether there's any more to read after consuming whitespace.
        while (stream && stream >> std::ws && stream.peek() != std::istream::traits_type::eof()) {
            std::string directive = nextWord("variable name, new_cut, histogram_ints, or histogram");
            if (directive == "new_cut") {
                finishCut();
            } else if (directive == "histogram_ints:") {
                std::string varName = nextWord("variable name");
                size_t varIndex = format.var(varName);
                cut.intHistograms.emplace_back(varName, varIndex);
            } else if (directive == "histogram:") {
                std::string varName = nextWord("variable name");
                size_t varIndex = format.var(varName);
                double min = std::atof(nextWord("min value for " + varName).c_str());
                double max = std::atof(nextWord("max value for " + varName).c_str());
                size_t numBins = std::atoi(nextWord("number of bins for " + varName).c_str());
                cut.binHistograms.emplace_back(varName, varIndex, min, max, numBins);
            } else if (directive == "histogram_custom:") {
                std::string varName = nextWord("variable name");
                size_t varIndex = format.var(varName);
                cut.binHistograms.emplace_back(varName, varIndex, spaceSeparatedDoublesToEndOfLine());
            } else {
                size_t varIndex = format.var(directive);
                double min = std::atof(nextWord("min value for " + directive).c_str());
                double max = std::atof(nextWord("max value for " + directive).c_str());
                cut.clauses.push_back({varIndex, min, max});
            }
        }

        finishCut();
    }
};

CutJetsResult getCutJets(const Format& format, const char* filename, const GetCutJetsSpec& spec);
