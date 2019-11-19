#include <cassert>
#include <iostream>
#include <sstream>

#include "Histogram.h"
#include "get_cuts.h"

template<typename T>
static bool vectorsEqual(const std::vector<T>& v1, const std::vector<T>& v2) {
    return std::equal(v1.begin(), v1.end(), v2.begin(), v2.end());
}

template<typename Fn>
static void assertThrows(const std::string& expected, Fn&& fn) {
    try {
        fn();
    } catch (const std::exception& e) {
        if (e.what() != expected) {
            throw std::runtime_error("Expected '" + expected + "', but threw '" + e.what() + "'");
        }
        return;
    }
    throw std::runtime_error("No error was thrown, expected'" + expected + "'");
}

static void testParseSpec() {
    Format format({
        "VAR_0", "VAR_1", "VAR_2",
        "VAR_WEIGHT", "Z_PX", "GLUON_FLAG_1", // required for all Formats
    });

    assertThrows("Expected 'takeNum:' in spec", [&]{ GetCutJetsSpec(format, ""); });
    assertThrows("Expected integer in spec", [&]{ GetCutJetsSpec(format, "takeNum:"); });
    assertThrows("Expected 'skipNum:' in spec", [&]{ GetCutJetsSpec(format, "takeNum: 1"); });
    assertThrows("Expected 'strict:' in spec", [&]{ GetCutJetsSpec(format, "takeNum: 1\nskipNum: 2"); });
    assertThrows("Expected boolean in spec", [&]{ GetCutJetsSpec(format, "takeNum: 1\nskipNum: 2\nstrict:"); });
    assertThrows("Expected 'maxWeight:' in spec", [&]{ GetCutJetsSpec(format, "takeNum: 1\nskipNum: 2\nstrict: true"); });
    assertThrows("Expected double in spec", [&]{ GetCutJetsSpec(format, "takeNum: 1\nskipNum: 2\nstrict: true\nmaxWeight:"); });

    {
        GetCutJetsSpec spec(format, "takeNum: 1 \n skipNum: 2 \n strict: true \n maxWeight: 1.5");
        assert(spec.takeNum == 1);
        assert(spec.skipNum == 2);
        assert(spec.strict == true);
        assert(spec.maxWeight == 1.5);
    }
    {
        GetCutJetsSpec spec(format, "takeNum: 20 \n skipNum: 30 \n strict: false \n maxWeight: inf");
        assert(spec.takeNum == 20);
        assert(spec.skipNum == 30);
        assert(spec.strict == false);
        assert(spec.maxWeight == INFINITY);
    }

    assertThrows("unrecognized variable VAR_3", [&]{
        GetCutJetsSpec(format, R"(
            takeNum: 1
            skipNum: 2
            strict: true
            maxWeight: 1.5

            new_cut
            VAR_3 0.1 2.5
        )");
    });

    GetCutJetsSpec spec(format, std::istringstream(R"(
        takeNum: 1
        skipNum: 2
        strict: true
        maxWeight: 1.5

        new_cut
        VAR_1 -0.2 5.1e6
        histogram_ints: VAR_2
        histogram: VAR_0 0.2 0.5 10
        histogram_ints: VAR_1
        histogram: VAR_1 1 2 3

        new_cut
        VAR_1 1 0x123
        VAR_2 -100 100
        VAR_0 6 6.1
        histogram: VAR_1 0 10.5 5
    )"));

    assert(spec.takeNum == 1);
    assert(spec.skipNum == 2);
    assert(spec.strict == true);
    assert(spec.cuts.size() == 2);

    assert(vectorsEqual(spec.cuts[0].clauses, {
        CutClause{.varIndex = 1, .min = -0.2, .max = 5.1e6},
    }));
    assert(spec.cuts[0].intHistograms.size() == 2);
    assert(spec.cuts[0].intHistograms[0].varIndex == 2);
    assert(spec.cuts[0].intHistograms[1].varIndex == 1);

    assert(spec.cuts[0].binHistograms.size() == 2);

    assert(spec.cuts[0].binHistograms[0].varIndex == 0);
    assert(spec.cuts[0].binHistograms[0].binEndpoints.front() == 0.2);
    assert(spec.cuts[0].binHistograms[0].binEndpoints.back() == 0.5);
    assert(spec.cuts[0].binHistograms[0].binSums.size() == 10);

    assert(spec.cuts[0].binHistograms[1].varIndex == 1);
    assert(spec.cuts[0].binHistograms[1].binEndpoints.front() == 1);
    assert(spec.cuts[0].binHistograms[1].binEndpoints.back() == 2);
    assert(spec.cuts[0].binHistograms[1].binSums.size() == 3);

    assert(vectorsEqual(spec.cuts[1].clauses, {
        CutClause{.varIndex = 1, .min = 1, .max = 0x123},
        CutClause{.varIndex = 2, .min = -100, .max = 100},
        CutClause{.varIndex = 0, .min = 6, .max = 6.1},
    }));
    assert(spec.cuts[1].binHistograms.size() == 1);

    assert(spec.cuts[1].binHistograms[0].varIndex == 1);
    assert(spec.cuts[1].binHistograms[0].binEndpoints.front() == 0);
    assert(spec.cuts[1].binHistograms[0].binEndpoints.back() == 10.5);
    assert(spec.cuts[1].binHistograms[0].binSums.size() == 5);
}

static void testIntHistogram() {
    IntHistogram h("foo", 1);

    h.add(0.5, {0, 1, 2});
    h.add(0.1, {0, 1, 1});
    h.add(2.0, {6, 4, 6});
    h.finish();

    assert(h.totalWeight == 0.5 + 0.1 + 2.0);
    assert(h.totalErr == 0.25 + 0.01 + 4.0);
    assert(h.binSums[1] == (0.5 + 0.1) / h.totalWeight);
    assert(h.binSums[4] == (2.0) / h.totalWeight);

    assert(h.binErrs[1] == std::sqrt(0.25 + 0.01) / h.totalWeight);
    assert(h.binErrs[4] == std::sqrt(4.0) / h.totalWeight);
}

static void testBinHistogram() {
    BinHistogram h("foo", 1, 2.0, 5.0, 6);

    assert(h.binWidth == 0.5);
    assert(vectorsEqual(h.binEndpoints, {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0}));

    // before first bin
    h.add(1, {0, std::nextafter(2.0, 0)});
    // [2, 2.5)
    h.add(0.1, {0, 2.0});
    h.add(0.2, {0, std::nextafter(2.5, 0)});
    // [2.5, 3)
    h.add(0.3, {0, 2.5});
    h.add(0.4, {0, std::nextafter(3.0, 0)});
    // [3, 3.5)
    h.add(0.5, {0, 3.0});
    h.add(0.6, {0, std::nextafter(3.5, 0)});
    // [3.5, 4)
    h.add(0.7, {0, 3.5});
    h.add(0.8, {0, std::nextafter(4.0, 0)});
    // [4, 4.5)
    h.add(0.9, {0, 4.0});
    h.add(1.0, {0, std::nextafter(4.5, 0)});
    // [4.5, 5]
    h.add(1.1, {0, 4.5});
    h.add(1.2, {0, 5.0});
    // beyond last bin
    h.add(1.3, {0, std::nextafter(5.0, 8)});

    h.finish();

    assert(h.totalWeight == 0.1 + 0.2 + 0.3 + 0.4 + 0.5 + 0.6 + 0.7 + 0.8 + 0.9 + 1.0 + 1.1 + 1.2);
    assert(h.totalErr == 0.01 + 0.04 + 0.09 + 0.16 + 0.25 + 0.36 + 0.49 + 0.64 + 0.81 + 1.0 + 1.21 + 1.44);
    assert(vectorsEqual(h.binSums, {
        (0.1 + 0.2) / h.binWidth / h.totalWeight,
        (0.3 + 0.4) / h.binWidth / h.totalWeight,
        (0.5 + 0.6) / h.binWidth / h.totalWeight,
        (0.7 + 0.8) / h.binWidth / h.totalWeight,
        (0.9 + 1.0) / h.binWidth / h.totalWeight,
        (1.1 + 1.2) / h.binWidth / h.totalWeight,
    }));
    auto exact_hypot = [](double a, double b) {
        return std::sqrt(a * a + b * b);
    };
    assert(vectorsEqual(h.binErrs, {
        exact_hypot(0.1, 0.2) / h.binWidth / h.totalWeight,
        exact_hypot(0.3, 0.4) / h.binWidth / h.totalWeight,
        exact_hypot(0.5, 0.6) / h.binWidth / h.totalWeight,
        exact_hypot(0.7, 0.8) / h.binWidth / h.totalWeight,
        exact_hypot(0.9, 1.0) / h.binWidth / h.totalWeight,
        exact_hypot(1.1, 1.2) / h.binWidth / h.totalWeight,
    }));
}

void runTests() {
    testParseSpec();
    testIntHistogram();
    testBinHistogram();
    std::cout << "All tests passed!" << std::endl;
}
