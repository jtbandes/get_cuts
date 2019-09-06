#include <cassert>
#include <iostream>

#include "Histogram.h"

static bool vectorsEqual(const std::vector<double>& v1, const std::vector<double>& v2) {
    return std::equal(v1.begin(), v1.end(), v2.begin(), v2.end());
}

static void testIntHistogram() {
    IntHistogram h(1);

    h.add(0.5, {0, 1, 2});
    h.add(0.1, {0, 1, 1});
    h.add(2.0, {6, 4, 6});

    assert(h.binSums[1] == 0.5 + 0.1);
    assert(h.binSums[4] == 2.0);
}

static void testBinHistogram() {
    BinHistogram h(1, 2.0, 7.0, 5);

    assert(vectorsEqual(h.binEndpoints, {2.0, 3.0, 4.0, 5.0, 6.0, 7.0}));

    // before first bin
    h.add(1, {0, std::nextafter(2.0, 0)});
    // [2, 3)
    h.add(0.1, {0, 2.0});
    h.add(0.2, {0, std::nextafter(3.0, 0)});
    // [3, 4)
    h.add(0.3, {0, 3.0});
    h.add(0.4, {0, std::nextafter(4.0, 0)});
    // [4, 5)
    h.add(0.5, {0, 4.0});
    h.add(0.6, {0, std::nextafter(5.0, 0)});
    // [5, 6)
    h.add(0.7, {0, 5.0});
    h.add(0.8, {0, std::nextafter(6.0, 0)});
    // [6, 7]
    h.add(0.9, {0, 6.0});
    h.add(1.0, {0, 7.0});
    // beyond last bin
    h.add(1.1, {0, std::nextafter(7.0, 8)});

    assert(vectorsEqual(h.binSums, {0.1 + 0.2, 0.3 + 0.4, 0.5 + 0.6, 0.7 + 0.8, 0.9 + 1.0}));
}

void runTests() {
    testIntHistogram();
    testBinHistogram();
    std::cout << "All tests passed!" << std::endl;
}
