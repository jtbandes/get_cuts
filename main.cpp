#include <cinttypes>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "get_cuts.h"
#include "test.h"

static Format NewerFormat({
    "VAR_NUM", "VAR_WEIGHT", "VAR_PT", "VAR_PSEUDORAP", "VAR_PHI", "VAR_M", "VAR_CONST", "VAR_RAP", "Z_PX", "Z_PY", "Z_PZ", "Z_E", "Z_RAP", "GLUON_FLAG_1", "GLUON_FLAG_2", "VAR_CONST_SD",
});

static Format NewFormat({
    "VAR_NUM", "VAR_WEIGHT", "VAR_PT", "VAR_PSEUDORAP", "VAR_PHI", "VAR_M", "VAR_CONST", "VAR_RAP", "Z_PX", "Z_PY", "Z_PZ", "Z_E", "Z_RAP", "GLUON_FLAG_1", "GLUON_FLAG_2", "VAR_C11", "VAR_C10", "VAR_ANG1", "VAR_ANG05", "VAR_CONST_SD", "VAR_C11_SD", "VAR_C10_SD", "VAR_ANG1_SD",
});

int main(int argc, char** argv) {
    std::vector<std::string> args(argv+1, argv+argc);
    if (args.size() > 0 && args[0] == "--test") {
        runTests();
        return 0;
    }

    if (args.size() != 2) {
        std::cerr << std::string(R"(
Usage: get_cuts [--new|--newer] input.txt < spec.txt
Spec file format:
  takeNum: 2
  skipNum: 2
  strict: true

  new_cut
  VAR_1 min1 max1
  VAR_2 min2 max2
  histogram_ints: VAR_3
  histogram: VAR_4 0.2 0.5 20

  new_cut
  VAR_1 min1 max1
  VAR_2 min2 max2
  VAR_3 min3 max3
  histogram: VAR_8 0 10.5 4
)").substr(1) << std::endl;
        return 1;
    }

    const Format* format;
    {
        const auto& formatArg = args[0];
        if (formatArg == "--new") {
            format = &NewFormat;
        } else if (formatArg == "--newer") {
            format = &NewerFormat;
        } else {
            throw std::runtime_error("Expected --new or --newer");
        }
    }

    const auto& filename = args[1];

    GetCutJetsSpec spec(*format, std::cin);
    CutJetsResult result = getCutJets(*format, filename.c_str(), spec);

    std::printf("%lf\n", result.csOnW);
    for (const auto& cutResult : result.cutResults) {
        std::printf("----------------------------------\n");
        for (const auto& hist : cutResult.intHistograms) {
            std::printf("%s\n", hist.varName.c_str());
            for (const auto& [k, v] : hist.binSums) {
                std::printf("%" PRIdMAX ": %lf err %lf\n", k, v, hist.binErrs.at(k));
            }
        }
        for (const auto& hist : cutResult.binHistograms) {
            std::printf("%s\n", hist.varName.c_str());
            for (size_t i = 0; i < hist.binSums.size(); i++) {
                std::printf("%lf-%lf: %lf err %lf\n",
                    hist.binEndpoints[i], hist.binEndpoints[i+1], hist.binSums[i], hist.binErrs[i]);
            }
        }
    }

    // Naive output as CSV
    // std::printf("%lf\n", result.csOnW);
    // for (size_t i = 0; i < result.jetsList.size(); i++) {
    //     // std::printf("# Cut %zu: %zu jets taken\n", i, result.jetsList[i].size());
    //     for (const auto& jet : result.jetsList[i]) {
    //         for (double val : jet) {
    //             std::printf("%lg,", val);
    //         }
    //         std::printf("\n");
    //     }
    //     std::printf("\n\n");
    // }

    return 0;
}
