#ifndef GENERATE_REAL_STATS_H
#define GENERATE_REAL_STATS_H
#include "thresholding/threshold_structs.h"
#include <string>

enum class AnalysisMethod { Threshold, HMM };

void generate_real_stats(ThresholdParameters params, const std::string &input_prefix,
                         const std::string &output_prefix, const std::string &final_merged_output,
                         int num_files = 3,
                         AnalysisMethod method = AnalysisMethod::Threshold // default only here
);

#endif