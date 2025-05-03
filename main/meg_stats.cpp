#include "meg/generate_real_stats.hpp"
#include "paths.hpp"
#include "thresholding/threshold_structs.h"
#include <string>

int main() {

    const std::string extension = ".csv";

    ThresholdParameters params;
    params.original_fs = 600.0;
    params.downsampled_rate = 100.0;
    params.percentile_threshold = 0.75;
    params.min_burst_duration_s = 0.05;
    params.merge_gap_s = 0.03;
    params.beta_low = 13.0;
    params.beta_high = 30.0;

    generate_real_stats(params);
}
