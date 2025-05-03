#ifndef THRESHOLD_STRUCTS_H
#define THRESHOLD_STRUCTS_H

#include <cstddef>


struct ThresholdParameters {
    //DEFAULTS
    double original_fs = 200.0;            // Sampling rate of input data (Hz)
    double downsampled_rate = 100.0;       // Downsample rate (Hz) - can be same as original
    double percentile_threshold = 0.75;    // e.g. 75th percentile
    double min_burst_duration_s = 0.05;    // Minimum duration in seconds
    double merge_gap_s = 0.03;             // Max gap for merging bursts (s)
    double beta_low = 13.0;
    double beta_high = 30.0;
};

struct BurstStats {
    double burst_rate;
    double mean_duration;
    double mean_amplitude;
};

#endif 
