#ifndef THRESHOLD_STRUCTS_H
#define THRESHOLD_STRUCTS_H

// Add this to the top or inside a relevant namespace
enum class ThresholdType { Percentile, ZScore };

// Update ThresholdParameters
struct ThresholdParameters {
    double original_fs;
    double downsampled_rate;
    double beta_low;
    double beta_high;
    double percentile_threshold = 0.95;
    double zscore_threshold = 2.0;
    ThresholdType threshold_type = ThresholdType::Percentile;
    double min_burst_duration_s = 0.05;
    double merge_gap_s = 0.03;
};

struct BurstStats {
    double burst_rate;
    double mean_duration;
    double mean_amplitude;
};

#endif
