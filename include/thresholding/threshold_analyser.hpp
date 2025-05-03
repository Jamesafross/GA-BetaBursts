#pragma once

#include "./threshold_structs.h"
#include <cmath>
#include <utility>
#include <vector>

// Signal processing and burst detection utilities
class SignalProcessor {
  public:
    explicit SignalProcessor(const ThresholdParameters &params);

    std::vector<double> downsample(const std::vector<double> &signal) const;

    static void normalize_zscore(std::vector<double> &data);

    std::vector<double> bandpass_filter(const std::vector<double> &signal) const;

    std::vector<double> lowpass_filter(const std::vector<double> &signal,
                                       double cutoff = 50.0) const;

    std::vector<double> compute_envelope(const std::vector<double> &signal) const;

    std::vector<std::pair<int, int>> detect_bursts(const std::vector<double> &envelope) const;

    std::vector<int> burst_mask_from_intervals(const std::vector<std::pair<int, int>> &bursts,
                                               std::size_t length) const;

    BurstStats analyze_current_trace(const std::vector<double> &signal) const;

  private:
    ThresholdParameters p;
};
