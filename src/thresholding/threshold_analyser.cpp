#include "thresholding/threshold_analyser.hpp"
#include "iir1/iir/Butterworth.h"
#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <numeric>
#include <stdexcept>

// === Constructor ===
SignalProcessor::SignalProcessor(const ThresholdParameters &params) : p(params) {
    if (p.percentile_threshold < 0.0 || p.percentile_threshold > 1.0) {
        throw std::invalid_argument("Percentile threshold must be between 0 and 1.");
    }
}

// === Downsampling ===
std::vector<double> SignalProcessor::downsample(const std::vector<double> &signal) const {
    std::vector<double> result;

    if (p.original_fs == p.downsampled_rate)
        return signal;
    if (p.original_fs < p.downsampled_rate) {
        throw std::invalid_argument("Downsample rate cannot exceed original sampling rate.");
    }

    size_t factor = static_cast<size_t>(p.original_fs / p.downsampled_rate);
    for (std::size_t i = 0; i < signal.size(); i += factor) {
        result.push_back(signal[i]);
    }

    return result;
}

void SignalProcessor::normalize_zscore(std::vector<double> &data) {
    if (data.empty())
        return;

    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double variance = sq_sum / data.size() - mean * mean;
    double stddev = std::sqrt(variance);

    if (stddev == 0.0)
        return; // Prevent division by zero

    for (auto &val : data) {
        val = (val - mean) / stddev;
    }
}

// === Bandpass Filtering ===
std::vector<double> SignalProcessor::bandpass_filter(const std::vector<double> &signal) const {
    double fs = p.original_fs;
    double low_cut = p.beta_low;
    double high_cut = p.beta_high;

    double center_freq = (low_cut + high_cut) / 2.0;
    double bandwidth = high_cut - low_cut;

    Iir::Butterworth::BandPass<4> bandpass;
    bandpass.setup(fs, center_freq, bandwidth);

    std::vector<double> filtered;
    filtered.reserve(signal.size());

    for (double x : signal) {
        filtered.push_back(bandpass.filter(x));
    }

    return filtered;
}

// === Lowpass Filtering ===
std::vector<double> SignalProcessor::lowpass_filter(const std::vector<double> &signal,
                                                    double cutoff) const {
    Iir::Butterworth::LowPass<4> lowpass;
    lowpass.setup(p.original_fs, cutoff);

    std::vector<double> filtered;
    filtered.reserve(signal.size());

    for (double x : signal) {
        filtered.push_back(lowpass.filter(x));
    }

    return filtered;
}

// === Envelope (Hilbert) ===
std::vector<double> SignalProcessor::compute_envelope(const std::vector<double> &signal) const {
    int N = signal.size();

    fftw_complex *freq = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    double *time = (double *)fftw_malloc(sizeof(double) * N);

    std::copy(signal.begin(), signal.end(), time);
    fftw_plan forward = fftw_plan_dft_r2c_1d(N, time, freq, FFTW_ESTIMATE);
    fftw_plan inverse = fftw_plan_dft_c2r_1d(N, freq, time, FFTW_ESTIMATE);

    fftw_execute(forward);

    // Apply Hilbert transform (analytic signal construction)
    for (int i = 0; i < N / 2 + 1; ++i) {
        if (i == 0 || (N % 2 == 0 && i == N / 2)) {
            freq[i][0] *= 1.0;
            freq[i][1] *= 1.0;
        } else {
            freq[i][0] *= 2.0;
            freq[i][1] *= 2.0;
        }
    }

    fftw_execute(inverse);

    std::vector<double> envelope(N);
    for (int i = 0; i < N; ++i) {
        double real = signal[i];
        double imag = time[i] / N;
        envelope[i] = std::sqrt(real * real + imag * imag);
    }

    fftw_destroy_plan(forward);
    fftw_destroy_plan(inverse);
    fftw_free(freq);
    fftw_free(time);

    return envelope;
}

// === Burst Detection ===
std::vector<std::pair<int, int>>
SignalProcessor::detect_bursts(const std::vector<double> &envelope) const {
    std::vector<std::pair<int, int>> bursts;
    bool in_burst = false;
    std::size_t burst_start = 0;

    double threshold;

    if (p.threshold_type == ThresholdType::Percentile) {
        std::vector<double> env_copy = envelope;
        size_t idx = static_cast<size_t>(p.percentile_threshold * env_copy.size());
        if (idx >= env_copy.size())
            idx = env_copy.size() - 1;
        std::nth_element(env_copy.begin(), env_copy.begin() + idx, env_copy.end());
        threshold = env_copy[idx];
    } else { // ZScore
        double mean = std::accumulate(envelope.begin(), envelope.end(), 0.0) / envelope.size();
        double sq_sum = std::inner_product(envelope.begin(), envelope.end(), envelope.begin(), 0.0);
        double stddev = std::sqrt(sq_sum / envelope.size() - mean * mean);
        threshold = mean + p.zscore_threshold * stddev;
    }

    for (std::size_t i = 0; i < envelope.size(); ++i) {
        if (!in_burst && envelope[i] > threshold) {
            in_burst = true;
            burst_start = i;
        } else if (in_burst && envelope[i] <= threshold) {
            in_burst = false;
            bursts.emplace_back(burst_start, i - 1);
        }
    }

    if (in_burst) {
        bursts.emplace_back(burst_start, envelope.size() - 1);
    }

    return bursts;
}

// === Burst Mask ===
std::vector<int>
SignalProcessor::burst_mask_from_intervals(const std::vector<std::pair<int, int>> &bursts,
                                           std::size_t length) const {
    std::vector<int> mask(length, 0);
    for (const auto &[start, end] : bursts) {
        std::size_t safe_start = static_cast<std::size_t>(start);
        std::size_t safe_end = static_cast<std::size_t>(end);
        for (std::size_t i = safe_start; i <= safe_end && i < length; ++i) {
            mask[i] = 1;
        }
    }
    return mask;
}

// === Burst Stats ===
BurstStats SignalProcessor::analyze_current_trace(const std::vector<double> &signal) const {
    auto filtered = bandpass_filter(signal);
    auto smoothed = lowpass_filter(filtered);
    auto downsampled = downsample(smoothed);
    auto envelope = compute_envelope(downsampled);
    normalize_zscore(envelope);
    auto raw_bursts = detect_bursts(envelope);

    std::vector<std::pair<int, int>> filtered_bursts;
    for (const auto &[start, end] : raw_bursts) {
        double dur = (end - start + 1) / p.downsampled_rate;
        if (dur >= p.min_burst_duration_s)
            filtered_bursts.emplace_back(start, end);
    }

    std::vector<std::pair<int, int>> bursts;
    int merge_gap = static_cast<int>(p.merge_gap_s * p.downsampled_rate);
    for (const auto &burst : filtered_bursts) {
        if (bursts.empty()) {
            bursts.push_back(burst);
        } else {
            auto &last = bursts.back();
            if (burst.first - last.second <= merge_gap) {
                last.second = std::max(last.second, burst.second);
            } else {
                bursts.push_back(burst);
            }
        }
    }

    double total_dur = 0.0, total_amp = 0.0;
    for (const auto &[start, end] : bursts) {
        double dur = (end - start + 1) / p.downsampled_rate;
        double amp = *std::max_element(envelope.begin() + start, envelope.begin() + end + 1);
        total_dur += dur;
        total_amp += amp;
    }

    int n_bursts = bursts.size();
    double total_time = envelope.size() / p.downsampled_rate;
    double burst_rate = n_bursts / total_time;
    double mean_dur = n_bursts > 0 ? total_dur / n_bursts : 0.0;
    double mean_amp = n_bursts > 0 ? total_amp / n_bursts : 0.0;

    return {burst_rate, mean_dur, mean_amp};
}
