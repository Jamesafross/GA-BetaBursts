#include "genetic_algorithm/fitness.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <json.hpp>
#include <numeric>
#include <sstream>
#include <stdexcept>

FitnessEvaluator::FitnessEvaluator() = default;

FitnessResult FitnessEvaluator::compute_fitness(
    const std::vector<double> &model_rate, const std::vector<double> &model_duration,
    const std::vector<double> &model_amplitude, const std::vector<double> &real_rate,
    const std::vector<double> &real_duration, const std::vector<double> &real_amplitude,
    double last_gen_mean_emd, double last_gen_mean_ks, double last_gen_mean_stat) const {

    // constexpr double eps = 1e-8;
    constexpr double min_ref = 1e-3; // Lower bound to avoid dividing by near-zero
    double eps = 1e-8;

    double weight_rate = 1.0 / (mean(real_rate) + eps);
    double weight_duration = 1.0 / (mean(real_duration) + eps);
    double weight_amplitude = 1.0 / (mean(real_amplitude) + eps);

    double raw_emd = wasserstein_distance(model_rate, real_rate) +
                     wasserstein_distance(model_duration, real_duration) +
                     wasserstein_distance(model_amplitude, real_amplitude);

    // --- Compute distances ---
    double emd_total = weight_rate * wasserstein_distance(model_rate, real_rate) +
                       weight_duration * wasserstein_distance(model_duration, real_duration) +
                       weight_amplitude * wasserstein_distance(model_amplitude, real_amplitude);

    double raw_stat = stat_diff(model_rate, real_rate) + stat_diff(model_duration, real_duration) +
                      stat_diff(model_amplitude, real_amplitude);

    double stat_total = weight_rate * stat_diff(model_rate, real_rate) +
                        weight_duration * stat_diff(model_duration, real_duration) +
                        weight_amplitude * stat_diff(model_amplitude, real_amplitude);

    double raw_ks = ks_distance(model_rate, real_rate) +
                    ks_distance(model_duration, real_duration) +
                    ks_distance(model_amplitude, real_amplitude);

    double ks_total = weight_rate * ks_distance(model_rate, real_rate) +
                      weight_duration * ks_distance(model_duration, real_duration) +
                      weight_duration * ks_distance(model_amplitude, real_amplitude);

    // --- Normalize relative to previous generation ---
    double norm_emd = emd_total / std::max(last_gen_mean_emd, min_ref);
    double norm_ks = ks_total / std::max(last_gen_mean_ks, min_ref);
    double norm_stat = stat_total / std::max(last_gen_mean_stat, min_ref);

    // --- Optional amplification ---
    double penalized_emd = std::pow(norm_emd, 1.5);
    double penalized_ks = std::pow(norm_ks, 1.2);

    double total_fitness = penalized_emd + penalized_ks + norm_stat;
    double raw_fitness = raw_emd + raw_stat + raw_ks;

    return FitnessResult{.total = total_fitness,
                         .emd_fitness = emd_total,
                         .ks_fitness = ks_total,
                         .stat_fitness = stat_total,
                         .raw_total = raw_fitness,
                         .raw_emd = raw_emd,
                         .raw_ks = raw_ks,
                         .raw_stat = raw_stat};
}

// === Distance calculation ===
double FitnessEvaluator::wasserstein_distance(std::vector<double> x, std::vector<double> y) const {
    if (x.empty() || y.empty())
        throw std::invalid_argument("Input distributions must not be empty.");

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    // Resample to the same size using ECDF quantiles
    size_t n = std::max(x.size(), y.size());
    std::vector<double> x_interp(n), y_interp(n);

    for (size_t i = 0; i < n; ++i) {
        double q = static_cast<double>(i) / (n - 1);
        x_interp[i] = interpolate_quantile(x, q);
        y_interp[i] = interpolate_quantile(y, q);
    }

    double dist = 0.0;
    for (size_t i = 0; i < n; ++i) {
        dist += std::abs(x_interp[i] - y_interp[i]);
    }

    return dist / n;
}

double FitnessEvaluator::ks_distance(std::vector<double> x, std::vector<double> y) const {
    if (x.empty() || y.empty())
        throw std::invalid_argument("Input distributions must not be empty.");

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    size_t i = 0, j = 0;
    double cdf_x = 0.0, cdf_y = 0.0;
    double step_x = 1.0 / x.size();
    double step_y = 1.0 / y.size();
    double max_diff = 0.0;

    while (i < x.size() && j < y.size()) {
        if (x[i] < y[j]) {
            cdf_x += step_x;
            ++i;
        } else if (y[j] < x[i]) {
            cdf_y += step_y;
            ++j;
        } else { // x[i] == y[j]
            cdf_x += step_x;
            cdf_y += step_y;
            ++i;
            ++j;
        }
        max_diff = std::max(max_diff, std::abs(cdf_x - cdf_y));
    }

    return max_diff;
}

double FitnessEvaluator::interpolate_quantile(const std::vector<double> &data, double q) const {
    if (data.empty())
        return 0.0;
    double pos = q * (data.size() - 1);
    size_t idx = static_cast<size_t>(pos);
    double frac = pos - idx;
    if (idx + 1 < data.size())
        return data[idx] * (1.0 - frac) + data[idx + 1] * frac;
    return data.back();
}

// === CSV Column Loader ===
std::unordered_map<std::string, std::vector<double>>
FitnessEvaluator::load_columns_by_name(const std::string &filename) const {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file: " + filename);

    std::unordered_map<std::string, std::vector<double>> columns;
    std::vector<std::string> headers;

    std::string line;
    if (!std::getline(file, line))
        throw std::runtime_error("Empty file: " + filename);
    std::istringstream header_stream(line);
    std::string header;
    while (std::getline(header_stream, header, ',')) {
        headers.push_back(header);
        columns[header] = {};
    }

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        size_t col = 0;
        while (std::getline(ss, token, ',') && col < headers.size()) {
            try {
                columns[headers[col]].push_back(std::stod(token));
            } catch (...) {
            }
            ++col;
        }
    }

    return columns;
}

// === Wrapper to compute fitness from two files ===
FitnessResult FitnessEvaluator::compute_fitness_from_csv(const std::string &model_file,
                                                         const std::string &real_file,
                                                         const std::string &summary_json) const {
    auto model = load_columns_by_name(model_file);
    auto real = load_columns_by_name(real_file);
    double last_gen_mean_emd = 1.0;  // default fallback
    double last_gen_mean_ks = 1.0;   // default fallback
    double last_gen_mean_stat = 1.0; // default fallback
    nlohmann::json summary;
    std::ifstream in(summary_json);
    if (in.is_open()) {
        try {
            in >> summary;
        } catch (const std::exception &e) {
            std::cerr << "Error parsing summary.json: " << e.what() << "\n";
        }
        in.close();
    } else {
        std::cerr << "Warning: summary.json not found. Using default EMD/KS weights.\n";
    }

    if (summary.is_array() && !summary.empty()) {
        const auto &last = summary.back();
        if (last.contains("mean_emd"))
            last_gen_mean_emd = last["mean_emd"].get<double>();
        if (last.contains("mean_ks"))
            last_gen_mean_ks = last["mean_ks"].get<double>();
        if (last.contains("mean_stat"))
            last_gen_mean_stat = last["mean_stat"].get<double>();
    }

    return compute_fitness(model["burst_rate"], model["mean_duration"], model["mean_amplitude"],
                           real["burst_rate"], real["mean_duration"], real["mean_amplitude"],
                           last_gen_mean_emd, last_gen_mean_ks, last_gen_mean_stat);
}

double FitnessEvaluator::mean(const std::vector<double> &v) const {
    if (v.empty())
        return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double FitnessEvaluator::stddev(const std::vector<double> &v) const {
    double m = mean(v);
    double accum = 0.0;
    for (double x : v)
        accum += (x - m) * (x - m);
    return std::sqrt(accum / v.size());
}

double FitnessEvaluator::skewness(const std::vector<double> &v) const {
    double m = mean(v), s = stddev(v);
    if (s == 0.0)
        return 0.0;
    double accum = 0.0;
    for (double x : v)
        accum += std::pow((x - m) / s, 3);
    return accum / v.size();
}

double FitnessEvaluator::kurtosis(const std::vector<double> &v) const {
    double m = mean(v), s = stddev(v);
    if (s == 0.0)
        return 0.0;
    double accum = 0.0;
    for (double x : v)
        accum += std::pow((x - m) / s, 4);
    return accum / v.size(); // Excess kurtosis
}

double FitnessEvaluator::stat_diff(const std::vector<double> &m,
                                   const std::vector<double> &r) const {
    return std::abs(mean(m) - mean(r)) + std::abs(stddev(m) - stddev(r));
}
