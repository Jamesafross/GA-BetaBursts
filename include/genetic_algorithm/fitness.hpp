#pragma once

#include <string>
#include <unordered_map>
#include <vector>

struct FitnessResult {
    double total;
    double emd_fitness;
    double ks_fitness;
    double stat_fitness;
};

class FitnessEvaluator {

  public:
    explicit FitnessEvaluator();
    // Compute fitness from raw vectors

    // Convenience method: compute fitness directly from two CSV files
    FitnessResult compute_fitness_from_csv(const std::string &model_file,
                                           const std::string &real_file,
                                           const std::string &summary_json) const;

  private:
    double stats_weight;
    // Core distance function
    double wasserstein_distance(std::vector<double> x, std::vector<double> y) const;

    double ks_distance(std::vector<double> x, std::vector<double> y) const;

    double mean(const std::vector<double> &v) const;

    double stddev(const std::vector<double> &v) const;

    double skewness(const std::vector<double> &v) const;

    double kurtosis(const std::vector<double> &v) const;

    double stat_diff(const std::vector<double> &m, const std::vector<double> &r) const;

    double interpolate_quantile(const std::vector<double> &data, double q) const;

    FitnessResult compute_fitness(
        const std::vector<double> &model_rate, const std::vector<double> &model_duration,
        const std::vector<double> &model_amplitude, const std::vector<double> &real_rate,
        const std::vector<double> &real_duration, const std::vector<double> &real_amplitude,
        double last_gen_mean_emd, double last_gen_mean_ks, double last_gen_mean_stat) const;

    // CSV utilities
    std::unordered_map<std::string, std::vector<double>>
    load_columns_by_name(const std::string &filename) const;
};
