#ifndef FITNESS_H
#define FITNESS_H

#include <string>
#include <unordered_map>
#include <vector>

struct FitnessResult {
    double total;
    double emd_fitness;
    double ks_fitness;
    double stat_fitness;
    double raw_total;
    double raw_emd;
    double raw_ks;
    double raw_stat;
};

class FitnessEvaluator {
  public:
    FitnessEvaluator(); // Just declare the default constructor

    // Compute fitness from two CSVs and a summary file
    FitnessResult compute_fitness_from_csv(const std::string &model_file,
                                           const std::string &real_file,
                                           const std::string &summary_json) const;

  private:
    // Distance metrics
    double wasserstein_distance(std::vector<double> x, std::vector<double> y) const;
    double ks_distance(std::vector<double> x, std::vector<double> y) const;

    // Summary statistics
    double mean(const std::vector<double> &v) const;
    double stddev(const std::vector<double> &v) const;
    double skewness(const std::vector<double> &v) const;
    double kurtosis(const std::vector<double> &v) const;
    double stat_diff(const std::vector<double> &m, const std::vector<double> &r) const;

    // Utilities
    double interpolate_quantile(const std::vector<double> &data, double q) const;
    std::unordered_map<std::string, std::vector<double>>
    load_columns_by_name(const std::string &filename) const;

    // Main internal fitness calculator
    FitnessResult compute_fitness(
        const std::vector<double> &model_rate, const std::vector<double> &model_duration,
        const std::vector<double> &model_amplitude, const std::vector<double> &real_rate,
        const std::vector<double> &real_duration, const std::vector<double> &real_amplitude,
        double last_gen_mean_emd, double last_gen_mean_ks, double last_gen_mean_stat) const;
};

#endif