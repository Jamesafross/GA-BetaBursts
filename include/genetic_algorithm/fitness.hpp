#pragma once

#include <vector>
#include <string>
#include <unordered_map>

class FitnessEvaluator {
public:
    FitnessEvaluator() = default;

    // Compute fitness from raw vectors
    double compute_fitness(const std::vector<double>& model_rate,
                           const std::vector<double>& model_duration,
                           const std::vector<double>& model_amplitude,
                           const std::vector<double>& real_rate,
                           const std::vector<double>& real_duration,
                           const std::vector<double>& real_amplitude) const;

    // Convenience method: compute fitness directly from two CSV files
    double compute_fitness_from_csv(const std::string& model_file,
                                    const std::string& real_file) const;

private:
    // Core distance function
    double wasserstein_distance(std::vector<double> x, std::vector<double> y) const;

    // CSV utilities
    std::unordered_map<std::string, std::vector<double>>
    load_columns_by_name(const std::string& filename) const;
};