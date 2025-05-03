#include "genetic_algorithm/fitness.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

// === Main fitness method ===
double FitnessEvaluator::compute_fitness(const std::vector<double>& model_rate,
                                         const std::vector<double>& model_duration,
                                         const std::vector<double>& model_amplitude,
                                         const std::vector<double>& real_rate,
                                         const std::vector<double>& real_duration,
                                         const std::vector<double>& real_amplitude) const
{
    double max_rate     = *std::max_element(real_rate.begin(), real_rate.end());
    double max_duration = *std::max_element(real_duration.begin(), real_duration.end());
    double max_amplitude= *std::max_element(real_amplitude.begin(), real_amplitude.end());

    double weight_rate     = (max_rate > 0.0)     ? 1.0 / max_rate     : 1.0;
    double weight_duration = (max_duration > 0.0) ? 1.0 / max_duration : 1.0;
    double weight_amplitude= (max_amplitude > 0.0)? 1.0 / max_amplitude: 1.0;

    return weight_rate     * wasserstein_distance(model_rate, real_rate) +
           weight_duration * wasserstein_distance(model_duration, real_duration) +
           weight_amplitude* wasserstein_distance(model_amplitude, real_amplitude);
}

// === Distance calculation ===
double FitnessEvaluator::wasserstein_distance(std::vector<double> x, std::vector<double> y) const {
    if (x.empty() || y.empty()) throw std::invalid_argument("Input distributions must not be empty.");

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    size_t n = x.size(), m = y.size();
    if (n == m) {
        double dist = 0.0;
        for (size_t i = 0; i < n; ++i) {
            dist += std::abs(x[i] - y[i]);
        }
        return dist / n;
    }

    size_t i = 0, j = 0;
    double fx = 0.0, fy = 0.0;
    double last_x = 0.0, last_y = 0.0;
    double step_x = 1.0 / n, step_y = 1.0 / m, area = 0.0;

    while (i < n || j < m) {
        double next_x = (i < n) ? x[i] : x.back();
        double next_y = (j < m) ? y[j] : y.back();

        if ((i < n && (j == m || next_x <= next_y))) {
            fx += step_x;
            last_x = next_x;
            ++i;
        } else {
            fy += step_y;
            last_y = next_y;
            ++j;
        }

        area += std::abs(fx - fy) * std::abs(last_x - last_y);
    }

    return area;
}

// === CSV Column Loader ===
std::unordered_map<std::string, std::vector<double>>
FitnessEvaluator::load_columns_by_name(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);

    std::unordered_map<std::string, std::vector<double>> columns;
    std::vector<std::string> headers;

    std::string line;
    if (!std::getline(file, line)) throw std::runtime_error("Empty file: " + filename);
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
            } catch (...) {}
            ++col;
        }
    }

    return columns;
}

// === Wrapper to compute fitness from two files ===
double FitnessEvaluator::compute_fitness_from_csv(const std::string& model_file,
                                                  const std::string& real_file) const {
    auto model = load_columns_by_name(model_file);
    auto real  = load_columns_by_name(real_file);

    return compute_fitness(model["burst_rate"],     model["mean_duration"],     model["mean_amplitude"],
                           real["burst_rate"],      real["mean_duration"],      real["mean_amplitude"]);
}
