#ifndef GENERATION_LOGGER_H
#define GENERATION_LOGGER_H

#include <json.hpp>
#include <string>
#include <tuple>
#include <vector>

class GenerationLogger {
  public:
    explicit GenerationLogger(const std::string &summary_file_path);

    void append_summary_from_directory(const std::string &directory, int generation);

  private:
    std::string summary_path;

    // Helper functions
    std::tuple<std::vector<std::pair<int, double>>, double, size_t, double, size_t, double, size_t>
    load_phenotype_data(const std::string &directory) const;

    std::tuple<double, double, int>
    compute_fitness_stats(const std::vector<std::pair<int, double>> &infos) const;

    void append_to_summary_log(const nlohmann::json &entry) const;
};

#endif
