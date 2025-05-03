#include "utils/python_plot.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <json.hpp>
#include <string>

void run_violin_plot(size_t generation, int phenotype_id) {
    std::string script_path = std::string(PYTHON_PLOTTING_DIR) + "/plot_violin.py";
    std::string cmd = "python3 " + script_path + " " + std::to_string(generation) + " " +
                      std::to_string(phenotype_id);

    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Python plotting script failed with code " << ret << std::endl;
    }
}

int get_best_phenotype_id(const std::string &summary_path) {
    std::ifstream infile(summary_path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open summary.json: " + summary_path);
    }

    nlohmann::json json_data;
    infile >> json_data;

    if (!json_data.is_array() || json_data.empty()) {
        throw std::runtime_error("summary.json is not a non-empty array.");
    }

    const auto &latest_entry = json_data.back();

    if (!latest_entry.contains("best_phenotype") ||
        !latest_entry["best_phenotype"].is_number_integer()) {
        throw std::runtime_error("Missing or invalid 'best_phenotype' in summary.json.");
    }

    return latest_entry["best_phenotype"];
}