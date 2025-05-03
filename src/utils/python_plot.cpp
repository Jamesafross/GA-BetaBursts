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

    if (!json_data.contains("phenotypes") || !json_data["phenotypes"].is_array()) {
        throw std::runtime_error("Missing or invalid 'phenotypes' array in summary.");
    }

    const auto &phenotypes = json_data["phenotypes"];

    int best_id = -1;
    double best_fitness = -std::numeric_limits<double>::infinity();

    for (const auto &pheno : phenotypes) {
        if (!pheno.is_object())
            continue;

        int id = pheno.value("id", -1);
        double fitness = pheno.value("fitness", -std::numeric_limits<double>::infinity());

        if (id >= 0 && fitness > best_fitness) {
            best_fitness = fitness;
            best_id = id;
        }
    }

    if (best_id == -1) {
        throw std::runtime_error("No valid phenotypes found in summary.");
    }

    return best_id;
}