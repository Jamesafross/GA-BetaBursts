#include "genetic_algorithm/fitness.hpp"
#include "genetic_algorithm/ga_parameters.h"

#include "meg/generate_real_stats.hpp"

#include "paths.hpp"
#include "thresholding/threshold_structs.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <json.hpp>

#include <string>

void ensure_output_directories(const std::string &base_dir, const std::string &stats_dir,
                               const std::string &simulated_dir,
                               const std::string &parameters_dir) {
    std::filesystem::create_directories(base_dir);
    std::filesystem::create_directories(stats_dir);
    std::filesystem::create_directories(simulated_dir);
    std::filesystem::create_directories(parameters_dir);
}

void cleanup_simulated_data(const std::string &simulated_dir) {
    try {
        if (std::filesystem::exists(simulated_dir)) {
            std::filesystem::remove_all(simulated_dir);
            std::cout << "Deleted simulated data directory: " << simulated_dir << "\n";
        }
    } catch (const std::filesystem::filesystem_error &e) {
        std::cerr << "Error deleting simulated data: " << e.what() << "\n";
    }
}

void clear_summary_file(const std::string &path) {
    std::ofstream out(path, std::ios::trunc); // truncate clears the file
    if (out.is_open()) {
        nlohmann::json empty_array = nlohmann::json::array();
        out << empty_array.dump(4); // pretty print
        out.close();
    }
}

const std::string meg_output_dir = OUTPUT_PATH + "/meg";

int main() {

    ThresholdParameters th_params;

    std::filesystem::create_directories(OUTPUT_PATH + "/meg");
    generate_real_stats(th_params, DATA_PATH + "/meg/meg_data_",
                        OUTPUT_PATH + "/meg/meg_burst_stats_",
                        OUTPUT_PATH + "/meg/meg_burst_stats_merged.csv", 3);

    return 0;
}
