#include "thresholding/threshold_analyser.hpp"
#include "thresholding/threshold_structs.h"
#include "thresholding/utils.hpp"
#include "meg/generate_real_stats.hpp"
#include "paths.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <iomanip>

const std::string input_prefix = DATA_PATH + "/meg/meg_data_";
const std::string output_prefix = OUTPUT_PATH + "/meg/meg_burst_stats_";
const std::string final_merged_output = OUTPUT_PATH + "/meg/meg_burst_stats_merged.csv";

void generate_real_stats(ThresholdParameters params) {
    const std::string extension = ".csv";
    const int num_files = 3;  // Adjust if needed


    SignalProcessor processor(params);

    // Process each MEG file
    for (int file_idx = 1; file_idx <= num_files; ++file_idx) {
        std::string input_file = input_prefix + std::to_string(file_idx) + extension;
        std::string output_file = output_prefix + std::to_string(file_idx) + extension;

        std::cout << "Loading " << input_file << "\n";
        auto meg_trials = load_region_signals(input_file);
        if (meg_trials.empty()) {
            std::cerr << "Error: failed to load " << input_file << "\n";
            continue;
        }

        std::ofstream out(output_file);
        out << "file,trial,burst_rate,mean_duration,mean_amplitude\n";

        std::size_t total_trials = meg_trials.size();
        for (std::size_t i = 0; i < total_trials; ++i) {
            auto stats = processor.analyze_current_trace(meg_trials[i]);
            out << file_idx << "," << i << "," << stats.burst_rate << ","
                << stats.mean_duration << "," << stats.mean_amplitude << "\n";

            // Progress bar
            int width = 50;
            float progress = static_cast<float>(i + 1) / total_trials;
            int pos = static_cast<int>(width * progress);

            std::cout << "\rFile " << file_idx << " progress [";
            for (int j = 0; j < width; ++j)
                std::cout << (j < pos ? "=" : (j == pos ? ">" : " "));
            std::cout << "] " << std::fixed << std::setprecision(1)
                      << (progress * 100.0) << "%";
            std::cout.flush();
        }

        std::cout << "\nSaved: " << output_file << "\n";
        out.close();
    }

    // Merge all output files
    std::ofstream merged(final_merged_output);
    merged << "file,trial,burst_rate,mean_duration,mean_amplitude\n";

    for (int file_idx = 1; file_idx <= num_files; ++file_idx) {
        std::string file = output_prefix + std::to_string(file_idx) + extension;
        std::ifstream in(file);
        if (!in.is_open()) {
            std::cerr << "Could not open " << file << " for merging.\n";
            continue;
        }

        std::string line;
        std::getline(in, line); // skip header
        while (std::getline(in, line)) {
            merged << line << "\n";
        }

        in.close();
    }

    merged.close();
    std::cout << "Merged output saved to: " << final_merged_output << "\n";
}
