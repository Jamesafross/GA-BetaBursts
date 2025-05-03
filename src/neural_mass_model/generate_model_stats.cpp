#include "thresholding/threshold_analyser.hpp"
#include "thresholding/utils.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <filesystem>

void generate_model_stats(const std::string& input_dir,
    const std::string& output_dir,
    const std::string& input_filename,
    const std::string& output_filename,
    const ThresholdParameters& params) {
    std::string input_path = input_dir + "/" + input_filename;
    std::string output_path = output_dir + "/" + output_filename;

    SignalProcessor processor(params);
    auto trials = load_region_signals(input_path);
    if (trials.empty()) {
    std::cerr << "Error: could not load simulation data from " << input_path << "\n";
    return;
    }

    std::ofstream out(output_path);
    if (!out.is_open()) {
    std::cerr << "Error: could not open output file " << output_path << "\n";
    return;
    }

    out << "trial,burst_rate,mean_duration,mean_amplitude\n";

    std::size_t total_trials = trials.size();
    for (std::size_t i = 0; i < total_trials; ++i) {
    auto stats = processor.analyze_current_trace(trials[i]);
    out << i << "," << stats.burst_rate << ","
    << stats.mean_duration << "," << stats.mean_amplitude << "\n";

    int width = 40;
    float progress = static_cast<float>(i + 1) / total_trials;
    int pos = static_cast<int>(width * progress);
    std::cout << "\rProcessing [" << std::string(pos, '=') << std::string(width - pos, ' ')
    << "] " << std::fixed << std::setprecision(1)
    << (progress * 100.0) << "%";
    std::cout.flush();
    }

    std::cout << "\nSaved model burst stats to: " << output_path << "\n";
    out.close();
}
