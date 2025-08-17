#pragma once

#include <string>

/**
 * @brief Runs the HMM burst detection analysis using MATLAB.
 *
 * This function starts a MATLAB engine session and calls the
 * `run_hmm_on_population` MATLAB function, passing in paths
 * to the JSON config, population current `.mat` file, and output directory.
 *
 * @param config_path Path to the JSON configuration file.
 * @param popcurrent_path Path to the .mat file containing population current data.
 * @param output_path Directory where the output (e.g., HMMStats.mat) will be saved.
 */
void run_hmm_analysis(const std::string &config_path, const std::string &popcurrent_path,
                      const std::string &output_path);