#include "HMM/run_hmm_analysis.hpp"
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include <filesystem>
#include <iostream>
#include <string>

void run_hmm_analysis(const std::string &config_path, const std::string &popcurrent_path,
                      const std::string &output_path) {
    namespace fs = std::filesystem;
    using namespace matlab::engine;

    try {
        // Start MATLAB Engine
        std::unique_ptr<MATLABEngine> matlab = startMATLAB();

        // Add MATLAB paths
        matlab->eval(u"addpath(genpath('../HMM-MAR'))");
        matlab->eval(u"addpath(genpath('../HMM-MAR/HMM-MAR'))");
        matlab->eval(u"addpath(genpath('../nutmeg'))");

        // Optional: ensure parallel pool is open
        matlab->eval(u"if isempty(gcp('nocreate')), parpool('Processes', 4); end");

        // Absolute file paths (to avoid relative path issues)
        std::string config_abs = fs::absolute(config_path);
        std::string popcurrent_abs = fs::absolute(popcurrent_path);
        std::string output_abs = fs::absolute(output_path);

        // Form MATLAB function call
        std::u16string cmd = u"run_hmm_on_population('" +
                             std::u16string(config_abs.begin(), config_abs.end()) + u"','" +
                             std::u16string(popcurrent_abs.begin(), popcurrent_abs.end()) + u"','" +
                             std::u16string(output_abs.begin(), output_abs.end()) + u"')";

        std::cout << "Running MATLAB command:\n"
                  << std::string(cmd.begin(), cmd.end()) << std::endl;

        matlab->eval(cmd);
    } catch (const matlab::engine::EngineException &ex) {
        std::cerr << "MATLAB Engine error: " << ex.what() << std::endl;
    }
}