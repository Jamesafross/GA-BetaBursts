#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

// ASCII-only helper; swap for a proper UTF-8→UTF-16 if needed
static std::u16string to_u16(const std::string &s) {
    return std::u16string(s.begin(), s.end());
}

void hmm_matlab_run(const std::string &popcurrent_path, const std::string &output_path,
                    double sampling_freq) {
    using namespace matlab::engine;
    using namespace matlab::data;

    try {
        // 1) Start MATLAB headless (faster, more stable on Manjaro)
        std::vector<std::u16string> opts{u"-nojvm", u"-nodisplay", u"-nosplash"};
        // setenv("LD_PRELOAD", "/usr/lib/libstdc++.so.6", /*overwrite=*/1);
        auto matlab = startMATLAB(opts);
        ArrayFactory factory;

        // 2) Add your code paths (HMM, HMM-MAR, nutmeg)
        const std::string matlabTestPath = std::string(EXTERNAL_DIR) + "/HMM/HMM";
        const std::string hmmPath = std::string(EXTERNAL_DIR) + "/HMM/HMM-MAR/";
        const std::string nutmegPath = std::string(EXTERNAL_DIR) + "/nutmeg/";

        matlab->setVariable(u"mtp", factory.createCharArray(to_u16(matlabTestPath)));
        matlab->setVariable(u"hp", factory.createCharArray(to_u16(hmmPath)));
        matlab->setVariable(u"nmp", factory.createCharArray(to_u16(nutmegPath)));

        matlab->eval(u"addpath(genpath(mtp))");
        matlab->eval(u"addpath(genpath(append(hp)))");
        matlab->eval(u"addpath(genpath(append(nmp)))");

        // Optional: show MATLAB finds your function
        matlab->eval(u"disp(which('run_hmm_on_population'))");

        // 3) Build arguments for run_hmm_on_population(...)
        std::vector<Array> args;
        args.emplace_back(factory.createCharArray(to_u16(popcurrent_path)));
        args.emplace_back(factory.createCharArray(to_u16(output_path)));

        args.emplace_back(factory.createScalar<double>(sampling_freq));

        // 4) Call your MATLAB function (no outputs expected)
        matlab->feval(u"run_hmm_on_population", 0, args);

        // If you want to ensure the pool is closed after (optional):
        // matlab->eval(u"pool = gcp('nocreate'); if ~isempty(pool), delete(pool); end");

    } catch (const matlab::engine::MATLABExecutionException &ex) {
        std::cerr << "MATLAB execution error: " << ex.what() << "\n";
    } catch (const matlab::engine::EngineException &ex) {
        std::cerr << "MATLAB Engine error: " << ex.what() << "\n";
    } catch (const std::exception &ex) {
        std::cerr << "Std error: " << ex.what() << "\n";
    }
}
