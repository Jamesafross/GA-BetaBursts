#include <complex> // ditto
#include <cstdint> // must come before Matlab headers
#include <iostream>

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

// naive ASCII->UTF16 (fine if your paths/messages are ASCII)
static std::u16string to_u16(const std::string &s) {
    return std::u16string(s.begin(), s.end());
}

void matlab_print_input(const std::string &message) {
    using namespace matlab::engine;
    using namespace matlab::data;

    try {
        auto matlab = startMATLAB();
        ArrayFactory factory;

        // Absolute path baked in via CMake: EXTERNAL_DIR="${PROJECT_SOURCE_DIR}/external"
        std::string matlabTestPath = std::string(EXTERNAL_DIR) + "/matlab_test";
        std::cout << "MATLAB_TEST_PATH = " << matlabTestPath << std::endl;

        // 1) Send path to MATLAB and add it (with subfolders)
        matlab->setVariable(u"p", factory.createCharArray(to_u16(matlabTestPath)));
        matlab->eval(u"addpath(genpath(p))");

        // 2) (Optional) Verify MATLAB sees the function
        matlab->eval(u"disp(which('print_input'))"); // should print a full path

        // 3) Call your MATLAB function
        matlab->setVariable(u"msg", factory.createCharArray(to_u16(message)));
        matlab->eval(u"print_input(msg)");

    } catch (const MATLABExecutionException &ex) {
        std::cerr << "MATLAB execution error: " << ex.what() << "\n";
    } catch (const EngineException &ex) {
        std::cerr << "MATLAB Engine error: " << ex.what() << "\n";
    }
}