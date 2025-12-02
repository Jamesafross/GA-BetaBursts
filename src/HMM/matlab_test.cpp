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
        std::vector<std::u16string> options{u"-nojvm", u"-nodisplay", u"-nosplash"};
        auto matlab = startMATLAB(options);

        ArrayFactory factory;

        std::string matlabTestPath = std::string(EXTERNAL_DIR) + "/matlab_test";
        std::cout << "MATLAB_TEST_PATH = " << matlabTestPath << std::endl;

        matlab->setVariable(u"p", factory.createCharArray(to_u16(matlabTestPath)));
        matlab->eval(u"addpath(genpath(p))");

        matlab->eval(u"disp(which('print_input'))");

        matlab->setVariable(u"msg", factory.createCharArray(to_u16(message)));
        matlab->eval(u"print_input(msg)");

    } catch (const MATLABExecutionException &ex) {
        std::cerr << "MATLAB execution error: " << ex.what() << "\n";
    } catch (const EngineException &ex) {
        std::cerr << "MATLAB Engine error: " << ex.what() << "\n";
    }
}
