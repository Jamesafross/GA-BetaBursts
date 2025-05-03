#include "thresholding/utils.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


std::vector<std::vector<double>> load_region_signals(const std::string& filename) {
    std::vector<std::vector<double>> region_signals;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return region_signals;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> region;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                region.push_back(std::stod(value));
            } catch (const std::invalid_argument&) {
                std::cerr << "Warning: invalid value '" << value << "' in line: " << line << std::endl;
            }
        }

        if (!region.empty()) {
            region_signals.push_back(region);
        }
    }

    return region_signals;
}