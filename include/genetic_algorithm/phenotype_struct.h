#include "neural_mass_model/parameters.h"

struct Individual {
    ModelParameters parameters;
    double fitness;
    std::string source_path;
};