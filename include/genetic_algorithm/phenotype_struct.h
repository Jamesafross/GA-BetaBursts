#ifndef PHENOTYPE_STRUCT_H
#define PHENOTYPE_STRUCT_H
#include "neural_mass_model/parameters.h"

struct Individual {
    ModelParameters parameters;
    double fitness;
    std::string source_path;
};
#endif