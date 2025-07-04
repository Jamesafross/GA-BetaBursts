#ifndef PARAMETER_UTILS_H
#define PARAMETER_UTILS_H

#include "genetic_algorithm/phenotype_struct.h"
#include "parameters.h"
#include <random>
#include <vector>

// Decode a normalized genome vector (values in [0,1]) into a ModelParameters struct
ModelParameters decode_parameters(const std::vector<double> &genome, const ParameterBounds &bounds);

// Generate one random parameter set by sampling a random genome and decoding it
ModelParameters sample_random_parameters(const ParameterBounds &bounds, std::mt19937 &rng);

// Generate a full population of parameter sets
std::vector<Individual> generate_parameter_population(size_t n, const ParameterBounds &bounds,
                                                      std::mt19937 &rng);

// Print a parameter set to std::cout
void print_parameters(const ModelParameters &params);

#endif // PARAMETER_UTILS_H
