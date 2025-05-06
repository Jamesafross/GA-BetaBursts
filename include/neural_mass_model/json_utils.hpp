#ifndef JSON_UTILS_HPP
#define JSON_UTILS_HPP

#include "json.hpp"
#include "neural_mass_model/parameters.h"
#include <string>

using json = nlohmann::json;

// Convert ModelParameters to JSON
void to_json(json &j, const ModelParameters &p);

// Optional: convert JSON to ModelParameters
void from_json(const json &j, ModelParameters &p);

// Save ModelParameters to JSON file
void save_parameters_to_json(const ModelParameters &p, const std::string &output_dir,
                             const std::string &filename);

// Save parameters and fitness to JSON
void save_parameters_with_fitness(const ModelParameters &p, double fitness, double emd, double ks,
                                  const std::string &output_dir, const std::string &filename);

#endif // JSON_UTILS_HPP
