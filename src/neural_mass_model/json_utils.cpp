#include "neural_mass_model/json_utils.hpp"
#include "genetic_algorithm/fitness.hpp"
#include "json.hpp"
#include "neural_mass_model/parameters.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using json = nlohmann::json;

// Convert ModelParameters to JSON
void to_json(json &j, const ModelParameters &p) {
    j = json{{"alphaee", p.alphaee},
             {"alphaei", p.alphaei},
             {"alphaie", p.alphaie},
             {"alphaii", p.alphaii},
             {"taue", p.taue},
             {"taui", p.taui},
             {"deltae", p.deltae},
             {"deltai", p.deltai},
             {"eta0e", p.eta0e},
             {"eta0i", p.eta0i},
             {"kappa_see", p.kappa_see},
             {"kappa_sei", p.kappa_sei},
             {"kappa_sie", p.kappa_sie},
             {"kappa_sii", p.kappa_sii},
             {"kappa_gapee", p.kappa_gapee},
             {"kappa_gapei", p.kappa_gapei},
             {"kappa_gapii", p.kappa_gapii},
             {"vsyn_ee", p.vsyn_ee},
             {"vsyn_ei", p.vsyn_ei},
             {"vsyn_ie", p.vsyn_ie},
             {"vsyn_ii", p.vsyn_ii},
             {"tauxe", p.tauxe},
             {"tauxi", p.tauxi},
             {"sigma_se", p.sigma_se},
             {"sigma_si", p.sigma_si}};
}

// Optional: convert JSON to ModelParameters
void from_json(const json &j, ModelParameters &p) {
    j.at("alphaee").get_to(p.alphaee);
    j.at("alphaei").get_to(p.alphaei);
    j.at("alphaie").get_to(p.alphaie);
    j.at("alphaii").get_to(p.alphaii);
    j.at("taue").get_to(p.taue);
    j.at("taui").get_to(p.taui);
    j.at("deltae").get_to(p.deltae);
    j.at("deltai").get_to(p.deltai);
    j.at("eta0e").get_to(p.eta0e);
    j.at("eta0i").get_to(p.eta0i);
    j.at("kappa_see").get_to(p.kappa_see);
    j.at("kappa_sei").get_to(p.kappa_sei);
    j.at("kappa_sie").get_to(p.kappa_sie);
    j.at("kappa_sii").get_to(p.kappa_sii);
    j.at("kappa_gapee").get_to(p.kappa_gapee);
    j.at("kappa_gapei").get_to(p.kappa_gapei);
    j.at("kappa_gapii").get_to(p.kappa_gapii);
    j.at("vsyn_ee").get_to(p.vsyn_ee);
    j.at("vsyn_ei").get_to(p.vsyn_ei);
    j.at("vsyn_ie").get_to(p.vsyn_ie);
    j.at("vsyn_ii").get_to(p.vsyn_ii);
    j.at("tauxe").get_to(p.tauxe);
    j.at("tauxi").get_to(p.tauxi);
    j.at("sigma_se").get_to(p.sigma_se);
    j.at("sigma_si").get_to(p.sigma_si);
}

// Save ModelParameters to JSON file
void save_parameters_to_json(const ModelParameters &p, const std::string &output_dir,
                             const std::string &filename) {
    json j = p;
    std::filesystem::path full_path = std::filesystem::path(output_dir) / filename;

    std::ofstream file(full_path);
    if (file.is_open()) {
        file << j.dump(4); // pretty print
        file.close();
        std::cout << "Saved parameters to " << full_path << std::endl;
    } else {
        std::cerr << "Could not open " << full_path << " for writing." << std::endl;
    }
}

// Save parameters and fitness to JSON

void save_parameters_with_fitness(const ModelParameters &p, FitnessResult result,
                                  const std::string &output_dir, const std::string &filename) {
    // Create JSON object and add subfields
    json j;
    j["parameters"] = p; // relies on a working to_json(ModelParameters, json&) function
    j["fitness"] = result.total;
    j["emd_fitness"] = result.emd_fitness;
    j["stat_fitness"] = result.stat_fitness;
    j["ks_fitness"] = result.ks_fitness;
    j["timestamp"] = std::time(nullptr); // optional: Unix timestamp
    j["generation"] = 0;                 // optional: default to 0 unless provided elsewhere

    std::filesystem::create_directories(output_dir);
    std::filesystem::path full_path = std::filesystem::path(output_dir) / filename;

    std::ofstream file(full_path);
    if (file.is_open()) {
        file << j.dump(4); // pretty print
        file.close();
        std::cout << "Saved parameters and fitness to " << full_path << std::endl;
    } else {
        std::cerr << "Could not open " << full_path << " for writing." << std::endl;
    }
}
