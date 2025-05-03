#pragma once

#include <vector>
#include <string>
#include "neural_mass_model/parameters.h"

struct Individual {
    ModelParameters parameters;
    double fitness;
    std::string source_path;
};

class Selection {
public:
    // Public interface
    std::vector<Individual> select_top_n_from_directory(const std::string& dir, size_t n) const;

    void save_selected_individuals(const std::vector<Individual>& selected,
        const std::string& output_dir) const;

private:
    // Internal helpers
    std::vector<Individual> load_generation(const std::string& dir) const;
    std::vector<Individual> select_top_n(const std::vector<Individual>& population, size_t n) const;
};