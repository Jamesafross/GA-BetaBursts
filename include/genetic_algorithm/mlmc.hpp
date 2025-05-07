#pragma once

#include "neural_mass_model/meanfield.h"
#include "neural_mass_model/parameters.h"
#include <string>
#include <vector>

struct MLMCLevel {
    int level_index;
    double sim_duration_sec;
    size_t num_trials;
};

class MLMCController {
  public:
    MLMCController(const std::vector<MLMCLevel> &levels, const std::string &output_dir);

    // Run MLMC estimation for a given parameter set
    void run_estimate(const ModelParameters &params);

  private:
    std::vector<MLMCLevel> levels;
    std::string output_dir;

    // Run simulations at a specific level
    void run_level(const ModelParameters &params, const MLMCLevel &level);

    // Possibly add methods for estimating Q_l, computing differences, and aggregating results
};
