#ifndef MEANFIELD_H
#define MEANFIELD_H

#include "parameters.h"
#include <random>
#include <string>
#include <vector>

// Forward declaration
struct SolverParameters;

// Mean-field dynamics class
class meanfield {
    const ModelParameters &p;

  public:
    meanfield(const ModelParameters &p_);

    void operator()(const std::vector<double> &state, std::vector<double> &dxdt, double t) const;
    void apply_noise(std::vector<double> &state, double dt, std::mt19937 &rng) const;
};

// Run multiple trials of the neural mass model
void run_neural_mass(const ModelParameters &p, const SolverParameters &opts,
                     const std::string &output_dir, const std::string &filename);

// Save a matrix of simulation output to CSV
void save_matrix(const std::vector<std::vector<double>> &data, const std::string &filename);

#endif