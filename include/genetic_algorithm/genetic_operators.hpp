#pragma once

#include "genetic_algorithm/selection.hpp"
#include "neural_mass_model/parameters.h"
#include <random>
#include <utility>

class GeneticOperators {
  public:
    GeneticOperators(double mutation_rate, double mutation_strength, size_t max_generations)
        : mutation_rate(mutation_rate), mutation_strength(mutation_strength),
          max_generations(max_generations) {}

    std::vector<ModelParameters> generate_next_population(const std::vector<Individual> &selected,
                                                          size_t total_size, size_t n_random,
                                                          std::mt19937 &rng,
                                                          const ParameterBounds &bounds,
                                                          size_t generation) const;

  private:
    double mutation_rate;
    double mutation_strength;
    size_t max_generations;

    void clamp_parameters_to_bounds(ModelParameters &p, const ParameterBounds &b) const;

    std::pair<ModelParameters, ModelParameters> crossover(const ModelParameters &parent1,
                                                          const ModelParameters &parent2,
                                                          std::mt19937 &rng,
                                                          size_t generation) const;

    void mutate(ModelParameters &individual, std::mt19937 &rng) const;
};
