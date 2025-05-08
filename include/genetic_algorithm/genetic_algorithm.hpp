#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "neural_mass_model/parameters.h"
#include "phenotype_struct.h"
#include <random>
#include <string>
#include <vector>
class GeneticAlgorithm {
  public:
    GeneticAlgorithm(ParameterBounds bounds, size_t total_size, double mutation_rate_init,
                     double mutation_strength_init, size_t max_generations);

    // Create a new generation by loading from disk and evolving
    std::vector<Individual>
    generate_next_generation_from_directory(const std::string &population_dir, double elite_frac,
                                            double random_frac, size_t tournament_size,
                                            size_t generation);

    // Save generation to directory (as JSON files)
    void save_generation(const std::vector<Individual> &individuals,
                         const std::string &output_dir) const;

  private:
    // Core GA parameters
    ParameterBounds bounds;
    size_t total_size;
    double mutation_rate_init;
    double mutation_strength_init;
    size_t max_generations;

    // Evolution steps
    std::vector<Individual> load_generation(const std::string &dir) const;
    std::vector<Individual> select_top_n(const std::vector<Individual> &, size_t n) const;
    std::vector<Individual> tournament_selection(const std::vector<Individual> &, size_t n,
                                                 size_t k) const;
    std::vector<Individual> generate_random_individuals(size_t count, std::mt19937 &rng) const;

    // Variation
    std::pair<ModelParameters, ModelParameters> crossover(const ModelParameters &p1,
                                                          const ModelParameters &p2,
                                                          std::mt19937 &rng,
                                                          size_t generation) const;

    void mutate(ModelParameters &individual, std::mt19937 &rng, double mutation_rate,
                double mutation_strength) const;

    void clamp(ModelParameters &individual) const;

    std::vector<ModelParameters> generate_variants(const std::vector<Individual> &parents,
                                                   size_t n_offspring, size_t n_random,
                                                   std::mt19937 &rng, size_t generation) const;

    double get_mutation_rate(size_t generation) const;

    double get_mutation_strength(size_t generation) const;
};

#endif
