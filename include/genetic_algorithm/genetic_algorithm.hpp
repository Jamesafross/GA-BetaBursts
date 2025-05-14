#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "genetic_algorithm/ga_parameters.h"
#include "neural_mass_model/parameters.h"
#include "phenotype_struct.h"
#include <random>
#include <string>
#include <vector>
class GeneticAlgorithm {
  public:
    GeneticAlgorithm(ParameterBounds bounds, GAParameters ga_params);

    // Create a new generation by loading from disk and evolving
    std::vector<Individual>
    generate_next_generation_from_directory(const std::string &population_dir, double elite_frac,
                                            double random_frac, size_t tournament_size,
                                            size_t generation);

    std::vector<Individual> generate_next_generation(std::vector<Individual> population,
                                                     double elite_frac, double random_frac,
                                                     size_t tournament_size, size_t generation);

    // Save generation to directory (as JSON files)
    void save_generation(const std::vector<Individual> &individuals,
                         const std::string &output_dir) const;

  private:
    // Core GA parameters
    ParameterBounds bounds;
    GAParameters ga_params;
    size_t total_size = ga_params.population_size;
    double mutation_rate_init = ga_params.mutation_rate_init;
    double mutation_strength_init = ga_params.mutation_strength_init;
    double crossover_alpha_init = ga_params.crossover_alpha_init;
    size_t max_generations = ga_params.num_generations;
    double current_mutation_strength;
    double current_mutation_rate;
    double crossover_alpha;
    double running_best_fitness = std::numeric_limits<double>::infinity();
    bool stagnation_flag = 0;
    double stagnation_threshold = 0.02;
    const double eta_c = 10.0;

    // Evolution steps
    std::vector<Individual> load_generation(const std::string &dir) const;
    std::vector<Individual> select_top_n(const std::vector<Individual> &, size_t n) const;
    Individual best_phenotype(const std::vector<Individual> &population) const;
    std::vector<Individual> tournament_selection(const std::vector<Individual> &, size_t n,
                                                 size_t k) const;
    std::vector<Individual> rank_selection(const std::vector<Individual> &population,
                                           size_t num_selected) const;
    std::vector<Individual> generate_random_individuals(size_t count, std::mt19937 &rng) const;

    // Variation
    std::pair<ModelParameters, ModelParameters> crossover(const ModelParameters &p1,
                                                          const ModelParameters &p2,
                                                          std::mt19937 &rng,
                                                          size_t generation) const;

    void mutate(ModelParameters &individual, std::mt19937 &rng, double mutation_rate,
                double mutation_strength) const;

    void update_adaptive_parameters(size_t generation, const Individual &best_phenotype);

    void clamp(ModelParameters &individual) const;
    void stagnation_detect(const Individual &best_phenotype);

    std::vector<ModelParameters> generate_variants(const std::vector<Individual> &parents,
                                                   size_t n_offspring, size_t n_random,
                                                   std::mt19937 &rng, size_t generation) const;
    void update_adaptive_parameters(size_t generation);
};

#endif
