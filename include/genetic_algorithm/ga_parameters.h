
#ifndef GA_PARAMETERS_H
#define GA_PARAMETERS_H

#include <cstddef>

struct GAParameters {
    size_t population_size = 10;    // Number of parameter sets per generation
    size_t num_generations = 500;    // Total number of generations to evolve
    double mutation_rate_init = 0.3; // Probability of mutation per gene
    double mutation_strength_init = 0.1;
    double crossover_rate = 0.4; // Probability of crossover between parents
    double crossover_alpha_init = 0.8;
    double elite_frac = 0.20;
    double random_frac = 0.10;

    unsigned int seed = 42; // RNG seed for reproducibility
};

#endif // GA_PARAMETERS_H