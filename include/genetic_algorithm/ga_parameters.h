
#ifndef GA_PARAMETERS_H
#define GA_PARAMETERS_H

#include <cstddef>

struct GAParameters {
    size_t population_size = 100;    // Number of parameter sets per generation
    size_t num_generations = 100;    // Total number of generations to evolve
    double mutation_rate_init = 0.5; // Probability of mutation per gene
    double mutation_strength_init = 1.0;
    double crossover_rate = 0.4; // Probability of crossover between parents
    double elite_fraction = 0.1; // Fraction of top performers to keep unaltered
    double crossover_alpha_init = 0.8;

    unsigned int seed = 42; // RNG seed for reproducibility
};

#endif // GA_PARAMETERS_H