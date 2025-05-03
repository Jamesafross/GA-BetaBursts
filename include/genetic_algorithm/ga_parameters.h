
#ifndef GA_PARAMETERS_H
#define GA_PARAMETERS_H

#include <cstddef>

struct GAParameters {
    size_t population_size = 100; // Number of parameter sets per generation
    size_t num_generations = 100; // Total number of generations to evolve
    double mutation_rate = 0.1;   // Probability of mutation per gene
    double crossover_rate = 0.8;  // Probability of crossover between parents
    double elite_fraction = 0.1;  // Fraction of top performers to keep unaltered

    unsigned int seed = 42; // RNG seed for reproducibility
};

#endif // GA_PARAMETERS_H