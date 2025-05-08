#include "genetic_algorithm/genetic_algorithm.hpp"
#include "genetic_algorithm/ga_parameters.h"
#include "neural_mass_model/json_utils.hpp"
#include "neural_mass_model/parameter_utils.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

#define FOR_EACH_PARAM(F)                                                                          \
    F(alphaee)                                                                                     \
    F(alphaei)                                                                                     \
    F(alphaie)                                                                                     \
    F(alphaii)                                                                                     \
    F(kappa_see)                                                                                   \
    F(kappa_sei)                                                                                   \
    F(kappa_sie)                                                                                   \
    F(kappa_sii)                                                                                   \
    F(vsyn_ee)                                                                                     \
    F(vsyn_ei)                                                                                     \
    F(vsyn_ie)                                                                                     \
    F(vsyn_ii)                                                                                     \
    F(tauxe)                                                                                       \
    F(tauxi)                                                                                       \
    F(taue)                                                                                        \
    F(taui)                                                                                        \
    F(deltae)                                                                                      \
    F(deltai) F(sigma_se) F(sigma_si) F(kappa_gapee) F(kappa_gapei) F(kappa_gapii) F(eta0e) F(eta0i)

// === Constructor ===
GeneticAlgorithm::GeneticAlgorithm(ParameterBounds bounds, GAParameters ga_params)
    : bounds(bounds), ga_params(ga_params) {}
// === Load individuals from a directory ===
std::vector<Individual> GeneticAlgorithm::load_generation(const std::string &dir) const {
    std::vector<Individual> generation;
    for (const auto &entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() != ".json")
            continue;

        std::ifstream in(entry.path());
        if (!in.is_open())
            continue;

        try {
            nlohmann::json j;
            in >> j;

            ModelParameters p = j.at("parameters").get<ModelParameters>();
            double fitness = j.value("fitness", -1.0);
            generation.push_back({p, fitness, entry.path().string()});
        } catch (const std::exception &e) {
            std::cerr << "Error parsing " << entry.path() << ": " << e.what() << "\n";
        }
    }
    return generation;
}

// === Save individuals to a directory ===
void GeneticAlgorithm::save_generation(const std::vector<Individual> &individuals,
                                       const std::string &output_dir) const {
    std::filesystem::create_directories(output_dir);
    for (const auto &ind : individuals) {
        nlohmann::json j;
        j["fitness"] = ind.fitness;
        j["parameters"] = ind.parameters;

        std::string filename = std::filesystem::path(ind.source_path).filename().string();
        if (filename.empty())
            filename = "individual_" + std::to_string(&ind - &individuals[0]) + ".json";

        std::ofstream out(std::filesystem::path(output_dir) / filename);
        if (out.is_open())
            out << j.dump(4);
    }
}

// === Select top N by fitness ===
std::vector<Individual> GeneticAlgorithm::select_top_n(const std::vector<Individual> &population,
                                                       size_t n) const {
    std::vector<Individual> sorted = population;
    std::sort(sorted.begin(), sorted.end(),
              [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });
    if (n > sorted.size())
        n = sorted.size();
    return {sorted.begin(), sorted.begin() + n};
}

Individual GeneticAlgorithm::best_phenotype(const std::vector<Individual> &population) const {
    std::vector<Individual> sorted = population;
    std::sort(sorted.begin(), sorted.end(),
              [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });
    return *sorted.begin(); // Dereference to get the actual object
}

// === Tournament selection ===
std::vector<Individual>
GeneticAlgorithm::tournament_selection(const std::vector<Individual> &population,
                                       size_t num_selected, size_t tournament_size) const {
    std::vector<Individual> selected;
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<> dis(0, population.size() - 1);

    for (size_t i = 0; i < num_selected; ++i) {
        std::vector<Individual> tournament;
        for (size_t j = 0; j < tournament_size; ++j)
            tournament.push_back(population[dis(rng)]);

        auto best = *std::min_element(
            tournament.begin(), tournament.end(),
            [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });
        selected.push_back(best);
    }

    return selected;
}

std::vector<Individual> GeneticAlgorithm::rank_selection(const std::vector<Individual> &population,
                                                         size_t num_selected) const {
    std::vector<Individual> sorted = population;
    std::sort(sorted.begin(), sorted.end(),
              [](const Individual &a, const Individual &b) { return a.fitness < b.fitness; });

    // Assign selection probabilities based on rank (linear scheme)
    const size_t N = sorted.size();
    std::vector<double> selection_probs(N);
    double sum_ranks = static_cast<double>(N * (N + 1)) / 2.0;

    for (size_t i = 0; i < N; ++i) {
        selection_probs[i] = static_cast<double>(N - i) / sum_ranks; // rank 1 gets highest prob
    }

    // Sample individuals based on rank probabilities
    std::discrete_distribution<size_t> dist(selection_probs.begin(), selection_probs.end());
    std::mt19937 rng(std::random_device{}());

    std::vector<Individual> selected;
    for (size_t i = 0; i < num_selected; ++i) {
        selected.push_back(sorted[dist(rng)]);
    }

    return selected;
}

// === Generate random individuals ===
std::vector<Individual> GeneticAlgorithm::generate_random_individuals(size_t count,
                                                                      std::mt19937 &rng) const {
    std::vector<Individual> randoms;
    for (size_t i = 0; i < count; ++i) {
        auto p = sample_random_parameters(bounds, rng);
        randoms.push_back({p, -1.0, ""});
    }
    return randoms;
}

// === Clamp to parameter bounds ===
void GeneticAlgorithm::clamp(ModelParameters &p) const {
    auto clamp = [](double val, double min, double max) { return std::clamp(val, min, max); };
    auto b = bounds;

#define CLAMP_FIELD(f) p.f = clamp(p.f, b.f##_min, b.f##_max);
    FOR_EACH_PARAM(CLAMP_FIELD)
#undef CLAMP_FIELD
}

// === Crossover ===
std::pair<ModelParameters, ModelParameters>
GeneticAlgorithm::crossover(const ModelParameters &parent1, const ModelParameters &parent2,
                            std::mt19937 &rng, size_t /*generation*/) const {
    ModelParameters child1, child2;
    const double alpha = 0.5;

#define CROSSOVER_FIELD(f)                                                                         \
    {                                                                                              \
        double min_val = std::min(parent1.f, parent2.f);                                           \
        double max_val = std::max(parent1.f, parent2.f);                                           \
        double range = max_val - min_val;                                                          \
        std::uniform_real_distribution<> dist(min_val - alpha * range, max_val + alpha * range);   \
        child1.f = dist(rng);                                                                      \
        child2.f = dist(rng);                                                                      \
    }

    FOR_EACH_PARAM(CROSSOVER_FIELD)
#undef CROSSOVER_FIELD

    return {child1, child2};
}
// === Mutation ===
void GeneticAlgorithm::mutate(ModelParameters &individual, std::mt19937 &rng, double mutation_rate,
                              double mutation_strength) const {
    std::uniform_real_distribution<> prob(0.0, 1.0);
    std::normal_distribution<> noise(0.0, mutation_strength);

    auto maybe_mutate = [&](double &val) {
        if (prob(rng) < mutation_rate)
            val += noise(rng);
    };

#define MUTATE_FIELD(f) maybe_mutate(individual.f);
    FOR_EACH_PARAM(MUTATE_FIELD)
#undef MUTATE_FIELD
}

// === Generate offspring + randoms ===
std::vector<ModelParameters>
GeneticAlgorithm::generate_variants(const std::vector<Individual> &parents, size_t n_offspring,
                                    size_t n_random, std::mt19937 &rng, size_t generation) const {
    std::vector<ModelParameters> children;

    std::uniform_int_distribution<size_t> dis(0, parents.size() - 1);
    while (children.size() < n_offspring - 1) {
        const auto &p1 = parents[dis(rng)].parameters;
        const auto &p2 = parents[dis(rng)].parameters;

        auto [c1, c2] = crossover(p1, p2, rng, generation);
        mutate(c1, rng, current_mutation_rate, current_mutation_strength);
        mutate(c2, rng, current_mutation_rate, current_mutation_strength);
        clamp(c1);
        clamp(c2);
        children.push_back(c1);
        if (children.size() < n_offspring)
            children.push_back(c2);
    }

    if (children.size() < n_offspring) {
        auto [c, _] = crossover(parents[0].parameters, parents[1].parameters, rng, generation);
        mutate(c, rng, current_mutation_rate, current_mutation_strength);
        clamp(c);
        children.push_back(c);
    }

    // Add random individuals
    for (size_t i = 0; i < n_random; ++i)
        children.push_back(sample_random_parameters(bounds, rng));

    return children;
}

// === Main generation method ===
std::vector<Individual> GeneticAlgorithm::generate_next_generation_from_directory(
    const std::string &population_dir, double elite_frac, double random_frac,
    size_t tournament_size, size_t generation) {
    std::vector<Individual> population = load_generation(population_dir);
    if (population.empty()) {
        std::cerr << "[GA] Warning: empty population — generating random fallback\n";
        std::mt19937 rng(std::random_device{}());
        return generate_random_individuals(total_size, rng);
    }

    size_t n_elite = static_cast<size_t>(elite_frac * total_size);
    size_t n_random = static_cast<size_t>(random_frac * total_size);
    size_t n_variants = total_size - n_elite;

    auto elites = select_top_n(population, n_elite);
    auto parents = tournament_selection(population, n_variants - n_random, tournament_size);

    std::mt19937 rng(std::random_device{}());
    auto children = generate_variants(parents, n_variants, n_random, rng, generation);

    std::vector<Individual> new_gen;
    for (const auto &e : elites)
        new_gen.push_back(e);
    for (const auto &c : children)
        new_gen.push_back({c, -1.0, ""});

    return new_gen;
}

void GeneticAlgorithm::update_adaptive_parameters(size_t generation,
                                                  const Individual &best_phenotype) {
    if (generation == 0) {
        current_mutation_rate = mutation_rate_init;
        current_mutation_strength = mutation_strength_init;
        crossover_alpha = crossover_alpha_init;
        return;
    }

    double delta = running_best_fitness - best_phenotype.fitness;
    running_best_fitness = best_phenotype.fitness;

    // Sensitivity thresholds
    if (delta > 0.1) { // big gain → exploit
        current_mutation_rate *= 0.1;
        current_mutation_strength *= 0.1;
        crossover_alpha *= 0.9;
    } else if (delta < 0.00) {
        current_mutation_rate *= 1.1;
        current_mutation_strength *= 1.1;
        crossover_alpha *= 1.05;
    } else { // moderate gain → gentle decay
        current_mutation_rate *= 0.95;
        current_mutation_strength *= 0.95;
        crossover_alpha *= 0.98;
    }

    // Clamp within safe bounds
    current_mutation_rate = std::clamp(current_mutation_rate, 1e-4, 1.0);
    current_mutation_strength = std::clamp(current_mutation_strength, 1e-4, 10.0);
    crossover_alpha = std::clamp(crossover_alpha, 0.05, 1.0);
}

void GeneticAlgorithm::stagnation_detect(const Individual &best_phenotype) {
    if (best_phenotype.fitness >= running_best_fitness ||
        running_best_fitness - best_phenotype.fitness < stagnation_threshold) {
        stagnation_flag = 1;
    } else {
        stagnation_flag = 0;
        running_best_fitness = best_phenotype.fitness;
    }
}

std::vector<Individual>
GeneticAlgorithm::generate_next_generation(std::vector<Individual> population, double elite_frac,
                                           double random_frac, size_t tournament_size,
                                           size_t generation) {

    if (population.empty()) {
        std::cerr << "[GA] Warning: empty population — generating random fallback\n";
        std::mt19937 rng(std::random_device{}());
        return generate_random_individuals(total_size, rng);
    }

    size_t n_elite = static_cast<size_t>(elite_frac * total_size);
    size_t n_random = static_cast<size_t>(random_frac * total_size);
    size_t n_variants = total_size - n_elite;
    auto best = best_phenotype(population);

    update_adaptive_parameters(generation, best);

    auto elites = select_top_n(population, n_elite);
    auto parents = rank_selection(population, n_variants - n_random);

    std::mt19937 rng(std::random_device{}());
    auto children = generate_variants(parents, n_variants, n_random, rng, generation);

    std::vector<Individual> new_gen;
    for (const auto &e : elites)
        new_gen.push_back(e);
    for (const auto &c : children)
        new_gen.push_back({c, -1.0, ""});

    return new_gen;
};