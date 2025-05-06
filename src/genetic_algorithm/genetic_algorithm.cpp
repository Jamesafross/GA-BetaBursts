#include "genetic_algorithm/genetic_algorithm.hpp"
#include "neural_mass_model/json_utils.hpp"

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
GeneticAlgorithm::GeneticAlgorithm(ParameterBounds bounds, size_t total_size, double mutation_rate,
                                   double mutation_strength, size_t max_generations)
    : bounds(bounds), total_size(total_size), mutation_rate(mutation_rate),
      mutation_strength(mutation_strength), max_generations(max_generations) {}

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
                            std::mt19937 &rng, size_t generation) const {
    ModelParameters child1 = parent1;
    ModelParameters child2 = parent2;

    // Adaptive linear crossover probability schedule
    double start_p = 0.5;
    double end_p = 0.1;
    double t = static_cast<double>(generation) / max_generations;
    double crossover_prob = start_p - t * (start_p - end_p);
    crossover_prob = std::clamp(crossover_prob, end_p, start_p);

    // Define macro to apply to each field
#define CROSSOVER_FIELD(f)                                                                         \
    if (std::bernoulli_distribution(crossover_prob)(rng))                                          \
        std::swap(child1.f, child2.f);

    FOR_EACH_PARAM(CROSSOVER_FIELD)
#undef CROSSOVER_FIELD

    return {child1, child2};
}

// === Mutation ===
void GeneticAlgorithm::mutate(ModelParameters &individual, std::mt19937 &rng) const {
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
        mutate(c1, rng);
        mutate(c2, rng);
        clamp(c1);
        clamp(c2);
        children.push_back(c1);
        if (children.size() < n_offspring)
            children.push_back(c2);
    }

    if (children.size() < n_offspring) {
        auto [c, _] = crossover(parents[0].parameters, parents[1].parameters, rng, generation);
        mutate(c, rng);
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
