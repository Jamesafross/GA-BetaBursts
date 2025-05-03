#include "genetic_algorithm/genetic_operators.hpp"
#include "neural_mass_model/parameter_utils.hpp"
#include <algorithm>
#include <random>
#include <vector>

void GeneticOperators::clamp_parameters_to_bounds(ModelParameters &p,
                                                  const ParameterBounds &b) const {
    auto clamp = [](double val, double min, double max) { return std::clamp(val, min, max); };

    p.alphaee = clamp(p.alphaee, b.alphaee_min, b.alphaee_max);
    p.alphaei = clamp(p.alphaei, b.alphaei_min, b.alphaei_max);
    p.alphaie = clamp(p.alphaie, b.alphaie_min, b.alphaie_max);
    p.alphaii = clamp(p.alphaii, b.alphaii_min, b.alphaii_max);

    p.kappa_see = clamp(p.kappa_see, b.kappa_see_min, b.kappa_see_max);
    p.kappa_sei = clamp(p.kappa_sei, b.kappa_sei_min, b.kappa_sei_max);
    p.kappa_sie = clamp(p.kappa_sie, b.kappa_sie_min, b.kappa_sie_max);
    p.kappa_sii = clamp(p.kappa_sii, b.kappa_sii_min, b.kappa_sii_max);

    p.vsyn_ee = clamp(p.vsyn_ee, b.vsyn_ee_min, b.vsyn_ee_max);
    p.vsyn_ei = clamp(p.vsyn_ei, b.vsyn_ei_min, b.vsyn_ei_max);
    p.vsyn_ie = clamp(p.vsyn_ie, b.vsyn_ie_min, b.vsyn_ie_max);
    p.vsyn_ii = clamp(p.vsyn_ii, b.vsyn_ii_min, b.vsyn_ii_max);

    p.tauxe = clamp(p.tauxe, b.tauxe_min, b.tauxe_max);
    p.tauxi = clamp(p.tauxi, b.tauxi_min, b.tauxi_max);
    p.taue = clamp(p.taue, b.taue_min, b.taue_max);
    p.taui = clamp(p.taui, b.taui_min, b.taui_max);

    p.deltae = clamp(p.deltae, b.deltae_min, b.deltae_max);
    p.deltai = clamp(p.deltai, b.deltai_min, b.deltai_max);
    p.sigma_se = clamp(p.sigma_se, b.sigma_se_min, b.sigma_se_max);
    p.sigma_si = clamp(p.sigma_si, b.sigma_si_min, b.sigma_si_max);

    p.kappa_gapee = clamp(p.kappa_gapee, b.kappa_gapee_min, b.kappa_gapee_max);
    p.kappa_gapei = clamp(p.kappa_gapei, b.kappa_gapei_min, b.kappa_gapei_max);
    p.kappa_gapii = clamp(p.kappa_gapii, b.kappa_gapii_min, b.kappa_gapii_max);

    p.eta0e = clamp(p.eta0e, b.eta0e_min, b.eta0e_max);
    p.eta0i = clamp(p.eta0i, b.eta0i_min, b.eta0i_max);
}

std::pair<ModelParameters, ModelParameters>
GeneticOperators::crossover(const ModelParameters &parent1, const ModelParameters &parent2,
                            std::mt19937 &rng, size_t generation) const {
    ModelParameters child1 = parent1;
    ModelParameters child2 = parent2;

    // Adaptive linear crossover probability schedule
    double start_p = 0.5;
    double end_p = 0.1;
    double t = static_cast<double>(generation) / max_generations;
    double crossover_prob = start_p - t * (start_p - end_p);
    crossover_prob = std::clamp(crossover_prob, end_p, start_p);

    // Updated macro with adaptive probability
#define CROSSOVER_FIELD(field, prob, rng)                                                          \
    if (std::bernoulli_distribution(prob)(rng))                                                    \
        std::swap(child1.field, child2.field);

    // Apply crossover to all model fields
    CROSSOVER_FIELD(alphaee, crossover_prob, rng)
    CROSSOVER_FIELD(alphaei, crossover_prob, rng)
    CROSSOVER_FIELD(alphaie, crossover_prob, rng)
    CROSSOVER_FIELD(alphaii, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_see, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_sei, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_sie, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_sii, crossover_prob, rng)
    CROSSOVER_FIELD(vsyn_ee, crossover_prob, rng)
    CROSSOVER_FIELD(vsyn_ei, crossover_prob, rng)
    CROSSOVER_FIELD(vsyn_ie, crossover_prob, rng)
    CROSSOVER_FIELD(vsyn_ii, crossover_prob, rng)
    CROSSOVER_FIELD(tauxe, crossover_prob, rng)
    CROSSOVER_FIELD(tauxi, crossover_prob, rng)
    CROSSOVER_FIELD(taue, crossover_prob, rng)
    CROSSOVER_FIELD(taui, crossover_prob, rng)
    CROSSOVER_FIELD(deltae, crossover_prob, rng)
    CROSSOVER_FIELD(deltai, crossover_prob, rng)
    CROSSOVER_FIELD(sigma_se, crossover_prob, rng)
    CROSSOVER_FIELD(sigma_si, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_gapee, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_gapei, crossover_prob, rng)
    CROSSOVER_FIELD(kappa_gapii, crossover_prob, rng)
    CROSSOVER_FIELD(eta0e, crossover_prob, rng)
    CROSSOVER_FIELD(eta0i, crossover_prob, rng)

    return {child1, child2};
}

void GeneticOperators::mutate(ModelParameters &individual, std::mt19937 &rng) const {
    std::uniform_real_distribution<> prob(0.0, 1.0);
    std::normal_distribution<> noise(0.0, mutation_strength);

    auto maybe_mutate = [&](double &val) {
        if (prob(rng) < mutation_rate)
            val += noise(rng);
    };

    maybe_mutate(individual.alphaee);
    maybe_mutate(individual.alphaei);
    maybe_mutate(individual.alphaie);
    maybe_mutate(individual.alphaii);
    maybe_mutate(individual.kappa_see);
    maybe_mutate(individual.kappa_sei);
    maybe_mutate(individual.kappa_sie);
    maybe_mutate(individual.kappa_sii);
    maybe_mutate(individual.vsyn_ee);
    maybe_mutate(individual.vsyn_ei);
    maybe_mutate(individual.vsyn_ie);
    maybe_mutate(individual.vsyn_ii);
    maybe_mutate(individual.tauxe);
    maybe_mutate(individual.tauxi);
    maybe_mutate(individual.taue);
    maybe_mutate(individual.taui);
    maybe_mutate(individual.deltae);
    maybe_mutate(individual.deltai);
    maybe_mutate(individual.sigma_se);
    maybe_mutate(individual.sigma_si);
    maybe_mutate(individual.kappa_gapee);
    maybe_mutate(individual.kappa_gapei);
    maybe_mutate(individual.kappa_gapii);
    maybe_mutate(individual.eta0e);
    maybe_mutate(individual.eta0i);
}

std::vector<ModelParameters>
GeneticOperators::generate_next_population(const std::vector<Individual> &selected,
                                           size_t total_size, size_t n_random, std::mt19937 &rng,
                                           const ParameterBounds &bounds, size_t generation) const {
    std::vector<ModelParameters> next_gen;

    // Step 1: Elites
    for (const auto &ind : selected) {
        next_gen.push_back(ind.parameters);
    }

    // Step 2: Children (fill up to total_size - n_random)
    std::uniform_int_distribution<size_t> dist(0, selected.size() - 1);
    while (next_gen.size() < total_size - n_random - 1) {
        const auto &parent1 = selected[dist(rng)].parameters;
        const auto &parent2 = selected[dist(rng)].parameters;

        auto [child1, child2] = crossover(parent1, parent2, rng, generation);
        mutate(child1, rng);
        mutate(child2, rng);
        clamp_parameters_to_bounds(child1, bounds);
        clamp_parameters_to_bounds(child2, bounds);

        next_gen.push_back(child1);
        if (next_gen.size() < total_size - n_random)
            next_gen.push_back(child2);
    }

    // Step 3: Optional one more child if needed
    if (next_gen.size() < total_size - n_random) {
        const auto &p1 = selected[dist(rng)].parameters;
        const auto &p2 = selected[dist(rng)].parameters;
        auto [child, _] = crossover(p1, p2, rng, generation);
        mutate(child, rng);
        clamp_parameters_to_bounds(child, bounds);
        next_gen.push_back(child);
    }

    // Step 4: Add n_random brand new individuals
    for (size_t i = 0; i < n_random; ++i) {
        next_gen.push_back(sample_random_parameters(bounds, rng));
    }

    return next_gen;
}