#include "neural_mass_model/parameter_utils.hpp"
#include <iomanip> // required for std::fixed and std::setprecision
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

ModelParameters decode_parameters(const std::vector<double> &genome,
                                  const ParameterBounds &bounds) {
    if (genome.size() != 25) {
        throw std::invalid_argument("Genome must have 25 elements.");
    }

    ModelParameters params;
    int i = 0;

    params.alphaee = bounds.alphaee_min + genome[i++] * (bounds.alphaee_max - bounds.alphaee_min);
    params.alphaei = bounds.alphaei_min + genome[i++] * (bounds.alphaei_max - bounds.alphaei_min);
    params.alphaie = bounds.alphaie_min + genome[i++] * (bounds.alphaie_max - bounds.alphaie_min);
    params.alphaii = bounds.alphaii_min + genome[i++] * (bounds.alphaii_max - bounds.alphaii_min);

    params.kappa_see =
        bounds.kappa_see_min + genome[i++] * (bounds.kappa_see_max - bounds.kappa_see_min);
    params.kappa_sei =
        bounds.kappa_sei_min + genome[i++] * (bounds.kappa_sei_max - bounds.kappa_sei_min);
    params.kappa_sie =
        bounds.kappa_sie_min + genome[i++] * (bounds.kappa_sie_max - bounds.kappa_sie_min);
    params.kappa_sii =
        bounds.kappa_sii_min + genome[i++] * (bounds.kappa_sii_max - bounds.kappa_sii_min);

    params.kappa_gapee =
        bounds.kappa_gapee_min + genome[i++] * (bounds.kappa_gapee_max - bounds.kappa_gapee_min);
    params.kappa_gapei =
        bounds.kappa_gapei_min + genome[i++] * (bounds.kappa_gapei_max - bounds.kappa_gapei_min);
    params.kappa_gapii =
        bounds.kappa_gapii_min + genome[i++] * (bounds.kappa_gapii_max - bounds.kappa_gapii_min);

    params.vsyn_ee = bounds.vsyn_ee_min + genome[i++] * (bounds.vsyn_ee_max - bounds.vsyn_ee_min);
    params.vsyn_ei = bounds.vsyn_ei_min + genome[i++] * (bounds.vsyn_ei_max - bounds.vsyn_ei_min);
    params.vsyn_ie = bounds.vsyn_ie_min + genome[i++] * (bounds.vsyn_ie_max - bounds.vsyn_ie_min);
    params.vsyn_ii = bounds.vsyn_ii_min + genome[i++] * (bounds.vsyn_ii_max - bounds.vsyn_ii_min);

    params.tauxe = bounds.tauxe_min + genome[i++] * (bounds.tauxe_max - bounds.tauxe_min);
    params.tauxi = bounds.tauxi_min + genome[i++] * (bounds.tauxi_max - bounds.tauxi_min);
    params.taue = bounds.taue_min + genome[i++] * (bounds.taue_max - bounds.taue_min);
    params.taui = bounds.taui_min + genome[i++] * (bounds.taui_max - bounds.taui_min);

    params.deltae = bounds.deltae_min + genome[i++] * (bounds.deltae_max - bounds.deltae_min);
    params.deltai = bounds.deltai_min + genome[i++] * (bounds.deltai_max - bounds.deltai_min);

    params.sigma_se =
        bounds.sigma_se_min + genome[i++] * (bounds.sigma_se_max - bounds.sigma_se_min);
    params.sigma_si =
        bounds.sigma_si_min + genome[i++] * (bounds.sigma_si_max - bounds.sigma_si_min);

    params.eta0e = bounds.eta0e_min + genome[i++] * (bounds.eta0e_max - bounds.eta0e_min);
    params.eta0i = bounds.eta0i_min + genome[i++] * (bounds.eta0i_max - bounds.eta0i_min);

    return params;
}

ModelParameters sample_random_parameters(const ParameterBounds &bounds, std::mt19937 &rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<double> genome(25);
    for (auto &gene : genome) {
        gene = dist(rng);
    }
    return decode_parameters(genome, bounds);
}

std::vector<Individual> generate_parameter_population(size_t n, const ParameterBounds &bounds,
                                                      std::mt19937 &rng) {
    std::vector<Individual> population;
    population.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        population[i].parameters = sample_random_parameters(bounds, rng);
    }
    return population;
}

void print_parameters(const ModelParameters &params) {
    std::cout << std::fixed << std::setprecision(3);

    std::cout << "alphaee: " << params.alphaee << ", "
              << "alphaei: " << params.alphaei << ", "
              << "alphaie: " << params.alphaie << ", "
              << "alphaii: " << params.alphaii << "\n"
              << "kappa_see: " << params.kappa_see << ", "
              << "kappa_sei: " << params.kappa_sei << ", "
              << "kappa_sie: " << params.kappa_sie << ", "
              << "kappa_sii: " << params.kappa_sii << "\n"
              << "vsyn_ee: " << params.vsyn_ee << ", "
              << "vsyn_ei: " << params.vsyn_ei << ", "
              << "vsyn_ie: " << params.vsyn_ie << ", "
              << "vsyn_ii: " << params.vsyn_ii << "\n"
              << "tauxe: " << params.tauxe << ", "
              << "tauxi: " << params.tauxi << ", "
              << "sigma_se: " << params.sigma_se << ", "
              << "sigma_si: " << params.sigma_si << "\n"
              << "kappa_gapee: " << params.kappa_gapee << ", "
              << "kappa_gapei: " << params.kappa_gapei << ", "
              << "kappa_gapii: " << params.kappa_gapii << "\n"
              << "eta0e: " << params.eta0e << ", "
              << "eta0i: " << params.eta0i << ", "
              << "taue: " << params.taue << ", "
              << "taui: " << params.taui << ", "
              << "deltae: " << params.deltae << ", "
              << "deltai: " << params.deltai << "\n"
              << "----------------------\n";
}
