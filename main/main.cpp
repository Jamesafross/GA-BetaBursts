#include "genetic_algorithm/fitness.hpp"
#include "genetic_algorithm/ga_parameters.h"
#include "genetic_algorithm/generation_logger.hpp"
#include "genetic_algorithm/genetic_operators.hpp"
#include "genetic_algorithm/selection.hpp"
#include "meg/generate_real_stats.hpp"
#include "neural_mass_model/generate_model_stats.hpp"
#include "neural_mass_model/json_utils.hpp"
#include "neural_mass_model/meanfield.h"
#include "neural_mass_model/parameter_utils.hpp"
#include "paths.hpp"
#include "thresholding/threshold_structs.h"
#include "utils/python_plot.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <random>

void ensure_output_directories(const std::string &base_dir, const std::string &stats_dir,
                               const std::string &simulated_dir,
                               const std::string &parameters_dir) {
    std::filesystem::create_directories(base_dir);
    std::filesystem::create_directories(stats_dir);
    std::filesystem::create_directories(simulated_dir);
    std::filesystem::create_directories(parameters_dir);
}

void cleanup_simulated_data(const std::string &simulated_dir) {
    try {
        if (std::filesystem::exists(simulated_dir)) {
            std::filesystem::remove_all(simulated_dir);
            std::cout << "Deleted simulated data directory: " << simulated_dir << "\n";
        }
    } catch (const std::filesystem::filesystem_error &e) {
        std::cerr << "Error deleting simulated data: " << e.what() << "\n";
    }
}

const std::string meg_output_dir = OUTPUT_PATH + "/meg";

int main() {
    GAParameters ga_p;
    ThresholdParameters th_params;
    FitnessEvaluator evaluator(0.1);
    Selection selection;
    GeneticOperators genetic_operators(ga_p.mutation_rate, 0.2, ga_p.num_generations);

    std::random_device rd;
    std::mt19937 rng(rd());

    // === Step 1: Prepare MEG stats ===
    th_params.original_fs = 600.0;
    th_params.downsampled_rate = 100.0;
    th_params.percentile_threshold = 0.75;
    th_params.min_burst_duration_s = 0.05;
    th_params.zscore_threshold = 2.0;
    th_params.merge_gap_s = 0.03;
    th_params.beta_low = 13.0;
    th_params.beta_high = 30.0;
    th_params.threshold_type = ThresholdType::Percentile;

    std::filesystem::create_directories(OUTPUT_PATH + "/meg");
    generate_real_stats(th_params);

    // logging
    std::string summary_path = OUTPUT_PATH + "/summary.json";

    // Delete existing summary if present
    if (std::filesystem::exists(summary_path)) {
        std::filesystem::remove(summary_path);
    }

    // clear existing files if present
    if (std::filesystem::exists(OUTPUT_PATH + "/model")) {
        std::filesystem::remove_all(OUTPUT_PATH + "/model");
    }

    GenerationLogger logger(summary_path);

    // === Step 2: Initial random population ===
    ParameterBounds bounds;
    auto population = generate_parameter_population(ga_p.population_size, bounds, rng);

    // === Step 3: Evolve over generations ===
    for (size_t gen = 0; gen < ga_p.num_generations; ++gen) {
        std::cout << "\n=== Generation " << gen << " ===\n";

        std::string gen_dir = OUTPUT_PATH + "/model/generation_" + std::to_string(gen);
        std::string stats_dir = gen_dir + "/stats";
        std::string sim_dir = gen_dir + "/simulated_data";
        std::string params_dir = gen_dir + "/parameters";
        std::string selected_dir = gen_dir + "/selected";

        ensure_output_directories(gen_dir, stats_dir, sim_dir, params_dir);

        // --- Run simulations and evaluate ---
        for (size_t i = 0; i < population.size(); ++i) {
            std::string base = "phenotype_" + std::to_string(i);
            std::string sim_csv = base + "_simulated.csv";
            std::string stats_csv = base + "_stats.csv";
            std::string json_file = base + ".json";

            SolverParameters opts;

            auto start = std::chrono::high_resolution_clock::now();

            run_neural_mass(population[i], opts, sim_dir, sim_csv);
            th_params.original_fs = 200.0;
            generate_model_stats(sim_dir, stats_dir, sim_csv, stats_csv, th_params);

            double fitness = evaluator.compute_fitness_from_csv(
                stats_dir + "/" + stats_csv, OUTPUT_PATH + "/meg/meg_burst_stats_merged.csv");

            save_parameters_with_fitness(population[i], fitness, params_dir, json_file);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            std::cout << "Sim " << i + 1 << "/" << population.size() << " (fitness=" << fitness
                      << ") in " << elapsed.count() << "s\n";
        }

        // or whatever generation index you're on
        logger.append_summary_from_directory(params_dir, gen);

        // plot violins
        auto best_id = get_best_phenotype_id(summary_path);
        run_violin_plot(gen, best_id);

        // --- Clean up large simulated files ---
        cleanup_simulated_data(sim_dir);

        // --- Selection & Next Gen ---
        auto top_n = selection.select_top_n_from_directory(
            params_dir, std::max<size_t>(1, ga_p.elite_fraction * ga_p.population_size));
        selection.save_selected_individuals(top_n, selected_dir);

        // --- Generate next population ---
        size_t n_random = std::max<size_t>(1, static_cast<size_t>(0.2 * ga_p.population_size));

        population = genetic_operators.generate_next_population(top_n, ga_p.population_size,
                                                                n_random, rng, bounds, gen);
    }

    return 0;
}
