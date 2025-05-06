#include "genetic_algorithm/fitness.hpp"
#include "genetic_algorithm/ga_parameters.h"
#include "genetic_algorithm/generation_logger.hpp"
#include "genetic_algorithm/genetic_algorithm.hpp"
#include "meg/generate_real_stats.hpp"
#include "neural_mass_model/generate_model_stats.hpp"
#include "neural_mass_model/json_utils.hpp"
#include "neural_mass_model/meanfield.h"
#include "neural_mass_model/parameter_utils.hpp"
#include "neural_mass_model/parameters.h"
#include "paths.hpp"
#include "thresholding/threshold_structs.h"
#include "utils/python_plot.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <json.hpp>
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

void clear_summary_file(const std::string &path) {
    std::ofstream out(path, std::ios::trunc); // truncate clears the file
    if (out.is_open()) {
        nlohmann::json empty_array = nlohmann::json::array();
        out << empty_array.dump(4); // pretty print
        out.close();
    }
}

const std::string meg_output_dir = OUTPUT_PATH + "/meg";

int main() {
    GAParameters ga_p;
    ThresholdParameters th_params;
    FitnessEvaluator evaluator;
    auto bounds = ParameterBounds();

    double elite_frac = 0.05;
    double random_frac = 0.10;
    size_t tournament_size = 3;

    GeneticAlgorithm ga(bounds, ga_p.population_size, ga_p.mutation_rate, 0.1,
                        ga_p.num_generations);

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

            auto result = evaluator.compute_fitness_from_csv(
                stats_dir + "/" + stats_csv, OUTPUT_PATH + "/meg/meg_burst_stats_merged.csv",
                summary_path);

            double fitness = result.total;

            save_parameters_with_fitness(population[i], result, params_dir, json_file);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            std::cout << "Sim " << i + 1 << "/" << population.size() << " (fitness=" << fitness
                      << ") in " << elapsed.count() << "s\n";
        }

        logger.append_summary_from_directory(params_dir, gen);

        // plot violins
        auto best_id = get_best_phenotype_id(summary_path);
        run_violin_plot(gen, best_id);

        // --- Clean up large simulated files ---
        cleanup_simulated_data(sim_dir);

        std::vector<Individual> population = ga.generate_next_generation_from_directory(
            params_dir, elite_frac, random_frac, tournament_size, gen);

        // std::string output_dir = base_dir + "/gen" + std::to_string(gen + 1);
        // ga.save_generation(next_gen, output_dir);
    }

    return 0;
}
