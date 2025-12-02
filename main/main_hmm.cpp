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
#include "HMM/hmm_test.hpp"
#include "paths.hpp"
#include "thresholding/threshold_structs.h"
#include "utils/python_plot.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <json.hpp>
#include <random>
#include <string>

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

    GeneticAlgorithm ga(bounds, ga_p);

    std::random_device rd;
    std::mt19937 rng(rd());

    std::filesystem::create_directories(OUTPUT_PATH + "/meg");

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
    ga_p.num_generations = 1;
    SolverParameters opts;
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
        for (size_t i = 0; i < ga_p.population_size; ++i) {
            std::cout << "\n=== Generation " << gen << " ===\n";

            std::string base = "phenotype_" + std::to_string(i);
            std::string sim_csv = base + "_simulated.csv";
            std::string stats_csv = base + "_stats.csv";
            std::string json_file = base + ".json";

            

            run_neural_mass(population[i].parameters, opts, sim_dir, sim_csv);

            
        
           
        };
        std::string popcurrent_path = sim_dir;
        std::string hmm_output_path = gen_dir + "/HMM";
        double sampling_freq = opts.fs;

        hmm_matlab_run(popcurrent_path, hmm_output_path, sampling_freq);
    };

    return 0;
}