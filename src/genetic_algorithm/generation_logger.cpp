#include "genetic_algorithm/generation_logger.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

GenerationLogger::GenerationLogger(const std::string &summary_file_path)
    : summary_path(summary_file_path) {}

#include <regex> // Add this at the top

void GenerationLogger::append_summary_from_directory(const std::string &directory, int generation) {
    std::vector<std::pair<int, double>> phenotype_infos;
    std::regex number_pattern(R"((\d+))");
    std::smatch match;

    double emd_sum = 0.0;
    double ks_sum = 0.0;
    size_t emd_count = 0;
    size_t ks_count = 0;

    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().extension() != ".json")
            continue;

        std::ifstream in(entry.path());
        if (!in.is_open()) {
            std::cerr << "Could not open " << entry.path() << "\n";
            continue;
        }

        try {
            nlohmann::json j;
            in >> j;

            if (j.contains("fitness")) {
                double fitness = j["fitness"].get<double>();

                // Extract phenotype number from filename
                std::string filename = entry.path().filename().string();
                int phenotype_number = -1;
                if (std::regex_search(filename, match, number_pattern)) {
                    phenotype_number = std::stoi(match[1]);
                }

                phenotype_infos.emplace_back(phenotype_number, fitness);
            }

            if (j.contains("emd_fitness")) {
                emd_sum += j["emd_fitness"].get<double>();
                ++emd_count;
            }

            if (j.contains("ks_fitness")) {
                ks_sum += j["ks_fitness"].get<double>();
                ++ks_count;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error reading " << entry.path() << ": " << e.what() << "\n";
        }
    }

    if (phenotype_infos.empty()) {
        std::cerr << "No valid fitness entries found in directory: " << directory << "\n";
        return;
    }

    // Compute best and mean fitness
    double best = std::numeric_limits<double>::infinity();
    double sum = 0.0;
    int best_phenotype = -1;

    for (const auto &[phenotype, fit] : phenotype_infos) {
        if (fit < best) {
            best = fit;
            best_phenotype = phenotype;
        }
        sum += fit;
    }

    double mean_fitness = sum / phenotype_infos.size();
    double mean_emd = (emd_count > 0) ? emd_sum / emd_count : -1.0;
    double mean_ks = (ks_count > 0) ? ks_sum / ks_count : -1.0;

    nlohmann::json summary_entry = {
        {"generation", generation},         {"best_fitness", best}, {"mean_fitness", mean_fitness},
        {"best_phenotype", best_phenotype}, {"mean_emd", mean_emd}, {"mean_ks", mean_ks}};

    nlohmann::json full_log;

    std::ifstream infile(summary_path);
    if (infile.is_open()) {
        try {
            infile >> full_log;
        } catch (...) {
            full_log = nlohmann::json::array();
        }
    }

    full_log.push_back(summary_entry);

    std::ofstream out(summary_path);
    out << full_log.dump(4);
}
