#include "genetic_algorithm/generation_logger.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

GenerationLogger::GenerationLogger(const std::string &summary_file_path)
    : summary_path(summary_file_path) {}

void GenerationLogger::append_summary_from_directory(const std::string &directory, int generation) {
    auto [phenotype_infos, emd_sum, emd_count, ks_sum, ks_count, stat_sum, stat_count] =
        load_phenotype_data(directory);

    if (phenotype_infos.empty()) {
        std::cerr << "No valid fitness entries found in directory: " << directory << "\n";
        return;
    }

    auto [best_fit, mean_fit, best_phenotype] = compute_fitness_stats(phenotype_infos);

    double mean_emd = (emd_count > 0) ? emd_sum / emd_count : -1.0;
    double mean_ks = (ks_count > 0) ? ks_sum / ks_count : -1.0;
    double mean_stat = (stat_count > 0) ? stat_sum / stat_count : -1.0;

    nlohmann::json summary_entry = {{"generation", generation}, {"best_fitness", best_fit},
                                    {"mean_fitness", mean_fit}, {"best_phenotype", best_phenotype},
                                    {"mean_emd", mean_emd},     {"mean_ks", mean_ks},
                                    {"mean_stat", mean_stat}};

    append_to_summary_log(summary_entry);
}

std::tuple<std::vector<std::pair<int, double>>, double, size_t, double, size_t, double, size_t>
GenerationLogger::load_phenotype_data(const std::string &directory) const {
    std::vector<std::pair<int, double>> infos;
    double emd_sum = 0.0, ks_sum = 0.0, stat_sum = 0.0;
    ;
    size_t emd_count = 0, ks_count = 0, stat_count = 0.0;
    ;

    std::regex number_pattern(R"((\d+))");
    std::smatch match;

    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().extension() != ".json")
            continue;

        std::ifstream in(entry.path());
        if (!in.is_open())
            continue;

        try {
            nlohmann::json j;
            in >> j;

            if (j.contains("fitness")) {
                double fitness = j["fitness"].get<double>();
                std::string filename = entry.path().filename().string();
                int num = -1;
                if (std::regex_search(filename, match, number_pattern))
                    num = std::stoi(match[1]);
                infos.emplace_back(num, fitness);
            }

            if (j.contains("emd_fitness")) {
                emd_sum += j["emd_fitness"].get<double>();
                ++emd_count;
            }

            if (j.contains("ks_fitness")) {
                ks_sum += j["ks_fitness"].get<double>();
                ++ks_count;
            }

            if (j.contains("stat_fitness")) {
                stat_sum += j["stat_fitness"].get<double>();
                ++stat_count;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error reading " << entry.path() << ": " << e.what() << "\n";
        }
    }

    return {infos, emd_sum, emd_count, ks_sum, ks_count, stat_sum, stat_count};
}

std::tuple<double, double, int>
GenerationLogger::compute_fitness_stats(const std::vector<std::pair<int, double>> &infos) const {
    double best = std::numeric_limits<double>::infinity();
    double sum = 0.0;
    int best_pheno = -1;

    for (const auto &[id, fit] : infos) {
        if (fit < best) {
            best = fit;
            best_pheno = id;
        }
        sum += fit;
    }

    return {best, sum / infos.size(), best_pheno};
}

void GenerationLogger::append_to_summary_log(const nlohmann::json &entry) const {
    nlohmann::json full_log;

    std::ifstream infile(summary_path);
    if (infile.is_open()) {
        try {
            infile >> full_log;
        } catch (...) {
            full_log = nlohmann::json::array();
        }
    }

    full_log.push_back(entry);

    std::ofstream out(summary_path);
    out << full_log.dump(4);
}