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

std::tuple<std::vector<std::pair<int, FitnessResult>>, double, size_t, double, size_t, double,
           size_t>
GenerationLogger::load_phenotype_data(const std::string &directory) const {
    std::vector<std::pair<int, FitnessResult>> infos;

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

            if (!j.contains("fitness"))
                continue;

            FitnessResult result{};
            result.total = j.value("fitness", 0.0);
            result.emd_fitness = j.value("emd_fitness", 0.0);
            result.ks_fitness = j.value("ks_fitness", 0.0);
            result.stat_fitness = j.value("stat_fitness", 0.0);
            result.raw_total = j.value("raw_total", 0.0);
            result.raw_emd = j.value("raw_emd", 0.0);
            result.raw_ks = j.value("raw_ks", 0.0);
            result.raw_stat = j.value("raw_stat", 0.0);

            std::string filename = entry.path().filename().string();
            int num = -1;
            if (std::regex_search(filename, match, number_pattern))
                num = std::stoi(match[1]);
            infos.emplace_back(num, result);

        } catch (const std::exception &e) {
            std::cerr << "Error reading " << entry.path() << ": " << e.what() << "\n";
        }
    }

    // Sort by total fitness (lower is better)
    std::sort(infos.begin(), infos.end(),
              [](const auto &a, const auto &b) { return a.second.total < b.second.total; });

    // Keep only top 10
    if (infos.size() > 10)
        infos.resize(10);

    // Compute means across top 10 (or fewer)
    double emd_sum = 0.0, ks_sum = 0.0, stat_sum = 0.0;
    size_t emd_count = 0, ks_count = 0, stat_count = 0;

    for (const auto &[_, res] : infos) {
        if (res.emd_fitness != 0.0) {
            emd_sum += res.emd_fitness;
            ++emd_count;
        }
        if (res.ks_fitness != 0.0) {
            ks_sum += res.ks_fitness;
            ++ks_count;
        }
        if (res.stat_fitness != 0.0) {
            stat_sum += res.stat_fitness;
            ++stat_count;
        }
    }

    return {infos, emd_sum, emd_count, ks_sum, ks_count, stat_sum, stat_count};
}

std::tuple<double, double, int> GenerationLogger::compute_fitness_stats(
    const std::vector<std::pair<int, FitnessResult>> &infos) const {
    double best = std::numeric_limits<double>::infinity();
    double sum = 0.0;
    int best_pheno = -1;

    for (const auto &[id, res] : infos) {
        if (res.total < best) {
            best = res.total;
            best_pheno = id;
        }
        sum += res.total;
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