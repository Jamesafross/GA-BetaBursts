#include "genetic_algorithm/selection.hpp"
#include "neural_mass_model/json_utils.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

std::vector<Individual> Selection::select_top_n_from_directory(const std::string &dir,
                                                               size_t n) const {

    auto generation = load_generation(dir);
    auto top_n = select_top_n(generation, n);
    return top_n;
}

std::vector<Individual> Selection::load_generation(const std::string &dir) const {
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

void Selection::save_selected_individuals(const std::vector<Individual> &selected,
                                          const std::string &output_dir) const {
    std::filesystem::create_directories(output_dir);
    for (const auto &ind : selected) {
        nlohmann::json j;
        j["fitness"] = ind.fitness;
        j["parameters"] = ind.parameters;

        std::filesystem::path input_path = ind.source_path;
        std::string filename = input_path.filename().string();
        std::filesystem::path output_path = std::filesystem::path(output_dir) / filename;

        std::ofstream out(output_path);
        if (out.is_open()) {
            out << j.dump(4);
            out.close();
        } else {
            std::cerr << "Failed to write to " << output_path << std::endl;
        }
    }
}

std::vector<Individual> Selection::select_top_n(const std::vector<Individual> &population,
                                                size_t n) const {
    std::vector<Individual> sorted = population;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto &a, const auto &b) { return a.fitness < b.fitness; });

    if (n > sorted.size())
        n = sorted.size();
    return std::vector<Individual>(sorted.begin(), sorted.begin() + n);
}
