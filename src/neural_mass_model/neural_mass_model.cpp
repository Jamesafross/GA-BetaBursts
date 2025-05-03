#include "neural_mass_model/meanfield.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <atomic>
#include <mutex>
#include <filesystem>
#include "neural_mass_model/euler_maruyama.hpp"

// === Save multiple trials to CSV ===
void save_matrix(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1)
                file << ",";
        }
        file << "\n";
    }
    file.close();
    std::cout << "Saved to " << filename << std::endl;
}

// === Mean-field model implementation ===
meanfield::meanfield(const ModelParameters& p_) : p(p_) {}

void meanfield::operator()(const std::vector<double>& state, std::vector<double>& dxdt, double /*t*/) const {
    const double ve = state[0], vi = state[1];
    const double re = state[2], ri = state[3];
    const double gee = state[4], phiee = state[5];
    const double gei = state[6], phiei = state[7];
    const double gie = state[8], phiie = state[9];
    const double gii = state[10], phiii = state[11];
    const double se = state[12], si = state[13];

    const double r_gap_e = -re * (p.kappa_gapee + p.kappa_gapei);
    const double r_gap_i = -ri * (p.kappa_gapii + p.kappa_gapei);
    const double v_gap_e = p.kappa_gapei * (vi - ve);
    const double v_gap_i = p.kappa_gapei * (ve - vi);

    const double syn_e = gee * (p.vsyn_ee - ve) + gei * (p.vsyn_ei - ve);
    const double syn_i = gii * (p.vsyn_ii - vi) + gie * (p.vsyn_ie - vi);

    const double pi_tau_e = M_PI * p.taue;
    const double pi_tau_i = M_PI * p.taui;

    const double pi_tau_e2 = pi_tau_e * re;
    const double pi_tau_i2 = pi_tau_i * ri;

    dxdt[0] = (p.eta0e + se + ve * ve - pi_tau_e2 * pi_tau_e2 + syn_e + v_gap_e) / p.taue;
    dxdt[1] = (p.eta0i + si + vi * vi - pi_tau_i2 * pi_tau_i2 + syn_i + v_gap_i) / p.taui;
    dxdt[2] = (r_gap_e + 2 * re * ve + (p.deltae / pi_tau_e)) / p.taue;
    dxdt[3] = (r_gap_i + 2 * ri * vi + (p.deltai / pi_tau_i)) / p.taui;
    dxdt[4] = p.alphaee * (phiee - gee);
    dxdt[5] = p.alphaee * (p.kappa_see * re - phiee);
    dxdt[6] = p.alphaei * (phiei - gei);
    dxdt[7] = p.alphaei * (p.kappa_sei * ri - phiei);
    dxdt[8] = p.alphaie * (phiie - gie);
    dxdt[9] = p.alphaie * (p.kappa_sie * re - phiie);
    dxdt[10] = p.alphaii * (phiii - gii);
    dxdt[11] = p.alphaii * (p.kappa_sii * ri - phiii);
    dxdt[12] = -se / p.tauxe;
    dxdt[13] = -si / p.tauxi;
}

void meanfield::apply_noise(std::vector<double>& state, double dt, std::mt19937& rng) const {
    std::normal_distribution<double> normal(0.0, 1.0);
    state[12] += p.sigma_se * normal(rng) * std::sqrt(dt);
    state[13] += p.sigma_si * normal(rng) * std::sqrt(dt);
}


void run_neural_mass(const ModelParameters& p,
                     const SolverParameters& opts,
                     const std::string& output_dir,
                     const std::string& filename)
{
    const size_t num_steps = opts.t_end();
    const size_t start_step = opts.save_start_step();
    const size_t save_interval = opts.save_every();
    const size_t steps_to_save = (num_steps > start_step)
        ? (num_steps - start_step + save_interval - 1) / save_interval
        : 0;

    std::vector<std::vector<double>> all_i_net(opts.n_trial, std::vector<double>(steps_to_save));
    std::atomic<size_t> completed_trials(0);
    std::mutex cout_mutex;

    auto print_progress = [&](size_t done, size_t total) {
        std::lock_guard<std::mutex> lock(cout_mutex);
        double percent = 100.0 * done / total;
        std::cout << "\rProgress: " << done << "/" << total << " (" << percent << "%)" << std::flush;
    };

    #pragma omp parallel for
    for (size_t trial = 0; trial < opts.n_trial; ++trial) {
        std::vector<double> state(14, 1.0);
        meanfield model(p);
        EulerMaruyamaStepper<meanfield> stepper(model, opts.dt);

        double t = 0.0;
        size_t save_counter = 0;

        for (size_t step = 0; step < num_steps; ++step) {
            if (step >= start_step && step % save_interval == 0 && save_counter < steps_to_save) {
                auto i_e = state[4] * (p.vsyn_ee - state[0]) + state[6] * (p.vsyn_ei - state[0]);
                auto i_i = state[8] * (p.vsyn_ee - state[1]) + state[10] * (p.vsyn_ei - state[1]);
                all_i_net[trial][save_counter] = i_e + i_i;
                save_counter++;
            }

            stepper.do_step(state, t);
            t += opts.dt;
        }

        size_t done = ++completed_trials;
        if (done % 1 == 0 || done == opts.n_trial) {
            print_progress(done, opts.n_trial);
        }
    }

    std::cout << std::endl;

    // Create directory if it doesn't exist (C++17)
    std::filesystem::create_directories(output_dir);

    std::string full_path = output_dir + "/" + filename; 
    save_matrix(all_i_net, full_path);
    std::cout << "Saved simulation to " << full_path << "\n";
    std::cout << num_steps << std::endl;
}
