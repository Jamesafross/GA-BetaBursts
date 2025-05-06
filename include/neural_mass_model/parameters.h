#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

constexpr unsigned int FIXED_SEED = 42;

struct ModelParameters {
    // Alpha parameters
    double alphaee;
    double alphaei;
    double alphaie;
    double alphaii;

    // Time constants and delays
    double taue;
    double taui;
    double deltae;
    double deltai;

    // Input drive
    double eta0e;
    double eta0i;

    // Synaptic couplings
    double kappa_see;
    double kappa_sei;
    double kappa_sie;
    double kappa_sii;

    // Gap junctions
    double kappa_gapee;
    double kappa_gapei;
    double kappa_gapii;

    // Reversal potentials
    double vsyn_ee;
    double vsyn_ei;
    double vsyn_ie;
    double vsyn_ii;

    // Noise timescales
    double tauxe;
    double tauxi;

    // Noise amplitudes
    double sigma_se;
    double sigma_si;
};

struct ParameterBounds {
    // === Alpha parameters ===
    double alphaee_min = 0.01, alphaee_max = 1.0;
    double alphaei_min = 0.01, alphaei_max = 1.0;
    double alphaie_min = 0.01, alphaie_max = 1.0;
    double alphaii_min = 0.01, alphaii_max = 1.0;

    // === Synaptic couplings ===
    double kappa_see_min = 0.01, kappa_see_max = 25.0;
    double kappa_sei_min = 0.01, kappa_sei_max = 25.0;
    double kappa_sie_min = 0.01, kappa_sie_max = 25.0;
    double kappa_sii_min = 0.01, kappa_sii_max = 25.0;

    // === Gap junction couplings ===
    double kappa_gapee_min = 0.0, kappa_gapee_max = 0.2;
    double kappa_gapei_min = 0.0, kappa_gapei_max = 0.5;
    double kappa_gapii_min = 0.0, kappa_gapii_max = 0.8;

    // === Reversal potentials ===
    double vsyn_ee_min = 1.0, vsyn_ee_max = 30.0;
    double vsyn_ei_min = 1.0, vsyn_ei_max = 30.0;
    double vsyn_ie_min = -100.0, vsyn_ie_max = -1.0;
    double vsyn_ii_min = -100.0, vsyn_ii_max = -1.0;

    // === Time constants ===
    double tauxe_min = 20.0, tauxe_max = 100.0;
    double tauxi_min = 20.0, tauxi_max = 100.0;
    double taue_min = 5.0, taue_max = 25.0;
    double taui_min = 5.0, taui_max = 25.0;

    // === HWHM ===
    double deltae_min = 0.05, deltae_max = 0.6;
    double deltai_min = 0.05, deltai_max = 0.6;

    // === Noise amplitudes ===
    double sigma_se_min = 0.001, sigma_se_max = 0.1;
    double sigma_si_min = 0.001, sigma_si_max = 0.1;

    // === Mean inputs ===
    double eta0e_min = 0.0, eta0e_max = 5.0;
    double eta0i_min = 0.0, eta0i_max = 5.0;
};

struct SolverParameters {
    double dt = 0.01;                 // time step in milliseconds (ms)
    double t_seconds = 100.0;         // total simulation time in seconds
    double save_start_seconds = 20.0; // start saving after this many seconds
    size_t n_trial = 200;

    double fs = 200.0; // desired saving frequency in Hz

    size_t steps_per_second() const {
        return static_cast<size_t>(1000.0 / dt); // ms to steps
    }

    size_t save_every() const {
        return static_cast<size_t>(1000.0 / (dt * fs)); // save every N steps
    }

    size_t t_end() const { return static_cast<size_t>(t_seconds * steps_per_second()); }

    size_t save_start_step() const {
        return static_cast<size_t>(save_start_seconds * steps_per_second());
    }
};

#endif // PARAMETERS_H