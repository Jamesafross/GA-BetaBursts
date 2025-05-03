#ifndef EULER_MARUYAMA_HPP
#define EULER_MARUYAMA_HPP

#include <vector>
#include <random>
#include <optional>

template <typename System>
class EulerMaruyamaStepper {
    System& system;
    double dt;
    std::mt19937 rng;
    std::normal_distribution<double> normal{0.0, 1.0};

public:
    // Constructor with optional seed
    EulerMaruyamaStepper(System& sys, double dt_, std::optional<unsigned> seed = std::nullopt)
        : system(sys), dt(dt_) {
        if (seed) {
            rng.seed(*seed);
        } else {
            std::random_device rd;
            rng.seed(rd());
        }
    }

    
    void do_step(std::vector<double>& state, double t) {
        std::vector<double> dxdt(state.size());
        system(state, dxdt, t);

        for (std::size_t i = 0; i < state.size(); ++i) {
            state[i] += dxdt[i] * dt;
        }

        system.apply_noise(state, dt, rng);
    }
};

#endif // EULER_MARUYAMA_HPP
