#include "HMM/hmm_test.hpp"
#include <string>

int main() {
    std::string popcurrent_path = std::string(DATA_DIR) + "/meg";
    std::string output_path = std::string(OUTPUT_DIR) + "/meg";
    int size_pop = 1;
    double sampling_freq = 600;
    int num_trials = 78;
    hmm_matlab_run(popcurrent_path, output_path, sampling_freq);

    return 0;
}
