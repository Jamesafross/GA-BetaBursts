#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <cstdlib>
#include <string>

#ifndef PYTHON_PLOTTING_DIR

#endif

void run_violin_plot(size_t generation, int phenotype_id);

int get_best_phenotype_id(const std::string &summary_path);
#endif // PLOTTING_HPP