#ifndef GENERATE_MODEL_STATS_H
#define GENERATE_MODEL_STATS_H
#include "thresholding/threshold_structs.h"
#include <string>

void generate_model_stats(const std::string &input_dir, const std::string &output_dir,
                          const std::string &input_filename, const std::string &output_filename,
                          const ThresholdParameters &params);

#endif