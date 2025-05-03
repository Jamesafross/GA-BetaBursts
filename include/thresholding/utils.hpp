#pragma once

#include <vector>
#include <string>
#include <utility>  // for std::pair

// Load multi-region MEG data from a CSV file.
// Assumes each row is one region's full time series.
std::vector<std::vector<double>> load_region_signals(const std::string& filename);