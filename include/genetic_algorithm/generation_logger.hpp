#pragma once

#include "selection.hpp"
#include <json.hpp>
#include <string>

class GenerationLogger {
  public:
    GenerationLogger(const std::string &summary_file_path);

    void append_summary_from_directory(const std::string &directory, int generation);

  private:
    std::string summary_path;
};