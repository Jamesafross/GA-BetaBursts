#pragma once

#include <string>

/**
 * @brief Calls the MATLAB function 'print_input' in the external/ folder.
 *
 * This uses the MATLAB Engine API for C++ to start MATLAB, add the external folder
 * to the MATLAB path, and call the function with the given message.
 *
 * @param message The message to send to MATLAB for printing.
 */
void matlab_print_input(const std::string &message);
