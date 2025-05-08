cmake_minimum_required(VERSION 3.14)
project(BetaBurstsGeneticAlgorithm)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -O2)

# === External Libraries ===

# FFTW
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

# OpenMP
find_package(OpenMP REQUIRED)

# NLopt
find_path(NLOPT_INCLUDE_DIR nlopt.hpp PATHS /usr/include /usr/local/include)
find_library(NLOPT_LIBRARY nlopt PATHS /usr/lib /usr/local/lib)

# === Global Include Paths ===
set(PROJECT_INCLUDES
    ${PROJECT_SOURCE_DIR}/include                 # your own headers
    ${PROJECT_SOURCE_DIR}/external
    ${PROJECT_SOURCE_DIR}/external/iir1
    /usr/include/eigen3                           # Eigen headers (header-only)
)

# === Global Definitions ===
set(DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(OUTPUT_DIR "${PROJECT_SOURCE_DIR}/output")
set(PYTHON_PLOTTING_DIR "${CMAKE_SOURCE_DIR}/python_plotting")
add_definitions(-DPYTHON_PLOTTING_DIR="${PYTHON_PLOTTING_DIR}")

# === Sources ===
file(GLOB_RECURSE SUPPORT_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# === Executable: ga_run ===
add_executable(ga_run
    ${PROJECT_SOURCE_DIR}/main/main.cpp
    ${SUPPORT_SOURCES}
)

target_include_directories(ga_run PRIVATE ${PROJECT_INCLUDES})
target_compile_definitions(ga_run PRIVATE
    DATA_DIR="${DATA_DIR}"
    OUTPUT_DIR="${OUTPUT_DIR}"
)
target_link_libraries(ga_run PRIVATE
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
    ${PROJECT_SOURCE_DIR}/external/iir1/libiir.so
)

# === Executable: meg_burst_stats ===
add_executable(meg_burst_stats
    ${PROJECT_SOURCE_DIR}/main/meg_stats.cpp
    ${SUPPORT_SOURCES}
)

target_include_directories(meg_burst_stats PRIVATE ${PROJECT_INCLUDES})
target_compile_definitions(meg_burst_stats PRIVATE
    DATA_DIR="${DATA_DIR}"
    OUTPUT_DIR="${OUTPUT_DIR}"
)
target_link_libraries(meg_burst_stats PRIVATE
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
    ${PROJECT_SOURCE_DIR}/external/iir1/libiir.so
)
