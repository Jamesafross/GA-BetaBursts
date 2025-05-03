cmake_minimum_required(VERSION 3.14)
project(BetaBurstsGeneticAlgorithm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Compile flags ===
add_compile_options(-Wall -Wextra -O2)

# === External libraries ===

# FFTW via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

# OpenMP
find_package(OpenMP REQUIRED)

# === Include directories ===
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external/iir1
    ${PROJECT_SOURCE_DIR}/external
)

# === Shared definitions ===
set(DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(OUTPUT_DIR "${PROJECT_SOURCE_DIR}/output")

# === Shared sources ===
file(GLOB_RECURSE SUPPORT_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# === Executable: ga_run ===
add_executable(ga_run
    ${PROJECT_SOURCE_DIR}/main/main.cpp
    ${SUPPORT_SOURCES}
)

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

target_compile_definitions(meg_burst_stats PRIVATE
    DATA_DIR="${DATA_DIR}"
    OUTPUT_DIR="${OUTPUT_DIR}"
)

target_link_libraries(meg_burst_stats PRIVATE
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
    ${PROJECT_SOURCE_DIR}/external/iir1/libiir.so
)
