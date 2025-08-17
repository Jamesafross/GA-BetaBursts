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

# NLopt (found but not linked here; keep if you plan to use later)
find_path(NLOPT_INCLUDE_DIR nlopt.hpp PATHS /usr/include /usr/local/include)
find_library(NLOPT_LIBRARY nlopt PATHS /usr/lib /usr/local/lib)

# === MATLAB C++ Engine (R2024b) ===
set(MATLAB_ROOT /usr/local/MATLAB/R2024b)

# Headers
include_directories(
    ${MATLAB_ROOT}/extern/include
)

# Lib dirs: C++ engine libs live in extern/bin; keep bin as well
link_directories(
    ${MATLAB_ROOT}/extern/bin/glnxa64
    ${MATLAB_ROOT}/bin/glnxa64
)

# C++ Engine libs (link per-target)
set(MATLAB_CPP_LIBS
    MatlabEngine
    MatlabDataArray
)



# === Global Include Paths ===
set(PROJECT_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external
    ${PROJECT_SOURCE_DIR}/external/iir1
    /usr/include/eigen3
)

# === Global Definitions ===
set(DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(OUTPUT_DIR "${PROJECT_SOURCE_DIR}/output")
set(EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external")
set(PYTHON_PLOTTING_DIR "${CMAKE_SOURCE_DIR}/python_plotting")
add_definitions(-DPYTHON_PLOTTING_DIR="${PYTHON_PLOTTING_DIR}")

# === Sources ===
file(GLOB_RECURSE SUPPORT_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# === Shared Libraries (non-MATLAB) ===
set(SHARED_LIBS
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
    ${PROJECT_SOURCE_DIR}/external/iir1/libiir.so
)

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
target_link_libraries(ga_run PRIVATE ${SHARED_LIBS})



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
target_link_libraries(meg_burst_stats PRIVATE ${SHARED_LIBS})

# === Executable: test_hmm ===
add_executable(test_hmm
    ${PROJECT_SOURCE_DIR}/main/test_hmm.cpp
    ${SUPPORT_SOURCES}  # keep if you rely on shared src
)
target_include_directories(test_hmm PRIVATE ${PROJECT_INCLUDES})
target_compile_definitions(test_hmm PRIVATE
    DATA_DIR="${DATA_DIR}"
    OUTPUT_DIR="${OUTPUT_DIR}"
)
target_link_libraries(test_hmm PRIVATE ${SHARED_LIBS} ${MATLAB_CPP_LIBS})
set_target_properties(test_hmm PROPERTIES
    BUILD_RPATH "${MATLAB_ROOT}/extern/bin/glnxa64;${MATLAB_ROOT}/bin/glnxa64"
)

# === Executable: test_matlab ===
add_executable(test_matlab
    ${PROJECT_SOURCE_DIR}/main/test_matlab.cpp
    ${PROJECT_SOURCE_DIR}/src/HMM/matlab_test.cpp
)
target_include_directories(test_matlab PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(test_matlab PRIVATE ${SHARED_LIBS} ${MATLAB_CPP_LIBS})
set_target_properties(test_matlab PROPERTIES
    BUILD_RPATH "${MATLAB_ROOT}/extern/bin/glnxa64;${MATLAB_ROOT}/bin/glnxa64"
)

target_compile_definitions(test_matlab PRIVATE
    EXTERNAL_DIR="${EXTERNAL_DIR}"
)

