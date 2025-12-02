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
set(MATLAB_INC ${MATLAB_ROOT}/extern/include)
set(MATLAB_BIN ${MATLAB_ROOT}/bin/glnxa64)
set(MATLAB_EXT ${MATLAB_ROOT}/extern/bin/glnxa64)

# Headers (global for simplicity)
include_directories(${MATLAB_INC})

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
add_compile_definitions(PYTHON_PLOTTING_DIR="${PYTHON_PLOTTING_DIR}")
# Make EXTERNAL_DIR a C++ string literal
add_compile_definitions(EXTERNAL_DIR=\"${EXTERNAL_DIR}\")
add_compile_definitions(DATA_DIR="${PROJECT_SOURCE_DIR}/data")
add_compile_definitions(OUTPUT_DIR="${PROJECT_SOURCE_DIR}/output")

# === Sources ===
file(GLOB_RECURSE SUPPORT_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# --- EXCLUDE MATLAB-DEPENDENT FILES FROM THE GENERAL POOL ---
list(FILTER SUPPORT_SOURCES EXCLUDE REGEX ".*/src/HMM/matlab_test\\.cpp$")
list(FILTER SUPPORT_SOURCES EXCLUDE REGEX ".*/src/HMM/run_hmm_analysis\\.cpp$")
list(FILTER SUPPORT_SOURCES EXCLUDE REGEX ".*/src/HMM/hmm_test\\.cpp$")

# === Shared Libraries (non-MATLAB) ===
set(SHARED_LIBS
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
    ${PROJECT_SOURCE_DIR}/external/iir1/libiir.so
)

# === Executable: ga_run (no MATLAB) ===
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

# === Executable: test_matlab (MATLAB-only files) ===
add_executable(test_matlab
    ${PROJECT_SOURCE_DIR}/main/test_matlab.cpp
    ${PROJECT_SOURCE_DIR}/src/HMM/matlab_test.cpp
)
target_include_directories(test_matlab PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(test_matlab PRIVATE
    ${SHARED_LIBS}
    ${MATLAB_EXT}/libMatlabDataArray.so
    ${MATLAB_EXT}/libMatlabEngine.so
    pthread
)
target_link_options(test_matlab PRIVATE
    -Wl,--no-as-needed
    -Wl,-rpath,${MATLAB_EXT}
    -Wl,-rpath,${MATLAB_BIN}
)

# === Executable: test_hmm (MATLAB-only files) ===
add_executable(test_hmm
    ${PROJECT_SOURCE_DIR}/main/test_hmm.cpp
    ${PROJECT_SOURCE_DIR}/src/HMM/hmm_test.cpp
)
target_include_directories(test_hmm PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(test_hmm PRIVATE
    ${SHARED_LIBS}
    ${MATLAB_EXT}/libMatlabDataArray.so
    ${MATLAB_EXT}/libMatlabEngine.so
    pthread
)
target_link_options(test_hmm PRIVATE
    -Wl,--no-as-needed
    -Wl,-rpath,${MATLAB_EXT}
    -Wl,-rpath,${MATLAB_BIN}
)


add_executable(ga_run_hmm
    ${PROJECT_SOURCE_DIR}/main/main_hmm.cpp
    ${SUPPORT_SOURCES}
    ${PROJECT_SOURCE_DIR}/src/HMM/hmm_test.cpp
)

target_include_directories(ga_run_hmm PRIVATE ${PROJECT_INCLUDES})
target_compile_definitions(ga_run_hmm PRIVATE
    DATA_DIR="${DATA_DIR}"
    OUTPUT_DIR="${OUTPUT_DIR}"
)
target_link_libraries(ga_run_hmm PRIVATE PRIVATE
    ${SHARED_LIBS}
    ${MATLAB_EXT}/libMatlabDataArray.so
    ${MATLAB_EXT}/libMatlabEngine.so
    pthread)

target_link_options(ga_run_hmm PRIVATE
    -Wl,--no-as-needed
    -Wl,-rpath,${MATLAB_EXT}
    -Wl,-rpath,${MATLAB_BIN}
)
