cmake_minimum_required(VERSION 3.16)
project(mvtdst
  VERSION     1.0.0
  LANGUAGES   CXX
  DESCRIPTION "Header-only C++ port of Alan Genz's MVTDST (multivariate normal/t CDF)"
)

# ── C++17 ────────────────────────────────────────────────────────────────────
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 17)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# ── Interface target — header-only C++ API ───────────────────────────────────
# Consumers that include mvtdst.hpp directly (C++ code) link this target.
add_library(mvtdst_headers INTERFACE)
add_library(mvtdst::headers ALIAS mvtdst_headers)

target_include_directories(mvtdst_headers INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(mvtdst_headers INTERFACE cxx_std_17)

# ── Optional: OpenMP for parallel QMC (strong scaling ~8× for n≥10) ──────────
# On macOS, Apple Clang ships without libomp; help CMake find Homebrew's copy.
if(APPLE AND NOT OpenMP_CXX_FOUND)
  find_program(_brew brew)
  if(_brew)
    execute_process(
      COMMAND ${_brew} --prefix libomp
      OUTPUT_VARIABLE _libomp_prefix
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
  endif()
  if(_libomp_prefix AND EXISTS "${_libomp_prefix}/include/omp.h")
    set(OpenMP_CXX_FLAGS         "-Xpreprocessor -fopenmp -I${_libomp_prefix}/include")
    set(OpenMP_CXX_LIB_NAMES     "omp")
    set(OpenMP_omp_LIBRARY        "${_libomp_prefix}/lib/libomp.dylib")
  endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(mvtdst_headers INTERFACE OpenMP::OpenMP_CXX)
  message(STATUS "mvtdst: OpenMP found — parallel integration enabled")
else()
  message(STATUS "mvtdst: OpenMP not found — single-threaded mode")
endif()

# ── Compiled C-compatible static library ─────────────────────────────────────
# Provides the pmvnorm / pmvnorm_complement C API declared in mvtnorm.h.
add_library(mvtnorm STATIC mvtnorm.cpp)
add_library(mvtdst::mvtnorm ALIAS mvtnorm)

# PUBLIC: downstream C++ code that links mvtnorm gets the include path + OMP.
target_link_libraries(mvtnorm PUBLIC mvtdst_headers)
target_include_directories(mvtnorm PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# ── Optional fast-math flags ─────────────────────────────────────────────────
# PRIVATE: flags apply only when compiling mvtnorm.cpp, not to downstream code.
# -ffp-contract=fast   FMA for Horner chains in mvphnv  (~30% speedup)
# -fno-math-errno      removes libm errno overhead; lets erfc inline
# -fno-trapping-math   removes FP-trap checks; safe for HPC use
option(MVTDST_FAST_MATH
  "Enable FMA + no-errno math flags for ~2× speedup (PRIVATE to the library)"
  ON
)
if(MVTDST_FAST_MATH)
  target_compile_options(mvtnorm PRIVATE
    -ffp-contract=fast -fno-math-errno -fno-trapping-math
  )
endif()

# ── Tests ────────────────────────────────────────────────────────────────────
# Built automatically when this is the top-level project; opt-in otherwise.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  option(MVTDST_BUILD_TESTS "Build correctness / equivalence tests" ON)
else()
  option(MVTDST_BUILD_TESTS "Build correctness / equivalence tests" OFF)
endif()

if(MVTDST_BUILD_TESTS)
  enable_testing()
  add_executable(test_mvtdst test_mvtdst.cpp)
  target_link_libraries(test_mvtdst PRIVATE mvtnorm)
  add_test(NAME mvtdst COMMAND test_mvtdst)
endif()

# ── Benchmarks ───────────────────────────────────────────────────────────────
# Off by default (requires gfortran).  Enable with -DMVTDST_BUILD_BENCHMARKS=ON.
option(MVTDST_BUILD_BENCHMARKS
  "Build Fortran comparison benchmarks (requires gfortran + C compiler)" OFF
)
if(MVTDST_BUILD_BENCHMARKS)
  enable_language(C Fortran)
  add_library(mvt_fortran_ref STATIC mvt.f randomF77.c)
  add_executable(benchmark_mvtdst benchmark_mvtdst.cpp)
  target_link_libraries(benchmark_mvtdst PRIVATE mvtnorm mvt_fortran_ref)
endif()
