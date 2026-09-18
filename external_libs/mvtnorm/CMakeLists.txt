# Header-only C++ port of Alan Genz's MVTDST, plus the compiled pmvnorm wrapper.
# mvt.f / randomF77.c are the Fortran reference implementation and are only used
# by the upstream benchmarks, so they are deliberately not built here.

add_library(regenie_mvtnorm STATIC mvtnorm.cpp)
add_library(regenie::mvtnorm ALIAS regenie_mvtnorm)

target_include_directories(regenie_mvtnorm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(regenie_mvtnorm PRIVATE cxx_std_17)
target_link_libraries(regenie_mvtnorm PRIVATE Eigen3::Eigen)

if(OpenMP_CXX_FOUND)
  target_link_libraries(regenie_mvtnorm PRIVATE OpenMP::OpenMP_CXX)
endif()

# FMA and no-errno math, as the hand-written Makefile used. -march=native is
# deliberately not carried over: it produced binaries that only ran on the
# build machine.
target_compile_options(regenie_mvtnorm PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
    -ffp-contract=fast -fno-math-errno -fno-trapping-math>)
target_compile_definitions(regenie_mvtnorm PRIVATE NDEBUG)

set_target_properties(regenie_mvtnorm PROPERTIES POSITION_INDEPENDENT_CODE ON)
