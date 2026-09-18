# LD-matrix / summary output for the REMETA meta-analysis format. Only built
# when HTSlib is available.

add_library(regenie_remeta STATIC
  bgz_writer.cpp
  regenie_ld_matrix_writer.cpp)
add_library(regenie::remeta ALIAS regenie_remeta)

target_include_directories(regenie_remeta PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(regenie_remeta PRIVATE cxx_std_17)
target_compile_definitions(regenie_remeta PUBLIC WITH_HTSLIB)
target_link_libraries(regenie_remeta PUBLIC regenie::htslib PRIVATE Eigen3::Eigen)

set_target_properties(regenie_remeta PROPERTIES POSITION_INDEPENDENT_CODE ON)
