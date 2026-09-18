# Pure-C translation of the QUADPACK routines regenie uses (adaptive integration
# with the 21-point Gauss-Kronrod rule). Produces SKAT-O p-values bit-identical
# to the Fortran original.

add_library(regenie_cquadpack STATIC
  dqags.c
  dqext.c
  dqsort.c
  gk21.c)
add_library(regenie::cquadpack ALIAS regenie_cquadpack)

target_include_directories(regenie_cquadpack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(regenie_cquadpack PRIVATE NDEBUG)

set_target_properties(regenie_cquadpack PROPERTIES POSITION_INDEPENDENT_CODE ON)
