# Original Fortran QUADPACK, kept as a fallback for USE_C_QUADPACK=OFF. Only
# added to the build in that case, so a default build needs no Fortran compiler.

enable_language(Fortran)

add_library(regenie_quadpack STATIC
  d1mach.f
  dqags.f
  dqagse.f
  dqelg.f
  dqk21.f
  dqpsrt.f)
add_library(regenie::quadpack ALIAS regenie_quadpack)

# Upstream is fixed-form FORTRAN 77 and trips modern gfortran's default checks.
target_compile_options(regenie_quadpack PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:-std=legacy>)

set_target_properties(regenie_quadpack PROPERTIES POSITION_INDEPENDENT_CODE ON)
