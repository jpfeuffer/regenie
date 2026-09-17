# regenie build.
#
# Dependencies are resolved in cmake/Dependencies.cmake: find_package first,
# with downloads only for header-only libraries that ship no CMake package.
# See CMakePresets.json for ready-made configurations.

cmake_minimum_required(VERSION 3.25)

# Prefer upstream BoostConfig.cmake over CMake's removed FindBoost module.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# Get Regenie version
file(STRINGS "VERSION" RG_VERSION)

project(regenie
  VERSION ${RG_VERSION}
  LANGUAGES C CXX
  DESCRIPTION "Whole genome regression for large genome-wide association studies"
  )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)

# Set output directories for build artifacts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Ensures -std=c++17

######################################
######## options

# Remote input (s3:// and http(s)://) through BGEN::s3; curl-based SigV4.
option(WITH_S3 "Enable reading s3:// and http(s):// inputs" ON)

# Use the pure-C cquadpack integrator instead of the Fortran quadpack library
# (removes the gfortran build/link dependency)
option(USE_C_QUADPACK "Use pure C cquadpack instead of Fortran quadpack" OFF)

option(HAS_BOOST_IOSTREAM "Read and write gzipped files via Boost.Iostreams" OFF)
option(BUILD_STATIC "Link libgcc/libstdc++ statically (Linux)" OFF)

# -ffast-math implies -ffinite-math-only, which lets the compiler assume no NaN
# reaches the isnan() guards that regenie uses for missing data. Historically
# always on; keep it that way by default but make it switchable so the effect on
# p-values can be measured.
option(REGENIE_FAST_MATH "Compile with -ffast-math" ON)

# Historical environment-variable spellings, still honoured.
if($ENV{STATIC})
  set(BUILD_STATIC ON)
endif()
if($ENV{HAS_BOOST_IOSTREAM})
  set(HAS_BOOST_IOSTREAM ON)
endif()

if(HAS_BOOST_IOSTREAM)
  set(RG_VERSION "${RG_VERSION}.gz")
endif()

set(Boost_USE_STATIC_LIBS    ${BUILD_STATIC})
set(Boost_USE_DEBUG_LIBS     OFF)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BLA_STATIC               ${BUILD_STATIC})

include(Dependencies)

######################################
######## regenie executable

# list each file specifically
add_executable(regenie
  ${CMAKE_SOURCE_DIR}/src/Data.cpp
  ${CMAKE_SOURCE_DIR}/src/Files.cpp
  ${CMAKE_SOURCE_DIR}/src/Geno.cpp
  ${CMAKE_SOURCE_DIR}/src/HLM.cpp
  ${CMAKE_SOURCE_DIR}/src/Interaction.cpp
  ${CMAKE_SOURCE_DIR}/src/Joint_Tests.cpp
  ${CMAKE_SOURCE_DIR}/src/Masks.cpp
  ${CMAKE_SOURCE_DIR}/src/NNLS.cpp
  ${CMAKE_SOURCE_DIR}/src/Pheno.cpp
  ${CMAKE_SOURCE_DIR}/src/Regenie.cpp
  ${CMAKE_SOURCE_DIR}/src/SKAT.cpp
  ${CMAKE_SOURCE_DIR}/src/Step1_Models.cpp
  ${CMAKE_SOURCE_DIR}/src/Step2_Models.cpp
  ${CMAKE_SOURCE_DIR}/src/MultiTrait_Tests.cpp
  ${CMAKE_SOURCE_DIR}/src/MCC.cpp
  ${CMAKE_SOURCE_DIR}/src/Ordinal.cpp
  ${CMAKE_SOURCE_DIR}/src/survival_data.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_ridge.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_score.cpp
  ${CMAKE_SOURCE_DIR}/src/cox_firth.cpp
  ${CMAKE_SOURCE_DIR}/src/S3_Utils.cpp
  )
target_include_directories(regenie PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_compile_definitions(regenie PRIVATE
  VERSION_NUMBER="${RG_VERSION}"
  $<$<BOOL:${HAS_BOOST_IOSTREAM}>:HAS_BOOST_IOSTREAM>
  $<$<BOOL:${REGENIE_HAVE_S3}>:WITH_S3>
  $<$<BOOL:${REGENIE_HAVE_HTSLIB}>:WITH_HTSLIB>
  $<$<BOOL:${USE_C_QUADPACK}>:USE_C_QUADPACK>)

# -O3 comes from the Release build type; only the warning selection and the
# fast-math opt-in are project policy.
target_compile_options(regenie PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
    -Wall -pedantic
    -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-long-long
    $<$<BOOL:${REGENIE_FAST_MATH}>:-ffast-math>>
  $<$<CXX_COMPILER_ID:AppleClang,Clang>:-Wno-c11-extensions>)

set_target_properties(regenie PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(regenie PRIVATE
  Eigen3::Eigen
  cxxopts::cxxopts
  LBFGSpp::LBFGSpp
  Boost::filesystem
  Boost::headers
  ZLIB::ZLIB
  Threads::Threads)

if(OpenMP_CXX_FOUND)
  target_link_libraries(regenie PRIVATE OpenMP::OpenMP_CXX)
endif()

if(BUILD_STATIC AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_options(regenie BEFORE PRIVATE -static-libgcc -static-libstdc++)
endif()

set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")
target_include_directories(regenie PRIVATE ${EXTERN_LIBS_PATH}/)

target_link_libraries(regenie PRIVATE BGEN::bgen)

if(REGENIE_HAVE_S3)
  target_link_libraries(regenie PRIVATE BGEN::s3)
endif()


# MVTNorm library
set(MVTN_PATH "${EXTERN_LIBS_PATH}/mvtnorm")
add_custom_target(
   libMvtnorm
   COMMAND ${CMAKE_MAKE_PROGRAM}
   WORKING_DIRECTORY ${MVTN_PATH}
)
target_link_libraries(regenie PRIVATE ${MVTN_PATH}/libMvtnorm.a)
add_dependencies(regenie libMvtnorm)

# QF library
set(QF_PATH "${EXTERN_LIBS_PATH}/qf")
add_custom_target(
   libqf
   COMMAND ${CMAKE_MAKE_PROGRAM}
   WORKING_DIRECTORY ${QF_PATH}
)
target_link_libraries(regenie PRIVATE ${QF_PATH}/qf.a)
add_dependencies(regenie libqf)

# Quadpack library (Fortran quadpack, or pure-C cquadpack when USE_C_QUADPACK=ON)
if(USE_C_QUADPACK)
  set(QUAD_PATH "${EXTERN_LIBS_PATH}/cquadpack")
  add_custom_target(
     libquad
     COMMAND ${CMAKE_MAKE_PROGRAM}
     WORKING_DIRECTORY ${QUAD_PATH}
  )
  target_link_libraries(regenie PRIVATE ${QUAD_PATH}/libcquadpack.a)
  target_include_directories(regenie PRIVATE ${QUAD_PATH})
else()
  set(QUAD_PATH "${EXTERN_LIBS_PATH}/quadpack")
  add_custom_target(
     libquad
     COMMAND ${CMAKE_MAKE_PROGRAM}
     WORKING_DIRECTORY ${QUAD_PATH}
  )
  target_link_libraries(regenie PRIVATE ${QUAD_PATH}/libquad.a)
  # gfortran is only needed to link the Fortran quadpack library
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(regenie PRIVATE gfortran)
  else()
    find_library(GFORTRAN_LIBRARY gfortran REQUIRED)
    target_link_libraries(regenie PRIVATE ${GFORTRAN_LIBRARY})
  endif()
endif()
add_dependencies(regenie libquad)

# PGEN library
set(PGEN_PATH "${EXTERN_LIBS_PATH}/pgenlib")
add_custom_target(
   pgenlib
   COMMAND ${CMAKE_MAKE_PROGRAM} CXX=${CMAKE_CXX_COMPILER}
   WORKING_DIRECTORY ${PGEN_PATH}
)
target_link_libraries(regenie PRIVATE ${PGEN_PATH}/pgenlib.a)
target_include_directories(regenie PRIVATE ${PGEN_PATH} ${PGEN_PATH}/simde/ ${PGEN_PATH}/include/)
add_dependencies(regenie pgenlib)

# REMETA library
if(REGENIE_HAVE_HTSLIB)
  set(REMETA_PATH "${EXTERN_LIBS_PATH}/remeta")
  add_custom_target(
     remeta
     COMMAND ${CMAKE_MAKE_PROGRAM} HTSLIB_PATH=${REGENIE_HTSLIB_PATH}
     WORKING_DIRECTORY ${REMETA_PATH}
  )
  target_link_libraries(regenie PRIVATE ${REMETA_PATH}/remeta.a regenie::htslib)
  target_include_directories(regenie PRIVATE ${REMETA_PATH})
  add_dependencies(regenie remeta)
endif()

# BLAS backend for Eigen (resolved in cmake/Dependencies.cmake)
if(REGENIE_BLAS_BACKEND STREQUAL "mkl")
  target_compile_definitions(regenie PRIVATE WITH_MKL EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
  target_link_libraries(regenie PRIVATE MKL::MKL)
elseif(REGENIE_BLAS_BACKEND STREQUAL "openblas")
  target_compile_definitions(regenie PRIVATE WITH_OPENBLAS EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
  target_link_libraries(regenie PRIVATE ${LAPACKE_LIBRARY} LAPACK::LAPACK BLAS::BLAS)
endif()

if(HAS_BOOST_IOSTREAM)
  find_package(Boost REQUIRED COMPONENTS iostreams)
  target_link_libraries(regenie PRIVATE Boost::iostreams)
endif()

install(TARGETS regenie RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
message(STATUS "REGENIE v${RG_VERSION}")

add_custom_target(full-clean
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${MVTN_PATH}" clean
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${QF_PATH}" clean
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${QUAD_PATH}" clean
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${PGEN_PATH}" clean
  )
