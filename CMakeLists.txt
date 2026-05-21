cmake_minimum_required(VERSION 3.25)

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 1)
project(regenie VERSION ${_version} LANGUAGES CXX C Fortran)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ─── Options ─────────────────────────────────────────────────────────────────
option(WITH_AWS_S3         "Build with AWS S3 support"  ON)
option(WITH_MKL            "Build with Intel MKL"       OFF)
option(WITH_OPENBLAS       "Build with OpenBLAS"        OFF)
option(WITH_BOOST_IOSTREAM "Build with Boost IOstreams" OFF)
option(WITH_HTSLIB         "Build with HTSlib/REMETA"   OFF)
option(BUILD_STATIC        "Link C/C++ runtime statically (GCC only)" OFF)

# ─── Dependencies ───────────────────────────────────────────────────────────
find_package(bgen REQUIRED)
find_package(ZLIB REQUIRED)
find_package(cxxopts REQUIRED)
find_package(OpenMP)

# Eigen 5.x has breaking API changes; pin to 3.4.x
find_package(Eigen3 3.4 NO_MODULE QUIET)
if(NOT Eigen3_FOUND OR Eigen3_VERSION VERSION_GREATER_EQUAL "5.0")
  message(STATUS "Fetching Eigen 3.4.0 (system version ${Eigen3_VERSION} incompatible)")
  include(FetchContent)
  FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(Eigen3)
endif()

if(WITH_AWS_S3)
  find_package(AWSSDK REQUIRED COMPONENTS s3)
endif()
if(WITH_MKL)
  find_package(MKL REQUIRED)
elseif(WITH_OPENBLAS)
  find_package(OpenBLAS REQUIRED)
endif()
if(WITH_BOOST_IOSTREAM)
  find_package(Boost REQUIRED COMPONENTS iostreams)
endif()

# ─── Regenie executable ─────────────────────────────────────────────────────
add_executable(regenie
  src/Data.cpp
  src/Files.cpp
  src/Geno.cpp
  src/HLM.cpp
  src/Interaction.cpp
  src/Joint_Tests.cpp
  src/Masks.cpp
  src/NNLS.cpp
  src/Pheno.cpp
  src/Regenie.cpp
  src/SKAT.cpp
  src/Step1_Models.cpp
  src/Step2_Models.cpp
  src/MultiTrait_Tests.cpp
  src/MCC.cpp
  src/Ordinal.cpp
  src/survival_data.cpp
  src/cox_ridge.cpp
  src/cox_score.cpp
  src/cox_firth.cpp
  src/S3_Utils.cpp
)

target_include_directories(regenie PRIVATE src)
target_compile_definitions(regenie PRIVATE VERSION_NUMBER="${PROJECT_VERSION}")
target_compile_options(regenie PRIVATE
  -O3 -Wall -pedantic -ffast-math
  -Wno-unused-local-typedefs -Wno-deprecated-declarations
  -Wno-long-long -Wno-c11-extensions
)

# -static-libgcc/-static-libstdc++ are GCC-only; macOS/AppleClang links its runtime from the OS
target_link_options(regenie BEFORE PRIVATE
  $<$<AND:$<BOOL:${BUILD_STATIC}>,$<CXX_COMPILER_ID:GNU>>:-static-libgcc;-static-libstdc++>
)

# ─── Link dependencies ──────────────────────────────────────────────────────
set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")

target_link_libraries(regenie PRIVATE
  bgen::bgen
  ZLIB::ZLIB
  Eigen3::Eigen
  cxxopts::cxxopts
  m
  ${CMAKE_DL_LIBS}
  $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# Link Fortran runtime libs from the detected gfortran install, skipping
# GCC-internal stubs (emutls_w, heapt_w) that don't exist as real libraries.
foreach(_flib IN LISTS CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
  find_library(_fortran_lib_${_flib} NAMES ${_flib}
    HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
    NO_DEFAULT_PATH)
  if(_fortran_lib_${_flib})
    target_link_libraries(regenie PRIVATE ${_fortran_lib_${_flib}})
  endif()
  unset(_fortran_lib_${_flib} CACHE)
endforeach()

target_include_directories(regenie PRIVATE
  ${EXTERN_LIBS_PATH}/
  ${EXTERN_LIBS_PATH}/LBFGSpp/include/
)

# ─── Bundled static libraries ────────────────────────────────────────────────
set(QF_PATH   "${EXTERN_LIBS_PATH}/qf")
set(QUAD_PATH "${EXTERN_LIBS_PATH}/quadpack")
set(PGEN_PATH "${EXTERN_LIBS_PATH}/pgenlib")

# mvtdst: proper CMake subdirectory — no separate make step needed
add_subdirectory("${EXTERN_LIBS_PATH}/mvtdst" "${CMAKE_BINARY_DIR}/extern/mvtdst")

add_custom_target(libqf   COMMAND make CXX=${CMAKE_CXX_COMPILER} WORKING_DIRECTORY ${QF_PATH})
add_custom_target(libquad COMMAND make FC=${CMAKE_Fortran_COMPILER} WORKING_DIRECTORY ${QUAD_PATH})

if(WITH_AWS_S3)
  add_custom_target(pgenlib
    COMMAND ${CMAKE_MAKE_PROGRAM} WITH_S3=1 AWS_INCLUDES=-I${AWSSDK_INCLUDE_DIRS} CXX=${CMAKE_CXX_COMPILER}
    WORKING_DIRECTORY ${PGEN_PATH})
else()
  add_custom_target(pgenlib COMMAND make WORKING_DIRECTORY ${PGEN_PATH})
endif()

add_dependencies(regenie libqf libquad pgenlib)
target_link_libraries(regenie PRIVATE
  mvtnorm
  ${QF_PATH}/qf.a
  ${QUAD_PATH}/libquad.a
  ${PGEN_PATH}/pgenlib.a
)
target_include_directories(regenie PRIVATE ${PGEN_PATH} ${PGEN_PATH}/simde/ ${PGEN_PATH}/include/)

# ─── Optional: HTSlib / REMETA ───────────────────────────────────────────────
if(WITH_HTSLIB)
  find_package(PkgConfig REQUIRED)
  if(DEFINED HTSLIB_PATH)
    set(ENV{PKG_CONFIG_PATH} "${HTSLIB_PATH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  endif()
  pkg_check_modules(htslib REQUIRED IMPORTED_TARGET htslib)

  set(REMETA_PATH "${EXTERN_LIBS_PATH}/remeta")
  add_custom_target(remeta
    COMMAND ${CMAKE_MAKE_PROGRAM} HTSLIB_PATH=1 "CPPFLAGS=-I${htslib_INCLUDE_DIRS}"
    WORKING_DIRECTORY ${REMETA_PATH})
  add_dependencies(regenie remeta)

  target_compile_definitions(regenie PRIVATE WITH_HTSLIB)
  target_include_directories(regenie PRIVATE ${REMETA_PATH})
  target_link_libraries(regenie PRIVATE ${REMETA_PATH}/remeta.a PkgConfig::htslib)
endif()

# ─── Optional: BLAS ──────────────────────────────────────────────────────────
if(WITH_MKL)
  target_compile_definitions(regenie PRIVATE WITH_MKL EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
  target_link_libraries(regenie PRIVATE MKL::MKL)
elseif(WITH_OPENBLAS)
  target_compile_definitions(regenie PRIVATE WITH_OPENBLAS EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
  target_link_libraries(regenie PRIVATE OpenBLAS::OpenBLAS)
endif()

# ─── Optional: Boost IOstreams ───────────────────────────────────────────────
if(WITH_BOOST_IOSTREAM)
  target_compile_definitions(regenie PRIVATE HAS_BOOST_IOSTREAM)
  target_link_libraries(regenie PRIVATE Boost::iostreams)
  message(STATUS "Boost Iostreams enabled")
endif()

# ─── Optional: AWS S3 ────────────────────────────────────────────────────────
if(WITH_AWS_S3)
  target_compile_definitions(regenie PRIVATE WITH_AWS_S3)
  target_link_libraries(regenie PRIVATE ${AWSSDK_LINK_LIBRARIES})
  message(STATUS "AWS S3 support enabled")
endif()

# ─── Install ─────────────────────────────────────────────────────────────────
install(TARGETS regenie RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

message(STATUS "REGENIE v${PROJECT_VERSION}")
