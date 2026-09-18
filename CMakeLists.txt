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

# Imported targets' include dirs are passed with -I, not -isystem. -isystem
# paths are searched after /usr/local/include, so a stale library installed
# there silently shadows the one find_package actually resolved.
set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)

# Get Regenie version
file(STRINGS "VERSION" RG_VERSION)

project(regenie
  VERSION ${RG_VERSION}
  LANGUAGES C CXX
  DESCRIPTION "Whole genome regression for large genome-wide association studies"
  )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)

# MSVC's <cmath>/<math.h> only define M_PI/M_SQRT1_2/etc. (BSD/POSIX
# extensions, not standard C/C++) when this is defined before their first
# inclusion. A directory-scope definition applies before any translation
# unit in this build is even opened, unlike an in-file #define, which would
# depend on include order in each of the several files that use them. A
# no-op on every other compiler/platform.
add_compile_definitions(_USE_MATH_DEFINES)

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

# Reads bgenix-generated .bgi index files, for compatibility with pipelines
# that already have one. bgen-limix's own metafile is preferred whenever both
# exist, and is what regenie itself creates; this only reads .bgi, never
# writes it. Off means such an input falls back to the metafile/scan path.
option(WITH_BGI "Read legacy .bgi index files via system SQLite3" ON)

# The pure-C cquadpack integrator produces bit-identical SKAT-O p-values to the
# Fortran quadpack it replaces, and removes the gfortran build/link dependency.
# Set OFF to fall back to Fortran quadpack.
option(USE_C_QUADPACK "Use pure C cquadpack instead of Fortran quadpack" ON)

option(HAS_BOOST_IOSTREAM "Read and write gzipped files via Boost.Iostreams" OFF)
option(BUILD_STATIC "Link libgcc/libstdc++ statically (Linux)" OFF)

# qfc.cpp keeps upstream's file-scope state; qf_mt() is reentrant.
option(WITH_THREADSAFE_QF "Use the thread-safe qf_mt() in SKAT-O" OFF)

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
  ${CMAKE_SOURCE_DIR}/src/BgenReader.cpp
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
  $<$<BOOL:${REGENIE_HAVE_BGI}>:WITH_BGI>
  $<$<BOOL:${REGENIE_HAVE_HTSLIB}>:WITH_HTSLIB>
  $<$<BOOL:${WITH_THREADSAFE_QF}>:WITH_THREADSAFE_QF>
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
target_include_directories(regenie PRIVATE
  ${EXTERN_LIBS_PATH}/
  ${EXTERN_LIBS_PATH}/mvtnorm/)

target_link_libraries(regenie PRIVATE BGEN::bgen)

if(REGENIE_HAVE_S3)
  target_link_libraries(regenie PRIVATE BGEN::s3)
endif()

if(REGENIE_HAVE_BGI)
  target_link_libraries(regenie PRIVATE SQLite::SQLite3)
endif()


# Vendored numerical kernels and readers. These have no upstream packaging, so
# they stay in-tree, but they are built as ordinary CMake targets rather than by
# shelling out to their Makefiles.
add_subdirectory(${EXTERN_LIBS_PATH}/mvtnorm)
add_subdirectory(${EXTERN_LIBS_PATH}/qf)
add_subdirectory(${EXTERN_LIBS_PATH}/pgenlib)

target_link_libraries(regenie PRIVATE
  regenie::mvtnorm
  regenie::qf
  regenie::pgenlib)

if(USE_C_QUADPACK)
  add_subdirectory(${EXTERN_LIBS_PATH}/cquadpack)
  target_link_libraries(regenie PRIVATE regenie::cquadpack)
else()
  add_subdirectory(${EXTERN_LIBS_PATH}/quadpack)
  target_link_libraries(regenie PRIVATE regenie::quadpack)
endif()

if(REGENIE_HAVE_HTSLIB)
  add_subdirectory(${EXTERN_LIBS_PATH}/remeta)
  target_link_libraries(regenie PRIVATE regenie::remeta)
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
