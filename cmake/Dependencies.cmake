# Central dependency resolution for regenie.
#
# Policy: every external dependency is looked up with find_package first, so an
# existing system / conda / homebrew copy is reused and can be a shared library.
# Only header-only dependencies that ship no CMake package are downloaded, and
# then into the build tree -- never vendored into the source tree.
#
# Configure with -DREGENIE_FETCH_DEPS=OFF to forbid downloads entirely, which is
# what distribution packagers want: the configure step then fails with the name
# of the missing package instead of silently reaching for the network.

include(FetchContent)

option(REGENIE_FETCH_DEPS "Download header-only dependencies that are not installed" ON)

# Exposes a header-only dependency that has no installed CMake package as an
# INTERFACE target. SOURCE_SUBDIR deliberately names a directory that does not
# exist, so FetchContent populates the sources without ever add_subdirectory()ing
# the upstream CMakeLists -- which would otherwise drag in that project's tests,
# docs and cache variables (Eigen is the notorious example).
function(_regenie_fetch_headers name target repo tag include_subdir)

  if(NOT REGENIE_FETCH_DEPS)
    message(FATAL_ERROR
      "${name} was not found and REGENIE_FETCH_DEPS=OFF.\n"
      "Install it and point CMAKE_PREFIX_PATH at it, or allow downloads with "
      "-DREGENIE_FETCH_DEPS=ON.")
  endif()

  message(STATUS "regenie: ${name} not installed -- fetching ${tag}")

  FetchContent_Declare(${name}
    GIT_REPOSITORY ${repo}
    GIT_TAG        ${tag}
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  headers-only-no-cmake-entry-point)
  FetchContent_MakeAvailable(${name})

  string(TOLOWER ${name} _lc)
  add_library(${target} INTERFACE IMPORTED GLOBAL)
  set_target_properties(${target} PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${${_lc}_SOURCE_DIR}/${include_subdir}")

endfunction()

######################################
# Linear algebra and numerics

find_package(Eigen3 3.4 CONFIG QUIET)
if(NOT Eigen3_FOUND)
  _regenie_fetch_headers(Eigen3 Eigen3::Eigen
    https://gitlab.com/libeigen/eigen.git 3.4.0 "")
endif()

######################################
# Header-only utilities

find_package(cxxopts CONFIG QUIET)
if(NOT cxxopts_FOUND)
  _regenie_fetch_headers(cxxopts cxxopts::cxxopts
    https://github.com/jarro2783/cxxopts.git v3.2.0 include)
endif()

# LBFGSpp ships no CMake package in any distribution, so this one is normally
# fetched. Headers live at <root>/include/{LBFGS.h,LBFGSpp/*}.
find_package(LBFGSpp CONFIG QUIET)
if(NOT LBFGSpp_FOUND)
  _regenie_fetch_headers(LBFGSpp LBFGSpp::LBFGSpp
    https://github.com/yixuan/LBFGSpp.git v0.4.0 include)
endif()

######################################
# Boost: math and exception are header-only. Filesystem is no longer needed --
# std::filesystem (C++17) replaced it -- which also drops the one compiled,
# non-header-only Boost component this project used, and with it a class of
# Windows link failures (DLL import-lib naming/availability for Boost's
# compiled filesystem library) that came with it.

find_package(Boost 1.71 REQUIRED)

######################################
# Compression and threading

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Apple Clang ships without libomp; point CMake at the homebrew copy so the
# find_package below succeeds on macOS the same way it does on Linux.
if(APPLE AND NOT DEFINED OpenMP_CXX_FLAGS)
  find_program(_brew brew)
  if(_brew)
    execute_process(COMMAND ${_brew} --prefix libomp
      OUTPUT_VARIABLE _libomp_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  endif()
  if(_libomp_prefix AND EXISTS "${_libomp_prefix}/include/omp.h")
    set(OpenMP_CXX_FLAGS     "-Xpreprocessor -fopenmp -I${_libomp_prefix}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY   "${_libomp_prefix}/lib/libomp.dylib")
  endif()
endif()

find_package(OpenMP COMPONENTS CXX)
if(NOT OpenMP_CXX_FOUND)
  message(WARNING "regenie: OpenMP not found -- building single-threaded")
endif()

######################################
# BGEN reader and remote input, both from jpfeuffer/bgen-limix.
#
# BGEN::bgen resolves its own zlib-ng/zstd dependencies, and BGEN::s3 does
# curl-based SigV4. Variant metadata itself comes from bgen-limix's own
# metafile format, which needs neither zstd nor SQLite from regenie; SQLite is
# pulled in separately, only for WITH_BGI, to read (not write) the legacy .bgi
# format that older pipelines already have lying around.

find_package(bgen 4.7.0 CONFIG REQUIRED)

set(REGENIE_HAVE_S3 OFF)
if(WITH_S3)
  if(TARGET BGEN::s3)
    set(REGENIE_HAVE_S3 ON)
    # pgenlib carries its own copy of the same curl-based reader.
    find_package(CURL REQUIRED)
    message(STATUS "regenie: remote input enabled (BGEN::s3)")
  else()
    message(WARNING
      "WITH_S3=ON but the bgen package does not provide BGEN::s3; building "
      "without remote input. Rebuild bgen-limix with BGEN_ENABLE_S3=ON, or "
      "configure regenie with -DWITH_S3=OFF.")
  endif()
endif()

set(REGENIE_HAVE_BGI OFF)
if(WITH_BGI)
  find_package(SQLite3 REQUIRED)
  set(REGENIE_HAVE_BGI ON)
  message(STATUS "regenie: .bgi reading enabled (SQLite3)")
endif()

######################################
# HTSlib, for the REMETA meta-analysis output. Optional.

set(REGENIE_HAVE_HTSLIB OFF)
set(REGENIE_HTSLIB_PATH "$ENV{HTSLIB_PATH}" CACHE PATH
  "Path to an HTSlib install (enables REMETA LD-matrix output)")

if(REGENIE_HTSLIB_PATH OR WITH_HTSLIB)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND AND NOT REGENIE_HTSLIB_PATH)
    pkg_check_modules(HTSLIB IMPORTED_TARGET htslib)
  endif()
  if(TARGET PkgConfig::HTSLIB)
    add_library(regenie::htslib ALIAS PkgConfig::HTSLIB)
    set(REGENIE_HAVE_HTSLIB ON)
  else()
    find_library(HTSLIB_LIBRARY NAMES hts HINTS "${REGENIE_HTSLIB_PATH}")
    if(HTSLIB_LIBRARY)
      add_library(regenie_htslib INTERFACE)
      add_library(regenie::htslib ALIAS regenie_htslib)
      target_link_libraries(regenie_htslib INTERFACE ${HTSLIB_LIBRARY})
      target_include_directories(regenie_htslib INTERFACE "${REGENIE_HTSLIB_PATH}/../include")
      set(REGENIE_HAVE_HTSLIB ON)
    endif()
  endif()
  if(NOT REGENIE_HAVE_HTSLIB)
    message(WARNING "regenie: HTSlib requested but not found -- REMETA output disabled")
  endif()
endif()

######################################
# BLAS / LAPACK backend for Eigen. Optional: without one, Eigen uses its own
# kernels, which is the default and is fine.

set(REGENIE_BLAS_BACKEND "none" CACHE STRING "BLAS backend: none, mkl or openblas")
set_property(CACHE REGENIE_BLAS_BACKEND PROPERTY STRINGS none mkl openblas)

# Preserve the historical environment-variable spellings.
if(REGENIE_BLAS_BACKEND STREQUAL "none")
  if(NOT "$ENV{MKLROOT}" STREQUAL "")
    set(REGENIE_BLAS_BACKEND "mkl")
  elseif(NOT "$ENV{OPENBLAS_ROOT}" STREQUAL "")
    set(REGENIE_BLAS_BACKEND "openblas")
  endif()
endif()

if(REGENIE_BLAS_BACKEND STREQUAL "mkl")
  find_package(MKL CONFIG REQUIRED)
  message(STATUS "regenie: BLAS backend = Intel MKL")
elseif(REGENIE_BLAS_BACKEND STREQUAL "openblas")
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
  message(STATUS "regenie: BLAS backend = OpenBLAS")
endif()
