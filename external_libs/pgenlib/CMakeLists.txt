# PGEN reader/writer, vendored from plink-ng. Intended to be replaced by a
# find_package(pgenlib) once plink-ng-s3 ships a CMake package.

add_library(regenie_pgenlib STATIC
  include/pgenlib_misc.cc
  include/pgenlib_read.cc
  include/pgenlib_write.cc
  include/plink2_base.cc
  include/plink2_bits.cc
  include/plink2_simd.cc
  pgenlib_ffi_support.cpp
  pgenlibr.cpp
  pvar_ffi_support.cc)
add_library(regenie::pgenlib ALIAS regenie_pgenlib)

target_include_directories(regenie_pgenlib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/simde)
target_compile_features(regenie_pgenlib PRIVATE cxx_std_17)

# S3/HTTP input. s3stream is a private copy of the same curl+SigV4 reader that
# bgen-limix carries: pgenlib stays independent of BGEN, so the sources are
# duplicated rather than shared.
if(REGENIE_HAVE_S3)
  target_sources(regenie_pgenlib PRIVATE
    plink2_s3.cc
    s3stream/s3stream.cc
    s3stream/s3stream_creds.cc
    s3stream/s3stream_http.cc)
  target_compile_definitions(regenie_pgenlib PRIVATE USE_S3 S3STREAM_ENABLE)
  target_link_libraries(regenie_pgenlib PRIVATE CURL::libcurl)
endif()

set_target_properties(regenie_pgenlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
