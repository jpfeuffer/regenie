cmake_minimum_required(VERSION 3.13)
project(regenie_s3_test)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(AWSSDK REQUIRED COMPONENTS s3)
add_definitions(-DWITH_AWS_S3)

# pgenlib path (bundled in the repo)
set(PGENLIB_PATH "${CMAKE_SOURCE_DIR}/../external_libs/pgenlib")

# Build pgenlib.a with S3 streaming support
add_custom_target(
  pgenlib_s3
  COMMAND ${CMAKE_MAKE_PROGRAM}
          WITH_S3=1
          AWS_INCLUDES=-I${AWSSDK_INCLUDE_DIRS}
          CXX=${CMAKE_CXX_COMPILER}
  WORKING_DIRECTORY ${PGENLIB_PATH}
)

add_executable(test_s3
  ${CMAKE_SOURCE_DIR}/../src/S3_Utils.cpp
  ${CMAKE_SOURCE_DIR}/test_s3.cpp
)
target_include_directories(test_s3 PRIVATE
  ${CMAKE_SOURCE_DIR}/../src
  ${PGENLIB_PATH}
  ${PGENLIB_PATH}/include
)
target_link_libraries(test_s3 PRIVATE
  ${AWSSDK_LINK_LIBRARIES}
  ${PGENLIB_PATH}/pgenlib.a
)
add_dependencies(test_s3 pgenlib_s3)
