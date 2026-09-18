# Davies' method for the distribution of a linear combination of chi-squares.
# qfc.cpp keeps upstream's file-scope state; qfc_mt.cpp is the thread-safe
# variant exposing qf_mt(). They define different symbols, so both are built and
# SKAT.cpp picks one via WITH_THREADSAFE_QF.

add_library(regenie_qf STATIC qfc.cpp qfc_mt.cpp)
add_library(regenie::qf ALIAS regenie_qf)

target_include_directories(regenie_qf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(regenie_qf PRIVATE cxx_std_17)
target_compile_definitions(regenie_qf PRIVATE NDEBUG)

set_target_properties(regenie_qf PROPERTIES POSITION_INDEPENDENT_CODE ON)
