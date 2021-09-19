# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered utilities subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(utilities C CXX Fortran)
if(NOT BUILDING_BENCHMARKS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(OomphGetGoogleBenchmark)

# -------------------------------[ BENCHMARKS ]---------------------------------

add_executable(benchmark_string benchmark_string.cc)
target_link_libraries(benchmark_string benchmark::benchmark)

message(VERBOSE "Leaving utilities subdirectory")
# ------------------------------------------------------------------------------
