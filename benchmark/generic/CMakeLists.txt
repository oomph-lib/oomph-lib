# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered generic subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(generic C CXX Fortran)
if(NOT BUILDING_BENCHMARKS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(OomphGetGoogleBenchmark)

# -------------------------------[ BENCHMARKS ]---------------------------------

add_executable(benchmark_vector benchmark_Vector.cc)
target_link_libraries(benchmark_vector benchmark::benchmark oomph::generic)

message(VERBOSE "Leaving generic subdirectory")
# ------------------------------------------------------------------------------
