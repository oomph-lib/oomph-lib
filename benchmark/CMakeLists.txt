# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered benchmark subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(benchmark C CXX Fortran)
if(NOT BUILDING_BENCHMARKS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(OomphGetGoogleBenchmark)

# -------------------------------[ BENCHMARKS ]---------------------------------

# The subdirectories to build
set(SUBDIRS generic utilities)

# Add each subdirectory target to the build
foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving benchmark subdirectory")
# ------------------------------------------------------------------------------
