# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(
  VERBOSE
  "Entered time_harmonic_fourier_decomposed_linear_elasticity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(time_harmonic_fourier_decomposed_linear_elasticity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS cylinder)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(
  VERBOSE
  "Leaving time_harmonic_fourier_decomposed_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
