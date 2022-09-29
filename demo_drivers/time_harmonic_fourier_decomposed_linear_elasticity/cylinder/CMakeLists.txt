# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered cylinder subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(cylinder C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 14)

oomph_add_executable(
  NAME adaptive_unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 14)

oomph_add_executable(
  NAME cylinder
  SOURCES cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
  CXX_STANDARD 14)

oomph_add_executable(
  NAME pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
  CXX_STANDARD 14)

oomph_add_executable(
  NAME unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME time_harmonic_fourier_decomposed_linear_elasticity.cylinder
  TARGET_DEPENDENCIES
    adaptive_pressure_loaded_cylinder adaptive_unstructured_cylinder cylinder
    pressure_loaded_cylinder unstructured_cylinder
  EXTRA_REQUIRES
  LABELS "time_harmonic_fourier_decomposed_linear_elasticity;cylinder")

message(VERBOSE "Leaving cylinder subdirectory")
# ------------------------------------------------------------------------------
