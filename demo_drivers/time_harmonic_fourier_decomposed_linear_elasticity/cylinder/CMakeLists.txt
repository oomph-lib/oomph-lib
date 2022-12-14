# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered cylinder subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(cylinder C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME adaptive_unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME cylinder
  SOURCES cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME time_harmonic_fourier_decomposed_linear_elasticity.cylinder
  TARGET_DEPENDENCIES
    adaptive_pressure_loaded_cylinder adaptive_unstructured_cylinder cylinder
    pressure_loaded_cylinder unstructured_cylinder
  EXTRA_REQUIRES
  LABELS time_harmonic_fourier_decomposed_linear_elasticity cylinder)

message(VERBOSE "Leaving cylinder subdirectory")
# ------------------------------------------------------------------------------
