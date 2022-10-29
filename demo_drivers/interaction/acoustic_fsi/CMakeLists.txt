# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered acoustic_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(acoustic_fsi C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME acoustic_fsi
  SOURCES acoustic_fsi.cc
  LIBRARIES oomph::helmholtz oomph::time_harmonic_linear_elasticity
            oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME unstructured_acoustic_fsi
  SOURCES unstructured_acoustic_fsi.cc
  LIBRARIES oomph::helmholtz oomph::time_harmonic_linear_elasticity
            oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.acoustic_fsi
  TARGET_DEPENDENCIES acoustic_fsi unstructured_acoustic_fsi
  EXTRA_REQUIRES
  LABELS interaction acoustic_fsi)

message(VERBOSE "Leaving acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
