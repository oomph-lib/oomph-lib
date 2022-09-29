# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oscillating_sphere subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(oscillating_sphere C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_oscillating_sphere
  SOURCES oscillating_sphere.cc
  LIBRARIES oomph::pml_fourier_decomposed_helmholtz
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 14)

oomph_add_executable(
  NAME oscillating_sphere
  SOURCES oscillating_sphere.cc
  LIBRARIES oomph::pml_fourier_decomposed_helmholtz
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME pml_fourier_decomposed_helmholtz.oscillating_sphere
  TARGET_DEPENDENCIES adaptive_oscillating_sphere oscillating_sphere
  EXTRA_REQUIRES
  LABELS "pml_fourier_decomposed_helmholtz;oscillating_sphere")

message(VERBOSE "Leaving oscillating_sphere subdirectory")
# ------------------------------------------------------------------------------
