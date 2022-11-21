# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_torus subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_torus C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME pressure_driven_torus
  SOURCES pressure_driven_torus.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::axisym_navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME unstructured_torus
  SOURCES unstructured_torus.cc
  LIBRARIES oomph::axisym_navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.unstructured_torus
  TARGET_DEPENDENCIES pressure_driven_torus unstructured_torus
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes unstructured_torus)

message(VERBOSE "Leaving unstructured_torus subdirectory")
# ------------------------------------------------------------------------------
