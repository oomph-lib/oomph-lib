# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_fluid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_fluid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_fluid
  SOURCES unstructured_two_d_fluid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_fluid
  TARGET_DEPENDENCIES unstructured_two_d_fluid
  EXTRA_REQUIRES fluid.fig fluid.fig.1.ele fluid.fig.1.node fluid.fig.1.poly
                 fluid.fig.poly
  LABELS navier_stokes unstructured_fluid)

message(VERBOSE "Leaving unstructured_fluid subdirectory")
# ------------------------------------------------------------------------------
