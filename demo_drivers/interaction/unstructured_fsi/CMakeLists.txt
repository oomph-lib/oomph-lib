# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_fsi C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_fsi
  SOURCES unstructured_two_d_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unstructured_fsi
  TARGET_DEPENDENCIES unstructured_two_d_fsi
  EXTRA_REQUIRES
    fluid.fig
    fluid.fig.1.ele
    fluid.fig.1.node
    fluid.fig.1.poly
    fluid.fig.poly
    solid.fig
    solid.fig.1.ele
    solid.fig.1.node
    solid.fig.1.poly
    solid.fig.poly
  LABELS "interaction;unstructured_fsi")

message(VERBOSE "Leaving unstructured_fsi subdirectory")
# ------------------------------------------------------------------------------
