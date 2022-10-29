# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_three_d_fluid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_three_d_fluid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_fluid
  SOURCES unstructured_three_d_fluid.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_three_d_fluid
  TARGET_DEPENDENCIES unstructured_three_d_fluid
  EXTRA_REQUIRES
    fsi_bifurcation_fluid.1.ele
    fsi_bifurcation_fluid.1.face
    fsi_bifurcation_fluid.1.node
    fsi_bifurcation_fluid.poly
    fsi_bifurcation_fluid_with_extended_tubes.1.ele
    fsi_bifurcation_fluid_with_extended_tubes.1.face
    fsi_bifurcation_fluid_with_extended_tubes.1.node
    fsi_bifurcation_fluid_with_extended_tubes.poly
    three_d_unstructured_fluid.mcr
  LABELS navier_stokes unstructured_three_d_fluid)

message(VERBOSE "Leaving unstructured_three_d_fluid subdirectory")
# ------------------------------------------------------------------------------
