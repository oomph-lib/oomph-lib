# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_three_d_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_three_d_fsi C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_fsi
  SOURCES unstructured_three_d_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unstructured_three_d_fsi
  TARGET_DEPENDENCIES unstructured_three_d_fsi
  EXTRA_REQUIRES
    fsi_bifurcation_fluid.1.ele
    fsi_bifurcation_fluid.1.face
    fsi_bifurcation_fluid.1.node
    fsi_bifurcation_fluid.poly
    fsi_bifurcation_solid.1.ele
    fsi_bifurcation_solid.1.face
    fsi_bifurcation_solid.1.node
    fsi_bifurcation_solid.poly
  LABELS interaction unstructured_three_d_fsi)

message(VERBOSE "Leaving unstructured_three_d_fsi subdirectory")
# ------------------------------------------------------------------------------
