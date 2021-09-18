# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_three_d_solid subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(unstructured_three_d_solid C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_solid
  SOURCES unstructured_three_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.unstructured_three_d_solid
  TARGET_DEPENDENCIES unstructured_three_d_solid
  EXTRA_REQUIRES
    fsi_bifurcation_solid.1.ele fsi_bifurcation_solid.1.face
    fsi_bifurcation_solid.1.node fsi_bifurcation_solid.poly
    three_d_unstructured_solid.mcr
  LABELS "solid;unstructured_three_d_solid")

message(VERBOSE "Leaving unstructured_three_d_solid subdirectory")
# ------------------------------------------------------------------------------
