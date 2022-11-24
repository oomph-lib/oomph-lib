# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_solid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_solid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_solid
  SOURCES unstructured_three_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive)

oomph_add_executable(
  NAME unstructured_two_d_solid
  SOURCES unstructured_two_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.unstructured_solid
  TARGET_DEPENDENCIES unstructured_three_d_solid unstructured_two_d_solid
  EXTRA_REQUIRES
    cube_hole.1.ele
    cube_hole.1.face
    cube_hole.1.node
    solid.fig
    solid.fig.1.ele
    solid.fig.1.node
    solid.fig.1.poly
    solid.fig.poly
  LABELS solid unstructured_solid)

message(VERBOSE "Leaving unstructured_solid subdirectory")
# ------------------------------------------------------------------------------
