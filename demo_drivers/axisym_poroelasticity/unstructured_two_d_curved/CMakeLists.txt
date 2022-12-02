# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_two_d_curved subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_two_d_curved C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_two_d_curved
  SOURCES unstructured_two_d_curved.cc
  LIBRARIES oomph::axisym_poroelasticity oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_curved
  SOURCES unstructured_two_d_curved.cc
  LIBRARIES oomph::axisym_poroelasticity oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_poroelasticity.unstructured_two_d_curved
  TARGET_DEPENDENCIES adaptive_unstructured_two_d_curved
                      unstructured_two_d_curved
  EXTRA_REQUIRES
  LABELS axisym_poroelasticity unstructured_two_d_curved)

message(VERBOSE "Leaving unstructured_two_d_curved subdirectory")
# ------------------------------------------------------------------------------
