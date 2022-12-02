# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_solid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_adaptive_solid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_adaptive_solid
  SOURCES unstructured_adaptive_solid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.unstructured_adaptive_solid
  TARGET_DEPENDENCIES unstructured_adaptive_solid
  EXTRA_REQUIRES
  LABELS solid unstructured_adaptive_solid adapt)

message(VERBOSE "Leaving unstructured_adaptive_solid subdirectory")
# ------------------------------------------------------------------------------
