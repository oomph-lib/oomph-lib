# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_fsi subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(unstructured_adaptive_fsi C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_adaptive_2d_fsi
  SOURCES unstructured_adaptive_2d_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unstructured_adaptive_fsi
  TARGET_DEPENDENCIES unstructured_adaptive_2d_fsi
  EXTRA_REQUIRES
  LABELS "interaction;unstructured_adaptive_fsi;adapt")

message(VERBOSE "Leaving unstructured_adaptive_fsi subdirectory")
# ------------------------------------------------------------------------------
