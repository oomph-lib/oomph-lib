# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_two_d_circle subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_two_d_circle C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_two_d_circle
  SOURCES unstructured_two_d_circle.cc
  LIBRARIES oomph::darcy
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_circle
  SOURCES unstructured_two_d_circle.cc
  LIBRARIES oomph::darcy)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME darcy.unstructured_two_d_circle
  TARGET_DEPENDENCIES adaptive_unstructured_two_d_circle
                      unstructured_two_d_circle
  EXTRA_REQUIRES
  LABELS darcy unstructured_two_d_circle)

message(VERBOSE "Leaving unstructured_two_d_circle subdirectory")
# ------------------------------------------------------------------------------
