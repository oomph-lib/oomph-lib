# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_variable_diff subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(two_d_variable_diff C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_variable_diff_adapt
  SOURCES two_d_variable_diff_adapt.cc
  LIBRARIES oomph::advection_diffusion
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_variable_diff
  TARGET_DEPENDENCIES two_d_variable_diff_adapt
  EXTRA_REQUIRES
  LABELS "advection_diffusion;two_d_variable_diff")

message(VERBOSE "Leaving two_d_variable_diff subdirectory")
# ------------------------------------------------------------------------------
