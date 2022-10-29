# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flow_past_oscillating_cylinder subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(flow_past_oscillating_cylinder C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME flow_past_oscillating_cylinder
  SOURCES flow_past_oscillating_cylinder.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.flow_past_oscillating_cylinder
  TARGET_DEPENDENCIES flow_past_oscillating_cylinder
  EXTRA_REQUIRES
  LABELS navier_stokes flow_past_oscillating_cylinder)

message(VERBOSE "Leaving flow_past_oscillating_cylinder subdirectory")
# ------------------------------------------------------------------------------
