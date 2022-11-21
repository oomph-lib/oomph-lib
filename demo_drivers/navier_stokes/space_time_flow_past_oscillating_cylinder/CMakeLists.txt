# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered space_time_flow_past_oscillating_cylinder subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(space_time_oscillating_cylinder C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME space_time_oscillating_cylinder
  SOURCES space_time_oscillating_cylinder.cc
  LIBRARIES oomph::space_time_navier_stokes
            oomph::space_time_block_preconditioner oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.space_time_oscillating_cylinder
  TARGET_DEPENDENCIES space_time_oscillating_cylinder
  EXTRA_REQUIRES
  LABELS navier_stokes space_time_oscillating_cylinder)

message(VERBOSE
        "Leaving space_time_flow_past_oscillating_cylinder subdirectory")
# ------------------------------------------------------------------------------
