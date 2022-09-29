# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered circular_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(circular_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME circular_driven_cavity
  SOURCES circular_driven_cavity.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME circular_driven_cavity2
  SOURCES circular_driven_cavity2.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.circular_driven_cavity
  TARGET_DEPENDENCIES circular_driven_cavity circular_driven_cavity2
  EXTRA_REQUIRES
  LABELS "navier_stokes;circular_driven_cavity")

message(VERBOSE "Leaving circular_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
