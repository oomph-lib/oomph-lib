# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(adaptive_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.adaptive_driven_cavity
  TARGET_DEPENDENCIES adaptive_driven_cavity
  EXTRA_REQUIRES
  LABELS "navier_stokes;adaptive_driven_cavity;adapt")

message(VERBOSE "Leaving adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
