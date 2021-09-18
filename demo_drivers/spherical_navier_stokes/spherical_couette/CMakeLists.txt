# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_couette subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(spherical_couette C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME spherical_couette
  SOURCES spherical_couette.cc
  LIBRARIES oomph::spherical_navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.spherical_couette
  TARGET_DEPENDENCIES spherical_couette
  EXTRA_REQUIRES
  LABELS "spherical_navier_stokes;spherical_couette")

message(VERBOSE "Leaving spherical_couette subdirectory")
# ------------------------------------------------------------------------------
