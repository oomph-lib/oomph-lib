# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flux_control subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(flux_control C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME flux_control
  SOURCES flux_control.cc
  LIBRARIES oomph::womersley oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.flux_control
  TARGET_DEPENDENCIES flux_control
  EXTRA_REQUIRES
  LABELS "navier_stokes;flux_control")

message(VERBOSE "Leaving flux_control subdirectory")
# ------------------------------------------------------------------------------
