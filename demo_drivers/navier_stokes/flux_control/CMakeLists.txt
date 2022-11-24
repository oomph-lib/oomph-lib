# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flux_control subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(flux_control C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME flux_control
  SOURCES flux_control.cc
  LIBRARIES oomph::womersley oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.flux_control
  TARGET_DEPENDENCIES flux_control
  EXTRA_REQUIRES
  LABELS navier_stokes flux_control)

message(VERBOSE "Leaving flux_control subdirectory")
# ------------------------------------------------------------------------------
