# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_instability subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(rayleigh_instability C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_instability
  SOURCES rayleigh_instability.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.rayleigh_instability
  TARGET_DEPENDENCIES rayleigh_instability
  EXTRA_REQUIRES
  LABELS "axisym_navier_stokes;rayleigh_instability")

message(VERBOSE "Leaving rayleigh_instability subdirectory")
# ------------------------------------------------------------------------------
