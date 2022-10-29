# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(navier_stokes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_channel
  SOURCES rayleigh_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME rayleigh_circular_driven_cavity_adapt
  SOURCES ray_circ_cavity_adapt.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.navier_stokes
  TARGET_DEPENDENCIES rayleigh_channel rayleigh_circular_driven_cavity_adapt
  EXTRA_REQUIRES
  LABELS optimisation navier_stokes)

message(VERBOSE "Leaving navier_stokes subdirectory")
# ------------------------------------------------------------------------------
