# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_traction_channel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(rayleigh_traction_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_traction_channel
  SOURCES rayleigh_traction_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.rayleigh_traction_channel
  TARGET_DEPENDENCIES rayleigh_traction_channel
  EXTRA_REQUIRES
  LABELS "navier_stokes;rayleigh_traction_channel")

message(VERBOSE "Leaving rayleigh_traction_channel subdirectory")
# ------------------------------------------------------------------------------
