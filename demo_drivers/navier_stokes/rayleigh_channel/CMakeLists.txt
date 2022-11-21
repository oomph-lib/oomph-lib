# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rayleigh_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_channel
  SOURCES rayleigh_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.rayleigh_channel
  TARGET_DEPENDENCIES rayleigh_channel
  EXTRA_REQUIRES
  LABELS navier_stokes rayleigh_channel)

message(VERBOSE "Leaving rayleigh_channel subdirectory")
# ------------------------------------------------------------------------------
