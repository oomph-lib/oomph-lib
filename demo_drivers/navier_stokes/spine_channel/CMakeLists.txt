# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spine_channel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(spine_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME simple_spine_channel
  SOURCES simple_spine_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME spine_channel
  SOURCES spine_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME spine_channel2
  SOURCES spine_channel2.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.spine_channel
  TARGET_DEPENDENCIES simple_spine_channel spine_channel spine_channel2
  EXTRA_REQUIRES
  LABELS "navier_stokes;spine_channel")

message(VERBOSE "Leaving spine_channel subdirectory")
# ------------------------------------------------------------------------------
