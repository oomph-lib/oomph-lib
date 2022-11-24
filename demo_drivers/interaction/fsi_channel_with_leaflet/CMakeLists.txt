# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_channel_with_leaflet C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_channel_with_leaflet
  SOURCES fsi_channel_with_leaflet.cc
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_channel_with_leaflet
  TARGET_DEPENDENCIES fsi_channel_with_leaflet
  EXTRA_REQUIRES
  LABELS interaction fsi_channel_with_leaflet)

message(VERBOSE "Leaving fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
