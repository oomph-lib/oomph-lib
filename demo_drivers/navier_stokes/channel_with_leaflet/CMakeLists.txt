# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(channel_with_leaflet C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME channel_with_leaflet
  SOURCES channel_with_leaflet.cc
  LIBRARIES oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.channel_with_leaflet
  TARGET_DEPENDENCIES channel_with_leaflet
  EXTRA_REQUIRES
  LABELS navier_stokes channel_with_leaflet)

message(VERBOSE "Leaving channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
