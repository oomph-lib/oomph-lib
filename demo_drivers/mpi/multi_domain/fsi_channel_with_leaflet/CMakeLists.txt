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
  LIBRARIES oomph::beam oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fsi_channel_with_leaflet
  TARGET_DEPENDENCIES fsi_channel_with_leaflet
  EXTRA_REQUIRES fsi_channel_with_leaflet_partition.dat
  LABELS mpi fsi_channel_with_leaflet
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.fsi_channel_with_leaflet PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
