# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_collapsible_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_collapsible_channel_adapt
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_adapt
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fsi_collapsible_channel
  TARGET_DEPENDENCIES fsi_collapsible_channel_adapt
                      fsi_collapsible_channel_macro_adapt
  EXTRA_REQUIRES fsi_collapsible_channel_partition.dat
  LABELS mpi fsi_collapsible_channel
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.fsi_collapsible_channel PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving fsi_collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
