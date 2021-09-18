# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_channel_seg_and_precond subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(fsi_channel_seg_and_precond C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_collapsible_channel_precond_driver
  SOURCES fsi_chan_precond_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::beam oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fsi_collapsible_channel_segregated_driver
  SOURCES fsi_chan_seg_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME simple_segregated_driver
  SOURCES simple_segregated_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_channel_seg_and_precond
  TARGET_DEPENDENCIES
    fsi_collapsible_channel_precond_driver
    fsi_collapsible_channel_segregated_driver simple_segregated_driver
  EXTRA_REQUIRES steady_precond.bash steady_seg.bash steady_seg_paper.bash
                 unsteady_precond.bash unsteady_seg.bash unsteady_seg_paper.bash
  LABELS "interaction;fsi_channel_seg_and_precond")

# This test needs more than the default allotted time to complete
set_tests_properties(interaction.fsi_channel_seg_and_precond PROPERTIES TIMEOUT
                                                                        10000)

message(VERBOSE "Leaving fsi_channel_seg_and_precond subdirectory")
# ------------------------------------------------------------------------------
