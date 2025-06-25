# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_channel_seg_and_precond subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_channel_seg_and_precond C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME fsi_collapsible_channel_precond_driver
  SOURCES fsi_chan_precond_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_collapsible_channel_segregated_driver
  SOURCES fsi_chan_seg_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME simple_segregated_driver
  SOURCES simple_segregated_driver.cc fsi_chan_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_channel_seg_and_precond
  DEPENDS_ON fsi_collapsible_channel_precond_driver
             fsi_collapsible_channel_segregated_driver simple_segregated_driver
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             steady_precond.bash
             steady_seg.bash
             steady_seg_paper.bash
             unsteady_precond.bash
             unsteady_seg.bash
             unsteady_seg_paper.bash)

# This test needs more than the default allotted time to complete
set_tests_properties(interaction.fsi_channel_seg_and_precond PROPERTIES TIMEOUT
                                                                        3600)

message(VERBOSE "Leaving fsi_channel_seg_and_precond subdirectory")
# ------------------------------------------------------------------------------
