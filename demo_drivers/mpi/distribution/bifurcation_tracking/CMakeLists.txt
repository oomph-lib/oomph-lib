# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered bifurcation_tracking subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(bifurcation_tracking C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME pitchfork
  SOURCES pitchfork.cc
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

oomph_add_executable(
  NAME track_pitch
  SOURCES track_pitch.cc
  LIBRARIES oomph::axisym_navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.bifurcation_tracking
  TARGET_DEPENDENCIES pitchfork track_pitch
  EXTRA_REQUIRES
  LABELS mpi bifurcation_tracking
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.bifurcation_tracking PROPERTIES RUN_SERIAL TRUE
                                                         DISABLED YES)

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
