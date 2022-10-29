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
  LIBRARIES oomph::generic MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME track_pitch
  SOURCES track_pitch.cc
  LIBRARIES oomph::axisym_navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.bifurcation_tracking
  TARGET_DEPENDENCIES pitchfork track_pitch
  EXTRA_REQUIRES
  LABELS mpi bifurcation_tracking
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.bifurcation_tracking
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

# TODO: Patch
set_tests_properties(mpi.bifurcation_tracking PROPERTIES DISABLED YES)

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
