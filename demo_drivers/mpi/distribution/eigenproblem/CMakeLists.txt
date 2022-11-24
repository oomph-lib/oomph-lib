# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblem subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigenproblem C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME harmonic
    SOURCES harmonic.cc
    LIBRARIES oomph::navier_stokes)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME mpi.eigenproblem
    TARGET_DEPENDENCIES harmonic
    EXTRA_REQUIRES
    LABELS mpi eigenproblem
    REQUIRES_MPI_RUN_COMMAND)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.eigenproblem PROPERTIES RUN_SERIAL TRUE)
endif()

message(VERBOSE "Leaving eigenproblem subdirectory")
# ------------------------------------------------------------------------------
