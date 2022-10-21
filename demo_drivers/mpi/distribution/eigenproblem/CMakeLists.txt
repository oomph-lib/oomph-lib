# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblem subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(eigenproblem C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME harmonic
    SOURCES harmonic.cc
    LIBRARIES oomph::navier_stokes MPI::MPI_CXX
    CXX_STANDARD 14)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME mpi.eigenproblem
    TARGET_DEPENDENCIES harmonic
    EXTRA_REQUIRES
    LABELS "mpi;eigenproblem"
    REQUIRES_MPI_RUN_COMMAND)

  # Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
  # We'll pass this to all other MPI-enabled self-tests to ensure that parallel
  # self-test do NOT run at the same time!
  set_tests_properties(mpi.eigenproblem PROPERTIES RESOURCE_LOCK
                                                   REQUIRES_SOLE_USE_OF_CPU_MPI)
endif()

message(VERBOSE "Leaving eigenproblem subdirectory")
# ------------------------------------------------------------------------------
