# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vector_split subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vector_split C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vector_split
  SOURCES vector_split.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME mpi.vector_split
    TARGET_DEPENDENCIES vector_split
    EXTRA_REQUIRES
    LABELS mpi vector
    REQUIRES_MPI_VARIABLENP_RUN_COMMAND)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.vector_split PROPERTIES PROCESSORS 4)
else()
  message(
    STATUS
      "Skipping test 'mpi.vector_split' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving vector_split subdirectory")
# ------------------------------------------------------------------------------
