# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered matrix_concatenation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(matrix_concatenation C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME matrix_concatenation
  SOURCES matrix_concatenation.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME mpi.matrix_concatenation
    TARGET_DEPENDENCIES matrix_concatenation
    EXTRA_REQUIRES
    LABELS mpi matrix
    REQUIRES_MPI_VARIABLENP_RUN_COMMAND)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.matrix_concatenation PROPERTIES RUN_SERIAL TRUE)
else()
  message(
    STATUS
      "Skipping test 'mpi.matrix_concatenation' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving matrix_concatenation subdirectory")
# ------------------------------------------------------------------------------
