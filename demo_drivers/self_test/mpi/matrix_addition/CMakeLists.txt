# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered matrix_addition subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(matrix_addition C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME matrix_addition
  SOURCES matrix_addition.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME self_test.mpi.matrix_addition
    DEPENDS_ON matrix_addition
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
    TEST_FILES validate.sh validata
    LABELS mpi matrix)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.matrix_addition PROPERTIES PROCESSORS 4)
else()
  message(
    STATUS
      "Skipping test 'mpi.matrix_addition' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving matrix_addition subdirectory")
# ------------------------------------------------------------------------------
