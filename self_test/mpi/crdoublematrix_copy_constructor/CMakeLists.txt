# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered crdoublematrix_copy_constructor subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(crdoublematrix_copy_constructor C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME crdoublematrix_copy_constructor
  SOURCES crdoublematrix_copy_constructor.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME mpi.crdoublematrix_copy_constructor
    TARGET_DEPENDENCIES crdoublematrix_copy_constructor
    VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
    EXTRA_REQUIRES
    LABELS mpi matrix)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.crdoublematrix_copy_constructor PROPERTIES PROCESSORS
                                                                      4)
else()
  message(
    STATUS
      "Skipping test 'mpi.crdoublematrix_copy_constructor' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving crdoublematrix_copy_constructor subdirectory")
# ------------------------------------------------------------------------------
