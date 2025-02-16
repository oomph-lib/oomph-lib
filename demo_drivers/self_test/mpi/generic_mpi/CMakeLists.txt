# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered generic_mpi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(generic_mpi C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME generic_mpi_test
  SOURCES generic_mpi_test.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.mpi.generic_mpi
  DEPENDS_ON generic_mpi_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata
  LABELS mpi matrix)

# Don't run this test alongside any other tests
set_tests_properties(self_test.mpi.generic_mpi PROPERTIES PROCESSORS
                                                          ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving generic_mpi subdirectory")
# ------------------------------------------------------------------------------
