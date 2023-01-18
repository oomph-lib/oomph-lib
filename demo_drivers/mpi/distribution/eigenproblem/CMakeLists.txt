# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblem subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigenproblem C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME harmonic
    SOURCES harmonic.cc
    LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME mpi.eigenproblem
    DEPENDS_ON harmonic
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
    TEST_FILES validate.sh validata
    LABELS mpi eigenproblem)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.eigenproblem PROPERTIES PROCESSORS
                                                   ${OOMPH_MPI_NUM_PROC})
endif()

message(VERBOSE "Leaving eigenproblem subdirectory")
# ------------------------------------------------------------------------------
