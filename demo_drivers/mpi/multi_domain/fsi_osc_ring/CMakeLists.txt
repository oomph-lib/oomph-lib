# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_osc_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_osc_ring C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_osc_ring
  SOURCES fsi_osc_ring.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fsi_osc_ring
  DEPENDS_ON fsi_osc_ring
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata fsi_osc_ring_partition.dat
  LABELS mpi fsi_osc_ring)

# Don't run this test alongside any other tests
set_tests_properties(mpi.fsi_osc_ring PROPERTIES PROCESSORS
                                                 ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving fsi_osc_ring subdirectory")
# ------------------------------------------------------------------------------
