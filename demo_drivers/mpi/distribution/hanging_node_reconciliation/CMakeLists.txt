# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hanging_node_reconciliation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(hanging_node_reconciliation C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.hanging_node_reconciliation
  DEPENDS_ON adaptive_driven_cavity hp_adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata nodes_cr.lay nodes.lay
             hp_adaptive_cavity_partition.dat
  LABELS mpi hanging_node_reconciliation)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.hanging_node_reconciliation
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving hanging_node_reconciliation subdirectory")
# ------------------------------------------------------------------------------
