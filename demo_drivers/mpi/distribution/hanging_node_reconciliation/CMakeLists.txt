# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hanging_node_reconciliation subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(hanging_node_reconciliation C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.hanging_node_reconciliation
  TARGET_DEPENDENCIES adaptive_driven_cavity hp_adaptive_driven_cavity
  EXTRA_REQUIRES nodes_cr.lay nodes.lay hp_adaptive_cavity_partition.dat
  LABELS mpi hanging_node_reconciliation
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.hanging_node_reconciliation
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

message(VERBOSE "Leaving hanging_node_reconciliation subdirectory")
# ------------------------------------------------------------------------------
