# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered restart subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(restart C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat_2adapt_load_balance
  SOURCES two_d_unsteady_heat_2adapt_load_balance.cc
  LIBRARIES oomph::unsteady_heat MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.restart
  TARGET_DEPENDENCIES two_d_unsteady_heat_2adapt_load_balance
  EXTRA_REQUIRES load_balance.mcr show_dist.mcr prune_and_balance.mcr
                 two_meshes.mcr partitioning.dat
  LABELS "mpi;restart"
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.restart PROPERTIES RESOURCE_LOCK
                                            REQUIRES_SOLE_USE_OF_CPU_MPI)

message(VERBOSE "Leaving restart subdirectory")
# ------------------------------------------------------------------------------
