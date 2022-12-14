# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered restart subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(restart C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat_2adapt_load_balance
  SOURCES two_d_unsteady_heat_2adapt_load_balance.cc
  LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.restart
  TARGET_DEPENDENCIES two_d_unsteady_heat_2adapt_load_balance
  EXTRA_REQUIRES load_balance.mcr show_dist.mcr prune_and_balance.mcr
                 two_meshes.mcr partitioning.dat
  LABELS mpi restart
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.restart PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving restart subdirectory")
# ------------------------------------------------------------------------------
