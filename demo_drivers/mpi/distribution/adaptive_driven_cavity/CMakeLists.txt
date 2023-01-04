# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(adaptive_driven_cavity C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_driven_cavity_load_balance
  SOURCES adaptive_driven_cavity_load_balance.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.adaptive_driven_cavity
  TARGET_DEPENDENCIES adaptive_driven_cavity adaptive_driven_cavity_load_balance
  EXTRA_REQUIRES adaptive_cavity_1_partition.dat adaptive_cavity_2_partition.dat
                 load_balance.mcr prune_and_balance.mcr
  LABELS mpi adaptive_driven_cavity adapt
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.adaptive_driven_cavity
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.adaptive_driven_cavity PROPERTIES DISABLED TRUE)
endif()

message(VERBOSE "Leaving adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
