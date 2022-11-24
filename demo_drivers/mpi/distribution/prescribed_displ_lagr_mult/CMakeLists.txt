# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(prescribed_displ_lagr_mult C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult
  SOURCES prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive)

oomph_add_executable(
  NAME resize_hanging_node_tester
  SOURCES resize_hanging_node_tester.cc
  LIBRARIES oomph::solid oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.prescribed_displ_lagr_mult
  TARGET_DEPENDENCIES prescribed_displ_lagr_mult resize_hanging_node_tester
  EXTRA_REQUIRES presc_displ_lagr_mult_partition.dat load_balance.mcr
  LABELS mpi prescribed_displ_lagr_mult
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.prescribed_displ_lagr_mult PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
