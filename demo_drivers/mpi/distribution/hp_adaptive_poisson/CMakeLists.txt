# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(hp_adaptive_poisson C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_hp_adapt
  SOURCES two_d_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.hp_adaptive_poisson
  TARGET_DEPENDENCIES two_d_poisson_hp_adapt
  EXTRA_REQUIRES hp_adaptive_poisson_partition.dat
  LABELS mpi hp_adaptive_poisson adapt
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.hp_adaptive_poisson PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving hp_adaptive_poisson subdirectory")
# ------------------------------------------------------------------------------
