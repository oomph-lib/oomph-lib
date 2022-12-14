# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fish_poisson C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fish_poisson
  SOURCES fish_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fish_poisson
  TARGET_DEPENDENCIES fish_poisson
  EXTRA_REQUIRES
    fish_fully_automatic_partition.dat fish_incremental_partition.dat
    fish_selective_partition.dat distributed_adapt.mcr
  LABELS mpi fish_poisson
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.fish_poisson PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving fish_poisson subdirectory")
# ------------------------------------------------------------------------------
