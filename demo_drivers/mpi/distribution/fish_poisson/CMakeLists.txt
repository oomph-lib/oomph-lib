# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fish_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fish_poisson
  SOURCES fish_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.fish_poisson
  DEPENDS_ON fish_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata fish_fully_automatic_partition.dat
             fish_incremental_partition.dat fish_selective_partition.dat
             distributed_adapt.mcr)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.fish_poisson
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving fish_poisson subdirectory")
# ------------------------------------------------------------------------------
