# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_prescribed_displ_lagr_mult C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_prescribed_displ_lagr_mult
  SOURCES three_d_prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.three_d_prescribed_displ_lagr_mult
  TARGET_DEPENDENCIES three_d_prescribed_displ_lagr_mult
  EXTRA_REQUIRES distr_block_prescr.mcr
  LABELS mpi three_d_prescribed_displ_lagr_mult
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.three_d_prescribed_displ_lagr_mult
                     PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving three_d_prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
