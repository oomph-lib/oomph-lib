# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_prescribed_displ_lagr_mult C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_prescribed_displ_lagr_mult
  SOURCES three_d_prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.three_d_prescribed_displ_lagr_mult
  TARGET_DEPENDENCIES three_d_prescribed_displ_lagr_mult
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  EXTRA_REQUIRES distr_block_prescr.mcr
  LABELS mpi three_d_prescribed_displ_lagr_mult)

# Don't run this test alongside any other tests
set_tests_properties(mpi.three_d_prescribed_displ_lagr_mult
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving three_d_prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
