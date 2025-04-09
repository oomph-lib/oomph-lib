# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(prescribed_displ_lagr_mult C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult
  SOURCES prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME resize_hanging_node_tester
  SOURCES resize_hanging_node_tester.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.prescribed_displ_lagr_mult
  DEPENDS_ON prescribed_displ_lagr_mult resize_hanging_node_tester
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata presc_displ_lagr_mult_partition.dat
             load_balance.mcr)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.prescribed_displ_lagr_mult
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
