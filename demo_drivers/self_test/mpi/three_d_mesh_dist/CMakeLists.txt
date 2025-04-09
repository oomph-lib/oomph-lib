# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_mesh_dist subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_mesh_dist C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_mesh_dist
  SOURCES three_d_mesh_dist.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.mpi.three_d_mesh_dist
  DEPENDS_ON three_d_mesh_dist
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata three_d_mesh_dist_partition.dat)

# Don't run this test alongside any other tests
set_tests_properties(self_test.mpi.three_d_mesh_dist
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving three_d_mesh_dist subdirectory")
# ------------------------------------------------------------------------------
