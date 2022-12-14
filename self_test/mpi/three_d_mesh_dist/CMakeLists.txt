# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_mesh_dist subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_mesh_dist C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_mesh_dist
  SOURCES three_d_mesh_dist.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.three_d_mesh_dist
  TARGET_DEPENDENCIES three_d_mesh_dist
  EXTRA_REQUIRES three_d_mesh_dist_partition.dat
  LABELS mpi three_d_mesh
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.three_d_mesh_dist PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving three_d_mesh_dist subdirectory")
# ------------------------------------------------------------------------------
