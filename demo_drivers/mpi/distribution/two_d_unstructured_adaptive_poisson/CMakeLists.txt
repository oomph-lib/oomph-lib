# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unstructured_adaptive_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_unstructured_adaptive_poisson C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_parallel_unstructured_adaptive_poisson
  SOURCES two_d_parallel_unstructured_adaptive_poisson.cc
  LIBRARIES oomph::poisson oomph::rigid_body
  CXX_STANDARD 14)

oomph_add_executable(
  NAME unstructured_adaptive_mesh_two_outer_boundaries
  SOURCES unstructured_adaptive_mesh_two_outer_boundaries.cc
  LIBRARIES oomph::poisson oomph::rigid_body
  CXX_STANDARD 14)

oomph_add_executable(
  NAME unstructured_adaptive_mesh_two_outer_boundaries_crossed
  SOURCES unstructured_adaptive_mesh_two_outer_boundaries_crossed.cc
  LIBRARIES oomph::poisson oomph::rigid_body
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.two_d_unstructured_adaptive_poisson
  TARGET_DEPENDENCIES
    two_d_parallel_unstructured_adaptive_poisson
    unstructured_adaptive_mesh_two_outer_boundaries
    unstructured_adaptive_mesh_two_outer_boundaries_crossed
  EXTRA_REQUIRES
  LABELS mpi two_d_unstructured_adaptive_poisson adapt
  REQUIRES_MPI_VARIABLENP_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.two_d_unstructured_adaptive_poisson
                     PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving two_d_unstructured_adaptive_poisson subdirectory")
# ------------------------------------------------------------------------------
