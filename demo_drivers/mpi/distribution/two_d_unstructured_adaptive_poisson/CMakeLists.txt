# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unstructured_adaptive_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_unstructured_adaptive_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_parallel_unstructured_adaptive_poisson
  SOURCES two_d_parallel_unstructured_adaptive_poisson.cc
  LIBRARIES oomph::poisson oomph::rigid_body oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_adaptive_mesh_two_outer_boundaries
  SOURCES unstructured_adaptive_mesh_two_outer_boundaries.cc
  LIBRARIES oomph::poisson oomph::rigid_body oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_adaptive_mesh_two_outer_boundaries_crossed
  SOURCES unstructured_adaptive_mesh_two_outer_boundaries_crossed.cc
  LIBRARIES oomph::poisson oomph::rigid_body oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.two_d_unstructured_adaptive_poisson
  DEPENDS_ON two_d_parallel_unstructured_adaptive_poisson
             unstructured_adaptive_mesh_two_outer_boundaries
             unstructured_adaptive_mesh_two_outer_boundaries_crossed
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
  TEST_FILES validate.sh validata
  LABELS mpi two_d_unstructured_adaptive_poisson adapt)

# Don't run this test alongside any other tests
set_tests_properties(mpi.two_d_unstructured_adaptive_poisson
                     PROPERTIES PROCESSORS 4)

message(VERBOSE "Leaving two_d_unstructured_adaptive_poisson subdirectory")
# ------------------------------------------------------------------------------
