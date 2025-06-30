# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unstructured_adaptive_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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

set(EXTRA_FLAG --non_ref_bin)
if(OOMPH_HAS_CGAL)
  set(EXTRA_FLAG --cgal)
endif()

oomph_add_test(
  TEST_NAME mpi.distribution.two_d_unstructured_adaptive_poisson
  DEPENDS_ON two_d_parallel_unstructured_adaptive_poisson
             unstructured_adaptive_mesh_two_outer_boundaries
             unstructured_adaptive_mesh_two_outer_boundaries_crossed
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
          ${EXTRA_FLAG}
  TEST_FILES validate.sh validata)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.two_d_unstructured_adaptive_poisson
                     PROPERTIES PROCESSORS 4)

message(VERBOSE "Leaving two_d_unstructured_adaptive_poisson subdirectory")
# ------------------------------------------------------------------------------
