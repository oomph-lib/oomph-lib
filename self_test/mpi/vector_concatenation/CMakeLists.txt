# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vector_concatenation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vector_concatenation C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vector_concatenation
  SOURCES vector_concatenation.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.vector_concatenation
  TARGET_DEPENDENCIES vector_concatenation
  EXTRA_REQUIRES
  LABELS mpi vector
  REQUIRES_MPI_VARIABLENP_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.vector_concatenation PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving vector_concatenation subdirectory")
# ------------------------------------------------------------------------------
