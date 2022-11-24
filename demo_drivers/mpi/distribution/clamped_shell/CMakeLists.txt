# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(clamped_shell C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_shell_with_arclength_cont
  SOURCES clamped_shell_with_arclength_cont.cc
  LIBRARIES oomph::shell)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.clamped_shell
  TARGET_DEPENDENCIES clamped_shell_with_arclength_cont
  EXTRA_REQUIRES
  LABELS mpi clamped_shell
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.clamped_shell PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving clamped_shell subdirectory")
# ------------------------------------------------------------------------------
