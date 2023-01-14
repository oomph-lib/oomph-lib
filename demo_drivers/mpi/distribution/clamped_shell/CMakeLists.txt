# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(clamped_shell C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_shell_with_arclength_cont
  SOURCES clamped_shell_with_arclength_cont.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.clamped_shell
  TARGET_DEPENDENCIES clamped_shell_with_arclength_cont
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  EXTRA_REQUIRES
  LABELS mpi clamped_shell)

# Don't run this test alongside any other tests
set_tests_properties(mpi.clamped_shell PROPERTIES PROCESSORS
                                                  ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving clamped_shell subdirectory")
# ------------------------------------------------------------------------------
