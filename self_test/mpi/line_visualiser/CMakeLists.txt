# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered line_visualiser subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(line_visualiser C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.line_visualiser
  TARGET_DEPENDENCIES adaptive_driven_cavity
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  EXTRA_REQUIRES adaptive_cavity_1_partition.dat adaptive_cavity_2_partition.dat
  LABELS mpi adapt driven_cavity)

# Don't run this test alongside any other tests
set_tests_properties(mpi.line_visualiser PROPERTIES PROCESSORS
                                                    ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving line_visualiser subdirectory")
# ------------------------------------------------------------------------------
