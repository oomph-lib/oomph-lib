# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(hp_adaptive_driven_cavity C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.hp_adaptive_driven_cavity
  DEPENDS_ON hp_adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata hp_adaptive_cavity_partition.dat
  LABELS mpi hp_adaptive_driven_cavity adapt)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.hp_adaptive_driven_cavity
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.distribution.hp_adaptive_driven_cavity
                       PROPERTIES DISABLED FALSE)
endif()

message(VERBOSE "Leaving hp_adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
