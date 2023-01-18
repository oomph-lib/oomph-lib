# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(hp_adaptive_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_hp_adapt
  SOURCES two_d_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.hp_adaptive_poisson
  DEPENDS_ON two_d_poisson_hp_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata hp_adaptive_poisson_partition.dat
  LABELS mpi hp_adaptive_poisson adapt)

# Don't run this test alongside any other tests
set_tests_properties(mpi.hp_adaptive_poisson PROPERTIES PROCESSORS
                                                        ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving hp_adaptive_poisson subdirectory")
# ------------------------------------------------------------------------------
