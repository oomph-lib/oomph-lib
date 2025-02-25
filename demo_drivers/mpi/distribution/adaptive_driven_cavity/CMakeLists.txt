# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(adaptive_driven_cavity C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_driven_cavity_load_balance
  SOURCES adaptive_driven_cavity_load_balance.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.adaptive_driven_cavity
  DEPENDS_ON adaptive_driven_cavity adaptive_driven_cavity_load_balance
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata adaptive_cavity_1_partition.dat
             adaptive_cavity_2_partition.dat load_balance.mcr
             prune_and_balance.mcr
  LABELS mpi adaptive_driven_cavity adapt)

# Tell CTest how many processors the test is likely to use
set_tests_properties(mpi.distribution.adaptive_driven_cavity
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
