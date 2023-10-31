# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered restart subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(restart C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat_2adapt_load_balance
  SOURCES two_d_unsteady_heat_2adapt_load_balance.cc
  LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.restart
  DEPENDS_ON two_d_unsteady_heat_2adapt_load_balance
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh
             validata
             load_balance.mcr
             show_dist.mcr
             prune_and_balance.mcr
             two_meshes.mcr
             partitioning.dat
  LABELS mpi restart)

# Don't run this test alongside any other tests
set_tests_properties(mpi.restart PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving restart subdirectory")
# ------------------------------------------------------------------------------
