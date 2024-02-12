# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_collapsible_channel C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_collapsible_channel_adapt
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_adapt
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.multi_domain.fsi_collapsible_channel
  DEPENDS_ON fsi_collapsible_channel_adapt fsi_collapsible_channel_macro_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata fsi_collapsible_channel_partition.dat
  LABELS mpi fsi_collapsible_channel)

# Don't run this test alongside any other tests
set_tests_properties(mpi.multi_domain.fsi_collapsible_channel
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving fsi_collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
