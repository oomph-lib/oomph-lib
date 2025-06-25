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
message(VERBOSE "Entered fsi_channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_channel_with_leaflet C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME fsi_channel_with_leaflet
  SOURCES fsi_channel_with_leaflet.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.multi_domain.fsi_channel_with_leaflet
  DEPENDS_ON fsi_channel_with_leaflet
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata fsi_channel_with_leaflet_partition.dat)

# Don't run this test alongside any other tests
set_tests_properties(mpi.multi_domain.fsi_channel_with_leaflet
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
