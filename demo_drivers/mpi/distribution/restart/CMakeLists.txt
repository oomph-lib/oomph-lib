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
message(VERBOSE "Entered restart subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
  TEST_NAME mpi.distribution.restart
  DEPENDS_ON two_d_unsteady_heat_2adapt_load_balance
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh
             validata
             load_balance.mcr
             show_dist.mcr
             prune_and_balance.mcr
             two_meshes.mcr
             partitioning.dat)

# Tell CTest how many processors the test is likely to use and don't let it run
# for longer than 10 minutes
set_tests_properties(mpi.distribution.restart
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC} TIMEOUT 600)

message(VERBOSE "Leaving restart subdirectory")
# ------------------------------------------------------------------------------
