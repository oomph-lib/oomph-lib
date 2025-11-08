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
message(VERBOSE "Entered two_d_poisson_flux_bc_adapt subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_poisson_flux_bc_adapt C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_d_poisson_flux_bc_adapt
  SOURCES two_d_poisson_flux_bc_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.two_d_poisson_flux_bc_adapt
  DEPENDS_ON two_d_poisson_flux_bc_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata two_d_poisson_flux_partition.dat)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.two_d_poisson_flux_bc_adapt
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving two_d_poisson_flux_bc_adapt subdirectory")
# ------------------------------------------------------------------------------
