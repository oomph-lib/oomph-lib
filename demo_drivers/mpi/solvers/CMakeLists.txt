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
message(VERBOSE "Entered solvers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(solvers C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME driven_cavity
  SOURCES driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME direct_solver_test
  SOURCES direct_solver_test.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME airy_cantilever
  SOURCES airy_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_channel_with_leaflet
  SOURCES fsi_channel_with_leaflet.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME TrilinosSolver_test
    SOURCES TrilinosSolver_test.cc
    LIBRARIES oomph::poisson oomph::meshes oomph::generic)

  oomph_add_executable(
    NAME two_d_multi_poisson
    SOURCES two_d_multi_poisson.cc
    LIBRARIES oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS driven_cavity direct_solver_test airy_cantilever
    fsi_channel_with_leaflet)

if(OOMPH_HAS_TRILINOS)
  list(APPEND TEST_TARGET_DEPS TrilinosSolver_test two_d_multi_poisson)
endif()

oomph_add_test(
  TEST_NAME mpi.solvers
  DEPENDS_ON ${TEST_TARGET_DEPS}
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata multi_poisson_partition.dat)

# Don't run this test alongside any other tests
set_tests_properties(mpi.solvers PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving solvers subdirectory")
# ------------------------------------------------------------------------------
