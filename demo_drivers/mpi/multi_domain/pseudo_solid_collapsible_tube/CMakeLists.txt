# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_collapsible_tube subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pseudo_solid_collapsible_tube C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME pseudo_solid_collapsible_tube
    SOURCES pseudo_solid_collapsible_tube.cc
    LIBRARIES oomph::multi_physics oomph::navier_stokes oomph::solid
              MPI::MPI_CXX
    CXX_STANDARD 14)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS)

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  list(APPEND TEST_TARGET_DEPS pseudo_solid_collapsible_tube)
endif()

oomph_add_test(
  TEST_NAME mpi.pseudo_solid_collapsible_tube
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES pseudo_solid_collapsible_tube.mcr make_dummy_solid_solns.bash
  LABELS "mpi;pseudo_solid_collapsible_tube"
  REQUIRES_MPI_RUN_COMMAND)

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  # Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
  # We'll pass this to all other MPI-enabled self-tests to ensure that parallel
  # self-test do NOT run at the same time!
  set_tests_properties(mpi.pseudo_solid_collapsible_tube
                       PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)
endif()

message(VERBOSE "Leaving pseudo_solid_collapsible_tube subdirectory")
# ------------------------------------------------------------------------------
