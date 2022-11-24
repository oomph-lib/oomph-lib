# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_collapsible_tube subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pseudo_solid_collapsible_tube C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME pseudo_solid_collapsible_tube
    SOURCES pseudo_solid_collapsible_tube.cc
    LIBRARIES oomph::multi_physics oomph::navier_stokes oomph::solid)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME mpi.pseudo_solid_collapsible_tube
    TARGET_DEPENDENCIES pseudo_solid_collapsible_tube
    EXTRA_REQUIRES pseudo_solid_collapsible_tube.mcr make_dummy_solid_solns.bash
    LABELS mpi pseudo_solid_collapsible_tube
    REQUIRES_MPI_RUN_COMMAND)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.pseudo_solid_collapsible_tube PROPERTIES RUN_SERIAL
                                                                    TRUE)
endif()

message(VERBOSE "Leaving pseudo_solid_collapsible_tube subdirectory")
# ------------------------------------------------------------------------------
