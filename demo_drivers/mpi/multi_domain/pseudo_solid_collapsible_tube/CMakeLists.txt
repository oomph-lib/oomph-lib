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
message(VERBOSE "Entered pseudo_solid_collapsible_tube subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pseudo_solid_collapsible_tube C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME pseudo_solid_collapsible_tube
    SOURCES pseudo_solid_collapsible_tube.cc
    LIBRARIES oomph::multi_physics oomph::navier_stokes oomph::solid
              oomph::constitutive oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_HYPRE AND OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME mpi.multi_domain.pseudo_solid_collapsible_tube
    DEPENDS_ON pseudo_solid_collapsible_tube
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
    TEST_FILES validate.sh validata pseudo_solid_collapsible_tube.mcr
               make_dummy_solid_solns.bash)

  # Don't run this test alongside any other tests
  set_tests_properties(mpi.multi_domain.pseudo_solid_collapsible_tube
                       PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
endif()

message(VERBOSE "Leaving pseudo_solid_collapsible_tube subdirectory")
# ------------------------------------------------------------------------------
