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
message(VERBOSE "Entered bifurcation_tracking subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(bifurcation_tracking C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME pitchfork
  SOURCES pitchfork.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME track_pitch
  SOURCES track_pitch.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.distribution.bifurcation_tracking
  DEPENDS_ON pitchfork track_pitch
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.bifurcation_tracking
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.distribution.bifurcation_tracking PROPERTIES DISABLED
                                                                        FALSE)
endif()

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
