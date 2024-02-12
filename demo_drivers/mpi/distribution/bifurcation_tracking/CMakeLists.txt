# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered bifurcation_tracking subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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
  TEST_FILES validate.sh validata
  LABELS mpi bifurcation_tracking)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.bifurcation_tracking
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.distribution.bifurcation_tracking PROPERTIES DISABLED
                                                                        FALSE)
endif()

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
