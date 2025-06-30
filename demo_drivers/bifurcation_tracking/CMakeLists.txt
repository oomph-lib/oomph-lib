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
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_hopf
  SOURCES adaptive_hopf.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_hopf_with_separate_meshes
  SOURCES adaptive_hopf_with_separate_meshes.cc
  LIBRARIES oomph::linearised_navier_stokes oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME adaptive_pitchfork
  SOURCES adaptive_pitchfork.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fold
  SOURCES fold.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME hopf
  SOURCES hopf.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME periodic_orbit
  SOURCES periodic_orbit.cc
  LIBRARIES oomph::meshes oomph::generic)

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
  TEST_NAME bifurcation_tracking.bifurcation_tracking
  DEPENDS_ON adaptive_hopf
             adaptive_hopf_with_separate_meshes
             adaptive_pitchfork
             fold
             hopf
             periodic_orbit
             pitchfork
             track_pitch
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata adapt_hopf_eigen.dat)

# This test needs more than the default allotted time to complete
set_tests_properties(bifurcation_tracking.bifurcation_tracking
                     PROPERTIES TIMEOUT 5000)

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
