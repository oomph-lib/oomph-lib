# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered bifurcation_tracking subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(bifurcation_tracking C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_hopf
  SOURCES adaptive_hopf.cc
  LIBRARIES oomph::navier_stokes)

oomph_add_executable(
  NAME adaptive_hopf_with_separate_meshes
  SOURCES adaptive_hopf_with_separate_meshes.cc
  LIBRARIES oomph::linearised_navier_stokes oomph::navier_stokes)

oomph_add_executable(
  NAME adaptive_pitchfork
  SOURCES adaptive_pitchfork.cc
  LIBRARIES oomph::navier_stokes)

oomph_add_executable(
  NAME fold
  SOURCES fold.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME hopf
  SOURCES hopf.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME periodic_orbit
  SOURCES periodic_orbit.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME pitchfork
  SOURCES pitchfork.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME track_pitch
  SOURCES track_pitch.cc
  LIBRARIES oomph::axisym_navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME bifurcation_tracking.bifurcation_tracking
  TARGET_DEPENDENCIES
    adaptive_hopf
    adaptive_hopf_with_separate_meshes
    adaptive_pitchfork
    fold
    hopf
    periodic_orbit
    pitchfork
    track_pitch
  EXTRA_REQUIRES adapt_hopf_eigen.dat
  LABELS bifurcation_tracking)

# This test needs more than the default allotted time to complete
set_tests_properties(bifurcation_tracking.bifurcation_tracking
                     PROPERTIES TIMEOUT 5000)

message(VERBOSE "Leaving bifurcation_tracking subdirectory")
# ------------------------------------------------------------------------------
