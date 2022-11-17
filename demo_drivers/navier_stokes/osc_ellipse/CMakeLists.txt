# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered osc_ellipse subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(osc_ellipse C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME osc_quarter_ellipse
  SOURCES osc_quarter_ellipse.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.osc_ellipse
  TARGET_DEPENDENCIES osc_quarter_ellipse
  EXTRA_REQUIRES
  LABELS navier_stokes osc_ellipse)

message(VERBOSE "Leaving osc_ellipse subdirectory")
# ------------------------------------------------------------------------------
