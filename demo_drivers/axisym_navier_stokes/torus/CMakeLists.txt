# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered torus subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(torus C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME torus
  SOURCES torus.cc
  LIBRARIES oomph::navier_stokes oomph::axisym_navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.torus
  TARGET_DEPENDENCIES torus
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes torus)

message(VERBOSE "Leaving torus subdirectory")
# ------------------------------------------------------------------------------
