# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spin_up subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spin_up C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME spin_up
  SOURCES spin_up.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.spin_up
  TARGET_DEPENDENCIES spin_up
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes spin_up)

message(VERBOSE "Leaving spin_up subdirectory")
# ------------------------------------------------------------------------------
