# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered falling_jet subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(falling_jet C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME falling_jet
  SOURCES falling_jet.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::constitutive
            oomph::solid
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.falling_jet
  TARGET_DEPENDENCIES falling_jet
  EXTRA_REQUIRES
  LABELS navier_stokes falling_jet)

message(VERBOSE "Leaving falling_jet subdirectory")
# ------------------------------------------------------------------------------
