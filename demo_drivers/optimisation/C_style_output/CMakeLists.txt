# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered C_style_output subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(C_style_output C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME c_style_output
  SOURCES c_style_output.cc
  LIBRARIES oomph::unsteady_heat oomph::navier_stokes oomph::poisson oomph::beam
            oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.C_style_output
  TARGET_DEPENDENCIES c_style_output
  EXTRA_REQUIRES
  LABELS "optimisation;C_style_output"
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving C_style_output subdirectory")
# ------------------------------------------------------------------------------
