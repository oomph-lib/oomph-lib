# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_entry_flow subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_entry_flow C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME full_tube
  SOURCES full_tube.cc
  LIBRARIES oomph::navier_stokes)

oomph_add_executable(
  NAME three_d_entry_flow
  SOURCES three_d_entry_flow.cc
  LIBRARIES oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.three_d_entry_flow
  TARGET_DEPENDENCIES full_tube three_d_entry_flow
  EXTRA_REQUIRES
  LABELS navier_stokes three_d_entry_flow)

message(VERBOSE "Leaving three_d_entry_flow subdirectory")
# ------------------------------------------------------------------------------
