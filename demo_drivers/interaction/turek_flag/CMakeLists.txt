# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(turek_flag C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag
  SOURCES turek_flag.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.turek_flag
  TARGET_DEPENDENCIES turek_flag
  EXTRA_REQUIRES
  LABELS "interaction;turek_flag")

message(VERBOSE "Leaving turek_flag subdirectory")
# ------------------------------------------------------------------------------
