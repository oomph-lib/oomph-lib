# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered plate subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(plate C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME plate
  SOURCES plate.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.plate
  TARGET_DEPENDENCIES plate
  EXTRA_REQUIRES
  LABELS shell plate)

# This test contains a huge number of subtests and as such needs a larger amount
# of time than the default to complete
set_tests_properties(shell.plate PROPERTIES TIMEOUT 10000)

message(VERBOSE "Leaving plate subdirectory")
# ------------------------------------------------------------------------------
