# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered tensioned_string subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(tensioned_string C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME tensioned_string
  SOURCES tensioned_string.cc
  LIBRARIES oomph::beam
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.tensioned_string
  TARGET_DEPENDENCIES tensioned_string
  EXTRA_REQUIRES
  LABELS "beam;tensioned_string")

message(VERBOSE "Leaving tensioned_string subdirectory")
# ------------------------------------------------------------------------------
