# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered static_fish subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(static_fish C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME static_fish
  SOURCES static_fish.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.static_fish
  TARGET_DEPENDENCIES static_fish
  EXTRA_REQUIRES
  LABELS "solid;static_fish")

message(VERBOSE "Leaving static_fish subdirectory")
# ------------------------------------------------------------------------------
