# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered cylinder subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(cylinder C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME cylinder
  SOURCES cylinder.cc
  LIBRARIES oomph::axisym_linear_elasticity
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_linear_elasticity.cylinder
  TARGET_DEPENDENCIES cylinder
  EXTRA_REQUIRES
  LABELS "axisym_linear_elasticity;cylinder")

message(VERBOSE "Leaving cylinder subdirectory")
# ------------------------------------------------------------------------------
