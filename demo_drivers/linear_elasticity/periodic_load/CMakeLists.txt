# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered periodic_load subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(periodic_load C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME periodic_load
  SOURCES periodic_load.cc
  LIBRARIES oomph::linear_elasticity
  CXX_STANDARD 14)

oomph_add_executable(
  NAME refineable_periodic_load
  SOURCES refineable_periodic_load.cc
  LIBRARIES oomph::linear_elasticity
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linear_elasticity.periodic_load
  TARGET_DEPENDENCIES periodic_load refineable_periodic_load
  EXTRA_REQUIRES
  LABELS "linear_elasticity;periodic_load")

message(VERBOSE "Leaving periodic_load subdirectory")
# ------------------------------------------------------------------------------
