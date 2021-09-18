# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(two_d_unsteady_heat C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat
  SOURCES two_d_unsteady_heat.cc
  LIBRARIES oomph::unsteady_heat
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_unsteady_heat_restarted
  SOURCES two_d_unsteady_heat_restarted.cc
  LIBRARIES oomph::unsteady_heat
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME unsteady_heat.two_d_unsteady_heat
  TARGET_DEPENDENCIES two_d_unsteady_heat two_d_unsteady_heat_restarted
  EXTRA_REQUIRES
  LABELS "unsteady_heat;two_d_unsteady_heat")

message(VERBOSE "Leaving two_d_unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
