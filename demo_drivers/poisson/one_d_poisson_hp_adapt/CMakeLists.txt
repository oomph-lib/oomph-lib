# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson_hp_adapt subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(one_d_poisson_hp_adapt C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_poisson_hp_adapt
  SOURCES one_d_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.one_d_poisson_hp_adapt
  TARGET_DEPENDENCIES one_d_poisson_hp_adapt
  EXTRA_REQUIRES
  LABELS "poisson;one_d_poisson_hp_adapt;adapt")

message(VERBOSE "Leaving one_d_poisson_hp_adapt subdirectory")
# ------------------------------------------------------------------------------
