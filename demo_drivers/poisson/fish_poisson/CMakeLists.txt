# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(fish_poisson C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fish_poisson
  SOURCES fish_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fish_poisson_no_adapt
  SOURCES fish_poisson_no_adapt.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.fish_poisson
  TARGET_DEPENDENCIES fish_poisson fish_poisson_no_adapt
  EXTRA_REQUIRES
  LABELS "poisson;fish_poisson")

message(VERBOSE "Leaving fish_poisson subdirectory")
# ------------------------------------------------------------------------------
