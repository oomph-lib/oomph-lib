# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(one_d_poisson C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_poisson
  SOURCES one_d_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.one_d_poisson
  TARGET_DEPENDENCIES one_d_poisson
  EXTRA_REQUIRES
  LABELS "poisson;one_d_poisson")

message(VERBOSE "Leaving one_d_poisson subdirectory")
# ------------------------------------------------------------------------------
