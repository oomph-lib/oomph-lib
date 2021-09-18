# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(two_d_poisson C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_poisson2
  SOURCES two_d_poisson2.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# Do we have the hsl sources?
set(CXX_DEFINITION)
if(HAVE_HSL_SOURCES)
  set(CXX_DEFINITION HAVE_HSL_SOURCES)
endif()

oomph_add_executable(
  NAME two_d_poisson_compare_solvers
  SOURCES two_d_poisson_compare_solvers.cc
  LIBRARIES oomph::poisson
  CXX_DEFINITIONS ${CXX_DEFINITION}
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson
  TARGET_DEPENDENCIES two_d_poisson two_d_poisson2 two_d_poisson_compare_solvers
  EXTRA_REQUIRES
  LABELS "poisson;two_d_poisson")

message(VERBOSE "Leaving two_d_poisson subdirectory")
# ------------------------------------------------------------------------------
