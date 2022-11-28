# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(one_d_poisson_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_poisson_adapt
  SOURCES one_d_poisson_adapt.cc
  LIBRARIES oomph::poisson)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.one_d_poisson_adapt
  TARGET_DEPENDENCIES one_d_poisson_adapt
  EXTRA_REQUIRES
  LABELS poisson one_d_poisson adapt)

message(VERBOSE "Leaving one_d_poisson_adapt subdirectory")
# ------------------------------------------------------------------------------
