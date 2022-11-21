# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(elastic_poisson C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_poisson
  SOURCES elastic_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.elastic_poisson
  TARGET_DEPENDENCIES elastic_poisson
  EXTRA_REQUIRES
  LABELS poisson elastic_poisson)

message(VERBOSE "Leaving elastic_poisson subdirectory")
# ------------------------------------------------------------------------------
