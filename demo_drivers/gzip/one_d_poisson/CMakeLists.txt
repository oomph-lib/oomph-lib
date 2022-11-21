# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(one_d_poisson C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_poisson_with_gzip_validation
  SOURCES one_d_poisson_with_gzip_validation.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME gzip.one_d_poisson
  TARGET_DEPENDENCIES one_d_poisson_with_gzip_validation
  EXTRA_REQUIRES
  LABELS gzip one_d_poisson
  NO_VALIDATE_SH)

message(VERBOSE "Leaving one_d_poisson subdirectory")
# ------------------------------------------------------------------------------
