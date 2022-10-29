# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_biharmonic subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_biharmonic C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_biharmonic
  SOURCES two_d_biharmonic.cc
  LIBRARIES oomph::biharmonic
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME biharmonic.two_d_biharmonic
  TARGET_DEPENDENCIES two_d_biharmonic
  EXTRA_REQUIRES
  LABELS biharmonic two_d_biharmonic)

message(VERBOSE "Leaving two_d_biharmonic subdirectory")
# ------------------------------------------------------------------------------
