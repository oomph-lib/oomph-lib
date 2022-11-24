# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered FAQ subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(FAQ C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME broken_this_demo
  SOURCES broken_this_demo.cc
  SILENCE_NO_LIBS_SUPPLIED_WARNING)

oomph_add_executable(
  NAME customising_output
  SOURCES customising_output.cc
  LIBRARIES oomph::poisson)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME FAQ.FAQ
  TARGET_DEPENDENCIES broken_this_demo customising_output
  EXTRA_REQUIRES
  LABELS FAQ
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving FAQ subdirectory")
# ------------------------------------------------------------------------------
