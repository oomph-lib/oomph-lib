# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disk_oscillation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(disk_oscillation C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME disk_oscillation
  SOURCES disk_oscillation.cc
  LIBRARIES oomph::solid oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.disk_oscillation
  TARGET_DEPENDENCIES disk_oscillation
  EXTRA_REQUIRES
  LABELS solid disk_oscillation)

message(VERBOSE "Leaving disk_oscillation subdirectory")
# ------------------------------------------------------------------------------
