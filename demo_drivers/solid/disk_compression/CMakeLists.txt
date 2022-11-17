# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disk_compression subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(disk_compression C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME disk_compression
  SOURCES disk_compression.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.disk_compression
  TARGET_DEPENDENCIES disk_compression
  EXTRA_REQUIRES
  LABELS solid disk_compression)

message(VERBOSE "Leaving disk_compression subdirectory")
# ------------------------------------------------------------------------------
