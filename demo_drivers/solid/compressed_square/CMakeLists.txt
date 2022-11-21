# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered compressed_square subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(compressed_square C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME compressed_square
  SOURCES compressed_square.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.compressed_square
  TARGET_DEPENDENCIES compressed_square
  EXTRA_REQUIRES
  LABELS solid compressed_square)

message(VERBOSE "Leaving compressed_square subdirectory")
# ------------------------------------------------------------------------------
