# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered sparse_assemble subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(sparse_assemble C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME sparse_assemble_test
  SOURCES sparse_assemble_test.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.sparse_assemble
  TARGET_DEPENDENCIES sparse_assemble_test
  EXTRA_REQUIRES
  LABELS optimisation sparse_assemble)

message(VERBOSE "Leaving sparse_assemble subdirectory")
# ------------------------------------------------------------------------------
