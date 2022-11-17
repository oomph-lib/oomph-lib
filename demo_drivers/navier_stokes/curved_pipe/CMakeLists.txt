# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered curved_pipe subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(curved_pipe C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME curved_pipe
  SOURCES curved_pipe.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME helical_pipe
  SOURCES helical_pipe.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.curved_pipe
  TARGET_DEPENDENCIES curved_pipe helical_pipe
  EXTRA_REQUIRES
  LABELS navier_stokes curved_pipe)

message(VERBOSE "Leaving curved_pipe subdirectory")
# ------------------------------------------------------------------------------
