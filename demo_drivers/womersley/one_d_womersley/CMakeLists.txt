# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_womersley subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(one_d_womersley C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_womersley
  SOURCES one_d_womersley.cc
  LIBRARIES oomph::womersley oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME womersley.one_d_womersley
  TARGET_DEPENDENCIES one_d_womersley
  EXTRA_REQUIRES
  LABELS womersley one_d_womersley)

message(VERBOSE "Leaving one_d_womersley subdirectory")
# ------------------------------------------------------------------------------
