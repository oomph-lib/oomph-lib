# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered point_source subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(point_source C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_helmholtz_point_source
  SOURCES helmholtz_point_source.cc
  LIBRARIES oomph::helmholtz
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 14)

oomph_add_executable(
  NAME helmholtz_point_source
  SOURCES helmholtz_point_source.cc
  LIBRARIES oomph::helmholtz
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME helmholtz.point_source
  TARGET_DEPENDENCIES adaptive_helmholtz_point_source helmholtz_point_source
  EXTRA_REQUIRES
  LABELS helmholtz point_source)

message(VERBOSE "Leaving point_source subdirectory")
# ------------------------------------------------------------------------------
