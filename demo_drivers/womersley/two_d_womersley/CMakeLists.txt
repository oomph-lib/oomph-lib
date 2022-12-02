# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_womersley subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_womersley C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_womersley
  SOURCES two_d_womersley.cc
  LIBRARIES oomph::womersley oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME womersley.two_d_womersley
  TARGET_DEPENDENCIES two_d_womersley
  EXTRA_REQUIRES
  LABELS womersley two_d_womersley)

message(VERBOSE "Leaving two_d_womersley subdirectory")
# ------------------------------------------------------------------------------
