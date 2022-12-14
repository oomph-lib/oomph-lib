# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rotating_ends subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rotating_ends C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME rotating_ends
  SOURCES rotating_ends.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.rotating_ends
  TARGET_DEPENDENCIES rotating_ends
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes rotating_ends)

message(VERBOSE "Leaving rotating_ends subdirectory")
# ------------------------------------------------------------------------------
