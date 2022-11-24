# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered free_surface_rotation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(free_surface_rotation C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME free_surface_rotation
  SOURCES free_surface_rotation.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.free_surface_rotation
  TARGET_DEPENDENCIES free_surface_rotation
  EXTRA_REQUIRES
  LABELS navier_stokes free_surface_rotation
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving free_surface_rotation subdirectory")
# ------------------------------------------------------------------------------
