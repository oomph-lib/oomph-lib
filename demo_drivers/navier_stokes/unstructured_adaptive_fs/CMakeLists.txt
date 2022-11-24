# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_fs subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_adaptive_fs C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_bubble_in_channel
  SOURCES adaptive_bubble_in_channel.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::fluid_interface)

oomph_add_executable(
  NAME adaptive_drop_in_channel
  SOURCES adaptive_drop_in_channel.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_adaptive_fs
  TARGET_DEPENDENCIES adaptive_bubble_in_channel adaptive_drop_in_channel
  EXTRA_REQUIRES
  LABELS "navier_stokes;unstructured_adaptive_fs;adapt")

message(VERBOSE "Leaving unstructured_adaptive_fs subdirectory")
# ------------------------------------------------------------------------------
