# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pipe subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pipe C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME pipe
  SOURCES pipe.cc
  LIBRARIES oomph::axisym_advection_diffusion
  CXX_STANDARD 14)

oomph_add_executable(
  NAME pipe_variable_diff
  SOURCES pipe_variable_diff.cc
  LIBRARIES oomph::axisym_advection_diffusion
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_advection_diffusion.pipe
  TARGET_DEPENDENCIES pipe pipe_variable_diff
  EXTRA_REQUIRES
  LABELS axisym_advection_diffusion pipe)

message(VERBOSE "Leaving pipe subdirectory")
# ------------------------------------------------------------------------------
